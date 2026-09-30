#!/usr/bin/env python
"""KPMP Atlas GraphQL client - per-gene cross-omics evidence. Open access, no DAR required.

Enrollment categories used by the API:
  hrt = Healthy Reference Tissue   dmr = Diabetes Mellitus with Renal disease (= DKD)
  ckd = Chronic Kidney Disease     aki = Acute Kidney Injury     all = pooled

IMPORTANT quirks discovered by probing the live endpoint (2026-08-24):
  * The resolver does NOT bind GraphQL variables. Arguments must be inlined into the
    query string, exactly as the Atlas frontend does.
  * enrollmentCategory must be a real value ("all", "dmr", ...). An empty string returns 0 rows.
  * getRTGeneExpressionByEnrollment always returns null - it appears dead. Use
    regional_transcriptomics_by_structure() instead.
  * getRPGeneExpressionByEnrollment returns a LIST (one entry per protein accession),
    each wrapping rpExpressionByEnrollmentCategory.all.
  * Regional transcriptomics (rt) and regional proteomics (rp) have ZERO dmr participants.
    They contrast Glomerulus vs Tubulointerstitium, not disease vs control. Only the
    single-cell / single-nucleus layers carry DKD-specific contrasts.
"""
import urllib.request, json, sys, time

ENDPOINT = "https://atlas.kpmp.org/graphql"


def _esc(s):
    """Escape for embedding inside a GraphQL double-quoted string literal."""
    return json.dumps(s or '')[1:-1]


def gql(query, retries=3):
    body = json.dumps({"query": query}).encode()
    req = urllib.request.Request(
        ENDPOINT, data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    last = None
    for _ in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                d = json.loads(r.read().decode('utf-8', 'replace'))
            if 'errors' in d:
                raise RuntimeError(d['errors'])
            return d['data']
        except Exception as e:
            last = e
            time.sleep(2)
    raise last


RT_FIELDS = ("id geneSymbol segment segmentName foldChange pVal adjPVal pValLog10 "
             "stdDev enrollmentCategory sampleCount")
RP_FIELDS = ("id geneSymbol accession description comparison region foldChange adjPVal "
             "enrollmentCategory sampleCount coveragePct numPeptides numUniquePeptides fdrConfidence")
SC_FIELDS = ("gene enrollmentCategory dataType cluster clusterName foldChange pVal pValAdj "
             "avgExp pct1 pct2 cellCount")


def data_types_for_gene(gene):
    """Which KPMP omics layers carry evidence for this gene: subset of sc / sn / rt / rp."""
    return gql('query { dataTypesForConcept2025(geneSymbol:"%s") }' % _esc(gene))["dataTypesForConcept2025"]


def single_cell(gene, data_type="sn", enrollment="all", cell_type=""):
    """Per-cluster differential expression. data_type: 'sn' (single-nucleus) or 'sc' (single-cell).
    enrollment: 'all' | 'dmr' | 'ckd' | 'aki' | 'hrt'."""
    q = ('query { geneExpressionSummary2025(dataType:"%s", geneSymbol:"%s", cellType:"%s", '
         'enrollmentCategory:"%s") { %s } }') % (
        _esc(data_type), _esc(gene), _esc(cell_type), _esc(enrollment), SC_FIELDS)
    return gql(q)["geneExpressionSummary2025"] or []


def regional_proteomics(gene):
    """Laser-microdissected regional PROTEOMICS. Returns flat rows.
    Contrast is Glomerulus vs Tubulointerstitium (comparison 'G.vs.TI.in.All')."""
    q = ('query { getRPGeneExpressionByEnrollment(geneSymbol:"%s") { accession '
         'rpExpressionByEnrollmentCategory { all{%s} } } }') % (_esc(gene), RP_FIELDS)
    out = []
    for entry in gql(q)["getRPGeneExpressionByEnrollment"] or []:
        out += (entry.get('rpExpressionByEnrollmentCategory') or {}).get('all') or []
    return out


def regional_transcriptomics_by_structure(structure="Glomerulus"):
    """Laser-microdissected regional TRANSCRIPTOMICS for a whole structure (all genes).
    Known structures include 'Glomerulus'. Returns ~26k rows - cache it."""
    q = 'query { getRTGeneExpressionByStructure(structure:"%s") { %s } }' % (_esc(structure), RT_FIELDS)
    return gql(q)["getRTGeneExpressionByStructure"] or []


def data_type_summary():
    """Participant counts per data type and enrollment category."""
    q = ('query { getDataTypeSummaryInformation2025 { omicsType dataType dataTypeShort '
         'akiCount ckdCount hrtCount dmrCount totalCount participantCount } }')
    return gql(q)["getDataTypeSummaryInformation2025"]


def gene_evidence(gene):
    """One-call cross-omics evidence bundle for a candidate biomarker."""
    ev = {"gene": gene, "layers": data_types_for_gene(gene)}
    ev["sn_dmr"] = single_cell(gene, "sn", "dmr")
    ev["sn_hrt"] = single_cell(gene, "sn", "hrt")
    ev["sc_dmr"] = single_cell(gene, "sc", "dmr")
    ev["rp"] = regional_proteomics(gene)
    return ev


if __name__ == "__main__":
    g = sys.argv[1] if len(sys.argv) > 1 else "CCL2"
    ev = gene_evidence(g)
    print("### %s ###  layers=%s" % (g, ev["layers"]))
    for key in ("sn_dmr", "sn_hrt", "sc_dmr"):
        rows = [r for r in ev[key] if r.get('pValAdj') is not None and r['pValAdj'] < 0.05]
        rows.sort(key=lambda x: -(x.get('foldChange') or 0))
        print("-- %s: %d clusters, %d with adjP<0.05" % (key, len(ev[key]), len(rows)))
        for r in rows[:3]:
            print("     %-46s FC=%7.3f adjP=%.2e cells=%s" % (
                (r['clusterName'] or '')[:46], r['foldChange'], r['pValAdj'], r['cellCount']))
    print("-- regional proteomics: %d rows" % len(ev["rp"]))
    for r in ev["rp"][:4]:
        print("     %-6s %-18s FC=%7.3f adjP=%.3g n=%s peptides=%s" % (
            r['region'], r['comparison'], r['foldChange'], r['adjPVal'], r['sampleCount'], r['numPeptides']))
