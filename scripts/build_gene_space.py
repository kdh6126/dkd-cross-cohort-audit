#!/usr/bin/env python
"""Build the common gene space across DKD cohorts.

Hub ID = NCBI Entrez Gene ID (chosen because GSE104948/GSE104954 use Brainarray ENTREZG
custom CDF, i.e. their probe IDs already ARE Entrez IDs).

Per-cohort mapping route:
  GSE142025  HGNC symbol      -> Entrez   via NCBI Homo_sapiens.gene_info (symbol, then synonym)
  GSE30528   Affy GPL571      -> Entrez   via GEO GPL571.annot "Gene ID" column
  GSE30529   Affy GPL571      -> Entrez   (same)
  GSE96804   Affy GPL17586    -> Entrez   via GPL17586 platform table "gene_assignment" field
  GSE104948  ENTREZG CDF      -> Entrez   strip the "_at" suffix
  GSE104954  ENTREZG CDF      -> Entrez   (same)

Probes/transcript-clusters that map to more than one distinct Entrez gene are DROPPED
(standard conservative practice - a cross-mapping probe cannot be attributed to one gene).
Writes the mapping tables and the intersected gene list to data/processed/.
"""
import gzip, os, sys, json, csv, collections

ANN = 'data/raw/annotation'
OUT = 'data/processed'
os.makedirs(OUT, exist_ok=True)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- symbol -> entrez
def load_gene_info():
    """Return (symbol->entrez, synonym->entrez, entrez->symbol). Only current human genes."""
    sym, syn, rev = {}, {}, {}
    with gzip.open(os.path.join(ANN, 'Homo_sapiens.gene_info.gz'), 'rt',
                   encoding='utf-8', errors='replace') as f:
        hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip('\n').split('\t')
            gid = p[ix['GeneID']]
            s = p[ix['Symbol']]
            sym[s] = gid
            rev[gid] = s
            for a in p[ix['Synonyms']].split('|'):
                if a and a != '-' and a not in syn:
                    syn[a] = gid
    return sym, syn, rev


# ---------------------------------------------------------------- GPL571
def map_gpl571():
    """probe -> entrez, dropping probes that hit >1 gene."""
    out, dropped = {}, 0
    with gzip.open(os.path.join(ANN, 'GPL571.annot.gz'), 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('ID\t'):
                hdr = line.rstrip('\n').split('\t')
                ix = {c: i for i, c in enumerate(hdr)}
                break
        else:
            raise RuntimeError('GPL571 header not found')
        col = ix['Gene ID']
        for line in f:
            if line.startswith('!'):
                break
            p = line.rstrip('\n').split('\t')
            if len(p) <= col:
                continue
            raw = p[col].strip()
            if not raw:
                continue
            ids = {g for g in raw.split('///') if g.strip().isdigit()}
            if len(ids) == 1:
                out[p[0]] = ids.pop()
            elif len(ids) > 1:
                dropped += 1
    log('  GPL571: %d probes -> entrez, %d dropped (multi-gene)' % (len(out), dropped))
    return out


# ---------------------------------------------------------------- GPL17586
def map_gpl17586():
    """transcript cluster -> entrez, parsed from the gene_assignment field."""
    path = os.path.join(ANN, 'GPL17586.platform.txt')
    out, dropped, nocall = {}, 0, 0
    with open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('ID\t'):
                hdr = line.rstrip('\n').split('\t')
                ix = {c: i for i, c in enumerate(hdr)}
                break
        else:
            raise RuntimeError('GPL17586 header not found')
        col = ix['gene_assignment']
        for line in f:
            if line.startswith('!'):
                break
            p = line.rstrip('\n').split('\t')
            if len(p) <= col:
                continue
            ga = p[col]
            if not ga or ga == '---':
                nocall += 1
                continue
            ids = set()
            for rec in ga.split('///'):
                fld = [x.strip() for x in rec.split('//')]
                # accession // symbol // description // cytoband // entrez
                if len(fld) >= 5 and fld[4].isdigit():
                    ids.add(fld[4])
            if len(ids) == 1:
                out[p[0]] = ids.pop()
            elif len(ids) > 1:
                dropped += 1
            else:
                nocall += 1
    log('  GPL17586: %d clusters -> entrez, %d dropped (multi-gene), %d no gene call'
        % (len(out), dropped, nocall))
    return out


# ---------------------------------------------------------------- matrix row ids
def matrix_ids(path):
    ids, started = [], False
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                started = True
                continue
            if not started:
                continue
            if line.startswith('!series_matrix_table_end'):
                break
            ids.append(line.split('\t')[0].strip('"'))
    return ids[1:]          # first row is the sample-id header


def tsv_ids(path):
    with open(path, 'rt', encoding='utf-8', errors='replace') as f:
        f.readline()
        return [l.split('\t')[0] for l in f]


# ---------------------------------------------------------------- per cohort
def main():
    log('loading annotation ...')
    sym2e, syn2e, e2sym = load_gene_info()
    log('  gene_info: %d symbols, %d synonyms' % (len(sym2e), len(syn2e)))
    g571 = map_gpl571()
    g17586 = map_gpl17586()

    cohorts = {}

    # --- GSE142025 : HGNC symbols
    ids = tsv_ids(os.path.join(OUT, 'GSE142025_expr_matrix.tsv'))
    m, viasyn, unmapped = {}, 0, 0
    for i in ids:
        if i in sym2e:
            m[i] = sym2e[i]
        elif i in syn2e:
            m[i] = syn2e[i]
            viasyn += 1
        else:
            unmapped += 1
    log('GSE142025: %d/%d mapped (%d via synonym, %d unmapped)' % (len(m), len(ids), viasyn, unmapped))
    cohorts['GSE142025'] = m

    # --- GPL571 cohorts
    for gse in ('GSE30528', 'GSE30529'):
        ids = matrix_ids('data/raw/geo/%s/%s_series_matrix.txt.gz' % (gse, gse))
        m = {i: g571[i] for i in ids if i in g571}
        log('%s: %d/%d probes mapped' % (gse, len(m), len(ids)))
        cohorts[gse] = m

    # --- GSE96804
    ids = matrix_ids('data/raw/geo/GSE96804/GSE96804_series_matrix.txt.gz')
    m = {i: g17586[i] for i in ids if i in g17586}
    log('GSE96804: %d/%d clusters mapped' % (len(m), len(ids)))
    cohorts['GSE96804'] = m

    # --- ENTREZG CDF cohorts (id is already the Entrez gene, "<gid>_at")
    for gse, plat in (('GSE104948', 'GPL22945'), ('GSE104954', 'GPL22945')):
        ids = matrix_ids('data/raw/geo/%s/%s-%s_series_matrix.txt.gz' % (gse, gse, plat))
        m = {}
        for i in ids:
            core = i[:-3] if i.endswith('_at') else i
            if core.isdigit():
                m[i] = core
        log('%s: %d/%d probes are Entrez-native' % (gse, len(m), len(ids)))
        cohorts[gse] = m

    # ---------------------------------------------------------------- write maps
    with open(os.path.join(OUT, 'probe2gene.tsv'), 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['cohort', 'probe_id', 'entrez_id', 'symbol'])
        for gse, m in cohorts.items():
            for p, g in m.items():
                w.writerow([gse, p, g, e2sym.get(g, '')])

    genesets = {gse: set(m.values()) for gse, m in cohorts.items()}
    log('\n=== distinct Entrez genes per cohort ===')
    for gse, s in genesets.items():
        log('  %-11s %6d genes  (from %d probes)' % (gse, len(s), len(cohorts[gse])))

    DISCOVERY = ['GSE142025', 'GSE30528', 'GSE96804', 'GSE104948']
    ALL6 = DISCOVERY + ['GSE30529', 'GSE104954']

    inter_disc = set.intersection(*[genesets[g] for g in DISCOVERY])
    inter_all = set.intersection(*[genesets[g] for g in ALL6])
    log('\n=== INTERSECTION ===')
    log('  4 discovery cohorts (%s): %d genes' % (', '.join(DISCOVERY), len(inter_disc)))
    log('  all 6 incl. compartment validation: %d genes' % len(inter_all))

    log('\n=== leave-one-out: genes lost by including each cohort ===')
    for g in DISCOVERY:
        others = [x for x in DISCOVERY if x != g]
        without = set.intersection(*[genesets[x] for x in others])
        log('  drop %-11s -> %6d genes (+%d)' % (g, len(without), len(without) - len(inter_disc)))

    for name, s in (('discovery4', inter_disc), ('all6', inter_all)):
        with open(os.path.join(OUT, 'gene_space_%s.tsv' % name), 'w', encoding='utf-8', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['entrez_id', 'symbol'])
            for g in sorted(s, key=int):
                w.writerow([g, e2sym.get(g, '')])
    log('\nwrote %s/gene_space_discovery4.tsv and gene_space_all6.tsv' % OUT)

    json.dump({g: sorted(s, key=int) for g, s in genesets.items()},
              open(os.path.join(OUT, 'genesets_per_cohort.json'), 'w'), indent=0)


if __name__ == '__main__':
    main()
