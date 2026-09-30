# KPMP — what is actually obtainable, and the one limitation that matters
Date: 2026-08-24

## 1. Two access routes, both open, no DAR required

### (a) Bulk atlas files — downloaded
The KPMP kidney atlas (Lake et al.) is mirrored on CELLxGENE Discover and downloadable directly.

| File | Cells | Size | Local path |
|---|---|---|---|
| scRNA-seq, adult human kidney v1.5 | 225,177 | 1.82 GB | `data/raw/kpmp/KPMP_v1.5_scRNA_human_kidney.h5ad` |
| snRNA-seq, adult human kidney v1.5 | 304,989 | 3.03 GB | `data/raw/kpmp/KPMP_v1.5_snRNA_human_kidney.h5ad` |

A v2.0 exists (preprint doi 10.1101/2025.09.26.678707): human snRNA 1,388,643 nuclei / 13.8 GB
and scRNA 348,984 / 4.4 GB. Not downloaded — v1.5 is the published, citable version and is
sufficient for candidate-gene corroboration. Fetch v2.0 only if cell-state resolution turns out
to matter.

### (b) Per-gene query API — decoded and wrapped
`https://atlas.kpmp.org/graphql`, wrapped in `scripts/kpmp_client.py`.

```python
from kpmp_client import gene_evidence
ev = gene_evidence("NPHS2")     # layers, sn/sc per-cluster DE, regional proteomics
```

Non-obvious quirks, all handled in the client:
- The resolver **does not bind GraphQL variables** — arguments must be inlined into the query
  string, exactly as the Atlas frontend does. Using variables silently returns `null`.
- `enrollmentCategory` must be a real value (`"all"`, `"dmr"`, …). An empty string returns 0 rows.
- `getRTGeneExpressionByEnrollment` always returns `null` — it appears dead. Use
  `getRTGeneExpressionByStructure(structure:"Glomerulus")` instead (~26k gene rows).
- `getRPGeneExpressionByEnrollment` returns a **list**, one entry per protein accession, each
  wrapping `rpExpressionByEnrollmentCategory.all`.

Verified end to end: `NPHS2` → Podocyte cluster *g*-equivalent FC = 6.63 (adjP ≈ 0) and regional
proteomics Glomerulus FC = 6.84 (adjP = 3.8e-102). Correct biology, so the client is sound.

## 2. The limitation that changes the plan

Participant counts per data type and enrollment category, straight from
`getDataTypeSummaryInformation2025`:

| Data type | AKI | CKD | Healthy ref | **DKD (`dmr`)** | Total |
|---|---|---|---|---|---|
| Single-cell RNA-seq (sc) | 19 | 51 | 40 | **3** | 113 |
| Single-nucleus RNA-seq (sn) | 33 | 72 | 40 | **11** | 156 |
| Regional transcriptomics (rt) | 5 | 22 | 9 | **0** | 36 |
| Regional proteomics (rp) | 12 | 14 | 5 | **0** | 31 |

And, separately: for the `dmr` slice the API returns **cell counts only — every expression
statistic (`foldChange`, `pVal`, `pValAdj`, `avgExp`) is null**. The `all` / `hrt` / `ckd` / `aki`
slices all carry full statistics; `dmr` does not.

**Consequences:**
1. KPMP **regional proteomics has zero diabetic participants**. There is no DKD-vs-control
   protein contrast to be had from KPMP. What regional proteomics *does* give is a
   Glomerulus-vs-Tubulointerstitium contrast (`comparison: "G.vs.TI.in.All"`) — which is still
   directly useful, just for a different question.
2. There is **no metabolomics layer** exposed in the Atlas at all.
3. DKD-specific differential expression cannot be *fetched*; it must be **computed from the
   downloaded h5ad** (which does carry per-cell disease annotation). This is the reason the
   bulk files were downloaded rather than relying on the API.

## 3. What KPMP can therefore honestly contribute

| Question | Available? | Route |
|---|---|---|
| Is candidate gene X expressed in a DKD-relevant cell type? | yes | API, `enrollmentCategory="all"` |
| Is X glomerular or tubulointerstitial — at RNA level? | yes | `getRTGeneExpressionByStructure` |
| Is X glomerular or tubulointerstitial — **at protein level**? | yes | regional proteomics, `G.vs.TI` |
| Is X differentially expressed in DKD vs healthy, per cell type? | **compute it** | from the v1.5 snRNA h5ad (11 DKD participants) |
| Is X differentially *abundant* in DKD? | partial | `dmr` cell counts only |
| Metabolite evidence for X's pathway | **no** | not in KPMP Atlas — use Metabolomics Workbench instead |

The compartment-attribution rows are the genuinely valuable part, because they answer the
glomerulus-vs-tubulointerstitium question **at the protein level** — independent evidence for
the compartment analysis, from a different assay and a different cohort.

Framing this as "cross-omics biological corroboration" rather than "multi-omics validation"
was the right call: with 11 DKD participants at the RNA level and 0 at the protein level, KPMP
cannot carry a validation claim, but it can corroborate compartment and cell-type attribution.
