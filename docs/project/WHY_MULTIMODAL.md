# Why multi-modal data — the three rungs, measured
Date: 2026-08-24 · all numbers from `results/`, reproducible from `scripts/`

The usual argument for multi-omics is that more modalities give more information. That is true
but unfalsifiable as stated. This project accidentally produced a sharper version of the
argument, because it went up the ladder one rung at a time and measured each one.

> **Within a single modality, reproducibility saturates and then becomes actively misleading.
> The missing ingredient is not more cohorts — it is a second modality on matched patients.
> Public data does not provide that.**

---

## Rung 1 — one dataset at a time

`scripts/per_dataset_features.py`. Eleven datasets, each run on its own with bootstrap
stability selection (B=60, four selectors, consensus). Nine could never join the cross-cohort
rotation; running them anyway is the point.

| Dataset | contrast | case/ctrl | within-dataset bootstrap Jaccard | note |
|---|---|---|---|---|
| GSE142025 | DKD vs control | 27 / 9 | 0.217 |
| GSE30528 | DKD vs control (glom) | 9 / 13 | 0.166 |
| GSE96804 | DKD vs control (glom) | 41 / 20 | 0.191 |
| GSE104948 | DKD vs control (glom) | 12 / 26 | 0.132 |
| GSE30529 | DKD vs control (tubulo) | 10 / 12 | 0.166 |
| GSE104954 | DKD vs control (tubulo) | 17 / 26 | 0.155 |
| GSE1009 | DN vs control | 3 / 3 | 0.121 | *only 7,006 of 9,900 genes map on GPL8300*
| GSE111154 | early DN vs control | 4 / 4 | 0.133 |
| GSE20602 | nephrosclerosis vs control | 13 / 5 | 0.124 |
| GSE142153 | DN vs healthy (PBMC) | 23 / 10 | 0.134 |
| GSE175759 | any kidney disease vs nephrectomy | 65 / 22 | 0.257 |

**Mean pairwise agreement of the top-50 sets across all eleven datasets: 0.049.**

Structure worth noting in `results/per_dataset/cross_dataset_jaccard.tsv`:
- The three tiny/odd datasets (GSE1009, GSE111154, GSE20602) agree with **each other** far more
  (0.190–0.220) than with any real cohort (0.000–0.020). Agreement is not evidence of biology.
- GSE142153 (blood) agrees with nothing renal (0.000–0.075). Different tissue, no transfer.
- GSE175759's own top features are FOS, DUSP1, NR4A1, EGR1, NR4A2, FOSB, ZFP36 — its contrast
  is literally biopsy-versus-nephrectomy, so single-dataset selection returns the procurement
  module in its purest form. A nice accidental positive control.

**Rung 1 verdict:** a single dataset yields a signature that does not transfer. Nothing here is
surprising, but it establishes the baseline the next rung has to beat.

---

## Rung 2 — many datasets, one modality

Two things happen, and only one of them is good.

### The good half: reproducibility genuinely improves

Cross-cohort stability selection (RBS) versus the seven baselines, K=50:

| | baselines (cross-fold Jaccard) | RBS | external AUROC |
|---|---|---|---|
| 3 glomerular cohorts | 0.032 – 0.050 | **0.192** | tie (0.851 vs 0.863 best baseline) |
| 4 cohorts | 0.136 – 0.196 | **0.303** | **best** (0.952, +1.33 sd over null) |

And the dominant factor is simply how many training cohorts each fold has: 2 → 3 moves
cross-fold Jaccard from ~0.18 to ~0.30 regardless of algorithm. So yes, more cohorts of the
same modality buys reproducibility.

### The bad half: reproducibility is not validity, and more cohorts make it worse

Take the naive version of "use many datasets": count how often a gene lands in a dataset's
top-50. Across eleven datasets, 20 genes reach ≥4 (chance expectation: 0.002 genes, so the
recurrence is real). Then check them against the confounder-matched contrast:

| Gene | in N datasets | g vs control | g vs **other CKD** | g other-CKD vs control | gate |
|---|---|---|---|---|---|
| ADA | 6 | +1.82 | **+0.19** | +1.07 | fail |
| DUSP1 | 4 | −2.21 | **+0.03** | −2.16 | procurement |
| ZFP36 | 4 | −2.02 | **−0.07** | −2.14 | procurement |
| MICA | 5 | +0.85 | +0.23 | +0.55 | fail |
| PPP4R1 | 4 | +0.71 | +0.08 | +0.60 | procurement |
| A2M | 5 | +1.60 | +0.55 | +0.95 | **pass** |
| CDH1 | 4 | +0.16 | +0.54 | −0.39 | **pass** |

**2 of 20 pass. 5 are procurement-flagged.** The pattern is unmistakable: strong against
controls, near-zero against other kidney disease, and the non-diabetic contrast carries almost
all of it. What recurs across datasets of one modality is what those datasets *share* — and what
they share is the case/control design (biopsy versus nephrectomy) and generic kidney damage,
not the disease.

### And the confounder replicates independently

The immediate-early module, non-diabetic CKD versus nephrectomy controls:

| Cohort | platform | compartment | module mean *g* | rest of transcriptome |
|---|---|---|---|---|
| ERCB (GSE104948/104954) | microarray | glom + tubulo | −1.23 | ~0 |
| **GSE175759** | **RNA-seq** | tubulo | **−2.26** | +0.69 |

Different centre, different platform, and the effect is *larger*
(Mann-Whitney module vs rest, p = 5.4e-14). This was the project's stated weakness — the
confounder-matched contrast existed only in ERCB — and it is now closed.

**Rung 2 verdict:** stacking cohorts of the same modality raises reproducibility and, at the
same time, reliably promotes a shared confounder to the top of the ranking. No amount of
cross-cohort stability detects it, because the confounder is equally reproducible. Escaping it
required a contrast that holds procurement constant — available in only two of eleven datasets.

---

## Rung 3 — a second modality, and why it could not be built here

What is actually available in public data, and what each layer delivers:

| Layer | Source | Status | What it gives | What it cannot give |
|---|---|---|---|---|
| Protein | KPMP regional proteomics | usable | 22/30 candidates have protein evidence; glomerulus-vs-tubulointerstitium attribution (FBN1 +2.03 adjP 8e-19, DPP6 +3.97 adjP 2.5e-13) | **0 diabetic participants** — no DKD-vs-control protein contrast exists |
| Single cell | KPMP v1.5 (4.85 GB), GSE131882 | usable | cell-type attribution: fibroblast for FMOD/LUM/MMP2/COL1A2/OLFML3, LYVE1+ macrophage for MS4A6A/FOLR2 | `dmr` slice returns cell counts with all statistics null; 11 DKD participants |
| Blood transcriptome | GSE142153 | usable | a genuinely different tissue, 23 DN vs 10 healthy | agreement with kidney datasets 0.000–0.075 — no feature transfer |
| Metabolomics | ST004483, ST004442 (2025, human DKD urine/plasma, n=40 each) | **raw only** | — | deposited as mzXML; no processed metabolite matrix, no data table via REST |
| Metabolomics | ST000691 (human, n=32) | **no contrast** | — | every sample is `Progression:Y` plus 2 pooled QC |
| GWAS | GCST90018612/832, GCST90179152, … | not downloaded | germline risk | not an expression contrast |

So of four public human DKD metabolomics studies located, **none yields an analysable
case/control matrix**, and the one modality with a clean protein-level assay has zero diabetic
subjects. The rung cannot be climbed with these data.

**That is the argument, and it is specific.** The gap is not "more omics would be nice". It is:

> No public resource provides a second modality measured on the *same* DKD patients with a
> case/control design. Every public modality either lacks diabetic subjects (KPMP proteomics),
> lacks processed data (2025 metabolomics), lacks a contrast (ST000691), or lacks feature
> transfer because the tissue differs (blood).

A confounder shared by every cohort of one modality is invisible to any single-modality method,
however stable. Two ways out exist: a confounder-matched contrast within the modality — which
worked here but was available in only two of eleven datasets and has no independent replication
for the *specificity* claim — or an orthogonal modality on matched patients, which does not
exist publicly. The second is the one that scales.

---

## What this implies for the study design

1. **Report reproducibility and validity separately.** They came apart here, in a measurable
   way. A stability metric without a confounder-matched contrast is not a validity claim.
2. **The confounder-matched contrast is the cheapest fix and should be a design requirement.**
   Non-diabetic kidney disease biopsies procured identically to the DKD biopsies did more work
   than any algorithm in this project. Two of eleven public datasets have them.
3. **For new data generation, matched modalities on the same patients beat more patients in one
   modality.** Rung 2 shows the returns to more same-modality cohorts are real but bounded, and
   that they come with a systematic failure mode.
4. **Procurement must be recorded and balanced.** In all eleven datasets cases are biopsies and
   controls are nephrectomy or donor tissue. That single design choice generated the strongest
   "reproducible" signal in the entire analysis.

## Reproduce

```
python scripts/per_dataset_features.py          # rung 1 + naive recurrence
python scripts/compare_methods.py               # rung 2, cross-cohort gain
python scripts/validate_gse175759.py            # rung 2, independent artifact replication
python scripts/kpmp_corroborate.py              # rung 3, protein + cell type
```
