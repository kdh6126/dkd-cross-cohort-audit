# Three extensions: single-cell, confounder-aware selection, germline
Date: 2026-08-24 · scripts: `kpmp_singlecell.py`, `dkd_deconfound.py`, `gwas_layer.py`

Two of the three worked. One is a clean negative. All three are reported.

---

## 1. Single-cell — the artifact confirmed with its own negative control

The KPMP Atlas API returns cell counts and nothing else for the diabetic slice: every
`foldChange`, `pVal` and `pValAdj` is null. The downloaded h5ad has the data and much better
metadata, and it supports contrasts the API cannot:

| group | donors | cells | definition |
|---|---|---|---|
| DKD | **27** | 86,430 | `disease_category` CKD **and** `diabetes_history` Yes |
| CKD non-diabetic | 10 | 34,551 | CKD and diabetes No — the **confounder-matched** contrast |
| Healthy reference | 12 | 40,386 | `Healthy_reference_tissue` |
| Tumour nephrectomy | 6 | 43,809 | `Tumor_nephrectomy` |

Statistics are on **donor-level pseudobulk** (817 donor × cell-type profiles with ≥25 cells),
not per cell. A per-cell test treats thousands of cells from one patient as independent and
returns p-values that are mostly a function of cell count; with 27 versus 10 donors the donor
is the unit of replication.

### The result

Hedges' *g* per (gene, cell type), donor-level:

| Contrast | what it compares | IEG module mean *g* | \|g\|>0.8 | IEG vs random |
|---|---|---|---|---|
| DKD vs healthy reference | **biopsy vs non-biopsy** | **−0.843** | **59%** | **p = 6.0e-76** |
| DKD vs non-diabetic CKD | **biopsy vs biopsy** | +0.082 | 4% | p = 0.21 (n.s.) |
| reference vs nephrectomy | non-biopsy vs non-biopsy | −0.133 | 11% | p = 0.89 (n.s.) |

In KPMP, DKD tissue arrives as **biopsy**; both control groups are **non-biopsy**. The module
separates biopsy from non-biopsy at p = 6e-76 and separates nothing else. The middle row is the
negative control the earlier analyses lacked: two groups that are both biopsies, both diseased,
differing only in diabetes — and the module is indistinguishable from random genes there.

The procurement artifact is now shown in **three independent resources** and **three
platforms**: ERCB microarray (g = −1.23), GSE175759 RNA-seq (−2.26), KPMP single-nucleus
(−0.843 with its own negative control).

### And an honest negative for the candidates

| Contrast | candidates mean *g* | random mean *g* |
|---|---|---|
| DKD vs healthy reference | +0.301 (27% \|g\|>0.8) | +0.101 (10%) |
| DKD vs non-diabetic CKD | **+0.145 (5%)** | **+0.127 (4%)** |

Against healthy controls the candidate set is clearly enriched over random. Against the
confounder-matched contrast **it is not distinguishable from random**. Individual hits survive
(DPP6 in IMM g=+1.03 p=0.0024 and PT g=+0.92 p=0.0034; TPD52 VSM/P +1.23; MMP2 PT +0.99;
HDAC9 POD +1.07), but with ~420 uncorrected tests roughly two are expected by chance and we see
about five. Marginal, and it should be described that way.

DPP6 is the only candidate appearing twice with a consistent direction in the matched contrast.

---

## 2. Confounder-aware selection — the correction moved inside the objective

Until now the procurement correction was a post-hoc gate needing non-diabetic kidney disease
biopsies, available in **2 of 11** public datasets. The immediate-early module gives a
per-sample handling proxy computable from the same matrix with no extra labels, so the
correction can move inside the selector. Two variants, both scored against the uncorrected arm:

- **orth** — residualise every gene on the handling score within each cohort, select on residuals
- **penalty** — keep the data, multiply each gene's RBS by (1 − |corr(gene, handling score)|)

Mean over LODO folds, four cohorts:

| K | arm | procurement-flagged | specificity pass | external AUROC | cross-fold Jaccard |
|---|---|---|---|---|---|
| 10 | none | **55%** | 18% | 0.932 | 0.228 |
| 10 | **orth** | **5%** | **42%** | 0.860 | 0.204 |
| 10 | penalty | 10% | 35% | 0.750 | 0.126 |
| 20 | none | 44% | 25% | 0.926 | 0.241 |
| 20 | **orth** | **8%** | **42%** | **0.938** | 0.208 |
| 20 | penalty | 9% | 42% | 0.726 | 0.183 |
| 50 | none | 28% | 32% | 0.952 | 0.303 |
| 50 | **orth** | **6%** | **42%** | 0.936 | 0.248 |
| 50 | penalty | 11% | 38% | 0.934 | 0.277 |

**Residualising works.** At K=10 it cuts procurement-flagged genes from 55% to 5% and more than
doubles the specificity pass-rate, 18% → 42%. At every K the pass-rate lands at 42%, against
18–32% uncorrected.

**What it costs, stated plainly.** External AUROC falls 0.932 → 0.860 at K=10 (at K=20 it
actually rises, 0.926 → 0.938; at K=50 it falls 0.952 → 0.936), and cross-fold Jaccard falls
0.303 → 0.248 at K=50. The handling score is correlated with case status because cases are
biopsies, so residualising on it necessarily removes some real disease signal along with the
confounder. This is a trade, not a free lunch.

**The soft variant is worse.** `penalty` helps less on procurement (10–11% vs 5–8%) and costs far
more AUROC at small K (0.750 and 0.726). Demoting genes by correlation is not the same as
removing the axis.

**Why this matters more than the numbers.** The post-hoc gate needs a comparator disease cohort
procured identically — 2 of 11 datasets have one. The handling score needs nothing but the
expression matrix. That generalises the fix from a special case to any cohort.

---

## 3. Germline layer — a clean null

GCST90179152, diabetic kidney disease in a Korean cohort (2,532 cases / 31,347 controls,
PMID 36627639). 4,976,559 variants; the file carries p-values and positions only, so this is a
positional gene-level test and nothing more.

Two traps, both handled: min-p is biased by how many variants sit in a window, so the null is
built from background genes **matched on SNP count**; and LD makes variants non-independent,
which a matched empirical null absorbs to first order. This is not MAGMA and is not presented
as such.

| | |
|---|---|
| candidates with empirical p < 0.05 | **1 / 30** (expected by chance 1.5) |
| candidates reaching genome-wide significance | 0 |
| binomial test on the 0.05 tail | p = 0.785 |
| KS test of empirical p against uniform | D = 0.120, p = 0.736 |

**No germline enrichment near the candidate genes.** The empirical p-values are uniform. Best
ranked are TNNT2 (0.049) and C1orf21 (0.050), both Tier 4; PRSS23 is 0.094; **OLFML3 is 0.926**.

This is unsurprising and should not be spun. 2,532 cases is small, and expression-derived
candidates have no obligation to carry common germline risk — most GWAS signal for kidney
traits is regulatory and often distal to the gene it acts on. The honest reading is that this
layer adds no support and no contradiction.

---

## What the three extensions change

**Strengthened:** the procurement artifact is now a three-resource, three-platform result with a
proper negative control, and the correction for it has been moved from a post-hoc gate that
needs a rare comparator cohort into the selection objective, where it needs nothing extra.

**Weakened:** the candidate set. It is enriched over random against healthy controls but
indistinguishable from random in the confounder-matched single-cell contrast, and it has no
germline support. Combined with the earlier sign flips for OLFML3 and PRSS23 in GSE175759, the
biomarker claim is thinner than the methodology claim by a wide margin.

**Implication for the write-up:** lead with the method. The artifact and the correction are
replicated, controlled and generalisable. The candidates are a worked example that shows the
pipeline behaving correctly — including refusing to promote its own outputs.

## Reproduce

```
python scripts/kpmp_singlecell.py     # caches the extracted gene block on first run
python scripts/dkd_deconfound.py --B 200
python scripts/gwas_layer.py
python scripts/make_figures.py --only singlecell,deconfound
```
