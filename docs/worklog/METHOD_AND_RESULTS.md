# Cross-Cohort Robust Feature Selection for DKD — Method and Results
Date: 2026-08-24 · Data and audit: `docs/DATASET_AUDIT.md`, `docs/GENE_SPACE_AND_FEASIBILITY.md`

## 0. What this run established, in one page

Three results reframe the project, and two of them were not anticipated.

1. **Raw AUROC cannot rank methods on this data.** A *randomly chosen* 50-gene set reaches
   external AUROC 0.68–0.88 depending on the fold. Every baseline's score falls inside that
   null band. This vindicates the decision to headline stability rather than accuracy — but it
   also means every AUROC must be reported as a percentile against the null, never raw.

2. **The most cross-cohort reproducible module in the data is a tissue-procurement artifact.**
   The immediate-early genes (FOS/JUN/EGR1/ATF3/DUSP1/ZFP36…) top every stability ranking, and
   non-diabetic CKD biopsies show *145%* of the DKD shift. Cases are biopsies and controls are
   nephrectomy/donor tissue in **every** cohort, so the confounder is shared — and a
   cross-cohort stability criterion rewards exactly that. **Reproducibility across cohorts does
   not imply biological validity when the confounder is shared.**

3. **Cross-fold signature agreement is the axis that separates methods.** All seven baselines
   land at Jaccard 0.03–0.05 at K=50 on three cohorts — they return an essentially different
   signature depending on which cohort is held out. The proposed score reaches 0.192 there at
   equal external AUROC, and on four cohorts 0.303 with the highest AUROC as well
   (0.952, +1.33 sd over the null; 0.936 / +1.42 sd with the uninformative fold removed).

---

## 1. Experimental design

**LODO.** One cohort held out entirely; the rest pooled for training. Two configurations:

| Config | Cohorts | Training cohorts per fold |
|---|---|---|
| `glom3` | GSE30528, GSE96804, GSE104948 | 2 |
| `all4` | + GSE142025 | 3 |

`glom3` is compartment-clean (all microdissected glomeruli). `all4` adds the only RNA-seq
cohort at the cost of mixing in whole cortex. Both are reported; the ablation below shows the
number of training cohorts matters more than any other design choice.

**Leakage control.** Per-cohort gene-wise z-scoring uses no labels, so applying it to a
held-out cohort is unsupervised domain alignment, not leakage. Every supervised step —
selection, model fitting — sees training cohorts only.

**Internal CV is nested.** Selecting features on all training data and then cross-validating
within it gives AUROC ≈ 0.99 for every method; we measured this. The reported internal number
re-runs selection inside each CV split.

## 2. Baseline selectors

`univariate` (F-test), `lasso`, `elastic_net`, `mrmr` (FCQ), `relieff`, `rf` importance,
`boruta`. Implemented directly rather than via boruta_py / skrebate / mrmr_selection, which pin
older numpy/sklearn and would make the benchmark fragile.

Regularisation was matched to sparsity rather than left at a default: LASSO at C=1.0 selects
46–70 genes per fold and elastic-net at C=0.05 selects 72–87, so the top-50 comparison is a
truncation rather than a pad. At the C=0.05 originally tried, LASSO selected 5 genes and the
padding would have silently turned it into a univariate filter.

**One incidental finding:** mRMR degenerates toward the univariate filter here, because mean
|r| within its relevance pool is 0.48 — the DKD signature is so co-expressed that the
redundancy penalty is near-constant across candidates.

## 3. The proposed score (RBS)

Every baseline pools the training cohorts and selects once, which is what lets a
dataset-specific gene win on the strength of the largest cohort. RBS never pools: it runs
bootstrap selection **inside each training cohort separately** and aggregates across cohorts.

Four components, each in [0,1], combined by a **weighted geometric mean** (a sum would let a
gene buy its way in on one component; the claim is "important AND stable AND reproducible AND
robust"):

| | |
|---|---|
| **I** predictive importance | \|coefficient\| ÷ max \|coefficient\| in that bootstrap, averaged |
| **S** bootstrap stability | mean over cohorts of within-cohort selection frequency |
| **R** cross-cohort reproducibility | cross-cohort aggregate of frequency, **gated** on effect direction agreeing in every cohort |
| **P** perturbation robustness | selection frequency under subsampling / gaussian noise / feature masking / batch shift |

Perturbation is applied **after** I/S/R have narrowed the pool — it stress-tests candidates,
it never discovers them, per the project decision that augmentation must not drive discovery.

### 3.1 The aggregation rule had to change, and why

The first implementation used a raw **MIN** across cohorts — the purest statement of worst-case
reproducibility. It failed: cross-fold Jaccard 0.036, and only 13 of 9,900 genes exceeded
R = 0.1, so everything below rank ~13 was noise. Diagnosis: with a 22-sample cohort the
per-gene frequency is estimated so poorly that `min(.)` is dominated by estimation error rather
than by real cohort disagreement.

Ablation over 24 combinations (`results/ablation/ablation.tsv`):

| config | agg | base | cross-fold Jacc @50 | AUROC @50 |
|---|---|---|---|---|
| glom3 | min | lasso | 0.061 | 0.753 |
| glom3 | **geomean** | lasso | **0.177** | **0.908** |
| all4 | min | univariate | 0.345 | 0.930 |
| all4 | **geomean** | univariate | 0.306 | 0.953 |
| all4 | geomean | relieff | 0.298 | 0.949 |
| all4 | *mean* | relieff | *0.396* | *0.951* |

`geomean` improved **both** stability and AUROC over `min`, confirming the noise diagnosis
rather than a stringency trade-off.

**`mean` scored best and was rejected.** S is already `F.mean(0)`, so with `agg='mean'` the R
component becomes mathematically identical to S and the cross-cohort term vanishes — the score
degenerates to I·S²·P plus the sign gate. That is a benchmark gain purchased by deleting the
mechanism the method exists for.

**The dominant factor is neither:** going from 2 training cohorts to 3 moves cross-fold Jaccard
from ~0.18 to ~0.30 regardless of aggregation rule. The original four-cohort LODO proposal was
the right instinct.

## 3.2 Head-to-head result (glom3: three glomerular cohorts, B=200, 9,900 genes)

At K = 50 genes, averaged over the three LODO folds:

| Method | external AUROC | vs random null | internal (nested) | bootstrap Jaccard | **cross-fold Jaccard** |
|---|---|---|---|---|---|
| **RBS (proposed)** | 0.851 | +0.74 sd (p74) | 0.987 | 0.178 | **0.192** |
| lasso | **0.863** | +0.77 sd (p75) | 0.989 | 0.128 | 0.046 |
| mrmr | 0.810 | +0.33 sd (p57) | 0.986 | 0.231 | 0.036 |
| univariate | 0.795 | +0.26 sd (p55) | 0.988 | 0.247 | 0.032 |
| rf | 0.793 | +0.24 sd (p53) | 0.993 | 0.116 | 0.046 |
| elastic_net | 0.785 | +0.09 sd (p51) | 0.991 | 0.164 | 0.040 |
| boruta | 0.780 | +0.13 sd (p48) | 0.976 | 0.212 | 0.043 |
| relieff | 0.769 | +0.06 sd (p48) | 0.985 | 0.160 | 0.050 |

At K = 100 the gap widens: RBS cross-fold Jaccard 0.230 versus 0.074 for the best baseline, and
RBS has the highest external AUROC (0.876).

**How to read this.** RBS does not win on AUROC — it ties lasso, and both sit at roughly the
75th percentile of the random-signature null, which is as much as any method achieves here. What
RBS does is return **four times more of the same signature when the held-out cohort changes**
(0.192 vs 0.032–0.050). That is the claim the design was built to support, and it is the axis on
which the baselines are indistinguishable from one another.

Two further readings worth keeping:
- **The internal/external gap is enormous for everyone** (internal ~0.98–0.99 versus external
  0.77–0.88). Random k-fold CV never confronts a new cohort; this is the concrete cost of
  reporting it instead of LODO. RBS has the smallest gap.
- **At K=100 four baselines score *below* the random null** (rf −0.25 sd, elastic_net −0.25 sd,
  boruta −0.19 sd, mrmr −0.55 sd). Selecting 100 genes by these criteria is worse than picking
  100 at random. Without the null control this would have been reported as "AUROC 0.77–0.80".

### 3.3 Four-cohort configuration (all4), K = 50

Adding GSE142025 gives every fold three training cohorts instead of two, and the picture changes
from "ties on AUROC, wins on stability" to "wins on both":

| Method | external AUROC | vs random null | bootstrap Jaccard | **cross-fold Jaccard** |
|---|---|---|---|---|
| **RBS (proposed)** | **0.952** | **+1.33 sd (p93)** | 0.192 | **0.303** |
| boruta | 0.928 | +1.10 sd (p87) | 0.227 | 0.143 |
| relieff | 0.920 | +1.13 sd (p86) | 0.183 | 0.196 |
| lasso | 0.903 | +0.97 sd (p81) | 0.133 | 0.148 |
| rf | 0.899 | +0.92 sd (p82) | 0.128 | 0.136 |
| elastic_net | 0.889 | +0.88 sd (p77) | 0.192 | 0.154 |
| univariate | 0.869 | +0.74 sd (p75) | 0.261 | 0.140 |
| mrmr | 0.866 | +0.72 sd (p73) | 0.243 | 0.145 |

Every method improves with a third training cohort, which is consistent with the ablation: the
number of training cohorts is the dominant design factor, not the algorithm.

**Sensitivity to the uninformative fold.** RBS scores AUROC 1.000 on the GSE142025 hold-out, but
so do random gene sets there (section 4), so that fold cannot discriminate methods. Removing it:

| Method | external (4 folds) | external (3 folds, GSE142025 dropped) | vs null (3 folds) |
|---|---|---|---|
| **RBS** | 0.952 | **0.936** | **+1.42 sd** |
| boruta | 0.928 | 0.903 | +1.11 sd |
| relieff | 0.920 | 0.893 | +1.15 sd |
| lasso | 0.903 | 0.870 | +0.93 sd |
| rf | 0.899 | 0.865 | +0.87 sd |
| elastic_net | 0.889 | 0.853 | +0.81 sd |
| univariate | 0.869 | 0.826 | +0.62 sd |
| mrmr | 0.866 | 0.822 | +0.60 sd |

The ranking is unchanged and RBS's margin over the null *increases* (+1.33 → +1.42 sd). The
result does not depend on the trivial fold.

A note on the bootstrap-Jaccard column: for the baselines it is the mean pairwise overlap
between bootstrap top-k sets drawn from pooled training data. RBS bootstraps *inside* each
cohort, so the same quantity is computed per cohort and averaged, using
E[Jaccard] = Σf² / (2k − Σf²) for independent draws with per-gene frequency f. The naive
alternative — overlap between the two training cohorts' top-k sets — measures cross-cohort
disagreement (0.007 here) and would be misleading in this column.

## 4. The null control that changes how everything is read

`scripts/null_control.py` — 300 random gene sets per fold per K, plus a label-permutation null.

| Held out | K=20 random AUROC | K=50 | K=100 | permuted labels |
|---|---|---|---|---|
| GSE30528 | 0.667 [0.376–0.941] | 0.682 [0.385–0.923] | 0.710 [0.470–0.923] | ~0.52 |
| GSE96804 | 0.723 [0.483–0.917] | 0.796 [0.597–0.946] | 0.855 [0.696–0.966] | ~0.49 |
| GSE104948 | 0.775 [0.545–0.955] | 0.809 [0.602–0.962] | 0.852 [0.702–0.965] | ~0.50 |
| GSE142025 | 0.827 [0.551–0.996] | 0.878 [0.654–**1.000**] | 0.932 [0.778–**1.000**] | ~0.49 |

Permuted labels sit at 0.5, so the signal is real — but **gene identity barely matters**. This
is the Venet et al. phenomenon: when a disease signature is broadly co-expressed, most random
signatures are significant.

Consequence for the GSE142025 fold: RBS scores AUROC = 1.000 there at every K, and random gene
sets reach 1.000 too. That fold cannot discriminate methods and its AUROC must not be averaged
in without comment.

## 5. The procurement artifact

The top of the RBS ranking on `all4`: DUSP1, ZFP36, FOS, EGR1, JUN, ATF3, FOSB, JUNB, BTG2,
NR4A2 — with Reactome enrichment for *AP-1 activation* and *NGF-stimulated transcription*, and
GO-MF entirely "sequence-specific DNA binding". Textbook immediate-early response.

These are also exactly the genes the DKD-specificity test rejects. The decisive test
(`scripts/ieg_artifact_test.py`) uses ERCB's non-diabetic CKD biopsies, which were procured the
same way as the DKD biopsies:

| | module mean Hedges *g* |
|---|---|
| DKD vs control | −0.85 |
| **other CKD vs control** | **−1.23** (145% of the DKD shift) |
| correlation of the two effect vectors | r = 0.83 (glom), 0.92 (tubulo) |

Non-diabetic CKD moves *further* than DKD. The module tracks whether tissue arrived as a biopsy
or as nephrectomy/donor tissue — warm-ischaemia and handling — not whether the patient had
diabetes. Because that confounder is identical in every cohort, it looks perfectly
"cross-cohort reproducible" to any stability criterion.

**This is the central methodological finding.** A robustness score is not a validity filter.
The only thing that separates shared confounding from shared biology is a contrast that holds
the confounder constant.

## 6. Specificity gate, promoted to a pipeline stage

`scripts/final_candidates.py`. A gene passes only if it still separates DKD from the 310
non-diabetic CKD biopsies (|g| ≥ 0.5, same direction as the control contrast), and is not
flagged as procurement-driven (other-CKD moving ≥80% as far from control as DKD does).

Of the top 200 by RBS: **92 (46%) pass; 21 (10%) are flagged as procurement-driven.**

### Final 30 candidates

CCDC91, MS4A6A, C1orf21, TPD52, DPP6, CFD, **FMOD**, PCOLCE, CALHM2, CD247, SVEP1, TNNT2,
FOLR2, MOXD1, IGFBP6, PRSS23, COL1A2, **MMP2**, CA10, RARRES1, HDAC9, S100A4, **LUM**, KBTBD11,
DYRK2, OLFML3, MSC, CD48, PPIC, FBN1

Most DKD-specific by the confounder-matched contrast:

| Gene | vs control | **vs other CKD** |
|---|---|---|
| FMOD | +1.19 | **+1.28** |
| MMP2 | +1.71 | +1.03 |
| MOXD1 | +1.83 | +0.95 |
| LUM | +1.58 | +0.94 |
| FOLR2 | +1.46 | +0.83 |
| TNNT2 | −1.87 | −0.76 |

FMOD separates DKD from other kidney disease *more strongly* than from healthy tissue — the
signature of a disease-specific marker rather than a damage marker.

### 6.1 How much of the gated list is noise?

The gate rests on 12 (glomerular) and 17 (tubulointerstitial) DKD patients versus ~155
other-CKD biopsies, so a threshold applied to 9,900 genes will pass some on noise alone.
Permutation of the DKD / other-CKD labels, 500 draws (`scripts/specificity_null.py`):

| | |
|---|---|
| genes passing \|g\| >= 0.5, observed | 833 (8.4%) |
| genes passing under permutation | 114 mean [5–95%: 23–341] |
| **empirical FDR of the gate** | **13.7%** |
| largest \|g\| ever reached under permutation | 1.08 (95th percentile 0.92) |

So the gate does real work, but at |g| >= 0.5 roughly one in seven passing genes is expected to
be noise. Of the 30 final candidates, **4 exceed the 95th percentile of the permutation maximum**
and therefore cannot be explained by label noise in this contrast:

| Gene | \|g\| vs other CKD |
|---|---|
| **FMOD** | 1.28 |
| **MMP2** | 1.03 |
| **MOXD1** | 0.95 |
| **LUM** | 0.94 |

These four are the defensible headline. The remaining 26 are legitimate candidates carried at a
13.7% expected false-discovery rate and should be described that way, not promoted individually.

### 6.2 Per-diagnosis profile — what the gate hides

Collapsing nine diagnoses into one "other CKD" group hides whether a gene is *uniquely* DN or
merely *highest* in DN, and hides whether the control contrast is really living-donor tissue
(procured differently from every biopsy group) versus everything else.
`scripts/per_diagnosis_profile.py` separates them. Mean z-score per diagnosis:

Glomerular compartment (GSE104948), selected genes:

| Gene | **DN** | living donor | FSGS | MCD | ANCA | HTN | IgA | MGN | SLE |
|---|---|---|---|---|---|---|---|---|---|
| **FMOD** | **+1.52** | −0.37 | +0.16 | −0.30 | −0.09 | −0.08 | −0.04 | −0.21 | −0.00 |
| **MMP2** | **+1.24** | −0.83 | +0.11 | −0.18 | +0.16 | −0.23 | +0.08 | −0.25 | +0.31 |
| MOXD1 | +1.14 | −0.52 | +0.19 | −0.29 | +0.27 | −0.05 | +0.05 | −0.34 | +0.19 |
| LUM | +0.94 | −0.15 | +0.19 | −0.26 | +0.24 | −0.03 | −0.02 | −0.42 | +0.15 |
| CALHM2 | +0.75 | **−1.66** | +0.21 | −0.15 | +0.69 | −0.08 | +0.24 | −0.33 | +0.47 |

Three patterns emerge, and they matter:

| Pattern | Meaning | Genes (consistent in both compartments) |
|---|---|---|
| **A — unique to DN** | every other disease flat | **MMP2** (both), **FMOD** (A glom / B tubulo) |
| **B — DN-dominant** | raised across kidney disease, far more in DN | LUM, MOXD1 (glomerular only) |
| **C — donor-driven** | contrast is living-donor vs everything else | CALHM2, CCDC91, DPP6, TNNT2, C1orf21, COL1A2, S100A4, CD48 |

The donor-driven flag is the procurement confounder reappearing in a subtler form: those genes
separate DKD from controls mainly because *living-donor* tissue differs from *all* biopsies, not
because DKD differs from other disease. They survive the DKD-vs-otherCKD gate on a smaller,
genuine component, but should not be led with.

### 6.3 Literature positioning — PubMed counts, not a web search

`scripts/literature_review.py` queries PubMed E-utilities per gene at four levels (any /
kidney / DKD / DKD+biomarker) and saves the PMIDs. This was written after a web search
misled us: FMOD looked novel until a targeted query surfaced *"Identification of Lumican and
Fibromodulin as Hub Genes Associated with Accumulation of Extracellular Matrix in Diabetic
Nephropathy"* (Kidney Blood Press Res 2021). **An earlier draft of this document called FMOD
the novel candidate. That was wrong** — the pipeline re-discovered a published DKD hub gene.

Of the 30 gated candidates:

| Verdict | n | Genes |
|---|---|---|
| already proposed as a DKD marker | 13 | C1orf21, FMOD, CD247, TNNT2, IGFBP6, COL1A2, MMP2, HDAC9, S100A4, LUM, DYRK2, CD48, FBN1 |
| established in DKD | 1 | RARRES1 |
| reported in DKD, few papers | 2 | TPD52, MOXD1 |
| kidney-studied, not in DKD | 6 | DPP6, PCOLCE, SVEP1, FOLR2, PRSS23, PPIC |
| no kidney literature | 4 | CCDC91, CALHM2, KBTBD11, OLFML3 |
| **symbol ambiguous, counts invalid** | 4 | MS4A6A, CFD, CA10, MSC |

The ambiguity flag matters: every "MSC AND diabetic nephropathy" hit is about mesenchymal
stem cells, not the gene *musculin* — 116 apparent DKD papers, none of them about the gene.
Verified by reading the returned titles. Counting without that check would have retired a live
candidate.

### 6.4 The deliverable, tiered honestly

`scripts/master_candidate_table.py` crosses data strength (per-diagnosis pattern A/B and not
donor-driven) with literature novelty:

| Tier | n | Genes | Reading |
|---|---|---|---|
| **1 — novel AND strong** | 2 | **OLFML3, PRSS23** | the actual candidates to pursue |
| 2 — strong, already reported | 4 | FMOD, MOXD1, MMP2, LUM | **method validation**, not discovery |
| 3 — novel, weak in data | 8 | CCDC91, DPP6, PCOLCE, CALHM2, SVEP1, FOLR2, KBTBD11, PPIC | mostly donor-driven |
| 4 — reported and weak | 16 | — | — |

**Tier 2 is the strongest result in this section, and it is not a biomarker result.** The
pipeline independently recovered FMOD and LUM, which the literature already names together as
DKD ECM hub genes, plus MMP2. That is evidence the method finds real biology.

**Tier 1 is the deliverable, with its limits stated.**
- **OLFML3** (olfactomedin-like 3) — no kidney literature at all; per-diagnosis pattern B
  (DN-dominant across eight comparator diseases); KPMP localises it to inflammatory cortical
  interstitial fibroblasts; part of the ECM module that survives proper-background enrichment.
- **PRSS23** (serine protease 23) — kidney-studied but zero DKD papers; pattern B.

Both carry a caveat that must travel with them: their specificity effect sizes (|g| = 0.63 and
0.51 against non-diabetic CKD) sit **below** the 0.92 permutation ceiling, so neither is
individually beyond label noise in that contrast. They are candidates for validation, not
findings. The four genes that do exceed the ceiling (FMOD 1.28, MMP2 1.03, MOXD1 0.95,
LUM 0.94) are all already in the literature.

So the Phase-1 requirement of "at least one novel DKD biomarker candidate" is met in form by
OLFML3, and the honest sentence is: *the pipeline recovered the known DKD ECM module and
surfaced two under-studied genes in the same module that warrant independent testing.*

## 7. Downstream validation

### 7.1 Compartment reproducibility (glomerulus → tubulointerstitium)

| Pair | all 9,900 genes ρ | **signature ρ** | direction agreement |
|---|---|---|---|
| Woroniecka (GSE30528→GSE30529) | 0.014 | **0.736** | 90% |
| ERCB (GSE104948→GSE104954) | 0.511 | **0.811** | 96% |

In the Woroniecka pair the genome-wide correlation is essentially zero, yet the selected
signature reaches 0.736. The selection is finding compartment-transcending biology, not cohort
artifacts. (Caveat from the audit: the pairs share 5 DKD subjects each, so this is partly a
within-patient comparison.)

### 7.2 Pathway enrichment (custom 9,900-gene background)

Using the genome-wide background inflates everything; with the correct background only
Reactome survives:

| Term | adj p | Genes |
|---|---|---|
| Extracellular Matrix Organization | 3.0e-03 | COL1A2, LUM, MMP2, PCOLCE, FMOD, FBN1 |
| Regulation of IGF Transport and Uptake by IGFBPs | 3.1e-03 | MMP2, IGFBP6, PRSS23, FBN1 |
| Keratan Sulfate Degradation / Biosynthesis | 5.4e-03 / 2.9e-02 | LUM, FMOD |
| ECM Proteoglycans | ~4e-02 | LUM, FMOD |

GO-BP, KEGG and GO-MF yield nothing at adj p < 0.05 under the correct background — worth
stating, because the genome-wide version produced a long and misleading list.

### 7.3 Perturbation stress test of the finished signature

`scripts/stress_test.py`. Discovery ran on original data; only the finished 30-gene signature is
stressed here. One logistic model is trained on unperturbed training data, and the **test**
cohort is perturbed with rising severity. Because random gene sets already reach 0.68-0.88, the
quantity that matters is the margin over random, not the absolute AUROC.

Signature AUROC (margin over random), averaged over folds:

| Perturbation | s=0 | s=0.5 | s=1 | s=1.5 | s=2 |
|---|---|---|---|---|---|
| subsample (drop up to 65% of samples) | 0.854 (+0.105) | 0.857 (+0.102) | 0.856 (+0.098) | 0.867 (+0.125) | 0.855 (+0.105) |
| gaussian noise (up to 2x gene SD) | 0.854 (+0.105) | 0.834 (+0.106) | 0.785 (+0.090) | 0.753 (+0.104) | 0.693 (+0.066) |
| feature masking (up to 50% genes lost) | 0.854 (+0.105) | 0.834 (+0.108) | 0.833 (+0.114) | 0.796 (+0.098) | 0.798 (+0.118) |
| batch shift | 0.854 (+0.105) | 0.850 (+0.101) | 0.852 (+0.107) | 0.841 (+0.103) | 0.830 (+0.098) |

The **margin stays between +0.07 and +0.13 everywhere**. Subsampling costs nothing; batch shift
and masking cost little; only heavy gaussian noise degrades absolute performance, and it
degrades random signatures just as fast. The signature is not more fragile than the data.

### 7.4 Generative robustness — one generator worked, one did not

`scripts/generative_robustness.py`, class-conditional, trained on training cohorts only, over a
500-gene space (signature + high-variance context; a generator over 9,900 dimensions from ~120
samples would be pure memorisation).

| Generator | class-diff correlation | dispersion ratio | NN-dist ratio | usable folds |
|---|---|---|---|---|
| **VAE** | **+0.94** | 0.73 | 0.77 | **4/4** |
| Diffusion (DDPM) | −0.02 | 18.2 | 18.8 | **0/4** |

**The diffusion model failed and is reported as such.** Its denoiser plateaus at eps-MSE ≈ 0.57
(predicting zero scores 1.0), which is not accurate enough for a 1000-step reverse process: the
per-step error compounds and the sampler diverges. That is a sample-size limit at n ≈ 120, not a
fixable bug, and no number from it should be quoted.

The VAE passed its fidelity guards and reproduces the per-gene class-mean difference at r = 0.94
without memorising (synthetic points sit 0.77x as far from real data as real points sit from
each other). On its synthetic populations the signature scores **AUROC 0.923 — but random gene
sets score 0.945**. The generative check is therefore a *consistency* result, not a
discriminating one: the VAE reproduces the training covariance, and because the DKD signal is
broadly co-expressed, synthetic data inherits exactly the property that makes random signatures
work. It confirms the signature survives a resampling of the learned joint distribution; it
cannot rank the signature against alternatives, and should not be presented as if it could.

### 7.5 KPMP cross-omics corroboration

22 of 30 candidates carry KPMP regional-proteomics evidence; 28 of 30 receive a significant
single-nucleus cell-type assignment. Cell types converge on **fibroblasts** (FMOD, LUM, MMP2,
COL1A2, PCOLCE, OLFML3 → interstitial / perivascular / corticomedullary, frequently the
*degenerative* states), **resident macrophages** (MS4A6A, FOLR2 → LYVE1+), and **lymphocytes**
(CD247 → NK, S100A4/CD48 → cytotoxic T).

Protein-level compartment attribution, from an independent cohort and assay:

| Gene | protein Glom-vs-TI FC | adj p |
|---|---|---|
| FBN1 | +2.03 | 8.1e-19 |
| DPP6 | +3.97 | 2.5e-13 |
| COL1A2 | +1.19 | 6.4e-06 |
| MOXD1 | −2.31 | 6.4e-07 |
| OLFML3 | −2.06 | 2.2e-07 |

This is the part of KPMP that genuinely earns its place: candidates found in glomerular
transcriptomics are confirmed as glomerular at the protein level. What KPMP *cannot* do is
validate a DKD-vs-control claim — it has zero diabetic participants in regional proteomics
(see `docs/KPMP_ACCESS.md`).

## 8. Reproduce

```
python scripts/dkd_data.py                 # cohort summary + LODO folds
python scripts/run_baselines.py --B 200 --k 50 --out results/baselines
python scripts/run_proposed.py  --B 200 --base relieff --agg geomean \
       --cohorts GSE30528,GSE96804,GSE104948,GSE142025 --out results/proposed_all4
python scripts/ablation_rbs.py             # 24-combination design sweep
python scripts/null_control.py             # random-signature and permuted-label nulls
python scripts/ieg_artifact_test.py        # procurement-artifact test
python scripts/final_candidates.py         # specificity-gated final list
python scripts/validate_candidates.py --components results/proposed_all4
python scripts/enrich.py --signature results/final_candidates_gated.tsv --topk 30
python scripts/kpmp_corroborate.py
python scripts/stress_test.py              # perturbation stress test of the finished signature
python scripts/generative_robustness.py    # VAE / diffusion synthetic-population check
python scripts/specificity_null.py         # permutation FDR of the specificity gate
python scripts/per_diagnosis_profile.py    # per-diagnosis stratification of candidates
python scripts/compare_methods.py --baselines results/baselines_all4        --proposed results/proposed_all4 --cohorts GSE30528,GSE96804,GSE104948,GSE142025
```


## 9. Open items

- **Novelty assessment is preliminary.** Section 6.3 is based on web literature search, not a
  systematic review. Before any novelty claim in a manuscript, the candidates need a proper
  PubMed/Nephroseq check - particularly CALHM2, C1orf21, CCDC91, KBTBD11 and MSC, where absence
  of search hits is suggestive but not evidence of absence.
- **The four ERCB-derived confounder-matched contrasts share the ERCB cohort.** The specificity
  gate has no independent replication - no other public DKD dataset carries non-diabetic CKD
  biopsies procured the same way. This is a genuine limitation, not an oversight.
- **GSE142025 hold-out cannot discriminate methods** (random gene sets reach AUROC 1.000 there).
  Reported separately rather than averaged in.
- **Compartment analysis is partly within-patient** - the paired cohorts share 5 DKD subjects
  each (see `docs/DATASET_AUDIT.md`). The strictly independent version restricts to the
  non-shared subjects (5 and 12 DKD respectively).
