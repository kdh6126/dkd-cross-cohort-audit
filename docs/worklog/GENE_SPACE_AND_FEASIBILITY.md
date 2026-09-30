# Common Gene Space + Cross-Cohort Feasibility
Date: 2026-08-24 · Precondition check for the stability-aware cross-cohort feature selection design

## 1. Hub ID and mapping routes

Hub = **NCBI Entrez Gene ID**, chosen because GSE104948/GSE104954 use the Brainarray ENTREZG
custom CDF — their probe IDs (`10000_at`) already *are* Entrez IDs, so this route adds no
mapping loss on the tightest cohort.

| Cohort | Source IDs | Route | Mapped |
|---|---|---|---|
| GSE142025 | HGNC symbol | NCBI `Homo_sapiens.gene_info` (symbol → synonym fallback) | 16,808 / 17,184 (925 via synonym) |
| GSE30528 / GSE30529 | Affy GPL571 probe | GEO `GPL571.annot` *Gene ID* | 19,931 / 22,277 |
| GSE96804 | Affy GPL17586 transcript cluster | GPL17586 `gene_assignment` field | 23,521 / 70,523 |
| GSE104948 / GSE104954 | ENTREZG CDF | strip `_at` | 12,074 / 12,074 |

Probes mapping to **more than one** distinct Entrez gene are dropped (GPL571: 1,223;
GPL17586: 1,572). GPL17586's low yield is expected — 45,660 of its 70,523 transcript clusters
carry no gene assignment at all (non-coding / unannotated clusters).

## 2. Collapse rule — fixed up front

**MAX-MEAN**: where several probes map to one gene, keep the single probe with the highest mean
expression in that cohort (WGCNA `collapseRows` "MaxMean"). Applied identically in every cohort
so all LODO folds share one feature definition. Changing this rule later invalidates every
stability count, so it is frozen here.

## 3. The common gene space

| Set | Genes |
|---|---|
| GSE142025 | 16,779 |
| GSE30528 / GSE30529 | 12,502 |
| GSE96804 | 22,203 |
| GSE104948 / GSE104954 | 12,074 |
| **Intersection, 4 discovery cohorts** | **9,900** |
| Intersection, all 6 (incl. compartment validation) | **9,900** (no further loss) |

Cost of each cohort, as genes lost from the intersection:

| Excluded | Intersection becomes | Gain |
|---|---|---|
| GSE142025 | 11,365 | +1,465 |
| GSE104948 | 10,326 | +426 |
| GSE96804 | 10,204 | +304 |
| GSE30528 | 10,008 | +108 |

GSE142025 is the most expensive member (its RNA-seq quantification simply does not report
1,465 genes the arrays share), but +1,465 genes is a fair price for the only RNA-seq cohort.
Adding GSE30529/GSE104954 costs **nothing** — they share platforms with their glomerular
counterparts — so compartment validation is free.

Artifacts: `data/processed/gene_space_discovery4.tsv`, `gene_space_all6.tsv`,
`probe2gene.tsv`, and harmonised matrices in `data/processed/harmonized/`.

## 4. Harmonised cohorts as built

| Cohort | genes × samples | DKD | control | other CKD | value range | compartment |
|---|---|---|---|---|---|---|
| GSE142025 | 9,900 × 36 | 27 | 9 | 0 | 0.40 … 19.29 | whole cortex |
| GSE30528 | 9,900 × 22 | 9 | 13 | 0 | −7.68 … 7.91 | glomerulus |
| GSE96804 | 9,900 × 61 | 41 | 20 | 0 | 2.18 … 19.37 | glomerulus |
| GSE104948 | 9,900 × 196 | 12 | 26 | 158 | 2.51 … 14.94 | glomerulus |
| GSE30529 | 9,900 × 22 | 10 | 12 | 0 | −7.59 … 7.43 | tubulointerstitium |
| GSE104954 | 9,900 × 195 | 17 | 26 | 152 | 2.70 … 14.35 | tubulointerstitium |

Two things to carry forward:
- **Value scales differ fundamentally.** GSE30528/GSE30529 are mean-centred (negative values);
  the others are positive log-scale. Any cross-cohort step must use scale-free quantities —
  per-cohort z-scores, ranks, or effect sizes. Never pool raw values.
- **GSE104948/GSE104954 carry 158/152 non-diabetic CKD samples** (SLE, IgA, FSGS, MCD, MGN,
  RPGN, hypertensive). These are not waste: they are a ready-made **specificity** test —
  is a candidate DKD-specific, or just a generic CKD/fibrosis marker? That is a strong
  addition to the paper and costs no new data.

## 5. Feasibility result — a reproducible cross-cohort signal exists

Per-gene effect sizes (Hedges' *g*, small-sample corrected) computed per cohort, DKD vs control.

**Direction concordance across the 4 discovery cohorts: 2,689 / 9,900 genes = 27.2%**
(chance = 12.5%). Of those, **827 genes** also have |g| > 0.5 in all four.

The premise of the whole design — that some genes are repeatedly selectable across cohorts
and platforms — holds. It is not guaranteed a priori, and it is now measured.

### Biological sanity check (top concordant genes, unprompted)

| Gene | mean *g* | Interpretation |
|---|---|---|
| NPHS1 (nephrin) | −1.95 | podocyte loss — the canonical DKD lesion |
| COL1A2, PCOLCE, FMOD | +1.4 … +2.0 | ECM / fibrosis |
| C1QA, MS4A6A | +2.0 … +2.2 | complement + macrophage infiltration |
| DUSP1, ZFP36, NR4A2, BTG2 | −1.9 … −2.9 | immediate-early stress-response genes, down |

The pipeline recovers textbook DKD biology without being told to. That is the strongest
available evidence that the ID mapping and collapse are correct.

## 6. Cross-cohort effect-size correlation (Spearman ρ)

|  | GSE142025 | GSE30528 | GSE96804 | GSE104948 | GSE30529 | GSE104954 |
|---|---|---|---|---|---|---|
| **GSE142025** (cortex) | – | 0.217 | 0.308 | 0.441 | 0.252 | 0.381 |
| **GSE30528** (glom) | 0.217 | – | 0.204 | 0.232 | **0.014** | 0.033 |
| **GSE96804** (glom) | 0.308 | 0.204 | – | 0.376 | 0.169 | 0.385 |
| **GSE104948** (glom) | 0.441 | 0.232 | 0.376 | – | 0.396 | **0.511** |
| **GSE30529** (tub) | 0.252 | 0.014 | 0.169 | 0.396 | – | 0.577 |
| **GSE104954** (tub) | 0.381 | 0.033 | 0.385 | 0.511 | 0.577 | – |

Three findings that should shape the experimental design:

**(a) GSE30528 is the weakest link.** It correlates 0.20–0.23 with everything and essentially
**zero (0.014)** with its own lab's tubular counterpart GSE30529. With 9 cases / 13 controls its
effect sizes are the noisiest in the set. Expect the LODO fold that holds out GSE30528 to be the
hardest, and say so rather than hiding it — a method that stays stable on the noisiest cohort
is exactly the claim being made.

**(b) The compartment effect is real but wildly uneven.** The same glom↔tub contrast gives
ρ = 0.014 in the Woroniecka data and ρ = 0.511 in ERCB. So "does a glomerular biomarker
reproduce in tubulointerstitium?" has no single answer — it depends on the cohort. That is a
more interesting result than a flat yes/no, and it argues for reporting both pairs separately
rather than pooling them.

**(c) The cortex/compartment warning is confirmed empirically.** GSE142025 (whole cortex)
correlates 0.381 with GSE104954 (tubulointerstitium) versus 0.217 with GSE30528 (glomeruli) —
consistent with bulk cortex being volume-dominated by tubulointerstitium. Treating GSE142025 as
interchangeable with microdissected glomerular cohorts is not supportable.

## 7. Recommended design, updated by these numbers

1. **Primary LODO** over the three glomerular cohorts (GSE30528, GSE96804, GSE104948) —
   one compartment, clean interpretation.
2. **GSE142025 as a cross-compartment generalisation test**, not a LODO member.
3. **GSE30529 / GSE104954** for compartment reproducibility, reported as two separate pairs,
   restricted to non-shared subjects for the headline figure (5 and 12 DKD respectively).
4. **The 310 non-diabetic CKD samples in ERCB** as a DKD-specificity filter on final candidates.
5. Per-cohort z-scoring or rank transform before any cross-cohort step. Platform
   (GPL22945 vs GPL24120) as a nested batch factor inside GSE104948/GSE104954.

## 8. Reproduce

```
python scripts/build_gene_space.py       # ID maps + intersection
python scripts/harmonize.py              # collapse + harmonised matrices + labels
python scripts/sanity_effect_sizes.py    # effect sizes, concordance, correlation matrix
```
