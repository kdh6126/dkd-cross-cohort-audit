# Post hoc analyses: outcomes against the plan

Plan: `docs/xtissue/POSTHOC_PLAN.md`, SHA-256 7b2bb5a4b16865317e34ed486770b396623c5d56b01e515856728b5afccc9035,
unchanged at the time these outcomes were recorded. Every analysis below is post hoc and is
reported separately from the pre-registered P1-P3.

## H4. Control-definition substitution in liver (GSE48452)

| Control group | n case / control | g_IEG [95% CI] | genes with abs(g) >= 0.5 | median abs(g) |
|---|---|---|---|---|
| oncological surgery | 17 / 12 | -0.38 [-1.45, +0.31] | 5,238 | 0.312 |
| bariatric normal liver | 17 / 16 | -0.02 [-0.84, +0.61] | 4,962 | 0.293 |

Top-50 genes shared between the two control definitions: 2 of 50. Spearman over all 18,714
genes: 0.48. Median abs(g) ratio, matched over naive: 0.94.

Outcome: **direction as expected, not resolved.** g_IEG is more negative against surgical
controls, but the intervals overlap. The top of the ranking is rebuilt (2 of 50 shared, as in
ERCB diabetic nephropathy with 1 of 50), while overall signal size barely shrinks, unlike ERCB.

## H5. Inflammation adjustment of the immediate-early case effect

| Arms | Cohorts | Negative, as measured | Negative, inflammation removed | Median coefficient, as measured | Median, removed |
|---|---|---|---|---|---|
| asymmetric (kidney 7, liver 3) | 10 | 9 | 10 | -0.94 | -1.14 |
| matched (liver 4, colon 10) | 14 | 3 | 5 | +0.89 | +0.30 |

Asymmetric versus matched, Mann-Whitney: p = 9.9e-05 as measured, 1.6e-04 after removal.
Colon median retained fraction of the case coefficient: 0.51.

Outcome: **partly held.** Asymmetric effects stay negative and become uniformly negative after
removal (kidney GSE294519 moves from +0.16 to -0.46). Colon shrinks by 49%, just short of the
stated "more than half".

## H6. Downsampling colon to kidney cohort sizes

| Replicate | Concordance | 3-cohort Jaccard difference | Negative |
|---|---|---|---|
| 0 | 0.638 | -0.013 to +0.143 | 1 |
| 1 | 0.613 | +0.034 to +0.255 | 0 |
| 2 | 0.604 | -0.001 to +0.225 | 1 |
| 3 | 0.642 | +0.000 to +0.239 | 0 |
| 4 | 0.647 | +0.004 to +0.193 | 0 |

Outcome: **failed as stated.** A sign change appeared in 2 of 5 replicates (threshold 3), both
within 0.013 of zero, and concordance stayed near 0.64 rather than falling towards the kidney
value of 0.22. Smaller cohorts alone do not reproduce the kidney-sized swings; the large
swings co-occur with low cross-cohort concordance, which downsampling a concordant tissue does
not create.

## S1. Liver expression filter

Relaxing the RNA-seq filter to a raw count of at least 1 in at least 20% of samples raises the
liver core gene space from 7,805 to 9,444 genes and immediate-early coverage from 11 to 15 of
19. Uncorrected stability top-50 immediate-early count stays at 0.6.

Outcome: **held.** The P3 ordering kidney 8.0 > liver 0.6 > colon 0.0 is unchanged.
