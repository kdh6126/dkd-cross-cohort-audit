# Cross-tissue extension: post hoc analyses, stated before they were run

Written 2026-09-14, after the pre-registered verdicts in `results/xtissue/verdicts.tsv` were
known (P1 held in liver and failed in colon, P2 failed, P3 held) and before any of the analyses
below were run. These are post hoc by construction and will be reported as such, separately
from P1-P3. The SHA-256 of this file is recorded in `results/xtissue/posthoc_plan.sha256`.

## H4. Control-definition substitution outside kidney

GSE48452 contains, within one study, NASH cases alongside two control groups: non-obese
controls obtained during major oncological liver surgery, and obese individuals with normal
liver histology sampled during the same bariatric operations as most cases. This is the
within-resource design the manuscript used in ERCB.

Expectation: g_IEG (NASH minus control) is more negative against oncological-surgery controls
than against bariatric healthy-obese controls, and the top-50 genes by |g| under the two
control definitions share a minority of genes.

## H5. Separating the inflammation response of the immediate-early module from procurement

Colon showed a large positive g_IEG in procurement-matched cohorts. A leukocyte and acute
inflammation score is defined from genes chosen before looking at data and excluding every
immediate-early gene: PTPRC, CD14, CD68, FCGR3A, S100A8, S100A9, LCN2, CXCL8, IL1B, TNF, SAA1,
MMP3, CHI3L1. Within each cohort the handling score is regressed on case status with and
without that inflammation score, both standardised.

Expectation: in colon the case coefficient on the handling score shrinks by more than half after
adjustment; in procurement-asymmetric kidney and liver cohorts it stays negative.

## H6. Cross-cohort effect concordance as the condition for verdict instability

Observed post hoc: median pairwise Spearman correlation of gene-wise g across cohorts is 0.66 in
colon, 0.22 in kidney, 0.15 in liver, and only colon showed no sign change in the three-cohort
RBS minus WGCNA Jaccard difference. If low concordance is the condition, reducing concordance
in colon should introduce sign changes.

Procedure: each colon core cohort is subsampled without replacement to the case and control
counts of a kidney discovery cohort (matched by rank of cohort size), five independent
replicates. For each replicate the concordance and the ten three-cohort Jaccard differences
are recomputed with the manuscript code.

Expectation: concordance falls towards the kidney value and at least one three-cohort Jaccard
difference becomes negative in at least three of five replicates.

## S1. Sensitivity of liver immediate-early coverage to the RNA-seq expression filter

The CPM >= 1 in >= 20% filter removed low-expressed immediate-early genes from the two liver
RNA-seq cohorts, leaving 11 of 19 in the liver core gene space. The filter is relaxed to a raw
count of at least 1 in at least 20% of samples and the liver full-set P3 count is recomputed.
Expectation: the P3 ordering (kidney > liver > colon) is unchanged.
