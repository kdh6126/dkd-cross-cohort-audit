# Cross-tissue extension: predictions fixed before any expression data were read

Written 2026-09-14, after reading only GEO series-level and sample-level *metadata* for the
candidate liver (MASLD/MASH) and colon (ulcerative colitis) cohorts, and before downloading or
inspecting any expression values. The SHA-256 of this file is recorded in
`results/xtissue/predictions.sha256` at the moment it was written, so any later edit is visible.

## Why these two tissues

The manuscript makes two claims that the kidney data alone cannot separate.

- **C1.** A cross-cohort *reproducibility* verdict between selection methods depends on which
  cohorts are in the benchmark, while the *predictive* verdict is more stable.
- **C2.** A pre-analytical confounder becomes the most reproducible signal when, and because,
  every cohort shares the same procurement asymmetry between cases and controls.

C1 is a statistical claim and should hold in any tissue. C2 is a design claim and should hold
only where procurement differs between arms, and only become *reproducible* where that
difference is shared across cohorts.

| Tissue | Procurement of cases vs controls across cohorts | Role |
|---|---|---|
| Kidney, DKD (existing) | asymmetric in every cohort, same direction (needle biopsy vs nephrectomy/donor) | C2 present and shared |
| Liver, MASH | mixed: asymmetric in some cohorts (needle biopsy vs living-donor surgery), matched in others (all needle biopsy, or all intraoperative at bariatric surgery) | C2 present but not shared |
| Colon, UC | matched in every cohort (endoscopic biopsy in both arms) | C2 absent |

## Predictions

**P1 (C1, every tissue).** In each new tissue with at least five independent cohorts, the paired
difference in cross-fold Jaccard between the stability selector (RBS) and the connectivity
criterion (WGCNA_hub) spans zero across three-cohort subsets. The paired difference in
external AUROC changes sign in fewer three-cohort subsets than the Jaccard difference does.

**P2 (C2 at the cohort level).** For each cohort, let g_IEG be Hedges' g of the
immediate-early handling score, case minus control. In procurement-matched cohorts
|g_IEG| < 0.5. In cohorts where cases are needle or endoscopic biopsies and controls are
surgical or donor tissue, g_IEG is negative and |g_IEG| is larger than in every matched
cohort of the same tissue.

**P3 (C2 becomes reproducible only when shared).** The mean number of immediate-early genes in
the uncorrected RBS top 50 across leave-one-dataset-out folds is ordered kidney > liver >
colon. In colon it is at most 1. The handling-score correction reduces that number
substantially only in kidney.

## What would count against the manuscript

- P1 failing in both new tissues would mean the composition instability is a property of
  these five kidney cohorts rather than of cross-cohort benchmarking.
- P2 failing, in particular a large |g_IEG| in procurement-matched colon or liver cohorts,
  would mean the module tracks disease or inflammation rather than procurement, and the
  kidney attribution would need to be weakened.
- P3 failing, in particular immediate-early genes entering the colon top 50, would mean
  reproducibility of the module does not require a shared design.

Every outcome, including failures, will be reported.

## Decisions fixed now, before data

- Case and control definitions. Liver: NASH versus histologically normal liver; cohorts whose
  controls are obese without NAFLD are analysed but flagged. Colon: actively inflamed UC colon
  versus normal control colon; one biopsy per subject; ileal and Crohn's samples excluded.
- Procurement coding is taken from the series text, sample source fields and the source
  publication, and recorded per arm with its evidence before g_IEG is computed.
- The immediate-early gene list is the 19-gene module already fixed in the manuscript. It is
  not re-derived in any new tissue.
- Signature size K = 50 and the four selection strategies are unchanged from the manuscript.
- Patient-level reuse is audited across series from the same group or platform before any
  cohort enters an analysis, by the same procedure as for kidney.
