# Open items

Moved out of the manuscript on 2026-08-25 so that `docs/manuscript/PAPER_DRAFT.md` contains
only text intended for a reviewer. This file tracks work; it is not part of the submission.

## Closed

- **Single entry point.** `run_all.py` covers 46 stages (41 on the reproduction path, 5
  exploratory) with declared outputs and `--list / --stage / --from / --all / --paper /
  --dry-run / --skip-network / --force`. Every script in `scripts/` is now reachable; before
  this, seven were orphaned, including `paper_figures.py`, which draws the manuscript figures.
- **Environment pinned.** `requirements.txt` with `==` pins for the nine packages actually
  imported. statsmodels and scanpy are not needed and were removed from the assumed set.
- **Stage output is logged.** Each stage's combined stdout/stderr goes to `logs/<stage>.log`.
  Previously `run_all.py` only streamed to the console and the `logs/` directory held stale
  files from manual runs.
- **Determinism verified empirically**, not just by inspection: `spec_null` (500 permutations)
  rerun with `--force` produces a byte-identical output file.
- **Gate FDR propagated.** Pass-rates are reported as fold-enrichment over the 8.3% random
  background with the 13.9% empirical gate FDR alongside, in all four result tables
  (`scripts/propagate_gate_fdr.py`).
- **Figures.** P1-P6 written by `scripts/paper_figures.py`; P1 (cohort sensitivity) is the
  lead figure.
- **Heart analysis placement.** Kept in the main text. It bounds the claim to procurement
  rather than to tissue injury generally; moving a claim-bounding negative control to
  supplementary would be the selective reporting the paper objects to.
- **Prior art for cohort-composition sensitivity.** Peer-reviewed anchors located and cited:
  Boulesteix 2013 (PMID 23637855), Jelizarow 2010 (PMID 20581402), Ullmann 2023
  (PMID 36608142), Weber 2019 (PMID 31221194). The arXiv-only "benchmark lottery" reference is
  no longer load-bearing and the phrase has been removed from the Introduction.
- **Unsourced pre-analytical claim.** The sentence asserting that the pre-analytical stage
  accounts for the majority of downstream analytical error has been deleted. It was traceable
  only to a review summary, and the argument does not depend on it.

- **Figures are cited.** The manuscript previously contained six figures and zero references
  to them. Callouts were added to S3.1-S3.5 and a Figure legends section written. Two legend
  claims that did not match the plotted data were corrected against the result files.
- **P2 and P4 are reproducible.** They had been copied by hand from the exploratory set, so a
  clean checkout produced four of six manuscript figures. `paper_figures.py` now copies them
  explicitly and `--paper` runs `figures` first.
- **Inventory corrected.** `data_inventory.py` labelled GSE175759, the KPMP h5ad and the GWAS
  catalog as unused or not downloaded when all three feed reported results, omitted GSE162830,
  GSE5406 and GSE294519 entirely, and derived `state` from a prefix whitelist that silently
  demoted any asset whose description was reworded. 14 assets are used, not 8.


## Open
- **LODO-as-reported-practice citation.** Refs 1-4 establish that benchmark composition
  matters and should be varied. What is still not located is an omics *biomarker-discovery*
  paper that reports the resulting interval. If it stays unlocated, §3.2 should say so
  explicitly rather than implying the issue has never been raised.
- **Positioning.** RBS must be presented as a case-study instrument, not a proposed method.
  Stability selection has close relatives and the cohort-sensitivity result argues against
  claiming a winner.
