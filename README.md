# Cross-cohort feature selection in diabetic kidney disease

ETRI, Phase 1 Year 1 - multi-omics database construction and biomarker-candidate discovery for
diabetic kidney disease (DKD).

The work asks whether a feature-selection method can be shown to produce *reproducible*
biomarkers across independent DKD cohorts. The answer turned out to be that the question, as
usually posed, is not well defined: which method wins depends on which cohorts are in the
benchmark, and the most cross-cohort-reproducible signal in this disease tracks how the tissue
was obtained rather than whether the patient had DKD.

**Scope: single-omics.** Every analysis here is bulk or single-nucleus transcriptomics. The
reasons no second modality was usable are recorded in
[docs/project/WHY_MULTIMODAL.md](docs/project/WHY_MULTIMODAL.md); in short, KPMP proteomics has
no diabetic participants and the available human DKD metabolomics depositions have no usable
comparison group.

**Deliverable status.** No defensible novel biomarker. Candidates that survived the pipeline
failed independent validation, and the paper reports that rather than the candidate list.

## Which manuscript is in this repository

`submission_bib/` is the submitted version (Briefings in Bioinformatics, Problem Solving
Protocol) with its supplementary notes, files and figures. An extended manuscript in another
journal's format exists in the authors' working tree and is not distributed; the pipeline stages
that audit that file are skipped automatically in this copy, and none of them produce results.

## Where things are

| path | role |
|---|---|
| [REPRODUCE.md](REPRODUCE.md) | how to run everything; one entry point, 133 stages, 104 on the paper path |
| [run_all.py](run_all.py) | that entry point |
| [requirements.txt](requirements.txt) | pinned environment |
| `docs/project/` | ETRI deliverables - hybrid DB architecture, the multi-omics rationale |
| `docs/worklog/` | how the work actually went, including dead ends. Not for submission. |
| `scripts/` | the stage scripts and library modules (Python; one R script for the reference WGCNA check) |
| `results/` | every number the manuscript quotes, as TSV |
| `db/` | schema of the cohort-audit database; the SQLite and Parquet stores are rebuilt by the pipeline, not shipped |

Manuscript: `submission_bib/bib-manuscript.tex` (PDF alongside).
Open items: [docs/worklog/OPEN_ITEMS.md](docs/worklog/OPEN_ITEMS.md).

## Applying the protocol to another disease

The audit is four stages of `run_all.py`; each reads the previous stage's output. To run it on
new cohorts, replace the inputs of stage 1 and run the four in order.

| step | what you supply | stage | what comes out |
|---|---|---|---|
| 1 audit the input | one expression matrix per cohort (genes x samples), a sample table with `case/control` and the **procurement route of each arm** (biopsy, nephrectomy, donor, endoscopic, ...) | `python run_all.py --stage audit` | patient-overlap table across accessions; a cohort table with procurement per arm. Series sharing patients are excluded here |
| 2 vary the cohorts | the retained cohorts | `python run_all.py --stage cohort_sensitivity` (add `--stage sens_selectors` to repeat it under other inner selectors) | `results/cohort_sensitivity/subsets.tsv`: the method difference on every admissible subset, and its sign changes |
| 3 vary the control | a second control group obtained the same way as the cases, if one exists | `python run_all.py --stage reordering` | genes at effect-size threshold under each control definition; overlap of the two top-50 lists |
| 4 benchmark the corrections | the confounder score from step 1 (here the immediate-early handling score) | `python run_all.py --stage ruv` | per correction: share of selected features flagged as confounder-driven, and external AUROC |

Read the result against the checklist in the manuscript (Table 4): each row states what to
report and what to do when the check fails. If step 3 has no same-procedure control, say so and
treat the candidate order as conditional; if the confounder score in step 4 also carries disease
biology (it did in colon), do not residualise on it blanket-fashion.

Cohort names, compartments and the immediate-early gene list live in `scripts/dkd_data.py` and
`scripts/dkd_deconfound.py`; those are the two files a new disease changes.

## Data

25 public assets catalogued in `results/data_inventory.tsv`; 14 used, 7 downloaded but unused,
3 excluded by the audit, 1 not retrieved.

Five discovery cohorts, 258 labelled samples after harmonisation into a frozen 9,900-gene
Entrez space:

| cohort | compartment | DKD | control |
|---|---|---:|---:|
| GSE30528 | glomerulus | 9 | 13 |
| GSE96804 | glomerulus | 41 | 20 |
| GSE104948 | glomerulus | 12 | 26 |
| GSE142025 | whole cortex | 27 | 9 |
| GSE294519 | tubulointerstitium | 23 | 13 |
| GSE30529 | tubulointerstitium | 10 | 12 |
| GSE104954 | tubulointerstitium | 17 | 26 |

Three further series were **excluded for patient reuse**, which the audit found and which the
database encodes as views: GSE30122 duplicates GSE30528+GSE30529, GSE47183's DN subjects are
identical to GSE104948, and GSE99340 is an ERCB superset. Reusing them would have inflated
every cross-cohort number in this study.

## What the analysis found

1. **Absolute AUROC is uninformative here.** Random 50-gene signatures reach 0.68-0.88 by fold.
   Any classifier result must be read against that null, not against 0.5.
2. **Method rankings depend on cohort composition.** Over 3-cohort subsets the proposed
   selector wins 9 of 10; over all five cohorts it loses. The difference ranges from -0.106 to
   +0.170 depending only on which cohorts are included.
3. **The most reproducible signal is procurement, not disease.** Cases are biopsies, controls
   are nephrectomies or donor kidneys, and an immediate-early gene module tracks that split in
   8 of 9 kidney datasets (Spearman -0.752, p < 0.001). A biopsy-vs-biopsy control in KPMP
   removes it (p = 0.21 vs p = 6.0e-76). It is weak in heart; outside kidney it followed
   procurement in liver, and in colon, where both arms are biopsies, it rose with inflammation,
   so the correction built on it is tissue-conditional (paper, Results and Limitations).
4. **Corrections that work.** Handling-score residualisation and SVA cut procurement-driven
   features from 27% to 5-6% of the selected set; fold-mean AUROC is unchanged to three
   decimals, but the paired-bootstrap interval on that cost is wide, so performance
   preservation is not established. ComBat returns the same values as no correction because
   the confounder is within-cohort, not between.

## Honesty notes

Several claims in earlier drafts were wrong and were retracted rather than quietly edited:
FMOD/LUM were not novel (Feng 2021, PMID 33887734); a reproducibility win for the proposed
method did not survive the five-cohort bootstrap; a within-group residualisation variant made
AUROC worse, not better. All three are documented in `docs/worklog/`.


---

## What is and is not in this repository

Everything derived here is included: the harmonised expression matrices, the frozen 9,900-gene
Entrez space and probe map, the sample-level phenotype tables with the procurement annotation,
every result table the paper quotes, the figures, the database schema, and all code.

Three categories are deliberately absent.

**Primary third-party data** (`data/raw/`, about 7 GB) is not redistributed. The KPMP
single-nucleus atlases carry their own data use terms, and the GEO series, GWAS catalog and
platform annotations all have canonical homes that should remain the source of truth. Fetch them
with:

```bash
python run_all.py --stage download
```

**Regenerable stores** are omitted so they cannot drift from the code that builds them: the
SQLite database and its Parquet mirror (`python run_all.py --stage database`, 3 minutes) and the
single-nucleus expression caches, which their own stages rebuild.

**The status deck** is an internal progress artefact and is not part of the published work.

As a result, a fresh clone reports five stages as outstanding: `download`, `database`, `deck`,
and the two exploratory stages that declare no output files. Every stage backing a number in the
manuscript is already satisfied by the files shipped here, so
`python run_all.py --paper --dry-run` will show them all as present.
