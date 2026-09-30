# Licence for the data in this repository

## What this covers

The derived data files distributed here -- the harmonised expression matrices under
`data/processed/harmonized/`, the sample-level phenotype and metadata tables under
`data/metadata/`, and the result tables under `results/` -- are released under the
**Creative Commons Attribution 4.0 International licence (CC BY 4.0)**.

Full text: https://creativecommons.org/licenses/by/4.0/legalcode
Summary:   https://creativecommons.org/licenses/by/4.0/

Attribute by citing the associated publication (see `CITATION.cff`).

## What this does not cover

These files are *derived* from primary data that we did not generate and do not redistribute.
The primary data remain under the terms of their original sources, and those terms govern any
use of them:

- **Gene Expression Omnibus series** (GSE30528, GSE30529, GSE96804, GSE104948, GSE104954,
  GSE142025, GSE294519, GSE162830, GSE175759, GSE142153, GSE1009, GSE111154, GSE20602,
  GSE5406). https://www.ncbi.nlm.nih.gov/geo/
- **Kidney Precision Medicine Project** single-nucleus atlas and Atlas API.
  https://atlas.kpmp.org/ -- subject to KPMP's own data use terms. Only the open-access,
  de-identified portion was used; no controlled-access data was requested or held. The `.h5ad`
  atlas files are **not** included in this repository; `run_all.py --stage download` retrieves
  them from KPMP. Results derived from KPMP are distributed as aggregates: the modality table
  is summarised by diagnosis and participant identifiers are not reproduced.
- **Metabolomics Workbench** studies ST003255, ST001411, ST000691, ST004483, ST004442 and
  ST002145. https://www.metabolomicsworkbench.org/ -- the measurement tables are **not**
  included here; only per-metabolite contrast results computed from them.
- **Mendeley Data 83k89shdx5** -- SOMAscan proteomics of human kidney cortex biopsies
  (23 diabetic kidney disease, 10 healthy), CC BY 4.0.
  https://data.mendeley.com/datasets/83k89shdx5 -- the measurement workbook is **not**
  included here; only the per-protein contrast computed from it.
- **PRIDE / ProteomeXchange PXD041884** -- LC-MS/MS of archival formalin-fixed human kidney
  cortex (5 diabetic, 7 non-diabetic). https://www.ebi.ac.uk/pride/archive/projects/PXD041884
  -- the protein rollup workbook is **not** included here; `run_all.py --stage pride_fetch`
  retrieves it. Only the per-protein contrast is distributed.
- **NHGRI-EBI GWAS Catalog**. https://www.ebi.ac.uk/gwas/
- **Platform annotation files** from GEO and NCBI `gene_info`.

The code in this repository is licensed separately, under the MIT licence in `LICENSE`.
