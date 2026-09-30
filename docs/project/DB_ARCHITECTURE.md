# Hybrid multi-omics database — architecture and reference implementation
Phase 1 deliverable: *다중 오믹스 데이터 통합 관리 구조 및 하이브리드 DB 아키텍처 설계·구축*

Built: `db/dkd.sqlite` (6.6 MB relational) + `db/columnar/` (48 MB Parquet), from `db/schema.sql`
via `scripts/build_db.py`.

## 1. Why hybrid, concretely

The obvious relational design puts one row per (gene, sample, value) in a measurement table.
For this project that is 9,900 × 222 = 2.2M rows for the labelled samples alone, and 5.9M with
the non-diabetic CKD samples. It is the wrong shape for the only access pattern that matters —
*fetch a gene × sample block and multiply it* — which becomes an index scan over millions of
narrow rows, and it costs roughly 3× the storage of the dense array.

So the design splits by access pattern, not by data type:

| Tier | Holds | Store | Why |
|---|---|---|---|
| **Relational** | studies, subjects, samples, audit verdicts, gene space, analysis runs, results, candidate evidence, external evidence | SQLite here / PostgreSQL in production | needs joins, constraints, referential integrity, ad-hoc queries |
| **Columnar** | expression, effect-size and selection-frequency matrices | Parquet + zstd, one file per matrix | dense numeric blocks, read whole-column, compresses ~10× |
| **Catalogue** | `matrix` table pointing at Parquet URIs with checksums | relational | the two tiers are useless if you cannot find and verify the files |

SQLite rather than a server so the deliverable is one portable file. `db/schema.sql` is written
so the same DDL loads into PostgreSQL; three things change and each is marked inline:
`INTEGER PRIMARY KEY` → `SERIAL`, `TEXT` holding JSON → `JSONB`, and `GROUP_CONCAT` →
`STRING_AGG` in the views.

Not used, and why: a document store (the flexible metadata is small enough to sit in one JSON
column), a graph store (no traversal query in the workload), TileDB/Zarr (justified when
matrices exceed memory; 9,900 × 222 does not).

## 2. What the schema encodes that a directory of files cannot

Three findings from this project are enforced by the schema rather than documented in prose,
because prose does not stop the next query from getting it wrong.

**Subject identity is separate from sample identity.** The audit found the same patients in
several accessions and in both compartments. With only a `sample` table that overlap is
invisible and inflates every bootstrap estimate. `subject` is populated by re-deriving patient
IDs from sample titles with the audit's own rules, and `v_subject_reuse` turns the finding into
a query:

```
DN1    ERCB   in 4 studies: GSE104948, GSE104954, GSE47183, GSE99340
DN901  ERCB   in 3 studies: GSE104948, GSE47183, GSE99340
```

**Procurement is a first-class column.** `sample.procurement` records biopsy vs
tumour-nephrectomy vs living-donor, and `v_procurement_balance` crosses it with case status:

```
GSE104948  glomerulus  label=0  living_donor       n=21
GSE104948  glomerulus  label=0  tumor_nephrectomy  n=5
GSE104948  glomerulus  label=1  biopsy             n=12
```

Every cohort has cases as biopsies and controls as nephrectomy or donor tissue. That is the
confounder that produced the immediate-early artifact, and any analysis that cannot see this
column will rediscover it.

**Audit verdicts are data.** `study_audit.usable` and `exclusion_reason` carry the three
exclusions (GSE30122, GSE47183, GSE99340), and `v_usable_sample` filters them by default. A
downstream user cannot silently pull a duplicated cohort.

## 3. Tables

```
source_repository ──< study ──< sample >── subject
                       │         │
                       │         └─ label, compartment, procurement, platform, stage
                       ├─ study_audit   (usable, exclusion_reason, supersedes)
                       └─< matrix ──> array_store_uri (Parquet)

gene ──< platform_probe            (probe -> Entrez, multimapping flag)
gene ──< gene_space_member >── gene_space   (collapse_rule, n_genes, member_studies)

analysis_run ──< method_result       (AUROC, nested CV, null z, jaccards)
             └─< candidate_evidence  (RBS components + all three specificity contrasts)
gene ──< external_evidence           (KPMP cell type / regional protein, PubMed counts)
```

14 tables and 3 views. Contents as built:

| | |
|---|---|
| studies | 19 (15 GEO transcriptomics, 4 Metabolomics Workbench) |
| samples | 1,208 |
| subjects | 485 resolved, with reuse detectable |
| genes | 9,900 · probe mappings 72,334 |
| gene spaces | 2 (`discovery4` and `all6`, both 9,900 genes), collapse rule `max_mean` frozen |
| matrices in the columnar tier | 22 |
| analysis runs | 4 · method results 280 |
| candidate evidence | 30 genes · external evidence 80 rows |

## 4. Design choices worth defending

**`gene_space` is a first-class object.** The collapse rule and the intersection are frozen and
versioned, because a stability count is only comparable within one gene space. Re-deriving the
intersection later with a different collapse rule would silently invalidate every bootstrap
frequency already stored.

**`candidate_evidence` stores all three specificity contrasts** (`g_dkd_vs_control`,
`g_dkd_vs_otherckd`, `g_otherckd_vs_control`) rather than a single effect size and a pass/fail
flag. The third one is what distinguishes disease biology from procurement, so storing only the
first two would make the artifact undetectable from the database.

**`external_evidence` is generic** (`source`, `evidence_type`, `value_text`, `value_num`,
`p_adj`, `detail` JSON) rather than one table per resource. KPMP cell types, KPMP regional
proteomics and PubMed counts have nothing in common structurally, and a table per source would
mean schema changes every time a resource is added.

**Checksums on every Parquet file.** The relational tier claims a matrix has 9,900 rows; without
a checksum nothing detects the file being replaced.

## 5. Query examples

```sql
-- labelled glomerular samples from usable studies only
SELECT study_accession, COUNT(*) FROM v_usable_sample
WHERE compartment='glomerulus' AND label IN (0,1) GROUP BY study_accession;

-- candidates that survive specificity but are flagged as procurement-driven
SELECT g.symbol, c.g_dkd_vs_control, c.g_otherckd_vs_control, c.per_diagnosis_pattern
FROM candidate_evidence c JOIN gene g USING (entrez_id)
WHERE c.per_diagnosis_pattern LIKE '%donor%' ORDER BY c.rbs_rank;

-- method comparison straight out of the results tier
SELECT r.name, m.method, ROUND(AVG(m.external_auroc),3) auroc,
       ROUND(AVG(m.cross_fold_jaccard),3) xfold
FROM method_result m JOIN analysis_run r USING (run_id)
WHERE m.k=50 GROUP BY r.name, m.method ORDER BY xfold DESC;

-- everything known about one candidate
SELECT e.source, e.evidence_type, e.value_text, e.value_num, e.p_adj
FROM external_evidence e JOIN gene g USING (entrez_id) WHERE g.symbol='FMOD';
```

## 6. Not yet in the database

- **Clinical variables** (`subject.egfr`, `age_years`, `sex`, `proteinuria`) — columns exist but
  are empty. GEO metadata for these cohorts does not carry them; they would come from ERCB or
  KPMP controlled tiers.
- **Proteomics and metabolomics measurements** — the four Metabolomics Workbench studies are
  catalogued but their mwTab data blocks are not parsed into the columnar tier.
- **GWAS summary statistics** — not downloaded.
- **Single-cell matrices** — the KPMP h5ad and GSE131882 files are on disk but not registered;
  a 1.4M-nucleus matrix is the case where Parquet stops being the right choice and Zarr/TileDB
  starts.

## 7. Rebuild

```
python scripts/data_inventory.py     # inventory feeds the study table
python scripts/build_db.py --reset   # DDL + load + verification queries
```
