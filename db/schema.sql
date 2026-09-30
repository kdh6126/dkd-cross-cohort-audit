-- =====================================================================================
-- DKD multi-omics hybrid database - relational tier
-- =====================================================================================
-- Written for PostgreSQL; the reference implementation in scripts/build_db.py loads the
-- same DDL into SQLite so the whole thing runs with no server. Differences are confined to
-- three things, each marked inline: SERIAL vs INTEGER PRIMARY KEY, JSONB vs TEXT, and
-- partial indexes.
--
-- Design decision: the relational tier stores everything EXCEPT the expression matrices.
-- A 9,900 x 222 float matrix as 2.2M rows in a measurement table is the obvious relational
-- design and the wrong one - it makes the common operation (fetch a gene x sample block)
-- an index scan over millions of rows, and it triples storage. Matrices live in the
-- columnar tier (Parquet) and are addressed from here by array_store_uri. That split is
-- what makes this a hybrid rather than a relational database with a slow spot.
-- =====================================================================================

-- ------------------------------------------------------------------ provenance / sources
CREATE TABLE source_repository (
    repository_id     INTEGER PRIMARY KEY,
    name              TEXT NOT NULL UNIQUE,      -- GEO, KPMP, MetabolomicsWorkbench, PRIDE, GWASCatalog
    base_url          TEXT,
    access_tier       TEXT CHECK (access_tier IN ('open', 'registered', 'controlled')),
    notes             TEXT
);

CREATE TABLE study (
    study_id          INTEGER PRIMARY KEY,
    repository_id     INTEGER NOT NULL REFERENCES source_repository(repository_id),
    accession         TEXT NOT NULL,             -- GSE104948, ST004483, PXD058790
    title             TEXT,
    omics_type        TEXT NOT NULL,             -- transcriptomics | proteomics | metabolomics | genomics
    assay             TEXT,                      -- microarray | rna-seq | snrna-seq | lc-ms | gc-ms
    organism          TEXT DEFAULT 'Homo sapiens',
    n_samples         INTEGER,
    publication_pmid  TEXT,
    downloaded_at     TEXT,                      -- ISO8601
    local_path        TEXT,
    checksum_sha256   TEXT,
    license           TEXT,
    UNIQUE (repository_id, accession)
);

-- Audit verdicts live in the database, not in a README, because every downstream query
-- must be able to exclude reused cohorts without a human remembering to.
CREATE TABLE study_audit (
    study_id          INTEGER PRIMARY KEY REFERENCES study(study_id),
    usable            INTEGER NOT NULL,          -- 0/1 (BOOLEAN in PostgreSQL)
    exclusion_reason  TEXT,
    supersedes        TEXT,                      -- accession this one duplicates, if any
    audited_at        TEXT,
    audit_script      TEXT
);

-- ------------------------------------------------------------------ subjects and samples
-- subject is separate from sample because the audit found the same patients appearing in
-- several accessions and in both compartments; without a subject table that overlap is
-- invisible and inflates every bootstrap estimate.
CREATE TABLE subject (
    subject_id        INTEGER PRIMARY KEY,
    external_id       TEXT NOT NULL,             -- e.g. ERCB 'DN901', Woroniecka 'Kidney 62'
    source_cohort     TEXT NOT NULL,             -- ERCB | Woroniecka | PekingUniv | KPMP
    diagnosis         TEXT,                      -- DN, FSGS, IgA, living_donor, tumor_nephrectomy, ...
    is_diabetic       INTEGER,
    ancestry          TEXT,
    sex               TEXT,
    age_years         REAL,
    egfr              REAL,
    proteinuria       REAL,
    UNIQUE (source_cohort, external_id)
);

CREATE TABLE sample (
    sample_id         INTEGER PRIMARY KEY,
    study_id          INTEGER NOT NULL REFERENCES study(study_id),
    subject_id        INTEGER REFERENCES subject(subject_id),
    accession         TEXT NOT NULL,             -- GSM757014
    title             TEXT,
    compartment       TEXT,                      -- glomerulus | tubulointerstitium | whole_cortex | blood
    procurement       TEXT,                      -- biopsy | tumor_nephrectomy | living_donor
    platform          TEXT,                      -- GPL571, GPL22945, Illumina HiSeq 4000
    label             INTEGER,                   -- 1 = DKD, 0 = control, -1 = other kidney disease
    disease_group     TEXT,
    stage             TEXT,                      -- Early_DN | Advanced_DN
    raw_characteristics TEXT,                    -- JSONB in PostgreSQL
    UNIQUE (study_id, accession)
);

-- procurement is stored explicitly because it turned out to be the dominant confounder:
-- cases are biopsies and controls are nephrectomy/donor tissue in every cohort, and the
-- immediate-early module tracks that rather than diabetes. Any analysis that cannot see
-- this column will rediscover the artifact.
CREATE INDEX idx_sample_study     ON sample(study_id);
CREATE INDEX idx_sample_subject   ON sample(subject_id);
CREATE INDEX idx_sample_label     ON sample(label);
CREATE INDEX idx_sample_compart   ON sample(compartment);

-- ------------------------------------------------------------------ feature space
CREATE TABLE gene (
    entrez_id         INTEGER PRIMARY KEY,
    symbol            TEXT,
    name              TEXT,
    chromosome        TEXT
);

CREATE TABLE platform_probe (
    probe_uid         INTEGER PRIMARY KEY,
    platform          TEXT NOT NULL,
    probe_id          TEXT NOT NULL,
    entrez_id         INTEGER REFERENCES gene(entrez_id),
    is_multimapping   INTEGER DEFAULT 0,
    UNIQUE (platform, probe_id)
);
CREATE INDEX idx_probe_gene ON platform_probe(entrez_id);

-- The analysis gene space is a first-class object: the collapse rule and the intersection
-- are frozen, and every stability count is only comparable within one gene_space_id.
CREATE TABLE gene_space (
    gene_space_id     INTEGER PRIMARY KEY,
    name              TEXT NOT NULL UNIQUE,      -- 'discovery4', 'all6'
    hub_id_type       TEXT NOT NULL DEFAULT 'entrez',
    collapse_rule     TEXT NOT NULL,             -- 'max_mean'
    n_genes           INTEGER NOT NULL,
    member_studies    TEXT,                      -- JSONB in PostgreSQL
    created_at        TEXT,
    build_script      TEXT
);

CREATE TABLE gene_space_member (
    gene_space_id     INTEGER NOT NULL REFERENCES gene_space(gene_space_id),
    entrez_id         INTEGER NOT NULL REFERENCES gene(entrez_id),
    PRIMARY KEY (gene_space_id, entrez_id)
);

-- ------------------------------------------------------------------ columnar tier pointers
-- One row per matrix. The bytes live in Parquet; this table is the catalogue that makes
-- them discoverable and reproducible.
CREATE TABLE matrix (
    matrix_id         INTEGER PRIMARY KEY,
    study_id          INTEGER REFERENCES study(study_id),
    gene_space_id     INTEGER REFERENCES gene_space(gene_space_id),
    kind              TEXT NOT NULL,             -- expression | effect_size | selection_frequency
    normalisation     TEXT,                      -- 'per-cohort gene-wise z-score'
    n_rows            INTEGER,
    n_cols            INTEGER,
    array_store_uri   TEXT NOT NULL,             -- parquet path
    checksum_sha256   TEXT,
    created_at        TEXT
);

-- ------------------------------------------------------------------ analysis layer
CREATE TABLE analysis_run (
    run_id            INTEGER PRIMARY KEY,
    name              TEXT NOT NULL,
    script            TEXT NOT NULL,
    git_or_file_hash  TEXT,
    parameters        TEXT,                      -- JSONB in PostgreSQL
    gene_space_id     INTEGER REFERENCES gene_space(gene_space_id),
    started_at        TEXT,
    finished_at       TEXT
);

CREATE TABLE method_result (
    result_id         INTEGER PRIMARY KEY,
    run_id            INTEGER NOT NULL REFERENCES analysis_run(run_id),
    method            TEXT NOT NULL,
    held_out_study    TEXT,
    k                 INTEGER,
    external_auroc    REAL,
    nested_cv_auroc   REAL,
    null_z            REAL,                      -- vs the random-signature null
    null_percentile   REAL,
    bootstrap_jaccard REAL,
    cross_fold_jaccard REAL
);
CREATE INDEX idx_result_run ON method_result(run_id);

-- Candidate-level evidence, one row per gene per run, so a candidate can never be quoted
-- without the contrast that supports it.
CREATE TABLE candidate_evidence (
    run_id            INTEGER NOT NULL REFERENCES analysis_run(run_id),
    entrez_id         INTEGER NOT NULL REFERENCES gene(entrez_id),
    rbs_rank          REAL,
    rbs_score         REAL,
    comp_importance   REAL,
    comp_stability    REAL,
    comp_reproducibility REAL,
    comp_perturbation REAL,
    g_dkd_vs_control  REAL,
    g_dkd_vs_otherckd REAL,
    g_otherckd_vs_control REAL,
    passes_specificity INTEGER,
    procurement_flag  INTEGER,
    per_diagnosis_pattern TEXT,                  -- A unique-to-DN | B DN-dominant | C donor-driven
    above_perm_ceiling INTEGER,
    PRIMARY KEY (run_id, entrez_id)
);

-- ------------------------------------------------------------------ cross-omics evidence
CREATE TABLE external_evidence (
    evidence_id       INTEGER PRIMARY KEY,
    entrez_id         INTEGER NOT NULL REFERENCES gene(entrez_id),
    source            TEXT NOT NULL,             -- KPMP_sn | KPMP_rp | PubMed | MetabolomicsWorkbench
    evidence_type     TEXT NOT NULL,             -- cell_type | compartment_protein | literature_count
    value_text        TEXT,
    value_num         REAL,
    p_adj             REAL,
    retrieved_at      TEXT,
    detail            TEXT                       -- JSONB in PostgreSQL
);
CREATE INDEX idx_evidence_gene ON external_evidence(entrez_id);

-- ------------------------------------------------------------------ convenience views
CREATE VIEW v_usable_sample AS
SELECT s.*, st.accession AS study_accession, st.omics_type, st.assay,
       sub.external_id AS subject_external_id, sub.diagnosis AS subject_diagnosis
FROM sample s
JOIN study st        ON st.study_id = s.study_id
LEFT JOIN study_audit a ON a.study_id = s.study_id
LEFT JOIN subject sub   ON sub.subject_id = s.subject_id
WHERE COALESCE(a.usable, 1) = 1;

-- Subjects appearing in more than one study: the query the audit had to run by hand.
CREATE VIEW v_subject_reuse AS
SELECT sub.subject_id, sub.external_id, sub.source_cohort,
       COUNT(DISTINCT s.study_id) AS n_studies,
       GROUP_CONCAT(DISTINCT st.accession) AS studies
FROM subject sub
JOIN sample s ON s.subject_id = sub.subject_id
JOIN study st ON st.study_id = s.study_id
GROUP BY sub.subject_id, sub.external_id, sub.source_cohort
HAVING COUNT(DISTINCT s.study_id) > 1;

-- Case/control counts crossed with procurement, so the confounder is one query away.
CREATE VIEW v_procurement_balance AS
SELECT st.accession, s.compartment, s.label, s.procurement, COUNT(*) AS n
FROM sample s JOIN study st ON st.study_id = s.study_id
WHERE s.label IN (0, 1)
GROUP BY st.accession, s.compartment, s.label, s.procurement;
