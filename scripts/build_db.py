#!/usr/bin/env python
"""Build the hybrid database - reference implementation.

Two tiers, one catalogue:

  RELATIONAL (SQLite here, PostgreSQL DDL in db/schema.sql)
      studies, subjects, samples, audit verdicts, gene space, analysis runs, results,
      candidate evidence, external cross-omics evidence.

  COLUMNAR (Parquet, one file per matrix)
      the expression / effect-size / selection-frequency matrices. Addressed from the
      relational tier via matrix.array_store_uri.

SQLite rather than a PostgreSQL server so the deliverable is a single portable file that runs
anywhere; the DDL is written so the same file loads into PostgreSQL with three noted changes.

The subject table is populated by re-deriving patient identifiers from sample titles using the
same rules as the audit, so v_subject_reuse reproduces the reused-cohort finding as a query
instead of as a claim in a document. That is the point of putting the audit in the schema.
"""
import os, re, sys, json, glob, hashlib, sqlite3, argparse
from datetime import datetime, timezone
import numpy as np
import pandas as pd

DB = 'db/dkd.sqlite'
PARQUET = 'db/columnar'
H = 'data/processed/harmonized'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def now():
    return datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def sha256(path, limit=64 << 20):
    h = hashlib.sha256()
    try:
        with open(path, 'rb') as f:
            while True:
                b = f.read(1 << 20)
                if not b:
                    break
                h.update(b)
                if f.tell() > limit:
                    break
    except OSError:
        return None
    return h.hexdigest()


# ---------------------------------------------------------------- subject derivation
ERCB_RE = re.compile(r'([A-Za-z0-9]+)-(Glom|Tub)-([A-Za-z&_]+?)_?(\d+)', re.I)
WORO_RE = re.compile(r'Human (?:Kidney|Glomeruli|Tubuli) (\d+)', re.I)


def derive_subject(study_acc, title, source):
    """Return (source_cohort, external_id) or (None, None).

    The ERCB and Woroniecka identifiers are the ones that revealed cohort reuse, so they are
    parsed here rather than left implicit. Compartment is deliberately stripped: the same
    patient contributes a glomerular and a tubulointerstitial sample and must resolve to ONE
    subject row, otherwise the overlap the audit found stays invisible.
    """
    m = ERCB_RE.search(title or '')
    if m:
        return 'ERCB', '%s%s' % (m.group(3).upper().strip('_&'), m.group(4))
    m = WORO_RE.search(title or '')
    if m:
        return 'Woroniecka', m.group(1)
    if study_acc == 'GSE142025':
        return 'PekingUniv', (title or '').strip()
    if study_acc == 'GSE96804':
        return 'GSE96804', (title or '').strip()
    return None, None


AUDIT = {
    'GSE30122': (0, 'duplicates GSE30528 + GSE30529 and re-runs 25 of the same controls',
                 'GSE30528;GSE30529'),
    'GSE47183': (0, 'DN subject IDs identical to GSE104948', 'GSE104948'),
    'GSE99340': (0, 'ERCB superset of GSE104948/GSE104954; ~51 samples are cell lines',
                 'GSE104948;GSE104954'),
}

STUDIES = [
    # accession, title, omics, assay, compartment-hint, path
    ('GSE142025', 'Early vs advanced DN kidney transcriptome', 'transcriptomics', 'rna-seq'),
    ('GSE30528', 'DKD glomeruli vs control glomeruli', 'transcriptomics', 'microarray'),
    ('GSE30529', 'DKD tubuli vs control tubuli', 'transcriptomics', 'microarray'),
    ('GSE96804', 'Exon-level expression profiling of diabetic nephropathy', 'transcriptomics', 'microarray'),
    ('GSE104948', 'ERCB glomerular transcriptome + living donors', 'transcriptomics', 'microarray'),
    ('GSE104954', 'ERCB tubulointerstitial transcriptome + living donors', 'transcriptomics', 'microarray'),
    ('GSE30122', 'Transcriptome analysis of human DKD (superseries)', 'transcriptomics', 'microarray'),
    ('GSE47183', 'In silico nano-dissection, ERCB glomeruli', 'transcriptomics', 'microarray'),
    ('GSE99340', 'Renal cell-type-specific hypoxia dysregulation', 'transcriptomics', 'microarray'),
    ('GSE131882', 'Single-cell landscape of early human diabetic nephropathy', 'transcriptomics', 'snrna-seq'),
    ('GSE175759', 'Tubulointerstitial RNA-seq', 'transcriptomics', 'rna-seq'),
    ('GSE142153', 'Human PBMC: healthy vs DN vs ESRD', 'transcriptomics', 'microarray'),
    ('GSE1009', 'Diabetic nephropathy', 'transcriptomics', 'microarray'),
    ('GSE111154', 'Early diabetic nephropathy', 'transcriptomics', 'microarray'),
    ('GSE20602', 'Human nephrosclerosis hypoxia glomerulopathy', 'transcriptomics', 'microarray'),
]

METAB = [('ST000691', 'Metabolomics of diabetic nephropathy', 'lc-ms'),
         ('ST004483', 'Urine GC-MS metabolomic signatures of DKD', 'gc-ms'),
         ('ST004442', 'Plasma GC-MS metabolomic signatures of DKD', 'gc-ms'),
         ('ST002145', 'ChREBP links mitochondrial lipidomes to DN progression', 'lc-ms')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--db', default=DB)
    ap.add_argument('--reset', action='store_true')
    args = ap.parse_args()

    os.makedirs('db', exist_ok=True)
    os.makedirs(PARQUET, exist_ok=True)
    if args.reset and os.path.exists(args.db):
        os.remove(args.db)

    con = sqlite3.connect(args.db)
    con.executescript(open('db/schema.sql', encoding='utf-8').read())
    cur = con.cursor()

    # ---------------------------------------------------------------- repositories
    repos = [(1, 'GEO', 'https://www.ncbi.nlm.nih.gov/geo/', 'open', 'NCBI Gene Expression Omnibus'),
             (2, 'KPMP', 'https://atlas.kpmp.org/', 'open',
              'open tier via GraphQL + CELLxGENE mirror; controlled tier needs a DAR'),
             (3, 'MetabolomicsWorkbench', 'https://www.metabolomicsworkbench.org/', 'open', ''),
             (4, 'PubMed', 'https://pubmed.ncbi.nlm.nih.gov/', 'open', 'literature counts'),
             (5, 'GWASCatalog', 'https://www.ebi.ac.uk/gwas/', 'open', 'not downloaded yet')]
    cur.executemany('INSERT INTO source_repository VALUES (?,?,?,?,?)', repos)

    # ---------------------------------------------------------------- studies
    inv = pd.read_csv('results/data_inventory.tsv', sep='\t')
    inv_path = dict(zip(inv['asset'], inv['path'].fillna('')))
    sid = {}
    n = 0
    for acc, title, omics, assay in STUDIES:
        n += 1
        sid[acc] = n
        path = 'data/raw/geo/%s' % acc
        cur.execute('INSERT INTO study (study_id,repository_id,accession,title,omics_type,assay,'
                    'organism,downloaded_at,local_path,license) VALUES (?,?,?,?,?,?,?,?,?,?)',
                    (n, 1, acc, title, omics, assay, 'Homo sapiens', now(), path, 'GEO open'))
        u, reason, sup = AUDIT.get(acc, (1, None, None))
        cur.execute('INSERT INTO study_audit VALUES (?,?,?,?,?,?)',
                    (n, u, reason, sup, now(), 'scripts/parse_series_matrix.py'))
    for acc, title, assay in METAB:
        n += 1
        sid[acc] = n
        cur.execute('INSERT INTO study (study_id,repository_id,accession,title,omics_type,assay,'
                    'organism,downloaded_at,local_path,license) VALUES (?,?,?,?,?,?,?,?,?,?)',
                    (n, 3, acc, title, 'metabolomics', assay, 'Homo sapiens', now(),
                     'data/raw/metabolomics/%s' % acc, 'CC BY 4.0'))
        cur.execute('INSERT INTO study_audit VALUES (?,?,?,?,?,?)',
                    (n, 1, 'downloaded, not yet analysed', None, now(), None))
    log('studies: %d' % n)

    # ---------------------------------------------------------------- samples + subjects
    sinv = pd.read_csv('data/metadata/sample_inventory.tsv', sep='\t')
    pheno = {}
    for f in glob.glob(os.path.join(H, '*_pheno.tsv')):
        acc = os.path.basename(f).split('_')[0]
        pheno[acc] = pd.read_csv(f, sep='\t', dtype={'sample': str}).set_index('sample')

    subj_id = {}
    ns = nsub = 0
    for _, r in sinv.iterrows():
        acc = r['gse']
        if acc not in sid:
            continue
        cohort, ext = derive_subject(acc, str(r['title']), None)
        sub_pk = None
        if cohort:
            key = (cohort, ext)
            if key not in subj_id:
                nsub += 1
                subj_id[key] = nsub
                diag = None
                m = re.search(r'diagnosis:\s*([^|]+)', str(r['characteristics']))
                if m:
                    diag = m.group(1).strip()
                elif 'LD' in str(r['title']).upper():
                    diag = 'living_donor'
                cur.execute('INSERT INTO subject (subject_id,external_id,source_cohort,diagnosis,'
                            'is_diabetic) VALUES (?,?,?,?,?)',
                            (nsub, ext, cohort, diag,
                             1 if (diag and 'iabet' in diag) or ext.startswith('DN') else None))
            sub_pk = subj_id[key]

        ph = pheno.get(acc)
        row = ph.loc[r['gsm']] if (ph is not None and r['gsm'] in ph.index) else None
        ns += 1
        cur.execute('INSERT INTO sample (sample_id,study_id,subject_id,accession,title,compartment,'
                    'procurement,platform,label,disease_group,stage,raw_characteristics) '
                    'VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
                    (ns, sid[acc], sub_pk, r['gsm'], r['title'],
                     (row['compartment'] if row is not None else None),
                     (row['control_type'] if row is not None and isinstance(row['control_type'], str)
                      and row['control_type'] else ('biopsy' if row is not None else None)),
                     r['platform'],
                     int(row['label']) if row is not None else None,
                     (row['group'] if row is not None else None),
                     (row['stage'] if row is not None and isinstance(row['stage'], str) else None),
                     json.dumps({'characteristics': r['characteristics'], 'source': r['source']})))
    log('samples: %d   subjects: %d' % (ns, nsub))

    # ---------------------------------------------------------------- genes + gene space
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    cur.executemany('INSERT OR IGNORE INTO gene (entrez_id,symbol) VALUES (?,?)',
                    [(int(g), s if isinstance(s, str) else None)
                     for g, s in zip(gs['entrez_id'], gs['symbol'])])
    p2g = pd.read_csv('data/processed/probe2gene.tsv', sep='\t')
    plat_map = {'GSE30528': 'GPL571', 'GSE30529': 'GPL571', 'GSE96804': 'GPL17586',
                'GSE104948': 'ENTREZG', 'GSE104954': 'ENTREZG', 'GSE142025': 'SYMBOL'}
    seen = set()
    rows = []
    k = 0
    for _, r in p2g.iterrows():
        plat = plat_map.get(r['cohort'], r['cohort'])
        key = (plat, str(r['probe_id']))
        if key in seen:
            continue
        seen.add(key)
        k += 1
        rows.append((k, plat, str(r['probe_id']), int(r['entrez_id']), 0))
    cur.executemany('INSERT INTO platform_probe VALUES (?,?,?,?,?)', rows)
    log('genes: %d   probe mappings: %d' % (len(gs), len(rows)))

    for i, (name, fn) in enumerate([('discovery4', 'gene_space_discovery4.tsv'),
                                    ('all6', 'gene_space_all6.tsv')], start=1):
        d = pd.read_csv('data/processed/%s' % fn, sep='\t', dtype=str)
        members = ({'GSE142025', 'GSE30528', 'GSE96804', 'GSE104948'} if name == 'discovery4'
                   else {'GSE142025', 'GSE30528', 'GSE96804', 'GSE104948', 'GSE30529', 'GSE104954'})
        cur.execute('INSERT INTO gene_space VALUES (?,?,?,?,?,?,?,?)',
                    (i, name, 'entrez', 'max_mean', len(d), json.dumps(sorted(members)),
                     now(), 'scripts/build_gene_space.py'))
        cur.executemany('INSERT INTO gene_space_member VALUES (?,?)',
                        [(i, int(g)) for g in d['entrez_id']])
    log('gene spaces: 2')

    # ---------------------------------------------------------------- columnar tier
    nm = 0
    for f in sorted(glob.glob(os.path.join(H, '*_expr.tsv'))):
        acc = os.path.basename(f).split('_')[0]
        df = pd.read_csv(f, sep='\t', index_col=0)
        out = os.path.join(PARQUET, '%s_expr.parquet' % acc)
        df.reset_index().to_parquet(out, index=False, compression='zstd')
        nm += 1
        cur.execute('INSERT INTO matrix (matrix_id,study_id,gene_space_id,kind,normalisation,'
                    'n_rows,n_cols,array_store_uri,checksum_sha256,created_at) '
                    'VALUES (?,?,?,?,?,?,?,?,?,?)',
                    (nm, sid.get(acc), 2, 'expression', 'log-scale as published (not z-scored)',
                     df.shape[0], df.shape[1], out, sha256(out), now()))
    for src, kind in [('data/processed/effect_sizes_hedges_g.tsv', 'effect_size'),
                      ('data/processed/effect_sizes_auc.tsv', 'effect_size')]:
        if not os.path.exists(src):
            continue
        df = pd.read_csv(src, sep='\t', index_col=0)
        out = os.path.join(PARQUET, os.path.basename(src).replace('.tsv', '.parquet'))
        df.reset_index().to_parquet(out, index=False, compression='zstd')
        nm += 1
        cur.execute('INSERT INTO matrix (matrix_id,gene_space_id,kind,n_rows,n_cols,'
                    'array_store_uri,checksum_sha256,created_at) VALUES (?,?,?,?,?,?,?,?)',
                    (nm, 2, kind, df.shape[0], df.shape[1], out, sha256(out), now()))
    for d in ('results/baselines', 'results/baselines_all4'):
        for f in glob.glob(os.path.join(d, 'freq_*.tsv')):
            df = pd.read_csv(f, sep='\t', index_col=0)
            out = os.path.join(PARQUET, '%s_%s.parquet'
                               % (os.path.basename(d), os.path.basename(f)[:-4]))
            df.reset_index().to_parquet(out, index=False, compression='zstd')
            nm += 1
            cur.execute('INSERT INTO matrix (matrix_id,gene_space_id,kind,n_rows,n_cols,'
                        'array_store_uri,checksum_sha256,created_at) VALUES (?,?,?,?,?,?,?,?)',
                        (nm, 2, 'selection_frequency', df.shape[0], df.shape[1], out,
                         sha256(out), now()))
    log('columnar matrices: %d' % nm)

    # ---------------------------------------------------------------- analysis runs + results
    runs = {}
    rid = 0
    for name, path, cohorts in [('baselines_glom3', 'results/baselines',
                                 'GSE30528,GSE96804,GSE104948'),
                                ('baselines_all4', 'results/baselines_all4',
                                 'GSE30528,GSE96804,GSE104948,GSE142025'),
                                ('rbs_glom3', 'results/proposed_glom',
                                 'GSE30528,GSE96804,GSE104948'),
                                ('rbs_all4', 'results/proposed_all4',
                                 'GSE30528,GSE96804,GSE104948,GSE142025')]:
        f = os.path.join(path, 'results.tsv')
        if not os.path.exists(f):
            continue
        rid += 1
        runs[name] = rid
        script = ('scripts/run_baselines.py' if 'baselines' in name else 'scripts/run_proposed.py')
        cur.execute('INSERT INTO analysis_run (run_id,name,script,parameters,gene_space_id,'
                    'started_at,finished_at) VALUES (?,?,?,?,?,?,?)',
                    (rid, name, script,
                     json.dumps({'B': 200, 'k': 50, 'cohorts': cohorts.split(','),
                                 'base': 'relieff', 'agg': 'geomean'}),
                     2, now(), now()))
        df = pd.read_csv(f, sep='\t')
        cmp_file = ('results/comparison_glom3.tsv' if 'glom3' in name
                    else 'results/comparison_all4.tsv')
        cmp_df = pd.read_csv(cmp_file, sep='\t') if os.path.exists(cmp_file) else None
        for _, r in df.iterrows():
            z = pct = xf = None
            if cmp_df is not None:
                m = cmp_df[(cmp_df['method'] == r['method']) & (cmp_df['K'] == r['K'])]
                if len(m):
                    z = float(m['null_z'].iloc[0]); pct = float(m['null_pct'].iloc[0])
                    xf = float(m['cross_fold_jaccard'].iloc[0])
            cur.execute('INSERT INTO method_result (run_id,method,held_out_study,k,external_auroc,'
                        'nested_cv_auroc,null_z,null_percentile,bootstrap_jaccard,'
                        'cross_fold_jaccard) VALUES (?,?,?,?,?,?,?,?,?,?)',
                        (rid, r['method'], r['held_out'], int(r['K']),
                         float(r['external_auroc']) if pd.notna(r['external_auroc']) else None,
                         float(r['nested_cv_auroc']) if pd.notna(r['nested_cv_auroc']) else None,
                         z, pct, float(r['jaccard']) if pd.notna(r['jaccard']) else None, xf))
    log('analysis runs: %d' % rid)

    # ---------------------------------------------------------------- candidate evidence
    sym2id = {s: int(g) for g, s in zip(gs['entrez_id'], gs['symbol']) if isinstance(s, str)}
    mt = pd.read_csv('results/master_candidate_table.tsv', sep='\t')
    run = runs.get('rbs_all4', rid)
    nc = 0
    for _, r in mt.iterrows():
        gid = sym2id.get(r['gene'])
        if gid is None:
            continue
        nc += 1
        cur.execute('INSERT OR REPLACE INTO candidate_evidence VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                    (run, gid, float(r['mean_rank']), None, None,
                     float(r['mean_S']) if pd.notna(r['mean_S']) else None,
                     float(r['mean_R']) if pd.notna(r['mean_R']) else None,
                     float(r['mean_P']) if pd.notna(r['mean_P']) else None,
                     float(r['g_DKD_vs_control']), float(r['g_DKD_vs_otherCKD']), None,
                     1, 0, str(r['pattern_glom']), int(bool(r['above_perm_ceiling']))))
    log('candidate evidence rows: %d' % nc)

    # ---------------------------------------------------------------- external evidence
    ne = 0
    kp = pd.read_csv('results/kpmp_corroboration.tsv', sep='\t')
    for _, r in kp.iterrows():
        gid = sym2id.get(r['symbol'])
        if gid is None:
            continue
        if isinstance(r.get('sn_top_celltype'), str) and r['sn_top_celltype']:
            ne += 1
            cur.execute('INSERT INTO external_evidence (evidence_id,entrez_id,source,evidence_type,'
                        'value_text,value_num,retrieved_at) VALUES (?,?,?,?,?,?,?)',
                        (ne, gid, 'KPMP_sn', 'cell_type', r['sn_top_celltype'],
                         float(r['sn_top_foldchange']) if pd.notna(r.get('sn_top_foldchange')) else None,
                         now()))
        if pd.notna(r.get('protein_glom_vs_TI_FC')):
            ne += 1
            cur.execute('INSERT INTO external_evidence (evidence_id,entrez_id,source,evidence_type,'
                        'value_num,p_adj,retrieved_at) VALUES (?,?,?,?,?,?,?)',
                        (ne, gid, 'KPMP_rp', 'compartment_protein',
                         float(r['protein_glom_vs_TI_FC']),
                         float(r['protein_adjP']) if pd.notna(r.get('protein_adjP')) else None, now()))
    lit = pd.read_csv('results/literature_review.tsv', sep='\t')
    for _, r in lit.iterrows():
        gid = sym2id.get(r['gene'])
        if gid is None:
            continue
        ne += 1
        cur.execute('INSERT INTO external_evidence (evidence_id,entrez_id,source,evidence_type,'
                    'value_text,value_num,retrieved_at,detail) VALUES (?,?,?,?,?,?,?,?)',
                    (ne, gid, 'PubMed', 'literature_count', r['verdict'], float(r['n_dkd']), now(),
                     json.dumps({'n_any': int(r['n_any']), 'n_kidney': int(r['n_kidney']),
                                 'n_dkd_biomarker': int(r['n_dkd_biomarker']),
                                 'pmids': str(r['dkd_pmids'])})))
    log('external evidence rows: %d' % ne)

    con.commit()

    # ---------------------------------------------------------------- verify
    log('\n=== verification queries ===')
    q = ('SELECT external_id, source_cohort, n_studies, studies FROM v_subject_reuse '
         "WHERE source_cohort='ERCB' AND external_id LIKE 'DN%' ORDER BY n_studies DESC LIMIT 8")
    log('\nsubjects reused across studies (the audit finding, as a query):')
    for row in cur.execute(q):
        log('  %-8s %-12s in %d studies: %s' % row)
    log('\nprocurement crossed with case status:')
    for row in cur.execute('SELECT accession, compartment, label, procurement, n '
                           'FROM v_procurement_balance ORDER BY accession, label LIMIT 12'):
        log('  %-11s %-20s label=%-3s %-20s n=%d' % row)
    log('\nbest cross-fold stability per run:')
    for row in cur.execute(
            'SELECT r.name, m.method, m.k, ROUND(AVG(m.external_auroc),3), '
            'ROUND(MAX(m.cross_fold_jaccard),3) FROM method_result m '
            'JOIN analysis_run r ON r.run_id=m.run_id WHERE m.k=50 '
            'GROUP BY r.name, m.method ORDER BY 5 DESC LIMIT 6'):
        log('  %-16s %-12s K=%d  AUROC=%.3f  cross-fold=%.3f' % row)

    sizes = cur.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    log('\ntables: %d   db size: %.2f MB   parquet: %.2f MB'
        % (len(sizes), os.path.getsize(args.db) / 1e6,
           sum(os.path.getsize(os.path.join(PARQUET, f))
               for f in os.listdir(PARQUET)) / 1e6))
    con.close()
    log('\nbuilt %s' % args.db)


if __name__ == '__main__':
    main()
