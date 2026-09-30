#!/usr/bin/env python
"""Why does the immediate-early module behave differently across datasets?

Current state of the evidence, which is not consistent:

    ERCB microarray        DKD vs other-CKD biopsies   IEG ~ 0            supports artifact
    KPMP single-nucleus    DKD vs non-diabetic CKD     p = 0.21           supports artifact
    GSE175759 RNA-seq      any CKD vs nephrectomy      g = -2.26          supports artifact
    GSE162830 RNA-seq      DKD vs ING (both biopsies)  p = 2.1e-06        REFUTES it

GSE162830 has the best design of the four - the comparator is procurement-matched AND
histology-matched - and it is the one that disagrees. Declaring the artifact replicated while
that stands would be wrong.

The direct test: if the module tracks procurement, then across every dataset the module level
should order by how the tissue was obtained (biopsy low, nephrectomy/donor high) and should NOT
order by disease. Put every group from every dataset on one axis and look.

Each dataset is normalised to its own reference group before comparison, because absolute
expression is not comparable across platforms.
"""
import os, sys, gzip, glob, argparse
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

IEG = ['FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'ATF3', 'DUSP1', 'ZFP36', 'NR4A2', 'BTG2',
       'EGR2', 'EGR3', 'IER2', 'KLF2', 'KLF4', 'SOCS3', 'DUSP2', 'JUND', 'NR4A1']
H = 'data/processed/harmonized'
OUT = 'results/ieg_heterogeneity.tsv'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def sym2entrez():
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    return {s: g for g, s in zip(gs['entrez_id'], gs['symbol']) if isinstance(s, str)}


def module_score(X, ieg_ids):
    """Mean z-score of the module per sample, z computed within the dataset."""
    Z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1).replace(0, np.nan)
    cols = [c for c in ieg_ids if c in Z.columns]
    return Z[cols].mean(axis=1), len(cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args()
    s2e = sym2entrez()
    ieg_ids = [s2e[s] for s in IEG if s in s2e]

    rows = []

    # ---------------------------------------------------------------- harmonised cohorts
    for f in sorted(glob.glob(os.path.join(H, '*_expr.tsv'))):
        acc = os.path.basename(f).split('_')[0]
        e = pd.read_csv(f, sep='\t', index_col=0)
        e.index = e.index.astype(str)
        p = pd.read_csv(os.path.join(H, '%s_pheno.tsv' % acc), sep='\t', dtype={'sample': str})
        X = e.T
        X.index = X.index.astype(str)
        p = p.set_index('sample').loc[X.index]
        sc, n = module_score(X, ieg_ids)
        for grp, sub in sc.groupby(p['group'].values):
            proc = p.loc[sub.index, 'control_type'].fillna('').replace('', 'biopsy').mode()
            rows.append(dict(dataset=acc, group=str(grp), n=len(sub),
                             procurement=(proc.iloc[0] if len(proc) else 'biopsy'),
                             ieg_z=float(sub.mean()), n_ieg=n))

    # ---------------------------------------------------------------- GSE162830
    d = pd.read_csv('data/raw/geo/GSE162830/GSE162830_ING_quantile_normalized_final.csv.gz',
                    index_col=0)
    lg = np.log2(d / d.sum(axis=0) * 1e6 + 1)
    lg = lg.loc[[g for g in lg.index if g in s2e]]
    lg.index = [s2e[g] for g in lg.index]
    lg = lg[~lg.index.duplicated()]
    X = lg.T
    sc, n = module_score(X, ieg_ids)
    grp = pd.Series(['DN' if c.startswith('DIA') else 'ING' if c.startswith('ING') else 'reference'
                     for c in X.index], index=X.index)
    proc = {'DN': 'biopsy', 'ING': 'biopsy', 'reference': 'tumor_nephrectomy'}
    for g_, sub in sc.groupby(grp.values):
        rows.append(dict(dataset='GSE162830', group=str(g_), n=len(sub),
                         procurement=proc[g_], ieg_z=float(sub.mean()), n_ieg=n))

    # ---------------------------------------------------------------- GSE175759
    try:
        from validate_gse175759 import ensembl_to_entrez, load_pheno, build_matrix
        Xr = build_matrix(ensembl_to_entrez())
        ph = load_pheno().set_index('gsm')
        Xr = Xr[[c for c in Xr.columns if c in ph.index]]
        ph = ph.loc[Xr.columns]
        keep = ~ph['outlier'].fillna(False)
        Xr, ph = Xr.loc[:, keep.values], ph[keep]
        XT = Xr.T
        sc, n = module_score(XT, ieg_ids)
        g2 = np.where(ph['diagnosis'].str.lower().str.contains('control'), 'reference',
                      np.where(ph['diagnosis'].str.lower().str.contains('diabetic'), 'DN', 'otherCKD'))
        pr = {'reference': 'tumor_nephrectomy', 'DN': 'biopsy', 'otherCKD': 'biopsy'}
        for g_, sub in sc.groupby(g2):
            rows.append(dict(dataset='GSE175759', group=str(g_), n=len(sub),
                             procurement=pr[g_], ieg_z=float(sub.mean()), n_ieg=n))
    except Exception as e:
        log('  ! GSE175759 skipped: %s' % str(e)[:90])

    # ---------------------------------------------------------------- KPMP pseudobulk
    try:
        z = np.load('results/kpmp_singlecell/expr_cache.npz', allow_pickle=True)
        Xk, gk = z['X'], list(z['genes'])
        import h5py
        f = h5py.File('data/raw/kpmp/KPMP_v1.5_snRNA_human_kidney.h5ad', 'r')

        def cat(k):
            gg = f['obs/' + k]
            c = np.array([x.decode() if isinstance(x, bytes) else x for x in gg['categories'][:]])
            return c[gg['codes'][:]]
        dcat, dm, donor = cat('disease_category'), cat('diabetes_history'), cat('donor_id')
        grp2 = np.where((dcat == 'CKD') & (dm == 'Yes'), 'DN',
                        np.where((dcat == 'CKD') & (dm == 'No'), 'otherCKD',
                                 np.where(dcat == 'Healthy_reference_tissue', 'reference',
                                          np.where(dcat == 'Tumor_nephrectomy', 'nephrectomy', 'other'))))
        dfk = pd.DataFrame(Xk, columns=gk)
        cols = [c for c in IEG if c in dfk.columns]
        pb = dfk.groupby(donor)[cols].mean()
        gper = pd.Series(grp2).groupby(donor).first()
        Zk = (pb - pb.mean()) / pb.std(ddof=1)
        sck = Zk.mean(axis=1)
        pr2 = {'DN': 'biopsy', 'otherCKD': 'biopsy', 'reference': 'donor/nephrectomy',
               'nephrectomy': 'tumor_nephrectomy'}
        for g_ in ('DN', 'otherCKD', 'reference', 'nephrectomy'):
            sel = sck[gper.reindex(sck.index).values == g_]
            if len(sel):
                rows.append(dict(dataset='KPMP_snRNA', group=g_, n=len(sel),
                                 procurement=pr2[g_], ieg_z=float(sel.mean()), n_ieg=len(cols)))
    except Exception as e:
        log('  ! KPMP skipped: %s' % str(e)[:90])

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)

    log('\n=== immediate-early module level per group, z-scored within each dataset ===')
    log('  %-12s %-18s %5s %-20s %s' % ('dataset', 'group', 'n', 'procurement', 'IEG z'))
    for ds, sub in df.groupby('dataset'):
        for _, r in sub.sort_values('ieg_z').iterrows():
            log('  %-12s %-18s %5d %-20s %+.2f'
                % (r['dataset'], r['group'][:18], r['n'], r['procurement'], r['ieg_z']))
        log('')

    log('=== does procurement or disease explain the module level? ===')
    d2 = df.copy()
    d2['is_biopsy'] = (d2['procurement'] == 'biopsy').astype(int)
    d2['is_dkd'] = d2['group'].str.contains('DN|DKD', case=False, regex=True).astype(int)
    for ds, sub in d2.groupby('dataset'):
        if sub['is_biopsy'].nunique() < 2:
            continue
        b = sub[sub['is_biopsy'] == 1]['ieg_z'].mean()
        nb = sub[sub['is_biopsy'] == 0]['ieg_z'].mean()
        log('  %-12s biopsy %+.2f   non-biopsy %+.2f   difference %+.2f'
            % (ds, b, nb, b - nb))
    log('')
    both = d2[d2.groupby('dataset')['is_biopsy'].transform('nunique') > 1]
    if len(both):
        r = stats.spearmanr(both['is_biopsy'], both['ieg_z'])
        log('  across all datasets with both procurement types:')
        log('    Spearman(is_biopsy, IEG z) = %+.3f  p = %.3f' % (r.statistic, r.pvalue))
    dkd_pairs = []
    for ds, sub in d2.groupby('dataset'):
        a = sub[(sub['is_dkd'] == 1)]['ieg_z']
        b = sub[(sub['is_dkd'] == 0) & (sub['is_biopsy'] == 1)]['ieg_z']
        if len(a) and len(b):
            dkd_pairs.append((ds, float(a.mean() - b.mean())))
    if dkd_pairs:
        log('\n  DKD minus other-disease BIOPSY (procurement held constant):')
        for ds, v in dkd_pairs:
            log('    %-12s %+.2f' % (ds, v))
        log('    -> if this is ~0 everywhere the module is procurement; if it is consistently')
        log('       negative the module also carries a disease component')
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
