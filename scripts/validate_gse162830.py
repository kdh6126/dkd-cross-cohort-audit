#!/usr/bin/env python
"""GSE162830: the cleanest confounder-matched contrast available anywhere in public data.

32 laser-microdissected samples:
    18  diabetic nephropathy          BIOPSY
     5  idiopathic nodular glomerulosclerosis (ING)   BIOPSY
     9  reference                      NEPHRECTOMY

ING is nodular mesangial sclerosis *without* diabetes - histologically the closest possible
phenocopy of diabetic nodular glomerulosclerosis. So this one dataset contains both contrasts:

    DN vs reference   biopsy vs nephrectomy   -> procurement is confounded with disease
    DN vs ING         biopsy vs biopsy        -> procurement matched, AND histology matched

That second contrast is stronger than the ERCB gate this project has relied on, where the
comparator was a mixture of glomerulonephritides with different histology. If the immediate-
early module behaves as claimed it should be large in the first contrast and absent in the
second, exactly as it was in KPMP. And a genuinely DKD-specific candidate should survive the
second contrast.

The ING group is only 5 samples, so effect sizes there are imprecise; the module-level test
(19 genes against a random background) is what carries weight, not any single gene.
"""
import os, sys, gzip, argparse
import numpy as np
import pandas as pd
from scipy import stats

EXPR = 'data/raw/geo/GSE162830/GSE162830_ING_quantile_normalized_final.csv.gz'
IEG = ['FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'ATF3', 'DUSP1', 'ZFP36', 'NR4A2', 'BTG2',
       'EGR2', 'EGR3', 'IER2', 'KLF2', 'KLF4', 'SOCS3', 'DUSP2', 'JUND', 'NR4A1']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    n1, n0 = len(a), len(b)
    if n1 < 2 or n0 < 2:
        return np.full(a.shape[1], np.nan)
    sp = np.sqrt(((n1 - 1) * a.var(0, ddof=1) + (n0 - 1) * b.var(0, ddof=1)) / (n1 + n0 - 2))
    sp[sp == 0] = np.nan
    J = 1 - 3 / (4 * (n1 + n0) - 9)
    return (a.mean(0) - b.mean(0)) / sp * J


def sym_maps():
    sym, syn = {}, {}
    with gzip.open('data/raw/annotation/Homo_sapiens.gene_info.gz', 'rt',
                   encoding='utf-8', errors='replace') as f:
        hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip('\n').split('\t')
            sym[p[ix['Symbol']]] = p[ix['GeneID']]
            for a in p[ix['Synonyms']].split('|'):
                if a and a != '-' and a not in syn:
                    syn[a] = p[ix['GeneID']]
    return sym, syn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', default='results/candidates_v2/candidates_gated.tsv')
    ap.add_argument('--out', default='results/gse162830_validation.tsv')
    args = ap.parse_args()

    d = pd.read_csv(EXPR, index_col=0)
    log('raw matrix: %d genes x %d samples' % d.shape)
    groups = {c: ('DN' if c.startswith('DIA') else 'ING' if c.startswith('ING') else 'REF')
              for c in d.columns}
    log('  groups: %s' % pd.Series(list(groups.values())).value_counts().to_dict())

    # counts -> CPM, log2, then z-score within the cohort as everywhere else
    cpm = d / d.sum(axis=0) * 1e6
    lg = np.log2(cpm + 1)
    lg = lg.loc[lg.var(axis=1) > 0]

    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    space = set(gs['entrez_id'])
    symname = dict(zip(gs['entrez_id'], gs['symbol']))
    sym, syn = sym_maps()
    ent = pd.Series([sym.get(i) or syn.get(i) for i in lg.index], index=lg.index)
    lg = lg[ent.notna().values]
    ent = ent[ent.notna()]
    lg = lg.assign(__e=ent.values)
    means = lg.drop(columns='__e').mean(axis=1)
    order = np.lexsort((-means.values, lg['__e'].values))
    lg = lg.iloc[order]
    lg = lg[~lg['__e'].duplicated()]
    X = lg.drop(columns='__e')
    X.index = lg['__e'].values
    X = X.loc[[g for g in X.index if g in space]]
    log('  mapped into the 9,900-gene space: %d genes' % len(X))

    Z = ((X.T - X.T.mean()) / X.T.std(ddof=1).replace(0, np.nan)).fillna(0.0)
    grp = np.array([groups[c] for c in Z.index])
    dn = Z[grp == 'DN'].values
    ref = Z[grp == 'REF'].values
    ing = Z[grp == 'ING'].values

    sym2id = {s: g for g, s in symname.items() if isinstance(s, str)}
    ieg_ids = [sym2id[s] for s in IEG if s in sym2id and sym2id[s] in Z.columns]
    log('  immediate-early genes present: %d' % len(ieg_ids))

    rows = []
    log('\n=== the two contrasts inside one dataset ===')
    for name, A, B, note in [('DN vs REF', dn, ref, 'biopsy vs nephrectomy - CONFOUNDED'),
                             ('DN vs ING', dn, ing, 'biopsy vs biopsy - MATCHED'),
                             ('ING vs REF', ing, ref, 'biopsy vs nephrectomy - control check')]:
        g = pd.Series(hedges_g(A, B), index=Z.columns)
        ieg = g[ieg_ids].dropna()
        rest = g.drop(index=ieg_ids).dropna()
        p = stats.mannwhitneyu(ieg, rest, alternative='two-sided').pvalue
        log('  %-11s (%2d vs %2d)  %-38s IEG mean g = %+.2f   rest %+.2f   p = %.2e'
            % (name, len(A), len(B), note, ieg.mean(), rest.mean(), p))
        rows.append(dict(contrast=name, n_a=len(A), n_b=len(B), ieg_mean_g=ieg.mean(),
                         rest_mean_g=rest.mean(), mannwhitney_p=p))
        for s in ['DUSP1', 'ZFP36', 'FOS', 'EGR1', 'JUN', 'ATF3']:
            gid = sym2id.get(s)
            if gid in g.index:
                rows.append(dict(contrast=name, gene=s, g=float(g[gid])))

    log('\n  read the first two rows together: if the module tracked diabetes it would separate')
    log('  DN from ING as well. If it tracks procurement it separates DN from REF only.')

    # ---------------------------------------------------------------- candidates
    cand = pd.read_csv(args.candidates, sep='\t')
    cands = [g for g in cand['symbol'].dropna().astype(str)]
    g_ref = pd.Series(hedges_g(dn, ref), index=Z.columns)
    g_ing = pd.Series(hedges_g(dn, ing), index=Z.columns)
    log('\n=== candidates: does the DKD signal survive the histology-matched contrast? ===')
    log('  %-10s %12s %12s' % ('gene', 'DN vs REF', 'DN vs ING'))
    keep = []
    for s in cands:
        gid = sym2id.get(s)
        if gid not in Z.columns:
            continue
        a, b = float(g_ref.get(gid, np.nan)), float(g_ing.get(gid, np.nan))
        keep.append((s, a, b))
        rows.append(dict(contrast='candidate', gene=s, g_dn_vs_ref=a, g_dn_vs_ing=b))
    for s, a, b in sorted(keep, key=lambda r: -abs(r[2]))[:18]:
        mark = ' *' if abs(b) > 0.8 and np.sign(a) == np.sign(b) else ''
        log('  %-10s %+12.2f %+12.2f%s' % (s, a, b, mark))
    arr = np.array([[a, b] for _, a, b in keep], float)
    ok = np.isfinite(arr).all(1)
    log('\n  candidates with |g| > 0.8 and consistent sign in the matched contrast: %d / %d'
        % (int(((np.abs(arr[ok, 1]) > 0.8) & (np.sign(arr[ok, 0]) == np.sign(arr[ok, 1]))).sum()),
           int(ok.sum())))
    rnd = np.random.default_rng(0).choice(Z.columns, size=min(500, Z.shape[1]), replace=False)
    ga, gb = g_ref[rnd].values, g_ing[rnd].values
    m = np.isfinite(ga) & np.isfinite(gb)
    log('  same criterion on 500 random genes: %.0f%%'
        % (100 * ((np.abs(gb[m]) > 0.8) & (np.sign(ga[m]) == np.sign(gb[m]))).mean()))

    pd.DataFrame(rows).to_csv(args.out, sep='\t', index=False)
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
