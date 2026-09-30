#!/usr/bin/env python
"""Feasibility check for cross-cohort feature selection.

For every cohort, compute a scale-free per-gene effect size for DKD vs control
(Hedges' g, plus rank-biserial AUC). Then correlate effect-size vectors across cohorts.

If the harmonisation is sound and a reproducible DKD signal exists, the glomerular
cohorts must correlate positively. This is the precondition for the whole
cross-cohort stability-selection design - worth establishing before any modelling.
"""
import os, sys, itertools
import numpy as np
import pandas as pd
from scipy import stats

H = 'data/processed/harmonized'
OUT = 'data/processed'
COHORTS = ['GSE142025', 'GSE30528', 'GSE96804', 'GSE104948', 'GSE30529', 'GSE104954']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load(c):
    e = pd.read_csv(os.path.join(H, '%s_expr.tsv' % c), sep='\t', index_col=0)
    p = pd.read_csv(os.path.join(H, '%s_pheno.tsv' % c), sep='\t')
    p = p[p['label'].isin([0, 1])]
    e = e[p['sample'].astype(str).values]
    return e, p


def hedges_g(x, y):
    """x = cases, y = controls; arrays are genes x samples."""
    n1, n2 = x.shape[1], y.shape[1]
    m1, m2 = np.nanmean(x, 1), np.nanmean(y, 1)
    v1, v2 = np.nanvar(x, 1, ddof=1), np.nanvar(y, 1, ddof=1)
    sp = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    sp[sp == 0] = np.nan
    d = (m1 - m2) / sp
    J = 1 - 3 / (4 * (n1 + n2) - 9)          # small-sample bias correction
    return d * J


def auc_rank(x, y):
    """Mann-Whitney AUC per gene, cases vs controls."""
    n1, n2 = x.shape[1], y.shape[1]
    both = np.hstack([x, y])
    r = np.apply_along_axis(stats.rankdata, 1, both)
    r1 = r[:, :n1].sum(1)
    return (r1 - n1 * (n1 + 1) / 2) / (n1 * n2)


def main():
    eff, aucs, meta = {}, {}, []
    for c in COHORTS:
        e, p = load(c)
        case = e.loc[:, (p['label'] == 1).values].values
        ctrl = e.loc[:, (p['label'] == 0).values].values
        g = hedges_g(case, ctrl)
        a = auc_rank(case, ctrl)
        eff[c] = pd.Series(g, index=e.index)
        aucs[c] = pd.Series(a, index=e.index)
        meta.append((c, case.shape[1], ctrl.shape[1]))
        log('%-11s cases=%-3d controls=%-3d  |g|>0.8: %5d genes  AUC>0.8: %5d genes'
            % (c, case.shape[1], ctrl.shape[1], int(np.nansum(np.abs(g) > 0.8)),
               int(np.nansum(np.abs(a - 0.5) > 0.3))))

    E = pd.DataFrame(eff)
    A = pd.DataFrame(aucs)
    E.to_csv(os.path.join(OUT, 'effect_sizes_hedges_g.tsv'), sep='\t')
    A.to_csv(os.path.join(OUT, 'effect_sizes_auc.tsv'), sep='\t')

    log('\n=== Spearman correlation of per-gene effect size (Hedges g) between cohorts ===')
    names = list(E.columns)
    print('%-11s' % '', ' '.join('%10s' % n[:10] for n in names))
    for a in names:
        cells = []
        for b in names:
            if a == b:
                cells.append('%10s' % '-')
            else:
                m = E[[a, b]].dropna()
                rho = stats.spearmanr(m[a], m[b]).statistic
                cells.append('%10.3f' % rho)
        print('%-11s' % a, ' '.join(cells))

    glom = ['GSE142025', 'GSE30528', 'GSE96804', 'GSE104948']
    log('\n=== consistency of direction across the 4 discovery cohorts ===')
    sub = E[glom].dropna()
    signs = np.sign(sub.values)
    agree = (np.abs(signs.sum(1)) == 4)
    log('  genes with the SAME direction in all 4: %d / %d (%.1f%%)  [chance = 12.5%%]'
        % (agree.sum(), len(sub), 100 * agree.mean()))
    strong = sub[(np.abs(sub) > 0.5).all(axis=1) & agree]
    log('  ... and |g| > 0.5 in all 4: %d genes' % len(strong))
    if len(strong):
        sym = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str).set_index('entrez_id')['symbol']
        top = strong.assign(mean_g=strong.mean(axis=1)).sort_values('mean_g', key=lambda s: -s.abs())
        log('\n  top 25 by |mean g| (concordant in all 4 discovery cohorts):')
        for gid, row in top.head(25).iterrows():
            log('    %-10s %-12s mean_g=%+.2f   ' % (gid, sym.get(str(gid), '?'), row['mean_g'])
                + ' '.join('%s=%+.2f' % (c[3:], row[c]) for c in glom))
        top.assign(symbol=[sym.get(str(i), '') for i in top.index]).to_csv(
            os.path.join(OUT, 'concordant_genes_discovery4.tsv'), sep='\t')
        log('\n  wrote %s/concordant_genes_discovery4.tsv (%d genes)' % (OUT, len(top)))


if __name__ == '__main__':
    main()
