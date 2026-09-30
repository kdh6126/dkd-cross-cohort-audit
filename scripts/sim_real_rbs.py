#!/usr/bin/env python
"""Repeat the block-size experiment with the ACTUAL RBS selector, not a stand-in.

The simulation so far used a simplified stability selector written inside simulation.py:
per-cohort bootstrap of a t-statistic, geometric-mean aggregation. That captures the idea, but
a reviewer is entitled to ask whether the conclusion holds for the selector the paper actually
proposes, which additionally has an importance term, a perturbation term, a sign gate and
ReliefF as its base learner.

So here the same generative model feeds rbs_fit() directly. If the size-bias result survives,
the earlier sweep stands. If it does not, the earlier sweep was an artefact of the stand-in.
"""
import os, sys, argparse
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simulation import make_cohort, zscore, sel_connectivity, sel_univariate
from dkd_rbs import rbs_fit

OUT = 'results/module_size'
CONFIGS = [(30, 30), (100, 30), (200, 20), (100, 100), (20, 200)]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def one_run_real(rng, n_cohorts, rho, sigma_d, sigma_c, p, n, n_dis, n_conf,
                 delta_d, delta_c, block_r, K, B, base):
    dis = np.arange(n_dis)
    conf = np.arange(n_dis, n_dis + n_conf)
    Xs, ys = [], []
    for _ in range(n_cohorts):
        X, y = make_cohort(rng, n, p, dis, conf, delta_d, delta_c, sigma_d, sigma_c,
                           rho, block_r)
        Xs.append(zscore(X)); ys.append(y)

    out = {}
    cohorts = {'c%d' % i: (Xs[i], ys[i]) for i in range(n_cohorts)}
    r = rbs_fit(cohorts, k=K, B=B, base=base, seed=int(rng.integers(1 << 30)),
                agg='geomean', pert_B=8, pert_pool=200)
    for name, sc in (('RBS_real', r['score']),
                     ('connectivity', sel_connectivity(Xs, ys, K, rng)),
                     ('effect_size', sel_univariate(Xs, ys, K, rng))):
        top = np.argsort(-np.nan_to_num(sc, nan=-np.inf))[:K]
        out[name] = dict(confounder_frac=float(np.isin(top, conf).mean()),
                         disease_precision=float(np.isin(top, dis).mean()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reps', type=int, default=25)
    ap.add_argument('--B', type=int, default=40)
    ap.add_argument('--base', default='relieff')
    ap.add_argument('--p', type=int, default=1500)
    ap.add_argument('--n', type=int, default=40)
    ap.add_argument('--K', type=int, default=50)
    ap.add_argument('--seed', type=int, default=20260825)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    rows = []
    for n_dis, n_conf in CONFIGS:
        log('  n_disease=%d n_confounder=%d ...' % (n_dis, n_conf))
        for r in range(args.reps):
            o = one_run_real(rng, 4, 0.95, 0.4, 0.05, args.p, args.n, n_dis, n_conf,
                             1.0, 1.0, 0.35, args.K, args.B, args.base)
            for m, v in o.items():
                rows.append(dict(n_dis=n_dis, n_conf=n_conf, method=m, rep=r, **v))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'real_rbs_sweep.tsv'), sep='\t', index=False)

    log('\n=== confounder fraction in the top-%d, real RBS vs the others ===' % args.K)
    log('  %-9s %-8s | %s' % ('n_disease', 'n_conf',
                              ' '.join('%-22s' % m for m in
                                       ('effect_size', 'RBS_real', 'connectivity'))))
    for n_dis, n_conf in CONFIGS:
        cells = []
        for m in ('effect_size', 'RBS_real', 'connectivity'):
            v = df[(df['n_dis'] == n_dis) & (df['n_conf'] == n_conf)
                   & (df['method'] == m)]['confounder_frac']
            se = v.std(ddof=1) / np.sqrt(len(v))
            cells.append('%.2f +-%.2f       ' % (v.mean(), 1.96 * se))
        log('  %-9d %-8d | %s' % (n_dis, n_conf, ' '.join(cells)))

    log('\n=== paired: connectivity minus RBS_real, CI of the mean ===')
    log('  %-9s %-8s %-30s %-10s %s' % ('n_disease', 'n_conf', 'difference [95% CI]',
                                        'Wilcoxon', 'verdict'))
    for n_dis, n_conf in CONFIGS:
        sub = df[(df['n_dis'] == n_dis) & (df['n_conf'] == n_conf)]
        piv = sub.pivot(index='rep', columns='method', values='confounder_frac')
        x = (piv['connectivity'] - piv['RBS_real']).dropna().values
        m_, se = x.mean(), x.std(ddof=1) / np.sqrt(len(x))
        lo, hi = m_ - 1.96 * se, m_ + 1.96 * se
        try:
            pv = stats.wilcoxon(x).pvalue
        except Exception:
            pv = np.nan
        v = ('connectivity WORSE' if lo > 0 else
             'connectivity BETTER' if hi < 0 else 'no difference')
        log('  %-9d %-8d %+.3f [%+.3f, %+.3f]          p=%-8.4f %s'
            % (n_dis, n_conf, m_, lo, hi, pv, v))
    log('\n  compare the sign pattern with results/module_size/block_size_paired.tsv;')
    log('  if it matches, the stand-in selector was an adequate proxy.')
    log('\nwrote %s/real_rbs_sweep.tsv' % OUT)


if __name__ == '__main__':
    main()
