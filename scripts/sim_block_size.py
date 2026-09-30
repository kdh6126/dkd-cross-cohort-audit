#!/usr/bin/env python
"""Definitive block-size experiment, with intervals.

Claim under test: intramodular connectivity (the WGCNA hub criterion) is biased by module
size, while stability- and effect-size-based selection are not. A confounder that forms a
SMALL coherent block is therefore filtered by connectivity and not by the others; a confounder
that forms a LARGE block reverses that, and connectivity becomes the worst of the three.

In the real data the immediate-early block has ~19 co-expressed members inside a 990-gene
disease module whose top-connected genes sit in blocks of ~705. That is the small-confounder
regime, and it is why WGCNA looked robust there. Nothing guarantees that regime in general.

Reported with percentile intervals over replicates, because the earlier 8-replicate version is
not enough to support a general claim.
"""
import os, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simulation import one_run

OUT = 'results/module_size'
CONFIGS = [(30, 30), (60, 30), (100, 30), (200, 20), (300, 20), (100, 100), (20, 200)]
METHODS = ('effect_size', 'stability', 'connectivity')


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reps', type=int, default=40)
    ap.add_argument('--p', type=int, default=2000)
    ap.add_argument('--n', type=int, default=40)
    ap.add_argument('--rho', type=float, default=0.95)
    ap.add_argument('--sigma-d', type=float, default=0.4)
    ap.add_argument('--sigma-c', type=float, default=0.05)
    ap.add_argument('--K', type=int, default=50)
    ap.add_argument('--seed', type=int, default=20260825)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    rows = []
    for n_dis, n_conf in CONFIGS:
        log('  n_disease=%d  n_confounder=%d ...' % (n_dis, n_conf))
        for r in range(args.reps):
            o = one_run(rng, 4, args.rho, args.sigma_d, args.sigma_c, args.p, args.n,
                        n_dis, n_conf, 1.0, 1.0, 0.35, args.K)
            for m, v in o.items():
                rows.append(dict(n_dis=n_dis, n_conf=n_conf, ratio=n_dis / n_conf,
                                 method=m, rep=r, **v))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'block_size_ci.tsv'), sep='\t', index=False)

    log('\n=== confounder genes in the top-%d  [95%% percentile interval over %d reps] ==='
        % (args.K, args.reps))
    log('  %-9s %-9s | %s' % ('n_disease', 'n_conf',
                              ' '.join('%-26s' % m for m in METHODS)))
    for n_dis, n_conf in CONFIGS:
        cells = []
        for m in METHODS:
            v = df[(df['n_dis'] == n_dis) & (df['n_conf'] == n_conf)
                   & (df['method'] == m)]['confounder_frac']
            cells.append('%.2f [%.2f, %.2f]     '
                         % (v.mean(), np.percentile(v, 2.5), np.percentile(v, 97.5)))
        log('  %-9d %-9d | %s' % (n_dis, n_conf, ' '.join(cells)))

    log('\n=== paired difference: connectivity minus stability, same replicate ===')
    log('  positive = connectivity picks MORE confounder genes than stability')
    log('  %-9s %-9s %-28s %s' % ('n_disease', 'n_conf', 'difference [95% CI]', 'verdict'))
    drows = []
    for n_dis, n_conf in CONFIGS:
        sub = df[(df['n_dis'] == n_dis) & (df['n_conf'] == n_conf)]
        piv = sub.pivot(index='rep', columns='method', values='confounder_frac')
        d = (piv['connectivity'] - piv['stability']).dropna()
        lo, hi = np.percentile(d, 2.5), np.percentile(d, 97.5)
        verdict = ('connectivity WORSE' if lo > 0 else
                   'connectivity BETTER' if hi < 0 else 'no difference')
        log('  %-9d %-9d %+.3f [%+.3f, %+.3f]        %s'
            % (n_dis, n_conf, d.mean(), lo, hi, verdict))
        drows.append(dict(n_dis=n_dis, n_conf=n_conf, diff=d.mean(), lo=lo, hi=hi,
                          verdict=verdict))
    pd.DataFrame(drows).to_csv(os.path.join(OUT, 'block_size_paired.tsv'), sep='\t', index=False)

    log('\n  the sign of the difference should flip with the block-size ratio if size bias is')
    log('  the mechanism; a constant sign would mean something else is going on.')
    log('\nwrote %s/block_size_ci.tsv and block_size_paired.tsv' % OUT)


if __name__ == '__main__':
    main()
