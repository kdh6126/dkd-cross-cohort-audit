#!/usr/bin/env python
"""When does stability-based feature selection prefer a confounder over the disease?

The real-data finding needs a mechanism, otherwise it is one observation about one disease.
The proposed mechanism:

    Stability selection rewards features that are chosen again and again across resamples and
    cohorts. What gets chosen consistently is whatever has the LOWEST cross-cohort
    heterogeneity - not whatever is most disease-relevant. A confounder introduced by a shared
    protocol (every study uses biopsies for cases and nephrectomies for controls) has almost no
    heterogeneity, while a true disease effect varies between cohorts through ancestry, stage,
    platform and compartment. So the confounder is *more* stable than the biology, and a
    stability criterion prefers it.

This predicts something specific and falsifiable: the failure should depend on the RATIO of
disease heterogeneity to confounder heterogeneity, not on effect size, and effect-size-based
selection (which ignores consistency) should be less affected.

Generative model, per cohort c and patient i:
    y_ci ~ Bernoulli(0.5)                                  case status
    C_ci = rho * y_ci + sqrt(1-rho^2) * N(0,1)             confounder, correlated with y by rho
    disease genes    x += (delta_d + N(0, sigma_d)) * y     effect varies across cohorts
    confounder genes x += (delta_c + N(0, sigma_c)) * C     effect nearly constant across cohorts
    both blocks carry within-block correlation, as real modules do

sigma_d and sigma_c are drawn ONCE PER COHORT, which is what makes them heterogeneity terms
rather than noise.
"""
import os, sys, argparse, itertools
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = 'results/simulation'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def make_cohort(rng, n, p, dis, conf, delta_d, delta_c, sigma_d, sigma_c, rho, block_r):
    """One cohort. dis/conf are index arrays for the two blocks."""
    y = rng.integers(0, 2, n)
    C = rho * y + np.sqrt(max(1 - rho ** 2, 0)) * rng.normal(0, 1, n)
    X = rng.normal(0, 1, (n, p))

    # within-block correlation: add a shared latent per block
    for blk in (dis, conf):
        lat = rng.normal(0, 1, n)
        X[:, blk] = np.sqrt(1 - block_r) * X[:, blk] + np.sqrt(block_r) * lat[:, None]

    d_eff = delta_d + rng.normal(0, sigma_d)        # cohort-specific disease effect
    c_eff = delta_c + rng.normal(0, sigma_c)        # cohort-specific confounder effect
    X[:, dis] += d_eff * y[:, None]
    X[:, conf] += c_eff * C[:, None]
    return X, y


def zscore(X):
    return (X - X.mean(0)) / (X.std(0, ddof=1) + 1e-12)


def sel_univariate(Xs, ys, k, rng):
    """Effect-size selection: pooled |t|."""
    X = np.vstack(Xs); y = np.concatenate(ys)
    x1, x0 = X[y == 1], X[y == 0]
    sp = np.sqrt(x1.var(0, ddof=1) / len(x1) + x0.var(0, ddof=1) / len(x0)) + 1e-12
    return np.abs(x1.mean(0) - x0.mean(0)) / sp


def sel_stability(Xs, ys, k, rng, B=60):
    """Stability selection: per-cohort bootstrap, geometric-mean aggregation across cohorts."""
    freqs = []
    for X, y in zip(Xs, ys):
        cnt = np.zeros(X.shape[1])
        for _ in range(B):
            idx = np.concatenate([rng.choice(np.flatnonzero(y == lab), (y == lab).sum(),
                                             replace=True) for lab in (0, 1)])
            Xb, yb = X[idx], y[idx]
            x1, x0 = Xb[yb == 1], Xb[yb == 0]
            sp = np.sqrt(x1.var(0, ddof=1) / max(len(x1), 1)
                         + x0.var(0, ddof=1) / max(len(x0), 1)) + 1e-12
            t = np.abs(x1.mean(0) - x0.mean(0)) / sp
            cnt[np.argsort(-t)[:k]] += 1
        freqs.append(cnt / B)
    F = np.vstack(freqs)
    return np.exp(np.log(F + 1e-3).mean(0))          # geometric mean, as in RBS


def sel_connectivity(Xs, ys, k, rng, power=6, module_frac=0.10):
    """WGCNA-style: disease-correlated module, rank by intramodular connectivity."""
    X = zscore(np.vstack(Xs)); y = np.concatenate(ys).astype(float)
    n = len(y)
    yz = (y - y.mean()) / (y.std(ddof=1) + 1e-12)
    gs = (X.T @ yz) / (n - 1)
    m = max(50, int(module_frac * X.shape[1]))
    mod = np.argsort(-np.abs(gs))[:m]
    S = X[:, mod]
    A = np.abs((S.T @ S) / (n - 1)) ** power
    np.fill_diagonal(A, 0.0)
    sc = np.zeros(X.shape[1])
    sc[mod] = A.sum(1) / (A.sum(1).max() + 1e-12)
    return sc


SELECTORS = {'effect_size': sel_univariate, 'stability': sel_stability,
             'connectivity': sel_connectivity}


def one_run(rng, n_cohorts, rho, sigma_d, sigma_c, p, n, n_dis, n_conf,
            delta_d, delta_c, block_r, K):
    dis = np.arange(n_dis)
    conf = np.arange(n_dis, n_dis + n_conf)
    Xs, ys = [], []
    for _ in range(n_cohorts):
        X, y = make_cohort(rng, n, p, dis, conf, delta_d, delta_c, sigma_d, sigma_c,
                           rho, block_r)
        Xs.append(zscore(X)); ys.append(y)

    out = {}
    for name, fn in SELECTORS.items():
        sc = fn(Xs, ys, K, rng)
        top = np.argsort(-sc)[:K]
        out[name] = dict(
            disease_precision=float(np.isin(top, dis).mean()),
            confounder_frac=float(np.isin(top, conf).mean()))
        # cross-fold agreement: leave one cohort out
        if n_cohorts >= 3:
            sets = []
            for h in range(n_cohorts):
                keep = [i for i in range(n_cohorts) if i != h]
                s2 = fn([Xs[i] for i in keep], [ys[i] for i in keep], K, rng)
                sets.append(set(np.argsort(-s2)[:K]))
            js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
            out[name]['cross_fold_jaccard'] = float(np.mean(js))
        else:
            out[name]['cross_fold_jaccard'] = np.nan
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reps', type=int, default=20)
    ap.add_argument('--p', type=int, default=2000)
    ap.add_argument('--n', type=int, default=40)
    ap.add_argument('--n-dis', type=int, default=30)
    ap.add_argument('--n-conf', type=int, default=30)
    ap.add_argument('--delta-d', type=float, default=0.8)
    ap.add_argument('--delta-c', type=float, default=0.8)
    ap.add_argument('--block-r', type=float, default=0.35)
    ap.add_argument('--K', type=int, default=50)
    ap.add_argument('--seed', type=int, default=20260825)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    rows = []

    # ---- sweep 1: heterogeneity ratio, the mechanism under test
    log('sweep 1: disease heterogeneity vs confounder heterogeneity (4 cohorts, rho=0.9)')
    for sd, sc_ in itertools.product([0.0, 0.2, 0.4, 0.6], [0.0, 0.2, 0.4, 0.6]):
        for r in range(args.reps):
            o = one_run(rng, 4, 0.9, sd, sc_, args.p, args.n, args.n_dis, args.n_conf,
                        args.delta_d, args.delta_c, args.block_r, args.K)
            for m, v in o.items():
                rows.append(dict(sweep='heterogeneity', method=m, sigma_d=sd, sigma_c=sc_,
                                 rho=0.9, n_cohorts=4, rep=r, **v))
    log('  done')

    # ---- sweep 2: confounding strength
    log('sweep 2: confounder-disease correlation rho (4 cohorts, sigma_d=0.4, sigma_c=0.1)')
    for rho in [0.0, 0.3, 0.6, 0.9, 0.99]:
        for r in range(args.reps):
            o = one_run(rng, 4, rho, 0.4, 0.1, args.p, args.n, args.n_dis, args.n_conf,
                        args.delta_d, args.delta_c, args.block_r, args.K)
            for m, v in o.items():
                rows.append(dict(sweep='rho', method=m, sigma_d=0.4, sigma_c=0.1,
                                 rho=rho, n_cohorts=4, rep=r, **v))
    log('  done')

    # ---- sweep 3: number of cohorts
    log('sweep 3: number of cohorts (rho=0.9, sigma_d=0.4, sigma_c=0.1)')
    for nc in [2, 3, 4, 6, 8]:
        for r in range(args.reps):
            o = one_run(rng, nc, 0.9, 0.4, 0.1, args.p, args.n, args.n_dis, args.n_conf,
                        args.delta_d, args.delta_c, args.block_r, args.K)
            for m, v in o.items():
                rows.append(dict(sweep='n_cohorts', method=m, sigma_d=0.4, sigma_c=0.1,
                                 rho=0.9, n_cohorts=nc, rep=r, **v))
    log('  done')

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'simulation.tsv'), sep='\t', index=False)

    # ---------------------------------------------------------------- report
    log('\n=== sweep 1: confounder genes in the top-%d, by heterogeneity ===' % args.K)
    log('  rows = disease heterogeneity sigma_d, cols = confounder heterogeneity sigma_c')
    for m in SELECTORS:
        s = df[(df['sweep'] == 'heterogeneity') & (df['method'] == m)]
        piv = s.pivot_table(index='sigma_d', columns='sigma_c', values='confounder_frac')
        log('\n  %s' % m)
        log('        ' + ' '.join('%7.1f' % c for c in piv.columns))
        for i, row in piv.iterrows():
            log('  %5.1f ' % i + ' '.join('%7.2f' % v for v in row.values))

    log('\n=== sweep 1: TRUE disease genes recovered (precision@%d) ===' % args.K)
    for m in SELECTORS:
        s = df[(df['sweep'] == 'heterogeneity') & (df['method'] == m)]
        piv = s.pivot_table(index='sigma_d', columns='sigma_c', values='disease_precision')
        log('\n  %s' % m)
        log('        ' + ' '.join('%7.1f' % c for c in piv.columns))
        for i, row in piv.iterrows():
            log('  %5.1f ' % i + ' '.join('%7.2f' % v for v in row.values))

    log('\n=== sweep 2: effect of confounding strength rho ===')
    log('  %-13s %6s %18s %18s' % ('method', 'rho', 'confounder in top-K', 'disease precision'))
    for m in SELECTORS:
        s = df[(df['sweep'] == 'rho') & (df['method'] == m)]
        for rho, g in s.groupby('rho'):
            log('  %-13s %6.2f %18.2f %18.2f'
                % (m, rho, g['confounder_frac'].mean(), g['disease_precision'].mean()))

    log('\n=== sweep 3: effect of cohort count ===')
    log('  %-13s %6s %18s %18s %s'
        % ('method', 'cohorts', 'confounder in top-K', 'disease precision', 'cross-fold J'))
    for m in SELECTORS:
        s = df[(df['sweep'] == 'n_cohorts') & (df['method'] == m)]
        for nc, g in s.groupby('n_cohorts'):
            log('  %-13s %6d %18.2f %18.2f %14.3f'
                % (m, nc, g['confounder_frac'].mean(), g['disease_precision'].mean(),
                   g['cross_fold_jaccard'].mean()))
    log('\nwrote %s/simulation.tsv' % OUT)


if __name__ == '__main__':
    main()
