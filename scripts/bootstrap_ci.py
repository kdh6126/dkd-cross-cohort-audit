#!/usr/bin/env python
"""Bootstrap confidence intervals on every method comparison.

Every number reported so far is a mean over 3-4 LODO folds with no measure of uncertainty.
"RBS 0.303 vs WGCNA 0.258" is six pairwise Jaccard values against six others; whether that
difference is real has not been established, and a methods reviewer will ask first.

Design:
  * The resampling unit is the PATIENT, stratified by (cohort, class). Resampling folds would
    give 3-4 units and resampling the pairwise Jaccards would treat six correlated numbers as
    independent - neither supports an interval.
  * The comparison is PAIRED: every method sees the same resampled patients in every replicate,
    so the interval is on the DIFFERENCE, which is what the claim is about and which has far
    less variance than two independent intervals.
  * The whole LODO evaluation is re-run inside each replicate. Selection, model fitting and
    scoring all move with the resample; nothing is held fixed that would narrow the interval
    artificially.

RBS runs with a reduced inner bootstrap (B_inner) to make this affordable, so the intervals
describe that configuration. Stated rather than hidden.
"""
import os, sys, argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM
from dkd_rbs import rbs_fit
from dkd_deconfound import handling_score, residualise, IEG
from conventional_pipeline import deg_meta, wgcna_hub

OUT = 'results/bootstrap_ci'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def resample(rng, data, cohorts):
    """Stratified patient bootstrap inside each (cohort, class)."""
    out = {}
    for c in cohorts:
        X, y = data[c][0].values, data[c][1]
        idx = []
        for lab in (0, 1):
            cell = np.flatnonzero(y == lab)
            if len(cell):
                idx.append(rng.choice(cell, size=len(cell), replace=True))
        idx = np.concatenate(idx)
        out[c] = (X[idx], y[idx])
    return out


def auc_of(Xtr, ytr, Xte, yte, cols):
    if len(np.unique(yte)) < 2 or len(np.unique(ytr)) < 2:
        return np.nan
    m = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000,
                           class_weight='balanced').fit(Xtr[:, cols], ytr)
    return roc_auc_score(yte, m.predict_proba(Xte[:, cols])[:, 1])


def run_all_methods(dat, cohorts, ieg_cols, K, B_inner, base, seed, power, return_orders=False):
    """One LODO evaluation of every method on one (possibly resampled) dataset."""
    orders = {m: {} for m in ('DEG_meta', 'WGCNA_hub', 'RBS', 'RBS_orth')}
    aucs = {m: [] for m in orders}
    for held in cohorts:
        tr = [c for c in cohorts if c != held]
        sub = {c: dat[c] for c in tr}
        Xtr = np.vstack([dat[c][0] for c in tr])
        ytr = np.concatenate([dat[c][1] for c in tr])
        Xte, yte = dat[held]

        sc = {}
        sc['DEG_meta'] = np.abs(deg_meta({c: (pd.DataFrame(sub[c][0]), sub[c][1]) for c in tr}, tr))
        sc['WGCNA_hub'] = wgcna_hub({c: (pd.DataFrame(sub[c][0]), sub[c][1]) for c in tr},
                                    tr, power=power)[0]
        r = rbs_fit(sub, k=K, B=B_inner, base=base, seed=seed, agg='geomean',
                    pert_B=8, pert_pool=200)
        sc['RBS'] = r['score']

        sub_o = {}
        for c in tr:
            Z, yc = sub[c][0].copy(), sub[c][1]
            sub_o[c] = (residualise(Z, handling_score(Z, ieg_cols)), yc)
        ro = rbs_fit(sub_o, k=K, B=B_inner, base=base, seed=seed, agg='geomean',
                     pert_B=8, pert_pool=200)
        sc['RBS_orth'] = ro['score']

        for m, s in sc.items():
            o = np.argsort(-np.nan_to_num(s, nan=-np.inf))
            orders[m][held] = o
            aucs[m].append(auc_of(Xtr, ytr, Xte, yte, o[:K]))

    res = {}
    for m in orders:
        sets = [set(orders[m][h][:K]) for h in cohorts]
        js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
        res[m] = dict(cross_fold_jaccard=float(np.mean(js)),
                      external_auroc=float(np.nanmean(aucs[m])))
    if return_orders:
        return res, orders
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-boot', type=int, default=40)
    ap.add_argument('--K', type=int, default=50)
    ap.add_argument('--B-inner', type=int, default=50)
    ap.add_argument('--base', default='relieff')
    ap.add_argument('--power', type=int, default=6)
    ap.add_argument('--cohorts', default=','.join(GLOM) + ',GSE142025')
    ap.add_argument('--seed', type=int, default=20260825)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])
    gpos = {g: i for i, g in enumerate(genes)}
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = {s: g for g, s in zip(gs['entrez_id'], gs['symbol']) if isinstance(s, str)}
    ieg_cols = [gpos[sym2id[s]] for s in IEG if s in sym2id and sym2id[s] in gpos]

    log('point estimate (no resampling) ...')
    plain = {c: (data[c][0].values, data[c][1]) for c in cohorts}
    point = run_all_methods(plain, cohorts, ieg_cols, args.K, args.B_inner,
                            args.base, args.seed, args.power)
    for m, v in point.items():
        log('  %-11s jaccard=%.3f  auroc=%.3f' % (m, v['cross_fold_jaccard'], v['external_auroc']))

    rng = np.random.default_rng(args.seed)
    rows = []
    for b in range(args.n_boot):
        dat = resample(rng, data, cohorts)
        try:
            r = run_all_methods(dat, cohorts, ieg_cols, args.K, args.B_inner,
                                args.base, args.seed + b, args.power)
        except Exception as e:
            log('  replicate %d failed: %s' % (b, str(e)[:80]))
            continue
        for m, v in r.items():
            rows.append(dict(rep=b, method=m, **v))
        if (b + 1) % 5 == 0:
            log('  %d / %d replicates' % (b + 1, args.n_boot))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'replicates.tsv'), sep='\t', index=False)

    methods = ['DEG_meta', 'WGCNA_hub', 'RBS', 'RBS_orth']
    log('\n=== percentile bootstrap intervals (%d replicates) ===' % df['rep'].nunique())
    log('  %-11s %-26s %s' % ('method', 'cross-fold Jaccard', 'external AUROC'))
    summ = []
    for m in methods:
        d = df[df['method'] == m]
        cj = d['cross_fold_jaccard'].dropna()
        ca = d['external_auroc'].dropna()
        log('  %-11s %.3f [%.3f, %.3f]      %.3f [%.3f, %.3f]'
            % (m, point[m]['cross_fold_jaccard'], np.percentile(cj, 2.5), np.percentile(cj, 97.5),
               point[m]['external_auroc'], np.percentile(ca, 2.5), np.percentile(ca, 97.5)))
        summ.append(dict(method=m, jaccard_point=point[m]['cross_fold_jaccard'],
                         jaccard_lo=np.percentile(cj, 2.5), jaccard_hi=np.percentile(cj, 97.5),
                         auroc_point=point[m]['external_auroc'],
                         auroc_lo=np.percentile(ca, 2.5), auroc_hi=np.percentile(ca, 97.5)))
    pd.DataFrame(summ).to_csv(os.path.join(OUT, 'intervals.tsv'), sep='\t', index=False)

    log('\n=== PAIRED differences vs the best conventional pipeline ===')
    log('  the question is whether the advantage survives resampling, so the interval is on')
    log('  the difference computed within each replicate, not on two separate intervals')
    piv_j = df.pivot(index='rep', columns='method', values='cross_fold_jaccard')
    piv_a = df.pivot(index='rep', columns='method', values='external_auroc')
    drows = []
    for a, b in [('RBS', 'WGCNA_hub'), ('RBS', 'DEG_meta'),
                 ('RBS_orth', 'WGCNA_hub'), ('RBS', 'RBS_orth')]:
        for name, piv in (('cross-fold Jaccard', piv_j), ('external AUROC', piv_a)):
            d = (piv[a] - piv[b]).dropna()
            lo, hi = np.percentile(d, 2.5), np.percentile(d, 97.5)
            p = 2 * min((d <= 0).mean(), (d >= 0).mean())
            sig = 'significant' if lo > 0 or hi < 0 else 'NOT significant'
            log('  %-10s - %-10s  %-19s  %+.3f [%+.3f, %+.3f]  p=%.3f  %s'
                % (a, b, name, d.mean(), lo, hi, p, sig))
            drows.append(dict(a=a, b=b, metric=name, diff=d.mean(), lo=lo, hi=hi,
                              p_two_sided=p, significant=(lo > 0 or hi < 0)))
    pd.DataFrame(drows).to_csv(os.path.join(OUT, 'paired_differences.tsv'), sep='\t', index=False)
    log('\nwrote %s/' % OUT)


if __name__ == '__main__':
    main()
