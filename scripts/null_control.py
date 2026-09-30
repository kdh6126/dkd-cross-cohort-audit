#!/usr/bin/env python
"""Null control: how well does a RANDOM gene set do on each LODO fold?

This control is not a formality here. The DKD signal is broadly co-expressed (mean |r| about
0.48 among the top relevance genes), so a randomly chosen gene set already carries much of
it - the phenomenon Venet et al. described for breast-cancer signatures. Without this control
an AUROC of 0.94 reads as a triumph when it may sit at the 95th percentile of chance.

Two nulls per fold and per signature size K:
  RANDOM genes     real labels, random gene subset  -> what any gene set achieves
  PERMUTED labels  shuffled training labels         -> confirms the labels are what matter

Raw draws are written next to the summary so downstream comparisons can express each
method's AUROC as a percentile / z-score against the matching null.
"""
import os, sys, argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def auc(Xtr, ytr, Xte, yte, cols):
    m = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000, class_weight='balanced')
    m.fit(Xtr[:, cols], ytr)
    if len(np.unique(yte)) < 2:
        return np.nan
    return roc_auc_score(yte, m.predict_proba(Xte[:, cols])[:, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cohorts', default=','.join(GLOM) + ',GSE142025')
    ap.add_argument('--Ks', default='5,10,20,50,100')
    ap.add_argument('--n-rep', type=int, default=300)
    ap.add_argument('--seed', type=int, default=20260824)
    ap.add_argument('--out', default='results/null_control.tsv')
    args = ap.parse_args()

    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    Ks = [int(x) for x in args.Ks.split(',')]
    data, genes = load_many(cohorts)
    p = len(genes)
    rng = np.random.default_rng(args.seed)
    rows, draws = [], {}

    log('Ks=%s, %d random gene sets per fold per K, %d genes\n' % (Ks, args.n_rep, p))
    log('%-11s %-4s %-25s %-25s %s'
        % ('held out', 'K', 'RANDOM genes AUROC', 'PERMUTED labels AUROC', 'n test'))
    for held in cohorts:
        tr = [c for c in cohorts if c != held]
        Xtr = np.vstack([data[c][0].values for c in tr])
        ytr = np.concatenate([data[c][1] for c in tr])
        Xte, yte = data[held][0].values, data[held][1]
        for K in Ks:
            rand = np.array([auc(Xtr, ytr, Xte, yte, rng.choice(p, K, replace=False))
                             for _ in range(args.n_rep)], float)
            perm = np.array([auc(Xtr, rng.permutation(ytr), Xte, yte,
                                 rng.choice(p, K, replace=False))
                             for _ in range(max(args.n_rep // 5, 20))], float)
            log('%-11s %-4d %.3f [%.3f-%.3f]        %.3f [%.3f-%.3f]         %d (%d case)'
                % (held, K,
                   np.nanmean(rand), np.nanpercentile(rand, 5), np.nanpercentile(rand, 95),
                   np.nanmean(perm), np.nanpercentile(perm, 5), np.nanpercentile(perm, 95),
                   len(yte), int(yte.sum())))
            rows.append(dict(held_out=held, K=K,
                             random_mean=np.nanmean(rand), random_sd=np.nanstd(rand, ddof=1),
                             random_p05=np.nanpercentile(rand, 5),
                             random_p95=np.nanpercentile(rand, 95),
                             random_max=np.nanmax(rand),
                             permuted_mean=np.nanmean(perm),
                             n_test=len(yte), n_case=int(yte.sum())))
            draws['%s|%d' % (held, K)] = rand

    pd.DataFrame(rows).to_csv(args.out, sep='\t', index=False)
    pd.DataFrame(draws).to_csv(args.out.replace('.tsv', '_draws.tsv'), sep='\t', index=False)
    log('\nwrote %s and %s' % (args.out, args.out.replace('.tsv', '_draws.tsv')))
    log('\nRead this as: on a fold where random genes already reach ~0.9, raw AUROC cannot\n'
        'separate a good signature from a bad one. Report percentile-vs-null, not raw AUROC.')


if __name__ == '__main__':
    main()
