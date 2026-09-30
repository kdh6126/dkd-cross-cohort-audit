#!/usr/bin/env python
"""Ablation over the RBS design choices.

The question this answers: the cross-cohort aggregation rule is the heart of the method, and
a raw MIN turned out to be dominated by estimation noise from the smallest cohort. Rather
than quietly swapping it for whatever scores best, sweep the choices and report the sweep.

Dimensions
  cohorts   glom3 (three glomerular cohorts) vs all4 (+ GSE142025, whole cortex)
            - with 3 cohorts each LODO fold trains on only 2, so any across-cohort
              aggregation is computed from two numbers. all4 gives three.
  agg       min | geomean | q25 | mean      how per-cohort frequencies are combined
  base      which selector runs inside each cohort's bootstrap

Reported: cross-fold Jaccard (does the signature survive changing the held-out cohort?)
and mean external LODO AUROC. Both matter - a rule that maximises stability by selecting
the same useless genes every time is not an improvement.
"""
import os, sys, itertools, argparse, time
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM
from dkd_rbs import rbs_fit

CONFIGS = {'glom3': GLOM, 'all4': GLOM + ['GSE142025']}
KS = (20, 50)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def ext_auc(Xtr, ytr, Xte, yte, cols):
    m = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000, class_weight='balanced')
    m.fit(Xtr[:, cols], ytr)
    if len(np.unique(yte)) < 2:
        return np.nan
    return roc_auc_score(yte, m.predict_proba(Xte[:, cols])[:, 1])


def run(cfg_name, agg, base, B, k, seed):
    cohorts = CONFIGS[cfg_name]
    data, genes = load_many(cohorts)
    orders, aucs = {}, {K: [] for K in KS}
    for held in cohorts:
        tr = [c for c in cohorts if c != held]
        train = {c: (data[c][0].values, data[c][1]) for c in tr}
        Xtr = np.vstack([train[c][0] for c in tr])
        ytr = np.concatenate([train[c][1] for c in tr])
        Xte, yte = data[held][0].values, data[held][1]
        res = rbs_fit(train, k=k, B=B, base=base, seed=seed, agg=agg,
                      pert_B=16, pert_pool=200)
        order = np.argsort(-res['score'])
        orders[held] = order
        for K in KS:
            aucs[K].append(ext_auc(Xtr, ytr, Xte, yte, order[:K]))
    out = dict(cohorts=cfg_name, agg=agg, base=base)
    for K in KS:
        sets = [set(orders[h][:K]) for h in cohorts]
        js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
        out['jaccard_K%d' % K] = float(np.mean(js))
        out['auroc_K%d' % K] = float(np.nanmean(aucs[K]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--B', type=int, default=100)
    ap.add_argument('--k', type=int, default=50)
    ap.add_argument('--seed', type=int, default=20260824)
    ap.add_argument('--out', default='results/ablation')
    ap.add_argument('--bases', default='lasso,univariate,relieff')
    ap.add_argument('--aggs', default='min,geomean,q25,mean')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    rows = []
    combos = list(itertools.product(CONFIGS, args.aggs.split(','), args.bases.split(',')))
    log('%d combinations, B=%d' % (len(combos), args.B))
    for cfg, agg, base in combos:
        t = time.time()
        r = run(cfg, agg, base, args.B, args.k, args.seed)
        r['seconds'] = round(time.time() - t, 1)
        rows.append(r)
        log('  %-6s %-8s %-11s  Jacc@20=%.3f Jacc@50=%.3f | AUROC@20=%.3f AUROC@50=%.3f  (%.0fs)'
            % (cfg, agg, base, r['jaccard_K20'], r['jaccard_K50'],
               r['auroc_K20'], r['auroc_K50'], r['seconds']))
        pd.DataFrame(rows).to_csv(os.path.join(args.out, 'ablation.tsv'), sep='\t', index=False)

    df = pd.DataFrame(rows)
    log('\n=== ranked by cross-fold Jaccard @ K=50 ===')
    for _, r in df.sort_values('jaccard_K50', ascending=False).head(12).iterrows():
        log('  %-6s %-8s %-11s Jacc=%.3f  AUROC=%.3f'
            % (r['cohorts'], r['agg'], r['base'], r['jaccard_K50'], r['auroc_K50']))
    log('\nwrote %s/ablation.tsv' % args.out)


if __name__ == '__main__':
    main()
