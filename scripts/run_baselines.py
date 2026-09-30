#!/usr/bin/env python
"""Baseline benchmark: LODO x {feature selectors} x bootstrap stability selection.

Produces the table the paper needs:
    method | #genes | internal CV AUROC | external (LODO) AUROC | within-cohort stability
plus the per-gene selection-frequency matrices that the cross-cohort stability step consumes.

Usage:
    python scripts/run_baselines.py [--B 200] [--k 50] [--methods a,b,c] [--out DIR]
"""
import os, sys, json, time, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, lodo_folds, GLOM
from dkd_selectors import SELECTORS
from dkd_stability import (run_fold, fit_and_score, nested_cv_auc,
                           mean_pairwise_jaccard, kuncheva)

TOPK = [5, 10, 20, 50, 100]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--B', type=int, default=200)
    ap.add_argument('--k', type=int, default=50)
    ap.add_argument('--methods', default=','.join(SELECTORS))
    ap.add_argument('--cohorts', default=','.join(GLOM))
    ap.add_argument('--out', default='results/baselines')
    ap.add_argument('--seed', type=int, default=20260824)
    args = ap.parse_args()

    methods = [m.strip() for m in args.methods.split(',') if m.strip()]
    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    os.makedirs(args.out, exist_ok=True)

    data, genes = load_many(cohorts)
    genes = np.array(genes)
    sym = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t',
                      dtype=str).set_index('entrez_id')['symbol'].to_dict()
    log('cohorts=%s  genes=%d  B=%d  k=%d' % (cohorts, len(genes), args.B, args.k))

    rows, freq_tables = [], {}
    for held, Xtr_df, ytr, ctr, Xte_df, yte in lodo_folds(data, cohorts):
        Xtr, Xte = Xtr_df.values, Xte_df.values
        log('\n=== hold out %s (train %d, test %d) ===' % (held, len(ytr), len(yte)))
        for mi, method in enumerate(methods):
            seed = args.seed + 1000 * mi + hash(held) % 997
            res = run_fold(method, Xtr, ytr, ctr, Xte, yte, genes,
                           B=args.B, k=args.k, seed=seed)
            freq, meanrank, sets = res['freq'], res['meanrank'], res['sets']
            freq_tables[(held, method)] = freq

            order = np.lexsort((meanrank, -freq))       # freq desc, then mean rank asc
            stab_j = mean_pairwise_jaccard(sets, rng=np.random.default_rng(seed))
            stab_k = kuncheva(sets, p=len(genes), rng=np.random.default_rng(seed))

            for K in TOPK:
                cols = order[:K]
                ext = fit_and_score(Xtr, ytr, Xte, yte, cols)
                icv = nested_cv_auc(Xtr, ytr, method, K=K, k=args.k, seed=seed)
                rows.append(dict(held_out=held, method=method, K=K,
                                 external_auroc=ext, nested_cv_auroc=icv,
                                 jaccard=stab_j, kuncheva=stab_k,
                                 seconds=res['seconds'], B=len(sets)))
            log('  %-12s %5.1fs  jaccard=%.3f kuncheva=%.3f  ext AUROC@[%s] = %s'
                % (method, res['seconds'], stab_j, stab_k,
                   ','.join(str(x) for x in TOPK),
                   ' '.join('%.3f' % r['external_auroc']
                            for r in rows[-len(TOPK):])))
            pd.DataFrame(rows).to_csv(os.path.join(args.out, 'results.tsv'),
                                      sep='\t', index=False)

    # ---- persist per-gene selection frequencies (input to the cross-cohort stability step)
    for method in methods:
        cols = {held: freq_tables[(held, method)] for held in cohorts if (held, method) in freq_tables}
        if not cols:
            continue
        df = pd.DataFrame(cols, index=genes)
        df.index.name = 'entrez_id'
        df.insert(0, 'symbol', [sym.get(g, '') for g in genes])
        df['mean_freq'] = df[list(cols)].mean(axis=1)
        df['min_freq'] = df[list(cols)].min(axis=1)
        df.sort_values('mean_freq', ascending=False).to_csv(
            os.path.join(args.out, 'freq_%s.tsv' % method), sep='\t')

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(args.out, 'results.tsv'), sep='\t', index=False)

    log('\n=== mean over LODO folds ===')
    piv = (res.groupby(['method', 'K'])
              .agg(external=('external_auroc', 'mean'),
                   internal=('nested_cv_auroc', 'mean'),
                   jaccard=('jaccard', 'mean'),
                   kuncheva=('kuncheva', 'mean'))
              .reset_index())
    piv.to_csv(os.path.join(args.out, 'summary.tsv'), sep='\t', index=False)
    for K in TOPK:
        log('\n  K=%d' % K)
        sub = piv[piv['K'] == K].sort_values('external', ascending=False)
        for _, r in sub.iterrows():
            log('    %-12s external=%.3f  internal=%.3f  jaccard=%.3f  kuncheva=%.3f'
                % (r['method'], r['external'], r['internal'], r['jaccard'], r['kuncheva']))
    log('\nwrote %s/{results,summary}.tsv and freq_<method>.tsv' % args.out)


if __name__ == '__main__':
    main()
