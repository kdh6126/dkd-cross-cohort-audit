#!/usr/bin/env python
"""Run the proposed Robust Biomarker Score across LODO folds and score it like a baseline.

Emits the same columns as run_baselines.py so the two can be concatenated into one table,
plus the per-gene component breakdown (I / S / R / P) for every fold.
"""
import os, sys, time, argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM
from dkd_rbs import rbs_fit
from dkd_stability import mean_pairwise_jaccard, kuncheva

TOPK = [5, 10, 20, 50, 100]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def fit_score(Xtr, ytr, Xte, yte, cols):
    m = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000, class_weight='balanced')
    m.fit(Xtr[:, cols], ytr)
    if len(np.unique(yte)) < 2:
        return np.nan
    return roc_auc_score(yte, m.predict_proba(Xte[:, cols])[:, 1])


def nested_cv(train_cohorts, k, B, base, seed, folds=5, agg='geomean'):
    """Nested internal CV: RBS is re-fit inside each CV split. The ranking is computed once
    per split and then read off at every K, which is why this stays affordable."""
    names = list(train_cohorts)
    Xall = np.vstack([train_cohorts[c][0] for c in names])
    yall = np.concatenate([train_cohorts[c][1] for c in names])
    call = np.concatenate([[c] * len(train_cohorts[c][1]) for c in names])

    n_min = min((yall == 0).sum(), (yall == 1).sum())
    folds = max(2, min(folds, n_min))
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    out = {K: [] for K in TOPK}
    for tr, te in skf.split(Xall, yall):
        if len(np.unique(yall[te])) < 2:
            continue
        sub = {}
        for c in names:
            m = (call[tr] == c)
            if m.sum() >= 4 and len(np.unique(yall[tr][m])) == 2:
                sub[c] = (Xall[tr][m], yall[tr][m])
        if len(sub) < 2:
            continue
        res = rbs_fit(sub, k=k, B=B, base=base, seed=seed, agg=agg, pert_B=8, pert_pool=200)
        order = np.argsort(-res['score'])
        for K in TOPK:
            out[K].append(fit_score(Xall[tr], yall[tr], Xall[te], yall[te], order[:K]))
    return {K: (float(np.nanmean(v)) if v else np.nan) for K, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--B', type=int, default=200)
    ap.add_argument('--k', type=int, default=50)
    ap.add_argument('--base', default='lasso')
    ap.add_argument('--cohorts', default=','.join(GLOM))
    ap.add_argument('--out', default='results/proposed')
    ap.add_argument('--seed', type=int, default=20260824)
    ap.add_argument('--inner-B', type=int, default=40)
    ap.add_argument('--tag', default='rbs')
    ap.add_argument('--agg', default='geomean', choices=['min','geomean','q25','mean'])
    args = ap.parse_args()

    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    os.makedirs(args.out, exist_ok=True)
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])
    sym = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t',
                      dtype=str).set_index('entrez_id')['symbol'].to_dict()
    log('cohorts=%s genes=%d B=%d k=%d base=%s agg=%s' % (cohorts, len(genes), args.B, args.k, args.base, args.agg))

    rows, orders, comps = [], {}, {}
    for held in cohorts:
        tr_names = [c for c in cohorts if c != held]
        train = {c: (data[c][0].values, data[c][1]) for c in tr_names}
        Xtr = np.vstack([train[c][0] for c in tr_names])
        ytr = np.concatenate([train[c][1] for c in tr_names])
        Xte, yte = data[held][0].values, data[held][1]

        t0 = time.time()
        res = rbs_fit(train, k=args.k, B=args.B, base=args.base, seed=args.seed, agg=args.agg)
        order = np.argsort(-res['score'])
        orders[held] = order
        secs = time.time() - t0

        icv = nested_cv(train, k=args.k, B=args.inner_B, base=args.base, seed=args.seed, agg=args.agg)

        # Within-cohort bootstrap stability, made COMPARABLE to the baselines.
        # The baselines report the mean pairwise Jaccard between B bootstrap top-k sets drawn
        # from the pooled training data. RBS draws its bootstraps inside each cohort, so the
        # naive analogue (overlap between the two cohorts' top-k) measures cross-cohort
        # disagreement instead and reads as ~0.007 - a different quantity entirely, and
        # misleading if placed in the same column.
        # For independent draws with per-gene selection frequency f, E|A n B| = sum f^2 and
        # |A|=|B|=k, so E[Jaccard] = sum f^2 / (2k - sum f^2). Averaged over training cohorts,
        # this is the same quantity the baselines report.
        js = []
        for c in tr_names:
            f = res['freq'][c]
            inter = float((f ** 2).sum())
            js.append(inter / max(2 * args.k - inter, 1e-9))
        stab_j = float(np.mean(js))
        # and the cross-cohort agreement, kept separately because it is genuinely informative
        sets = [np.argsort(-res['freq'][c])[:args.k] for c in tr_names]
        stab_k = mean_pairwise_jaccard(sets) if len(sets) > 1 else np.nan

        for K in TOPK:
            cols = order[:K]
            rows.append(dict(held_out=held, method=args.tag, K=K,
                             external_auroc=fit_score(Xtr, ytr, Xte, yte, cols),
                             nested_cv_auroc=icv[K],
                             jaccard=stab_j, cross_cohort_jaccard=stab_k, kuncheva=np.nan,
                             seconds=secs, B=args.B))
        log('  hold out %-11s %5.1fs  concordant=%d  ext AUROC = %s'
            % (held, secs, int(res['concordant'].sum()),
               ' '.join('%.3f' % r['external_auroc'] for r in rows[-len(TOPK):])))

        comps[held] = pd.DataFrame({
            'entrez_id': genes, 'symbol': [sym.get(g, '') for g in genes],
            'RBS': res['score'], 'I': res['I'], 'S': res['S'], 'R': res['R'], 'P': res['P'],
            'concordant': res['concordant'].astype(int),
        }).sort_values('RBS', ascending=False)
        comps[held].to_csv(os.path.join(args.out, 'components_holdout_%s.tsv' % held),
                           sep='\t', index=False)

    res_df = pd.DataFrame(rows)
    res_df.to_csv(os.path.join(args.out, 'results.tsv'), sep='\t', index=False)

    # ---- cross-fold stability: does the signature survive changing which cohort is held out?
    log('\n=== cross-fold signature agreement (Jaccard of top-K between LODO folds) ===')
    xf = []
    for K in TOPK:
        sets = [set(orders[h][:K]) for h in cohorts]
        js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
        xf.append(dict(method=args.tag, K=K, cross_fold_jaccard=float(np.mean(js))))
        log('  K=%-4d %.3f' % (K, np.mean(js)))
    pd.DataFrame(xf).to_csv(os.path.join(args.out, 'cross_fold_stability.tsv'),
                            sep='\t', index=False)
    log('\nwrote %s/' % args.out)


if __name__ == '__main__':
    main()
