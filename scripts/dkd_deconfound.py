#!/usr/bin/env python
"""Confounder-aware selection: move the procurement correction into the objective.

At present the procurement artifact is removed *after* selection, by a gate that needs
non-diabetic kidney disease biopsies - available in only 2 of 11 public datasets. If the
correction can be moved inside the selector, it works in any cohort.

The handle: the immediate-early module is a per-sample, measurable proxy for how the tissue was
handled. Its mean z-score is a scalar per sample, computable from the same matrix, needing no
extra labels. Two ways to use it:

  ORTH     residualise every gene on the handling score within each cohort, then select on the
           residuals. Aggressive: the score is correlated with case status (cases are biopsies),
           so this also removes real disease signal.
  PENALTY  keep the data, and multiply each gene's RBS by (1 - |corr(gene, handling score)|).
           Softer: it demotes genes that track handling without deleting the axis.

Which is better is an empirical question with a real chance of a negative answer, so all three
arms (none / orth / penalty) are scored on the same three axes:

  specificity pass-rate   fraction of the selected top-K that survive the DKD-vs-other-CKD gate
                          -> did the correction do what it was for?
  external AUROC          -> did it cost predictive signal?
  cross-fold Jaccard      -> did it cost reproducibility?

A correction that raises the pass-rate while destroying AUROC has not solved anything.
"""
import os, sys, argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM
from dkd_rbs import rbs_fit

IEG = ['FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'ATF3', 'DUSP1', 'ZFP36', 'NR4A2', 'BTG2',
       'EGR2', 'EGR3', 'IER2', 'KLF2', 'KLF4', 'SOCS3', 'DUSP2', 'JUND', 'NR4A1']
TOPK = [10, 20, 50]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def handling_score(Z, ieg_cols):
    """Per-sample mean z-score of the immediate-early module, standardised."""
    s = Z[:, ieg_cols].mean(axis=1)
    sd = s.std(ddof=1)
    return (s - s.mean()) / (sd if sd > 0 else 1.0)


def residualise(Z, s, y=None, lam=1.0):
    """Remove the handling axis from every gene.

    Two estimators for the per-gene loading on the handling score:

      TOTAL slope (y=None) - ordinary OLS over all samples. Because cases are biopsies and
      controls are not, s is nearly collinear with case status, so this slope absorbs the
      case/control difference itself and OVER-corrects: it deletes real disease signal along
      with the confounder.

      WITHIN-GROUP slope (y given) - centre s and Z inside each class before estimating the
      slope. The group means, which carry the confounded between-group difference, are removed
      first, so the slope estimates the handling loading from variation that is not contaminated
      by disease status. Subtracting s * beta_within then removes exactly as much of the group
      difference as the handling axis explains *at the rate observed within groups*. This is the
      standard covariate adjustment for Simpson-type confounding.

    Note the limit of both: with procurement perfectly collinear with case status, no
    within-cohort method can fully separate them. The within-group slope helps by avoiding
    over-correction, not by identifying the confounder from data that cannot identify it.

    lam scales the correction (1.0 = full removal) for a shrinkage sweep.
    """
    if y is None:
        sw, Zw = s, Z
    else:
        sw, Zw = s.copy(), Z.copy()
        for lab in np.unique(y):
            m = (y == lab)
            if m.sum() < 2:
                continue
            sw[m] = s[m] - s[m].mean()
            Zw[m] = Z[m] - Z[m].mean(axis=0)
    denom = float(sw @ sw)
    if denom == 0:
        return Z.copy()
    beta = (sw @ Zw) / denom
    return Z - lam * np.outer(s - s.mean(), beta)


def ext_auc(Xtr, ytr, Xte, yte, cols):
    m = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000,
                           class_weight='balanced').fit(Xtr[:, cols], ytr)
    if len(np.unique(yte)) < 2:
        return np.nan
    return roc_auc_score(yte, m.predict_proba(Xte[:, cols])[:, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--B', type=int, default=200)
    ap.add_argument('--k', type=int, default=50)
    ap.add_argument('--base', default='relieff')
    ap.add_argument('--cohorts', default=','.join(GLOM) + ',GSE142025')
    ap.add_argument('--seed', type=int, default=20260824)
    ap.add_argument('--out', default='results/deconfound.tsv')
    args = ap.parse_args()

    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])
    gpos = {g: i for i, g in enumerate(genes)}

    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym = dict(zip(gs['entrez_id'], gs['symbol']))
    sym2id = {s: g for g, s in sym.items() if isinstance(s, str)}
    ieg_cols = [gpos[sym2id[s]] for s in IEG if s in sym2id and sym2id[s] in gpos]
    log('handling score from %d immediate-early genes present in the space' % len(ieg_cols))

    spec = pd.read_csv('results/final_candidates.tsv', sep='\t')
    spec = spec.rename(columns={spec.columns[0]: 'entrez_id'})
    spec['entrez_id'] = spec['entrez_id'].astype(str)
    passes = dict(zip(spec['entrez_id'], spec['passes_specificity'].astype(bool)))
    proc = dict(zip(spec['entrez_id'], spec['procurement_flag'].astype(bool)))

    rows, orders = [], {}
    ARMS = ('none', 'orth', 'orth_within', 'orth_within_0.5', 'penalty')
    for arm in ARMS:
        log('\n=== arm: %s ===' % arm)
        per_fold_orders = {}
        for held in cohorts:
            tr = [c for c in cohorts if c != held]
            train = {}
            for c in tr:
                Z = data[c][0].values.copy()
                s = handling_score(Z, ieg_cols)
                yc = data[c][1]
                if arm == 'orth':
                    Z = residualise(Z, s)                          # total slope
                elif arm == 'orth_within':
                    Z = residualise(Z, s, y=yc)                    # within-group slope
                elif arm == 'orth_within_0.5':
                    Z = residualise(Z, s, y=yc, lam=0.5)           # shrunk within-group slope
                train[c] = (Z, yc)
            res = rbs_fit(train, k=args.k, B=args.B, base=args.base,
                          seed=args.seed, agg='geomean')
            score = res['score'].copy()

            if arm == 'penalty':
                # correlation of each gene with the handling score, averaged over cohorts
                cors = []
                for c in tr:
                    Z = data[c][0].values
                    s = handling_score(Z, ieg_cols)
                    zc = (Z - Z.mean(0)) / (Z.std(0, ddof=1) + 1e-12)
                    cors.append(np.abs(zc.T @ s) / (len(s) - 1))
                cor = np.clip(np.mean(cors, axis=0), 0, 1)
                score = score * (1.0 - cor)

            order = np.argsort(-score)
            per_fold_orders[held] = order

            Xtr = np.vstack([data[c][0].values for c in tr])   # evaluation always on RAW data
            ytr = np.concatenate([data[c][1] for c in tr])
            Xte, yte = data[held][0].values, data[held][1]
            for K in TOPK:
                cols = order[:K]
                gid = genes[cols]
                npass = sum(1 for g in gid if passes.get(g, False))
                nproc = sum(1 for g in gid if proc.get(g, False))
                rows.append(dict(arm=arm, held_out=held, K=K,
                                 external_auroc=ext_auc(Xtr, ytr, Xte, yte, cols),
                                 spec_pass_frac=npass / K, procurement_frac=nproc / K))
            log('  hold out %-11s top10: %s' % (held, ', '.join(
                sym.get(g, g) for g in genes[order[:10]])))
        orders[arm] = per_fold_orders

    res = pd.DataFrame(rows)
    for arm in orders:
        for K in TOPK:
            sets = [set(orders[arm][h][:K]) for h in cohorts]
            js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
            res.loc[(res['arm'] == arm) & (res['K'] == K), 'cross_fold_jaccard'] = float(np.mean(js))
    res.to_csv(args.out, sep='\t', index=False)

    log('\n=== summary (mean over LODO folds) ===')
    log('  %-9s %-4s %-16s %-18s %-18s %s'
        % ('arm', 'K', 'external AUROC', 'specificity pass', 'procurement-flagged', 'cross-fold J'))
    agg = (res.groupby(['arm', 'K'])
              .agg(auc=('external_auroc', 'mean'), sp=('spec_pass_frac', 'mean'),
                   pr=('procurement_frac', 'mean'), xf=('cross_fold_jaccard', 'mean'))
              .reset_index())
    for _, r in agg.sort_values(['K', 'arm']).iterrows():
        log('  %-9s %-4d %-16.3f %-18.0f%% %-17.0f%% %.3f'
            % (r['arm'], r['K'], r['auc'], 100 * r['sp'], 100 * r['pr'], r['xf']))
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
