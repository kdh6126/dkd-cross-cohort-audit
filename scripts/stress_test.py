#!/usr/bin/env python
"""Standalone perturbation stress test of a FINISHED signature.

This is the project's stage 5: discovery runs on original data, and only then are the selected
biomarkers stress-tested. So this script never selects anything - it takes a fixed gene list,
trains one plain logistic model on unperturbed training data, and measures how the held-out
AUROC degrades as the TEST cohort is perturbed with increasing severity.

Perturbations (severity is the sweep axis):
  subsample        drop a fraction of test samples, stratified
  noise            add gaussian noise at s x the per-gene SD
  mask             zero a fraction of the signature's genes (measurement dropout)
  batch            shift a random half of the test samples by a per-gene offset

Every signature is compared against random gene sets of the same size, because on this data a
random signature already achieves AUROC 0.68-0.88 (scripts/null_control.py). Absolute retention
means little; retention *relative to random* is the question.
"""
import os, sys, argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM

SEVERITIES = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def perturb_test(X, y, kind, s, rng):
    """Perturb the TEST matrix only. s = severity (0 means untouched)."""
    if s == 0:
        return X, y, np.ones(X.shape[1], bool)
    avail = np.ones(X.shape[1], bool)
    if kind == 'subsample':
        frac = max(0.15, 1.0 - 0.4 * s)
        keep = []
        for lab in (0, 1):
            cell = np.flatnonzero(y == lab)
            m = max(2, int(round(frac * len(cell))))
            keep.append(rng.choice(cell, m, replace=False))
        keep = np.concatenate(keep)
        return X[keep], y[keep], avail
    if kind == 'noise':
        sd = X.std(0, ddof=1)
        return X + rng.normal(0, s * sd, X.shape), y, avail
    if kind == 'mask':
        frac = min(0.8, 0.2 * s + 0.1)
        avail = rng.random(X.shape[1]) > frac
        Xm = X.copy()
        Xm[:, ~avail] = 0.0
        return Xm, y, avail
    if kind == 'batch':
        half = rng.random(len(X)) < 0.5
        Xb = X.copy()
        Xb[half] += rng.normal(0, 0.5 * s, X.shape[1])
        return Xb, y, avail
    raise ValueError(kind)


def auc_of(model, Xte, yte, cols):
    if len(np.unique(yte)) < 2:
        return np.nan
    return roc_auc_score(yte, model.predict_proba(Xte[:, cols])[:, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', default='results/final_candidates_gated.tsv')
    ap.add_argument('--topn', type=int, default=30)
    ap.add_argument('--cohorts', default=','.join(GLOM) + ',GSE142025')
    ap.add_argument('--n-rep', type=int, default=40)
    ap.add_argument('--n-random', type=int, default=40)
    ap.add_argument('--seed', type=int, default=20260824)
    ap.add_argument('--out', default='results/stress_test.tsv')
    args = ap.parse_args()

    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])
    gidx = {g: i for i, g in enumerate(genes)}

    sym2id = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = sym2id.dropna(subset=['symbol']).set_index('symbol')['entrez_id'].to_dict()
    cand = pd.read_csv(args.candidates, sep='\t').head(args.topn)
    sig = [gidx[sym2id[s]] for s in cand['symbol'].dropna().astype(str)
           if s in sym2id and sym2id[s] in gidx]
    log('signature: %d genes on the shared space\n' % len(sig))

    rng = np.random.default_rng(args.seed)
    rows = []
    for held in cohorts:
        tr = [c for c in cohorts if c != held]
        Xtr = np.vstack([data[c][0].values for c in tr])
        ytr = np.concatenate([data[c][1] for c in tr])
        Xte, yte = data[held][0].values, data[held][1]

        model = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000,
                                   class_weight='balanced').fit(Xtr[:, sig], ytr)
        rand_sets = [rng.choice(len(genes), len(sig), replace=False) for _ in range(args.n_random)]
        rand_models = [LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000,
                                          class_weight='balanced').fit(Xtr[:, r], ytr)
                       for r in rand_sets]

        log('=== held out %s ===' % held)
        for kind in ('subsample', 'noise', 'mask', 'batch'):
            line = []
            for s in SEVERITIES:
                a_sig, a_rnd = [], []
                reps = 1 if s == 0 else args.n_rep
                for _ in range(reps):
                    Xp, yp, _ = perturb_test(Xte, yte, kind, s, rng)
                    a_sig.append(auc_of(model, Xp, yp, sig))
                for _ in range(1 if s == 0 else max(args.n_rep // 4, 5)):
                    Xp, yp, _ = perturb_test(Xte, yte, kind, s, rng)
                    a_rnd += [auc_of(m, Xp, yp, r) for m, r in zip(rand_models, rand_sets)]
                ms, mr = float(np.nanmean(a_sig)), float(np.nanmean(a_rnd))
                rows.append(dict(held_out=held, perturbation=kind, severity=s,
                                 signature_auroc=ms, random_auroc=mr, margin=ms - mr))
                line.append('%.3f/%.3f' % (ms, mr))
            log('  %-10s ' % kind + '  '.join('s=%-4g %s' % (s, v) for s, v in zip(SEVERITIES, line)))
        log('')

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)

    log('=== retention, averaged over folds (signature AUROC / margin over random) ===')
    log('  %-11s %s' % ('perturbation', '  '.join('s=%-11g' % s for s in SEVERITIES)))
    for kind in ('subsample', 'noise', 'mask', 'batch'):
        sub = df[df['perturbation'] == kind].groupby('severity')
        a = sub['signature_auroc'].mean()
        m = sub['margin'].mean()
        log('  %-11s %s' % (kind, '  '.join('%.3f (%+.3f)' % (a[s], m[s]) for s in SEVERITIES)))
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
