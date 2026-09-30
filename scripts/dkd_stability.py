#!/usr/bin/env python
"""Bootstrap stability selection + Leave-One-Dataset-Out evaluation.

Protocol
--------
For each LODO fold (one cohort held out entirely):
  1. Draw B bootstrap resamples of the TRAINING rows, stratified by (cohort, label) so no
     replicate collapses to one class or loses a cohort.
  2. Run each selector on each replicate, keep its top-k set.
  3. selection frequency f(g) = fraction of replicates in which g was selected.
  4. Final signature at size K = top-K genes by f, ties broken by mean rank.
  5. Fit a plain L2 logistic model on the FULL training data restricted to those K genes and
     score the held-out cohort (AUROC). The classifier is deliberately simple - this measures
     the features, not the model.

Nothing supervised ever sees the held-out cohort. Per-cohort z-scoring is unsupervised
(see dkd_data) so applying it to the test cohort is domain alignment, not leakage.

Stability metrics
-----------------
  within-cohort  mean pairwise Jaccard between the B bootstrap top-k sets
  Kuncheva index chance-corrected agreement, comparable across different k and p
"""
import os, sys, json, time, itertools
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_selectors import SELECTORS


def stratified_bootstrap(rng, y, cohort):
    """Resample with replacement inside each (cohort, label) cell."""
    idx = []
    for c in np.unique(cohort):
        for lab in (0, 1):
            cell = np.flatnonzero((cohort == c) & (y == lab))
            if len(cell):
                idx.append(rng.choice(cell, size=len(cell), replace=True))
    return np.concatenate(idx)


def jaccard(a, b):
    a, b = set(a), set(b)
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def mean_pairwise_jaccard(sets, max_pairs=2000, rng=None):
    n = len(sets)
    pairs = list(itertools.combinations(range(n), 2))
    if len(pairs) > max_pairs and rng is not None:
        pairs = [pairs[i] for i in rng.choice(len(pairs), max_pairs, replace=False)]
    return float(np.mean([jaccard(sets[i], sets[j]) for i, j in pairs])) if pairs else 0.0


def kuncheva(sets, p, max_pairs=2000, rng=None):
    """Kuncheva consistency index: chance-corrected overlap of equal-size subsets."""
    n = len(sets)
    pairs = list(itertools.combinations(range(n), 2))
    if len(pairs) > max_pairs and rng is not None:
        pairs = [pairs[i] for i in rng.choice(len(pairs), max_pairs, replace=False)]
    vals = []
    for i, j in pairs:
        a, b = set(sets[i]), set(sets[j])
        k = (len(a) + len(b)) / 2
        if k <= 0 or k >= p:
            continue
        r = len(a & b)
        exp = k * k / p
        denom = k - exp
        if denom != 0:
            vals.append((r - exp) / denom)
    return float(np.mean(vals)) if vals else 0.0


def fit_and_score(Xtr, ytr, Xte, yte, cols):
    """Plain L2 logistic model on the selected genes; AUROC on the held-out cohort."""
    m = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000, class_weight='balanced')
    m.fit(Xtr[:, cols], ytr)
    if len(np.unique(yte)) < 2:
        return np.nan
    return roc_auc_score(yte, m.predict_proba(Xte[:, cols])[:, 1])


def nested_cv_auc(X, y, method, K, k, seed=0, folds=5):
    """Honest internal estimate: feature selection is re-run INSIDE each CV fold.

    Without this nesting the internal number is meaningless - selecting on all of the
    training data and then cross-validating within it gives AUROC ~0.99 for every method
    (we measured it). Selection here is single-shot rather than bootstrap-stabilised, purely
    for cost; it is still a nested, unbiased estimate of the select-then-classify pipeline.
    """
    from dkd_selectors import SELECTORS
    fn = SELECTORS[method]
    y = np.asarray(y)
    n_min = min((y == 0).sum(), (y == 1).sum())
    folds = max(2, min(folds, n_min))
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    rng = np.random.default_rng(seed)
    aucs = []
    for tr, te in skf.split(X, y):
        if len(np.unique(y[te])) < 2 or len(np.unique(y[tr])) < 2:
            continue
        idx, _ = fn(X[tr], y[tr], k=k, rng=rng)
        cols = np.asarray(idx, dtype=int)[:K]
        m = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000, class_weight='balanced')
        m.fit(X[tr][:, cols], y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(X[te][:, cols])[:, 1]))
    return float(np.mean(aucs)) if aucs else np.nan


def run_fold(method, Xtr, ytr, ctr, Xte, yte, genes, B, k, seed, log=print):
    """Bootstrap stability selection for one method on one LODO fold."""
    fn = SELECTORS[method]
    rng = np.random.default_rng(seed)
    p = Xtr.shape[1]
    counts = np.zeros(p)
    ranksum = np.zeros(p)
    seen = np.zeros(p)
    sets = []
    t0 = time.time()
    for b in range(B):
        bi = stratified_bootstrap(rng, ytr, ctr)
        Xb, yb = Xtr[bi], ytr[bi]
        if len(np.unique(yb)) < 2:
            continue
        idx, _ = fn(Xb, yb, k=k, rng=rng)
        idx = np.asarray(idx, dtype=int)
        counts[idx] += 1
        ranksum[idx] += np.arange(1, len(idx) + 1)
        seen[idx] += 1
        sets.append(idx)
    freq = counts / max(len(sets), 1)
    meanrank = np.where(seen > 0, ranksum / np.maximum(seen, 1), np.inf)
    return dict(freq=freq, meanrank=meanrank, sets=sets, seconds=time.time() - t0)
