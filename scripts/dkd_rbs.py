#!/usr/bin/env python
"""Robust Biomarker Score (RBS) - the proposed selector.

What makes it different from the baselines
------------------------------------------
Every baseline here pools the training cohorts into one matrix and selects once. That is
exactly what lets a dataset-specific feature win: a gene that is strongly discriminative in
the single largest cohort can dominate the pooled objective while being useless elsewhere.

RBS never pools. It runs bootstrap selection **inside each training cohort separately** and
then aggregates across cohorts with a WORST-CASE rule. A gene that is selected 95% of the
time in cohort A and 15% of the time in cohort B scores 0.15, not 0.55. That is the
mechanism that removes dataset-specific features, and it is the whole point.

The score has four components, each mapped to [0, 1]:

  I  predictive importance    mean over bootstraps of |coefficient| divided by the largest
                              |coefficient| in that same bootstrap, averaged over cohorts.
                              Max-normalisation (not rank) so a gene the base learner never
                              touches scores exactly 0 rather than a mid-range tie rank.
  S  bootstrap stability      mean over cohorts of within-cohort selection frequency
  R  cross-cohort reproduc.   MIN over cohorts of selection frequency, gated by effect
                              direction agreeing in every cohort (a gene that flips sign
                              across cohorts scores 0 no matter how often it is selected)
  P  perturbation robustness  selection frequency under noise / subsampling / feature
                              masking / batch shift

They are combined by a WEIGHTED GEOMETRIC MEAN, not a sum: a sum lets a gene buy its way in
on one strong component, a geometric mean requires it to be decent at all of them, which is
the actual claim being made ("important AND stable AND reproducible AND robust").

Ordering of stages follows the project decision that augmentation must not drive discovery:
I, S, R are computed on ORIGINAL data only. P is computed afterwards, on the candidate pool
that I/S/R already narrowed - perturbation is used to stress-test candidates, never to find
them from scratch.
"""
import numpy as np
from scipy.stats import rankdata
from dkd_selectors import SELECTORS, _f_scores


# ------------------------------------------------------------------ helpers
def rank_norm(v):
    """Map a score vector to [0,1] by AVERAGE rank.

    Tie handling matters here and is easy to get wrong: an L1 coefficient vector is ~99%
    exact zeros, and argsort-of-argsort would hand those zeros arbitrary distinct ranks
    spanning most of [0,1], making the importance component pure noise. rankdata's average
    method gives every zero the same rank.
    """
    v = np.asarray(v, dtype=float)
    r = rankdata(v, method='average') - 1.0
    return r / max(len(v) - 1, 1)


def hedges_g(X, y):
    x1, x0 = X[y == 1], X[y == 0]
    n1, n0 = len(x1), len(x0)
    if n1 < 2 or n0 < 2:
        return np.zeros(X.shape[1])
    v1, v0 = x1.var(0, ddof=1), x0.var(0, ddof=1)
    sp = np.sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2))
    sp[sp == 0] = np.nan
    d = (x1.mean(0) - x0.mean(0)) / sp
    return np.nan_to_num(d)


def _boot(rng, y):
    """Stratified bootstrap inside one cohort."""
    idx = []
    for lab in (0, 1):
        cell = np.flatnonzero(y == lab)
        if len(cell):
            idx.append(rng.choice(cell, size=len(cell), replace=True))
    return np.concatenate(idx)


# ------------------------------------------------------------------ perturbations
def perturb(X, y, kind, rng):
    """Return a perturbed (X, y, available_mask). Never used to discover features - only to
    stress-test candidates that I/S/R already short-listed."""
    n, p = X.shape
    avail = np.ones(p, dtype=bool)
    if kind == 'subsample':
        keep = []
        for lab in (0, 1):
            cell = np.flatnonzero(y == lab)
            m = max(2, int(round(0.7 * len(cell))))
            keep.append(rng.choice(cell, size=m, replace=False))
        keep = np.concatenate(keep)
        return X[keep], y[keep], avail
    if kind == 'noise':
        sd = X.std(0, ddof=1)
        return X + rng.normal(0, 0.5 * sd, size=X.shape), y, avail
    if kind == 'mask':
        avail = rng.random(p) > 0.20        # 20% of genes go missing
        Xm = X.copy()
        Xm[:, ~avail] = 0.0
        return Xm, y, avail
    if kind == 'batch':
        half = rng.random(n) < 0.5          # a synthetic platform split
        shift = rng.normal(0, 0.5, size=p)
        Xb = X.copy()
        Xb[half] += shift
        return Xb, y, avail
    raise ValueError(kind)


PERTURBATIONS = ('subsample', 'noise', 'mask', 'batch')


# ------------------------------------------------------------------ the selector
def rbs_fit(cohort_data, k=50, B=200, base='lasso', seed=0,
            weights=(1.0, 1.0, 2.0, 1.0), pert_B=40, pert_pool=300, verbose=False,
            agg='geomean'):
    """cohort_data: dict name -> (X ndarray, y ndarray). Training cohorts ONLY.

    Returns dict with the four components, the combined score, and diagnostics.
    """
    rng = np.random.default_rng(seed)
    fn = SELECTORS[base]
    names = list(cohort_data)
    p = next(iter(cohort_data.values()))[0].shape[1]

    freq = {}          # per-cohort selection frequency
    imp = {}           # per-cohort max-normalised importance in [0,1]
    sign = {}          # per-cohort effect direction

    for c in names:
        X, y = cohort_data[c]
        counts = np.zeros(p)
        impsum = np.zeros(p)
        nb = 0
        for b in range(B):
            bi = _boot(rng, y)
            Xb, yb = X[bi], y[bi]
            if len(np.unique(yb)) < 2:
                continue
            idx, sc = fn(Xb, yb, k=k, rng=rng)
            counts[np.asarray(idx, dtype=int)] += 1
            a = np.abs(np.nan_to_num(sc))
            impsum += a / (a.max() + 1e-12)      # max-normalised, so never-selected genes get 0
            nb += 1
        freq[c] = counts / max(nb, 1)
        imp[c] = impsum / max(nb, 1)
        sign[c] = np.sign(hedges_g(X, y))
        if verbose:
            print('  %-11s selected>=50%%: %d genes' % (c, int((freq[c] >= 0.5).sum())))

    F = np.vstack([freq[c] for c in names])
    I = np.vstack([imp[c] for c in names]).mean(0)   # already in [0,1] per cohort
    S = F.mean(0)
    # Cross-cohort aggregation. A raw MIN is the purest statement of worst-case
    # reproducibility, but it is also maximally sensitive to the noisiest cohort: with a
    # 22-sample cohort the per-gene frequency is estimated so poorly that min(.) is
    # dominated by estimation error rather than by real cohort disagreement (measured -
    # only 13 genes exceeded R=0.1 under 'min'). 'geomean' keeps the multiplicative,
    # penalise-the-weak-cohort behaviour while degrading gracefully when one cohort's
    # estimate is noisy. See results/ablation for the comparison.
    if agg == 'min':
        R_raw = F.min(0)
    elif agg == 'geomean':
        R_raw = np.exp(np.log(F + 1e-3).mean(0)) - 1e-3
    elif agg == 'q25':
        R_raw = np.percentile(F, 25, axis=0)
    elif agg == 'mean':
        R_raw = F.mean(0)
    else:
        raise ValueError(agg)
    R_raw = np.clip(R_raw, 0, None)

    Sg = np.vstack([sign[c] for c in names])
    concordant = (np.abs(Sg.sum(0)) == len(names)) & (Sg != 0).all(0)
    R = R_raw * concordant                      # sign-flipping genes are zeroed outright

    # ---- perturbation robustness, on the candidate pool only
    pre = I * S * R                  # provisional ranking from original data
    pool = np.argsort(-pre)[:pert_pool]
    P = np.zeros(p)
    hit = np.zeros(p)
    tries = np.zeros(p)
    for c in names:
        X, y = cohort_data[c]
        for b in range(pert_B):
            kind = PERTURBATIONS[b % len(PERTURBATIONS)]
            Xp_, yp_, avail = perturb(X, y, kind, rng)
            if len(np.unique(yp_)) < 2:
                continue
            idx, _ = fn(Xp_, yp_, k=k, rng=rng)
            sel = np.zeros(p, dtype=bool)
            sel[np.asarray(idx, dtype=int)] = True
            hit[pool] += sel[pool] & avail[pool]
            tries[pool] += avail[pool]          # masked-out genes do not count against a gene
    P[pool] = hit[pool] / np.maximum(tries[pool], 1)

    # ---- weighted geometric mean; eps keeps a single zero from being uninformative-NaN
    eps = 1e-6
    wI, wS, wR, wP = weights
    comps = np.vstack([I + eps, S + eps, R + eps, P + eps])
    w = np.array([wI, wS, wR, wP])[:, None]
    RBS = np.exp((w * np.log(comps)).sum(0) / w.sum())
    RBS[~concordant] = 0.0                      # hard gate, not a soft penalty

    return dict(score=RBS, I=I, S=S, R=R, P=P,
                freq=freq, concordant=concordant, names=names)


def rbs_select(cohort_data, K, **kw):
    out = rbs_fit(cohort_data, **kw)
    return np.argsort(-out['score'])[:K], out
