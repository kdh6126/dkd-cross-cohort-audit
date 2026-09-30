#!/usr/bin/env python
"""Baseline feature selectors for the DKD cross-cohort benchmark.

Every selector has the same signature:

    select(X, y, k, rng, cohort=None) -> (selected_idx, scores)

  X        samples x genes, already per-cohort z-scored (numpy array)
  y        binary labels
  k        how many features to return (the intrinsic selectors are truncated/padded to k
           so that all methods are compared at matched sparsity)
  rng      numpy Generator, so every bootstrap replicate is reproducible
  scores   per-gene score, higher = more important (used for ranking / diagnostics)

Implemented here rather than pulled from boruta_py / skrebate / mrmr_selection deliberately:
those packages pin older numpy/sklearn versions and would make the whole benchmark fragile.
The algorithms are small enough that a direct implementation is the safer dependency.
"""
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier


# ------------------------------------------------------------------ helpers
def _f_scores(X, y):
    """Univariate ANOVA F statistic per gene (two-group case). Vectorised."""
    x1, x0 = X[y == 1], X[y == 0]
    n1, n0 = len(x1), len(x0)
    m1, m0 = x1.mean(0), x0.mean(0)
    m = X.mean(0)
    ssb = n1 * (m1 - m) ** 2 + n0 * (m0 - m) ** 2
    ssw = ((x1 - m1) ** 2).sum(0) + ((x0 - m0) ** 2).sum(0)
    dfw = max(n1 + n0 - 2, 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        f = (ssb / 1.0) / (ssw / dfw)
    return np.nan_to_num(f, nan=0.0, posinf=0.0)


def _topk(scores, k):
    k = min(k, len(scores))
    idx = np.argpartition(-scores, k - 1)[:k]
    return idx[np.argsort(-scores[idx])]


def _pad_to_k(chosen, scores, k):
    """Intrinsic selectors return a variable number of features; pad by score to reach k
    (and truncate if they overshoot) so sparsity is matched across methods."""
    chosen = list(dict.fromkeys(chosen))
    if len(chosen) >= k:
        s = np.array(chosen)
        return s[np.argsort(-scores[s])][:k]
    rest = [i for i in _topk(scores, len(scores)) if i not in set(chosen)]
    return np.array(chosen + rest[:k - len(chosen)], dtype=int)


# ------------------------------------------------------------------ selectors
def univariate(X, y, k, rng=None, cohort=None):
    """Plain F-test filter. The 'do nothing clever' baseline."""
    s = _f_scores(X, y)
    return _topk(s, k), s


def lasso(X, y, k, rng=None, cohort=None, C=1.0):
    # C chosen so the intrinsic selection lands near k=50 on these folds (46-70 nonzero),
    # so padding to k almost never fires and LASSO is not silently turned into a filter.
    m = LogisticRegression(solver='liblinear', l1_ratio=1.0, C=C, max_iter=3000)
    m.fit(X, y)
    coef = np.abs(m.coef_.ravel())
    return _pad_to_k(np.flatnonzero(coef > 0), coef, k), coef


def elastic_net(X, y, k, rng=None, cohort=None, C=0.05, l1_ratio=0.5):
    m = LogisticRegression(solver='saga', C=C, l1_ratio=l1_ratio,
                           max_iter=2000, tol=1e-3)
    m.fit(X, y)
    coef = np.abs(m.coef_.ravel())
    return _pad_to_k(np.flatnonzero(coef > 0), coef, k), coef


def mrmr(X, y, k, rng=None, cohort=None, pool=800):
    """Minimum Redundancy Maximum Relevance, FCQ variant.

    relevance = F statistic, redundancy = mean |Pearson r| with the already-selected set.
    Restricted to the top `pool` genes by relevance - the full 9,900 x 9,900 correlation
    matrix is neither affordable per bootstrap nor informative (mRMR never reaches down
    that far anyway).
    """
    rel = _f_scores(X, y)
    pool_idx = _topk(rel, min(pool, X.shape[1]))
    Xp = X[:, pool_idx]
    Xp = (Xp - Xp.mean(0)) / (Xp.std(0, ddof=1) + 1e-12)
    n = Xp.shape[0]
    R = np.abs((Xp.T @ Xp) / (n - 1))          # |correlation| within the pool
    np.fill_diagonal(R, 0.0)

    relp = rel[pool_idx]
    chosen = [int(np.argmax(relp))]
    red = R[chosen[0]].copy()
    for _ in range(min(k, len(pool_idx)) - 1):
        q = relp / (red / len(chosen) + 1e-12)
        q[chosen] = -np.inf
        nxt = int(np.argmax(q))
        chosen.append(nxt)
        red += R[nxt]
    order = pool_idx[np.array(chosen)]
    scores = np.zeros(X.shape[1])
    scores[order] = np.linspace(len(order), 1, len(order))   # rank-as-score
    return order, scores


def relieff(X, y, k, rng=None, cohort=None, n_neighbors=10, pool=2000):
    """ReliefF. Restricted to a relevance pool for tractability; distances are computed on
    the pooled features, which is what makes it affordable inside a bootstrap loop."""
    rel = _f_scores(X, y)
    pool_idx = _topk(rel, min(pool, X.shape[1]))
    Xp = X[:, pool_idx]
    n, p = Xp.shape
    span = Xp.max(0) - Xp.min(0)
    span[span == 0] = 1.0

    d = ((Xp[:, None, :] - Xp[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d, np.inf)

    W = np.zeros(p)
    kk = min(n_neighbors, max(1, min((y == 0).sum(), (y == 1).sum()) - 1))
    for i in range(n):
        same = np.flatnonzero(y == y[i])
        diff = np.flatnonzero(y != y[i])
        same = same[same != i]
        if len(same) == 0 or len(diff) == 0:
            continue
        hits = same[np.argsort(d[i, same])[:kk]]
        misses = diff[np.argsort(d[i, diff])[:kk]]
        W -= np.abs(Xp[i] - Xp[hits]).mean(0) / span
        W += np.abs(Xp[i] - Xp[misses]).mean(0) / span
    W /= n
    scores = np.full(X.shape[1], -np.inf)
    scores[pool_idx] = W
    return _topk(scores, k), np.nan_to_num(scores, neginf=0.0)


def rf_importance(X, y, k, rng=None, cohort=None, n_trees=300):
    seed = int(rng.integers(2 ** 31 - 1)) if rng is not None else 0
    m = RandomForestClassifier(n_estimators=n_trees, random_state=seed, n_jobs=-1,
                               class_weight='balanced_subsample')
    m.fit(X, y)
    imp = m.feature_importances_
    return _topk(imp, k), imp


def boruta(X, y, k, rng=None, cohort=None, max_iter=25, n_trees=150, alpha=0.05, pool=2000):
    """Boruta: compare each real feature against the best of a shadow (permuted) copy,
    accumulate hits, and test them against Binomial(iter, 0.5).

    Restricted to a relevance pool - Boruta on the full 9,900 genes means a 19,800-feature
    random forest per iteration, which is not affordable inside a bootstrap loop and mostly
    burns time on genes that never approach the shadow threshold.
    """
    from scipy.stats import binomtest
    rel = _f_scores(X, y)
    pool_idx = _topk(rel, min(pool, X.shape[1]))
    Xp = X[:, pool_idx]
    n, p = Xp.shape
    hits = np.zeros(p, dtype=int)
    seed0 = int(rng.integers(2 ** 31 - 1)) if rng is not None else 0

    for it in range(max_iter):
        r = np.random.default_rng(seed0 + it)
        shadow = np.column_stack([r.permutation(Xp[:, j]) for j in range(p)])
        Z = np.hstack([Xp, shadow])
        m = RandomForestClassifier(n_estimators=n_trees, random_state=seed0 + it,
                                   n_jobs=-1, class_weight='balanced_subsample')
        m.fit(Z, y)
        imp = m.feature_importances_
        hits += (imp[:p] > imp[p:].max()).astype(int)

    pv = np.array([binomtest(int(h), max_iter, 0.5, alternative='greater').pvalue for h in hits])
    confirmed = np.flatnonzero(pv < alpha / p)          # Bonferroni, as in the original method
    scores = np.full(X.shape[1], -np.inf)
    scores[pool_idx] = hits.astype(float)
    scores = np.nan_to_num(scores, neginf=0.0)
    return _pad_to_k(list(pool_idx[confirmed]), scores, k), scores


SELECTORS = {
    'univariate': univariate,
    'lasso': lasso,
    'elastic_net': elastic_net,
    'mrmr': mrmr,
    'relieff': relieff,
    'rf': rf_importance,
    'boruta': boruta,
}
