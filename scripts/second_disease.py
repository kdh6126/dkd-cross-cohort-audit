#!/usr/bin/env python
"""External replication in a second disease and a second organ: heart failure (GSE5406).

Everything so far rests on one disease. GSE5406 has the same structural features that made the
DKD analysis possible, in a completely different organ:

    194  explanted left ventricle at cardiac transplantation   (108 ischaemic, 86 idiopathic)
     16  unused donor heart with normal LV function

so it contains both contrasts:
    disease vs donor        explant vs donor heart  -> procurement differs
    ischaemic vs idiopathic explant vs explant      -> procurement matched

Four things are tested, in the order they were established in kidney:

  1. Does the immediate-early module separate by PROCUREMENT (explant vs donor)?
  2. Does it stay silent in the procurement-MATCHED contrast, as it did in KPMP?
  3. What is the block size of the IEG module relative to the top-connectivity block?
  4. Given that size ratio, does connectivity avoid the IEG module while stability does not -
     which is what the simulation predicts?

If all four behave as in kidney, the framework is not a DKD peculiarity. If they do not, that is
the single most important thing to report.
"""
import os, sys, gzip, argparse
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_deconfound import IEG
from per_dataset_features import annot_map, read_series_matrix

MATRIX = 'data/raw/geo/GSE5406/GSE5406_series_matrix.txt.gz'
OUT = 'results/second_disease'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    n1, n0 = len(a), len(b)
    if n1 < 2 or n0 < 2:
        return np.full(a.shape[1], np.nan)
    sp = np.sqrt(((n1 - 1) * a.var(0, ddof=1) + (n0 - 1) * b.var(0, ddof=1)) / (n1 + n0 - 2))
    sp[sp == 0] = np.nan
    return (a.mean(0) - b.mean(0)) / sp * (1 - 3 / (4 * (n1 + n0) - 9))


def collapse(df, pmap):
    keep = [p for p in df.index if p in pmap]
    sub = df.loc[keep]
    g = np.array([pmap[p] for p in keep])
    means = sub.mean(axis=1, skipna=True).values
    order = np.lexsort((-means, g))
    seen, sel = set(), []
    for i in order:
        if g[i] not in seen:
            seen.add(g[i]); sel.append(i)
    out = sub.iloc[sel]
    out.index = g[sel]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--K', type=int, default=50)
    ap.add_argument('--power', type=int, default=6)
    ap.add_argument('--module-frac', type=float, default=0.10)
    ap.add_argument('--B', type=int, default=100)
    ap.add_argument('--seed', type=int, default=20260825)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    df, hdr = read_series_matrix(MATRIX)
    log('loaded %d probes x %d samples' % df.shape)
    chars = hdr['Sample_characteristics_ch1'][0]
    src = hdr['Sample_source_name_ch1'][0]
    grp = np.array(['donor' if 'donor' in s else
                    ('ischaemic' if 'ischemic' in c else 'idiopathic')
                    for c, s in zip(chars, src)])
    proc = np.array(['donor_heart' if 'donor' in s else 'explant' for s in src])
    log('  groups: %s' % pd.Series(grp).value_counts().to_dict())

    X = collapse(df, annot_map('GPL96')).T.astype(float)
    X = X.loc[:, ~X.isna().all(axis=0)]
    log('  collapsed to %d Entrez genes' % X.shape[1])
    Z = ((X - X.mean()) / X.std(ddof=1).replace(0, np.nan)).fillna(0.0)

    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym = dict(zip(gs['entrez_id'], gs['symbol']))
    sym2id = {s: g for g, s in sym.items() if isinstance(s, str)}
    cols = np.array([str(c) for c in Z.columns])
    ieg_ids = [sym2id[s] for s in IEG if s in sym2id and sym2id[s] in set(cols)]
    ieg_pos = np.array([int(np.where(cols == g)[0][0]) for g in ieg_ids])
    log('  immediate-early genes present: %d' % len(ieg_pos))

    V = Z.values
    rows = []

    # ---------------------------------------------------------------- 1 & 2
    log('\n=== 1-2. the two contrasts ===')
    for name, a, b, note in [
            ('disease vs donor', V[grp != 'donor'], V[grp == 'donor'], 'explant vs donor - CONFOUNDED'),
            ('ischaemic vs idiopathic', V[grp == 'ischaemic'], V[grp == 'idiopathic'],
             'explant vs explant - MATCHED')]:
        g = pd.Series(hedges_g(a, b), index=cols)
        ieg = g.iloc[ieg_pos].dropna()
        rest = g.drop(index=cols[ieg_pos]).dropna()
        p = stats.mannwhitneyu(ieg, rest, alternative='two-sided').pvalue
        log('  %-24s (%3d vs %3d)  %-32s IEG mean g = %+.2f  rest %+.2f  p = %.2e'
            % (name, len(a), len(b), note, ieg.mean(), rest.mean(), p))
        rows.append(dict(test=name, n_a=len(a), n_b=len(b), ieg_g=ieg.mean(),
                         rest_g=rest.mean(), p=p))

    # ---------------------------------------------------------------- 3. block size
    log('\n=== 3. block sizes inside the disease-correlated module ===')
    y = (grp != 'donor').astype(float)
    yz = (y - y.mean()) / (y.std(ddof=1) + 1e-12)
    n = len(y)
    gsig = (V.T @ yz) / (n - 1)
    m = max(50, int(args.module_frac * V.shape[1]))
    mod = np.argsort(-np.abs(gsig))[:m]
    S = V[:, mod]
    R = np.abs((S.T @ S) / (n - 1))
    A = R ** args.power
    np.fill_diagonal(A, 0.0)
    kIN = A.sum(1)
    blk = (R > 0.5).sum(1)
    pos_in_mod = {g: i for i, g in enumerate(mod)}
    ieg_in = [pos_in_mod[p_] for p_ in ieg_pos if p_ in pos_in_mod]
    top_conn = np.argsort(-kIN)[:args.K]
    log('  module: %d genes;  IEG inside the module: %d / %d'
        % (m, len(ieg_in), len(ieg_pos)))
    if ieg_in:
        log('  block size (module genes with |r| > 0.5)')
        log('    immediate-early : median %d' % int(np.median(blk[ieg_in])))
        log('    top-%d by kIN   : median %d' % (args.K, int(np.median(blk[top_conn]))))
        log('    whole module    : median %d' % int(np.median(blk)))
        ratio = np.median(blk[top_conn]) / max(np.median(blk[ieg_in]), 1)
        log('  size ratio (disease block / IEG block) = %.1f' % ratio)
        rows.append(dict(test='block_size', ieg_block=float(np.median(blk[ieg_in])),
                         top_block=float(np.median(blk[top_conn])), ratio=float(ratio)))

    # ---------------------------------------------------------------- 4. selectors
    log('\n=== 4. do the selectors behave as the size ratio predicts? ===')
    rng = np.random.default_rng(args.seed)
    x1, x0 = V[y == 1], V[y == 0]
    sp = np.sqrt(x1.var(0, ddof=1) / len(x1) + x0.var(0, ddof=1) / len(x0)) + 1e-12
    eff = np.abs(x1.mean(0) - x0.mean(0)) / sp
    conn = np.zeros(V.shape[1]); conn[mod] = kIN / kIN.max()
    cnt = np.zeros(V.shape[1])
    for _ in range(args.B):
        idx = np.concatenate([rng.choice(np.flatnonzero(y == lab), int((y == lab).sum()),
                                         replace=True) for lab in (0, 1)])
        a1, a0 = V[idx][y[idx] == 1], V[idx][y[idx] == 0]
        s2 = np.sqrt(a1.var(0, ddof=1) / max(len(a1), 1)
                     + a0.var(0, ddof=1) / max(len(a0), 1)) + 1e-12
        cnt[np.argsort(-np.abs(a1.mean(0) - a0.mean(0)) / s2)[:args.K]] += 1
    stab = cnt / args.B

    iegset = set(cols[ieg_pos])
    log('  %-14s %14s %s' % ('selector', 'IEG in top-%d' % args.K, 'top 10 genes'))
    for nm, sc in (('effect_size', eff), ('stability', stab), ('connectivity', conn)):
        top = np.argsort(-sc)[:args.K]
        n_ieg = sum(1 for i in top if cols[i] in iegset)
        log('  %-14s %14d  %s' % (nm, n_ieg,
                                  ', '.join(sym.get(cols[i], cols[i]) for i in top[:10])))
        rows.append(dict(test='selector', selector=nm, n_ieg_top=int(n_ieg)))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'gse5406.tsv'), sep='\t', index=False)
    log('\nwrote %s/gse5406.tsv' % OUT)


if __name__ == '__main__':
    main()
