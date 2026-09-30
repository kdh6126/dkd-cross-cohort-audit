#!/usr/bin/env python
"""Single-cell analysis from the KPMP h5ad - what the Atlas API could not give.

The API returns cell counts and nothing else for the `dmr` (diabetic) slice: every
foldChange / pVal / pValAdj is null. The downloaded file has the actual data and much richer
metadata, and it turns out to support three contrasts the API cannot:

  donors  cells    group
      27  86,430   DKD                (disease_category CKD + diabetes_history Yes)
      10  34,551   CKD non-diabetic   -> the CONFOUNDER-MATCHED contrast, independent of ERCB
      12  40,386   Healthy reference
       6  43,809   Tumor nephrectomy  -> the PROCUREMENT axis, inside one resource

So:
  1. DKD vs healthy reference, per cell type
  2. DKD vs non-diabetic CKD, per cell type  - the specificity gate replicated in a third
     resource at single-cell resolution
  3. Healthy reference vs tumour nephrectomy for the immediate-early module - the procurement
     axis measured between two groups that are *both* controls

Statistics are computed on **donor-level pseudobulk**, not per cell. Per-cell tests treat
thousands of cells from one patient as independent and produce p-values that are essentially
a function of cell count; with 27 versus 10 donors the donor is the unit of replication.

Extraction is a single chunked pass over the CSR matrix, keeping only the target columns -
reading all 669M nonzeros to build a dense 305k x ~200 block.
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd
import h5py
from scipy import stats

H5 = 'data/raw/kpmp/KPMP_v1.5_snRNA_human_kidney.h5ad'
OUT = 'results/kpmp_singlecell'
IEG = ['FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'ATF3', 'DUSP1', 'ZFP36', 'NR4A2', 'BTG2',
       'EGR2', 'EGR3', 'IER2', 'KLF2', 'KLF4', 'SOCS3', 'DUSP2', 'JUND', 'NR4A1']
N_RANDOM = 150


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    n1, n0 = len(a), len(b)
    if n1 < 3 or n0 < 3:
        return np.nan, np.nan
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / (n1 + n0 - 2))
    if sp == 0:
        return np.nan, np.nan
    g = (a.mean() - b.mean()) / sp * (1 - 3 / (4 * (n1 + n0) - 9))
    p = stats.mannwhitneyu(a, b, alternative='two-sided').pvalue
    return g, p


def read_obs(f):
    def cat(k):
        g = f['obs/' + k]
        c = np.array([x.decode() if isinstance(x, bytes) else x for x in g['categories'][:]])
        return c[g['codes'][:]]
    return pd.DataFrame({'donor': cat('donor_id'), 'dcat': cat('disease_category'),
                         'dm': cat('diabetes_history'), 'ct': cat('subclass.l1'),
                         'ct2': cat('cell_type'), 'egfr': cat('eGFR'),
                         'region': cat('region')})


def group_of(dcat, dm):
    if dcat == 'CKD':
        return 'DKD' if dm == 'Yes' else ('CKD_nondiabetic' if dm == 'No' else 'CKD_unknown')
    return {'Healthy_reference_tissue': 'reference', 'Tumor_nephrectomy': 'nephrectomy',
            'AKI': 'AKI'}.get(dcat, dcat)


def extract(f, want_idx, chunk=20000):
    """One sequential pass over CSR, keeping only the wanted gene columns."""
    n_cells, n_genes = f['X'].attrs['shape']
    lut = np.full(n_genes, -1, dtype=np.int32)
    lut[want_idx] = np.arange(len(want_idx), dtype=np.int32)
    out = np.zeros((n_cells, len(want_idx)), dtype=np.float32)
    indptr = f['X/indptr'][:]
    data_ds, idx_ds = f['X/data'], f['X/indices']
    done = 0
    for a in range(0, n_cells, chunk):
        b = min(a + chunk, n_cells)
        lo, hi = int(indptr[a]), int(indptr[b])
        idx = idx_ds[lo:hi]
        col = lut[idx]
        keep = col >= 0
        if keep.any():
            dat = data_ds[lo:hi][keep]
            rows = np.repeat(np.arange(a, b), np.diff(indptr[a:b + 1]))[keep]
            out[rows, col[keep]] = dat
        done += 1
        if done % 4 == 0:
            log('    ... %d / %d cells' % (b, n_cells))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-cells', type=int, default=25,
                    help='minimum cells per (donor, cell type) to form a pseudobulk value')
    ap.add_argument('--min-donors', type=int, default=5)
    ap.add_argument('--seed', type=int, default=20260824)
    ap.add_argument('--candidates', default='results/final_candidates_gated.tsv')
    ap.add_argument('--out', default='results/kpmp_singlecell')
    args = ap.parse_args()
    OUT_D = args.out
    os.makedirs(OUT_D, exist_ok=True)

    cand = pd.read_csv(args.candidates, sep='\t')
    candidates = [g for g in cand['symbol'].dropna().astype(str) if g]

    f = h5py.File(H5, 'r')
    var = f['var/feature_name']
    names = np.array([x.decode() if isinstance(x, bytes) else x
                      for x in var['categories'][:]])[var['codes'][:]]
    pos = {}
    for i, n in enumerate(names):
        pos.setdefault(n, i)

    rng = np.random.default_rng(args.seed)
    targets = [g for g in candidates if g in pos] + [g for g in IEG if g in pos]
    pool = [n for n in pos if n not in set(targets)]
    randoms = list(rng.choice(pool, size=min(N_RANDOM, len(pool)), replace=False))
    targets += randoms
    idx = np.array([pos[g] for g in targets])
    log('extracting %d genes (%d candidates, %d immediate-early, %d random context)'
        % (len(targets), sum(g in pos for g in candidates), sum(g in pos for g in IEG), len(randoms)))

    cache = os.path.join(OUT_D, 'expr_cache.npz')
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        X, targets = z['X'], list(z['genes'])
        log('  loaded cached matrix %s' % (X.shape,))
    else:
        X = extract(f, idx)
        np.savez_compressed(cache, X=X, genes=np.array(targets, dtype=object))
        log('  extracted %s, cached' % (X.shape,))

    obs = read_obs(f)
    obs['group'] = [group_of(a, b) for a, b in zip(obs['dcat'], obs['dm'])]
    log('\ncells per group: %s' % obs['group'].value_counts().to_dict())

    # ---------------------------------------------------------------- pseudobulk
    log('\nbuilding donor x cell-type pseudobulk (mean of normalised expression) ...')
    key = obs['donor'] + '||' + obs['ct']
    dfX = pd.DataFrame(X, columns=targets)
    dfX['__k'] = key.values
    counts = dfX.groupby('__k').size()
    pb = dfX.groupby('__k').mean()
    pb = pb.loc[counts[counts >= args.min_cells].index]
    meta = pd.DataFrame({'donor': [k.split('||')[0] for k in pb.index],
                         'ct': [k.split('||')[1] for k in pb.index]}, index=pb.index)
    dgroup = obs.groupby('donor')['group'].first()
    meta['group'] = meta['donor'].map(dgroup)
    log('  %d (donor, cell type) pseudobulk profiles with >=%d cells' % (len(pb), args.min_cells))

    rows = []
    contrasts = [('DKD_vs_reference', 'DKD', 'reference'),
                 ('DKD_vs_CKDnondiabetic', 'DKD', 'CKD_nondiabetic'),
                 ('reference_vs_nephrectomy', 'reference', 'nephrectomy')]
    cts = [c for c in meta['ct'].unique()
           if (meta['ct'] == c).sum() >= 2 * args.min_donors]
    for cname, ga, gb in contrasts:
        for ct in cts:
            m = meta[(meta['ct'] == ct)]
            A = pb.loc[m.index[m['group'] == ga]]
            B = pb.loc[m.index[m['group'] == gb]]
            if len(A) < args.min_donors or len(B) < args.min_donors:
                continue
            for g in targets:
                gg, pp = hedges_g(A[g].values, B[g].values)
                rows.append(dict(contrast=cname, cell_type=ct, gene=g, g=gg, p=pp,
                                 n_a=len(A), n_b=len(B),
                                 kind=('candidate' if g in candidates else
                                       'IEG' if g in IEG else 'random')))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT_D, 'pseudobulk_de.tsv'), sep='\t', index=False)

    # ---------------------------------------------------------------- report
    for cname, _, _ in contrasts:
        sub = res[res['contrast'] == cname]
        if sub.empty:
            continue
        log('\n=== %s ===' % cname)
        for kind in ('IEG', 'candidate', 'random'):
            s = sub[sub['kind'] == kind]['g'].dropna()
            log('  %-10s n=%-5d mean g = %+.3f   |g|>0.8 in %.0f%% of (gene, cell type) pairs'
                % (kind, len(s), s.mean(), 100 * (s.abs() > 0.8).mean()))
        rnd = sub[sub['kind'] == 'random']['g'].dropna()
        ieg = sub[sub['kind'] == 'IEG']['g'].dropna()
        if len(rnd) > 10 and len(ieg) > 10:
            log('  IEG vs random: Mann-Whitney p = %.2e'
                % stats.mannwhitneyu(ieg, rnd, alternative='two-sided').pvalue)
        top = sub[(sub['kind'] == 'candidate') & sub['p'].notna()].nsmallest(8, 'p')
        for _, r in top.iterrows():
            log('     %-9s %-6s g=%+.2f p=%.3g  (%d vs %d donors)'
                % (r['gene'], r['cell_type'], r['g'], r['p'], r['n_a'], r['n_b']))

    # candidate summary across cell types
    log('\n=== candidates: strongest cell type per contrast ===')
    for g in candidates:
        s = res[(res['gene'] == g) & res['g'].notna()]
        if s.empty:
            continue
        line = []
        for cname, _, _ in contrasts:
            t = s[s['contrast'] == cname]
            if t.empty:
                line.append('%-22s' % '-'); continue
            b = t.loc[t['g'].abs().idxmax()]
            line.append('%-8s %+5.2f%s' % (b['cell_type'][:8], b['g'],
                                           '*' if b['p'] < 0.05 else ' '))
        log('  %-9s | %s' % (g, ' | '.join(line)))
    log('\n  columns: DKD vs reference | DKD vs non-diabetic CKD | reference vs nephrectomy')
    log('  * = p < 0.05 at donor level, uncorrected')
    log('\nwrote %s/' % OUT)


if __name__ == '__main__':
    main()
