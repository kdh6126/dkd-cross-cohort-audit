#!/usr/bin/env python
"""Does the Section 3.4 conclusion depend on how faithfully WGCNA was reimplemented?

Section 3.4 argues that a co-expression connectivity criterion misses the procurement
confounder not because it is robust but because that confounder forms a small block, and
intramodular connectivity is a *sum* over the block. The comparator used there
(`conventional_pipeline.wgcna_hub`) is a reimplementation, and it departs from the canonical
WGCNA pipeline in three ways that a reviewer will notice:

    soft power      fixed at beta = 6 rather than chosen by the scale-free topology criterion.
    module          top decile by |gene-trait correlation|, rather than hierarchical clustering
                    on TOM dissimilarity followed by a dynamic tree cut.
    hub rank        intramodular connectivity only. Canonical practice often ranks hubs by kME,
                    the correlation with the module eigengene, instead.

If the conclusion holds only under our particular shortcuts, it is an artifact of the
reimplementation and Section 3.4 collapses. This script therefore rebuilds the comparator four
ways along each axis and asks the same question of every variant: how many immediate-early
genes reach the top 50, and where do they rank?

Nothing here tries to argue the reimplementation *is* canonical WGCNA. The claim under test is
narrower and is the only one the manuscript needs: that the ranking behaviour Section 3.4
describes is a property of connectivity-style hub criteria in general, not of our shortcuts.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM
from dkd_deconfound import IEG

OUT = 'results/wgcna_validation'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ----------------------------------------------------------------- network construction

def adjacency(R, power, signed):
    """Canonical WGCNA adjacency. Unsigned |r|^b; signed ((1+r)/2)^b."""
    A = (((1.0 + R) / 2.0) ** power) if signed else (np.abs(R) ** power)
    np.fill_diagonal(A, 0.0)
    return A


def tom(A):
    """Topological overlap matrix, Zhang & Horvath's definition.

        TOM_ij = (sum_u a_iu a_uj + a_ij) / (min(k_i, k_j) + 1 - a_ij)
    """
    k = A.sum(1)
    num = A @ A + A
    denom = np.minimum.outer(k, k) + 1.0 - A
    T = num / np.maximum(denom, 1e-12)
    np.fill_diagonal(T, 1.0)
    return T


def scale_free_fit(A, nbins=20):
    """R^2 of log10(p(k)) on log10(k), the criterion `pickSoftThreshold` maximises."""
    k = A.sum(1)
    k = k[k > 0]
    if len(k) < nbins * 2:
        return np.nan
    counts, edges = np.histogram(k, bins=nbins)
    centres = (edges[:-1] + edges[1:]) / 2.0
    keep = counts > 0
    if keep.sum() < 4:
        return np.nan
    x = np.log10(centres[keep])
    yv = np.log10(counts[keep] / counts.sum())
    r = np.corrcoef(x, yv)[0, 1]
    return r ** 2


# ----------------------------------------------------------------- module definition

def module_by_gs(gs, frac):
    """Our shortcut: the top decile by |gene-trait correlation|."""
    k = max(50, int(frac * len(gs)))
    return np.argsort(-np.abs(gs))[:k]


def module_by_tom(Xz, gs, frac, power, signed, seed=0):
    """Closer to canonical: cluster on TOM dissimilarity, keep the most trait-correlated module.

    Dynamic tree cut is replaced by average-linkage clustering with a flat cut, which is the
    part of WGCNA hardest to reproduce faithfully outside R. The substitution is stated in the
    output rather than hidden: what matters for the test is that modules come from network
    topology instead of from the trait, which is the actual difference of interest.
    """
    n_gene = Xz.shape[1]
    # TOM on the full 9,900-gene space is 780 MB and slow; restrict to a variance-ranked
    # subspace the way WGCNA's blockwise modules do, keeping the trait-correlated genes.
    pool = np.argsort(-np.abs(gs))[:max(2000, int(frac * n_gene) * 4)]
    S = Xz[:, pool]
    R = np.corrcoef(S, rowvar=False)
    R = np.nan_to_num(R)
    A = adjacency(R, power, signed)
    T = tom(A)
    D = 1.0 - T
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2.0
    Z = linkage(squareform(D, checks=False), method='average')
    target = max(50, int(frac * n_gene))
    # take the cluster whose mean |gs| is highest among clusters large enough to be a module
    for n_clust in (4, 6, 8, 12, 20, 30):
        lab = fcluster(Z, n_clust, criterion='maxclust')
        cand = [(np.abs(gs[pool[lab == c]]).mean(), c, (lab == c).sum())
                for c in np.unique(lab) if (lab == c).sum() >= target // 2]
        if cand:
            cand.sort(reverse=True)
            _, best, size = cand[0]
            return pool[lab == best], n_clust, int(size)
    return pool[:target], -1, int(target)


# ----------------------------------------------------------------- hub ranking

def rank_kin(Xz, mod, power, signed):
    """Intramodular connectivity: the summed adjacency inside the module."""
    S = Xz[:, mod]
    R = np.nan_to_num(np.corrcoef(S, rowvar=False))
    A = adjacency(R, power, signed)
    return A.sum(1)


def rank_kme(Xz, mod):
    """Module membership: correlation of each gene with the first module eigengene."""
    S = Xz[:, mod]
    S = S - S.mean(0)
    # first principal component of the module = WGCNA's module eigengene
    _, _, Vt = np.linalg.svd(S, full_matrices=False)
    eg = S @ Vt[0]
    eg = (eg - eg.mean()) / (eg.std(ddof=1) + 1e-12)
    Sz = S / (S.std(0, ddof=1) + 1e-12)
    return np.abs((Sz.T @ eg) / (len(eg) - 1))


# ----------------------------------------------------------------- driver

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cohorts', default=','.join(GLOM) + ',GSE142025')
    ap.add_argument('--topk', type=int, default=50)
    ap.add_argument('--frac', type=float, default=0.10)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])

    gs_tab = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym = dict(zip(gs_tab['entrez_id'], gs_tab['symbol']))
    sym2id = {s: g for g, s in sym.items() if isinstance(s, str)}
    ieg_ids = set(sym2id[s] for s in IEG if s in sym2id)
    ieg_mask = np.array([g in ieg_ids for g in genes])
    log('immediate-early genes present in the gene space: %d of %d' % (ieg_mask.sum(), len(IEG)))

    X = np.vstack([data[c][0].values for c in cohorts])
    y = np.concatenate([data[c][1] for c in cohorts]).astype(float)
    Xz = (X - X.mean(0)) / (X.std(0, ddof=1) + 1e-12)
    n = len(y)
    yz = (y - y.mean()) / (y.std(ddof=1) + 1e-12)
    gs = (Xz.T @ yz) / (n - 1)

    # ---- 1. is beta = 6 defensible on the scale-free criterion?
    log('')
    log('=== soft-threshold power, scale-free topology fit ===')
    mod6 = module_by_gs(gs, args.frac)
    R6 = np.nan_to_num(np.corrcoef(Xz[:, mod6], rowvar=False))
    sft = []
    for b in (2, 4, 6, 8, 10, 12, 14):
        r2 = scale_free_fit(adjacency(R6, b, signed=False))
        sft.append(dict(power=b, scale_free_r2=r2))
        log('  beta = %-3d  R^2 = %.3f%s' % (b, r2, '   <- used' if b == 6 else ''))
    pd.DataFrame(sft).to_csv(os.path.join(OUT, 'soft_threshold.tsv'), sep='\t', index=False)

    # ---- 2. the grid
    log('')
    log('=== does the conclusion survive the variants? ===')
    log('  %-14s %-7s %-8s %-6s %6s %10s %10s'
        % ('module', 'rank', 'adjacency', 'beta', 'IEG@50', 'medianRank', 'moduleSize'))
    rows = []
    for mod_def in ('gs_top', 'tom_cluster'):
        for power in (4, 6, 8, 12):
            for signed in (False, True):
                if mod_def == 'gs_top':
                    mod = module_by_gs(gs, args.frac)
                    info = 'top-decile |gs|'
                else:
                    mod, nclust, _ = module_by_tom(Xz, gs, args.frac, power, signed)
                    info = 'TOM cluster (k=%d)' % nclust
                for rank_by in ('kIN', 'kME'):
                    v = (rank_kin(Xz, mod, power, signed) if rank_by == 'kIN'
                         else rank_kme(Xz, mod))
                    order = mod[np.argsort(-v)]
                    top = order[:args.topk]
                    n_ieg_top = int(ieg_mask[top].sum())
                    in_mod = ieg_mask[mod]
                    if in_mod.any():
                        ranks = np.argsort(np.argsort(-v))[in_mod] + 1
                        med = float(np.median(ranks))
                    else:
                        med = float('nan')
                    rows.append(dict(module=mod_def, module_detail=info, rank_by=rank_by,
                                     adjacency='signed' if signed else 'unsigned', power=power,
                                     module_size=len(mod), ieg_in_module=int(in_mod.sum()),
                                     ieg_in_top50=n_ieg_top, ieg_median_rank=med))
                    log('  %-14s %-7s %-8s %-6d %6d %10s %10d'
                        % (mod_def, rank_by, 'signed' if signed else 'unsigned', power,
                           n_ieg_top, ('%.0f' % med) if med == med else 'n/a', len(mod)))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(OUT, 'variants.tsv'), sep='\t', index=False)

    # ---- 3. verdict
    log('')
    log('=== verdict ===')
    tot = len(d)
    zero = int((d['ieg_in_top50'] == 0).sum())
    log('  variants tested                        : %d' % tot)
    log('  variants placing NO immediate-early    : %d (%.0f%%)' % (zero, 100 * zero / tot))
    log('  worst case (most IEG in top 50)        : %d' % d['ieg_in_top50'].max())
    med_all = d['ieg_median_rank'].dropna()
    if len(med_all):
        log('  median IEG rank within module, range   : %.0f to %.0f'
            % (med_all.min(), med_all.max()))
    if zero == tot:
        log('')
        log('  Every variant misses the confounder. Section 3.4 does not depend on the')
        log('  reimplementation shortcuts; it is a property of connectivity-style ranking.')
    elif d['ieg_in_top50'].max() <= 2:
        log('')
        log('  No variant promotes the confounder to any meaningful degree. The conclusion')
        log('  holds, with the qualification that it is near-total rather than total.')
    else:
        log('')
        log('  AT LEAST ONE VARIANT PROMOTES THE CONFOUNDER. Section 3.4 as written is too')
        log('  strong and must be qualified by module definition, not stated in general.')

    log('')
    log('  wrote %s/{soft_threshold,variants}.tsv' % OUT)


if __name__ == '__main__':
    main()
