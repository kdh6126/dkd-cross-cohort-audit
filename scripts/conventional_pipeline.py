#!/usr/bin/env python
"""Would the conventional DKD bioinformatics pipeline have promoted the artifact?

This is the experiment that decides whether the method here has novelty or merely re-labels
existing practice. The published work that named FMOD and LUM as DKD hub genes
(PMID 33887734) used WGCNA on a single dataset (GSE99339, 179 microdissected glomeruli),
selecting hub genes by intramodular connectivity. That is the standard recipe in this
literature: differential expression, co-expression modules, hub genes by network degree.

Three pipelines are run on the same cohorts and the same 9,900-gene space, and scored on the
same axes:

  DEG_meta   per-cohort t-statistic, combined across cohorts (Stouffer). The plain approach.
  WGCNA_hub  correlation network raised to a soft-thresholding power, module assignment by
             correlation with the disease trait, hub rank by intramodular connectivity. A
             connectivity-based hub ranking in the manner of WGCNA, NOT a reproduction of it:
             see wgcna_validation.py for the three ways it departs from the canonical pipeline
             and for the 32-variant grid showing the conclusion does not depend on them.
  RBS_corr   this project's confounder-aware cross-cohort selector.

The question is not which gets the higher AUROC. It is which of them fills its top ranks with
the tissue-procurement module - because that module is the most reproducible signal in the data
and any criterion that rewards reproducibility or connectivity without a confounder-matched
contrast should rank it first.
"""
import os, sys, argparse
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM
from dkd_deconfound import IEG

OUT = 'results/conventional_pipeline'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def deg_meta(data, cohorts):
    """Per-cohort two-sample t, combined by Stouffer's method weighted by sqrt(n)."""
    zs, ws = [], []
    for c in cohorts:
        Z, y = data[c][0].values, data[c][1]
        t, p = stats.ttest_ind(Z[y == 1], Z[y == 0], axis=0, equal_var=False)
        z = np.sign(t) * stats.norm.isf(np.clip(p, 1e-300, 1) / 2)
        zs.append(np.nan_to_num(z))
        ws.append(np.sqrt(len(y)))
    Zs = np.vstack(zs)
    w = np.array(ws)[:, None]
    return (w * Zs).sum(0) / np.sqrt((w ** 2).sum())


def wgcna_hub(data, cohorts, power=6, module_frac=0.10):
    """Intramodular connectivity as a hub criterion, in the manner of WGCNA.

    Soft-threshold the correlation matrix, take the module most correlated with the disease
    trait, rank genes inside it by summed adjacency. This is not canonical WGCNA and does not
    claim to be. It departs in three respects:

        power       fixed at 6, not chosen by the scale-free topology criterion. On these data
                    beta=6 gives a scale-free fit of R^2=0.62, below WGCNA's usual 0.80
                    threshold; pickSoftThreshold would have chosen 10 or 12.
        module      top decile by |gene-trait correlation|, not dynamic tree cut on the TOM.
        hub rank    intramodular connectivity only, not module membership (kME).

    wgcna_validation.py rebuilds this across all three axes -- 32 variants -- and finds that
    every one of them keeps the immediate-early module out of the top 50, which is the only
    property the Section 3.4 argument needs.
    """
    X = np.vstack([data[c][0].values for c in cohorts])
    y = np.concatenate([data[c][1] for c in cohorts]).astype(float)
    Xz = (X - X.mean(0)) / (X.std(0, ddof=1) + 1e-12)
    n = len(y)
    yz = (y - y.mean()) / (y.std(ddof=1) + 1e-12)
    gs = (Xz.T @ yz) / (n - 1)                        # gene significance
    k = max(50, int(module_frac * X.shape[1]))
    mod = np.argsort(-np.abs(gs))[:k]                 # the disease-correlated module
    S = Xz[:, mod]
    R = (S.T @ S) / (n - 1)
    A = np.abs(R) ** power                            # soft threshold
    np.fill_diagonal(A, 0.0)
    kIN = A.sum(1)                                    # intramodular connectivity
    score = np.zeros(X.shape[1])
    score[mod] = kIN / kIN.max()
    return score, mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cohorts', default=','.join(GLOM) + ',GSE142025')
    ap.add_argument('--topk', type=int, default=50)
    ap.add_argument('--power', type=int, default=6)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])
    gs_tab = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym = dict(zip(gs_tab['entrez_id'], gs_tab['symbol']))
    sym2id = {s: g for g, s in sym.items() if isinstance(s, str)}
    ieg_ids = set(sym2id[s] for s in IEG if s in sym2id)

    spec = pd.read_csv('results/final_candidates.tsv', sep='\t')
    spec = spec.rename(columns={spec.columns[0]: 'entrez_id'})
    spec['entrez_id'] = spec['entrez_id'].astype(str)
    spec = spec.set_index('entrez_id')

    pipelines = {}
    log('running DEG meta-analysis ...')
    pipelines['DEG_meta'] = np.abs(deg_meta(data, cohorts))
    log('running WGCNA-style hub ranking (soft power %d) ...' % args.power)
    pipelines['WGCNA_hub'], mod = wgcna_hub(data, cohorts, power=args.power)
    log('  module size %d genes' % len(mod))

    v1 = pd.read_csv('results/final_candidates.tsv', sep='\t')
    v1 = v1.rename(columns={v1.columns[0]: 'entrez_id'})
    v1['entrez_id'] = v1['entrez_id'].astype(str).values
    v2 = pd.read_csv('results/candidates_v2/candidates_all.tsv', sep='\t')
    v2 = v2.rename(columns={v2.columns[0]: 'entrez_id'})
    v2['entrez_id'] = v2['entrez_id'].astype(str).values
    rank_to_score = lambda d: (-d.set_index('entrez_id')['mean_rank'].reindex(genes).values)
    pipelines['RBS_uncorrected'] = rank_to_score(v1)
    pipelines['RBS_corrected'] = rank_to_score(v2)

    rows = []
    log('\n%-18s %6s %14s %16s %s'
        % ('pipeline', 'top-K', 'IEG in top-K', 'procurement %', 'specificity pass %'))
    for name, score in pipelines.items():
        order = np.argsort(-np.nan_to_num(score, nan=-np.inf))
        for K in (10, 20, 50, 100):
            top = genes[order[:K]]
            n_ieg = sum(1 for g in top if g in ieg_ids)
            pr = float(spec['procurement_flag'].reindex(top).fillna(False).mean())
            sp = float(spec['passes_specificity'].reindex(top).fillna(False).mean())
            rows.append(dict(pipeline=name, K=K, n_ieg=n_ieg,
                             procurement_frac=pr, spec_pass_frac=sp))
            if K == args.topk:
                log('%-18s %6d %14d %15.0f%% %17.0f%%' % (name, K, n_ieg, 100 * pr, 100 * sp))
        if name in ('DEG_meta', 'WGCNA_hub'):
            log('    top 15: %s' % ', '.join(sym.get(g, g) for g in genes[order[:15]]))

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, 'comparison.tsv'), sep='\t', index=False)

    log('\n=== IEG genes in the top 50, by pipeline ===')
    for name, score in pipelines.items():
        order = np.argsort(-np.nan_to_num(score, nan=-np.inf))
        top = [sym.get(g, g) for g in genes[order[:50]] if g in ieg_ids]
        log('  %-18s %d  %s' % (name, len(top), ', '.join(top)))

    log('\n=== where the published hub genes rank in each pipeline ===')
    published = ['FMOD', 'LUM', 'MMP2', 'COL1A2', 'THBS2', 'VCAN']
    log('  %-9s %s' % ('gene', ' '.join('%16s' % n for n in pipelines)))
    for s in published:
        g = sym2id.get(s)
        if g is None:
            continue
        i = int(np.where(genes == g)[0][0])
        cells = []
        for name, score in pipelines.items():
            order = np.argsort(-np.nan_to_num(score, nan=-np.inf))
            cells.append('%16d' % (int(np.where(order == i)[0][0]) + 1))
        log('  %-9s %s' % (s, ' '.join(cells)))

    # ---------------------------------------------------------------- LODO comparison
    # The specificity axis above favours the conventional pipelines. The axis this project
    # actually claimed was cross-cohort reproducibility, so it has to be measured for them too,
    # under the same leave-one-dataset-out protocol.
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    log('\n=== LODO: cross-fold reproducibility and external AUROC ===')
    lodo = {}
    for name in ('DEG_meta', 'WGCNA_hub'):
        orders, aucs = {}, {K: [] for K in (10, 20, 50)}
        for held in cohorts:
            tr = [c for c in cohorts if c != held]
            sub = {c: data[c] for c in tr}
            sc = (np.abs(deg_meta(sub, tr)) if name == 'DEG_meta'
                  else wgcna_hub(sub, tr, power=args.power)[0])
            o = np.argsort(-np.nan_to_num(sc, nan=-np.inf))
            orders[held] = o
            Xtr = np.vstack([data[c][0].values for c in tr])
            ytr = np.concatenate([data[c][1] for c in tr])
            Xte, yte = data[held][0].values, data[held][1]
            for K in aucs:
                m = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000,
                                       class_weight='balanced').fit(Xtr[:, o[:K]], ytr)
                aucs[K].append(roc_auc_score(yte, m.predict_proba(Xte[:, o[:K]])[:, 1]))
        lodo[name] = (orders, aucs)

    log('  %-16s %5s %20s %s' % ('pipeline', 'K', 'cross-fold Jaccard', 'external AUROC'))
    xrows = []
    for name, (orders, aucs) in lodo.items():
        for K in (10, 20, 50):
            sets = [set(orders[h][:K]) for h in cohorts]
            js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
            log('  %-16s %5d %20.3f %14.3f'
                % (name, K, float(np.mean(js)), float(np.mean(aucs[K]))))
            xrows.append(dict(pipeline=name, K=K, cross_fold_jaccard=float(np.mean(js)),
                              external_auroc=float(np.mean(aucs[K]))))
    dec = pd.read_csv('results/deconfound.tsv', sep='\t')
    ag = dec.groupby(['arm', 'K'])[['cross_fold_jaccard', 'external_auroc']].mean()
    for arm in ('none', 'orth'):
        for K in (10, 20, 50):
            log('  %-16s %5d %20.3f %14.3f'
                % ('RBS_' + arm, K, ag.loc[(arm, K), 'cross_fold_jaccard'],
                   ag.loc[(arm, K), 'external_auroc']))
            xrows.append(dict(pipeline='RBS_' + arm, K=K,
                              cross_fold_jaccard=float(ag.loc[(arm, K), 'cross_fold_jaccard']),
                              external_auroc=float(ag.loc[(arm, K), 'external_auroc'])))
    pd.DataFrame(xrows).to_csv(os.path.join(OUT, 'lodo_comparison.tsv'), sep='\t', index=False)

    log('\nwrote %s/comparison.tsv and lodo_comparison.tsv' % OUT)


if __name__ == '__main__':
    main()
