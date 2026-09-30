#!/usr/bin/env python
"""Re-derive the candidate list with the confounder-aware selector.

The existing 30 candidates were selected WITHOUT the correction and then filtered by a
post-hoc specificity gate. That ordering has a defect: the selector spends its capacity on the
procurement axis, and the gate can only throw things away afterwards - it cannot promote genes
the selector never ranked. Selecting on residualised data changes which genes are ranked at all,
so the candidate set should be rebuilt from scratch rather than re-filtered.

Emits the same column layout as final_candidates_gated.tsv so the downstream scripts
(per-diagnosis profile, PubMed review, KPMP corroboration) run against it unchanged.
"""
import os, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM
from dkd_rbs import rbs_fit
from dkd_deconfound import handling_score, residualise, IEG


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--arm', default='orth_within',
                    choices=['none', 'orth', 'orth_within', 'orth_within_0.5'])
    ap.add_argument('--B', type=int, default=200)
    ap.add_argument('--k', type=int, default=50)
    ap.add_argument('--base', default='relieff')
    ap.add_argument('--topn', type=int, default=30)
    ap.add_argument('--cohorts', default=','.join(GLOM) + ',GSE142025')
    ap.add_argument('--seed', type=int, default=20260824)
    ap.add_argument('--out', default='results/candidates_v2')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])
    gpos = {g: i for i, g in enumerate(genes)}

    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym = dict(zip(gs['entrez_id'], gs['symbol']))
    sym2id = {s: g for g, s in sym.items() if isinstance(s, str)}
    ieg_cols = [gpos[sym2id[s]] for s in IEG if s in sym2id and sym2id[s] in gpos]

    log('arm=%s  cohorts=%s  B=%d  k=%d' % (args.arm, cohorts, args.B, args.k))

    ranks, comps = {}, {}
    for held in cohorts:
        tr = [c for c in cohorts if c != held]
        train = {}
        for c in tr:
            Z = data[c][0].values.copy()
            s = handling_score(Z, ieg_cols)
            yc = data[c][1]
            if args.arm == 'orth':
                Z = residualise(Z, s)
            elif args.arm == 'orth_within':
                Z = residualise(Z, s, y=yc)
            elif args.arm == 'orth_within_0.5':
                Z = residualise(Z, s, y=yc, lam=0.5)
            train[c] = (Z, yc)
        res = rbs_fit(train, k=args.k, B=args.B, base=args.base, seed=args.seed, agg='geomean')
        order = np.argsort(-res['score'])
        r = np.empty(len(genes), float)
        r[order] = np.arange(1, len(genes) + 1)
        ranks[held] = r
        comps[held] = res
        log('  hold out %-11s top10: %s'
            % (held, ', '.join(sym.get(g, g) for g in genes[order[:10]])))

    R = pd.DataFrame(ranks, index=genes)
    S = pd.DataFrame({h: comps[h]['S'] for h in comps}, index=genes)
    Rr = pd.DataFrame({h: comps[h]['R'] for h in comps}, index=genes)
    P = pd.DataFrame({h: comps[h]['P'] for h in comps}, index=genes)

    spec = pd.read_csv('results/final_candidates.tsv', sep='\t')
    spec = spec.rename(columns={spec.columns[0]: 'entrez_id'})
    spec['entrez_id'] = spec['entrez_id'].astype(str)
    spec = spec.set_index('entrez_id')

    df = pd.DataFrame({
        'symbol': [sym.get(g, '') for g in genes],
        'mean_rank': R.mean(axis=1).values,
        'worst_rank': R.max(axis=1).values,
        'mean_RBS': np.nan,
        'mean_S': S.mean(axis=1).values,
        'mean_R': Rr.mean(axis=1).values,
        'mean_P': P.mean(axis=1).values,
    }, index=genes)
    for c in ('g_DKD_vs_control', 'g_DKD_vs_otherCKD', 'g_otherCKD_vs_control',
              'passes_specificity', 'procurement_flag'):
        df[c] = spec[c].reindex(df.index).values
    df.index.name = 'entrez_id'
    df = df.sort_values('mean_rank')
    df.to_csv(os.path.join(args.out, 'candidates_all.tsv'), sep='\t')

    top200 = df.head(200)
    log('\n=== of the top 200 by the %s selector ===' % args.arm)
    log('  pass DKD-specificity gate    : %d (%.0f%%)'
        % (top200['passes_specificity'].sum(), 100 * top200['passes_specificity'].mean()))
    log('  flagged as procurement-driven: %d (%.0f%%)'
        % (top200['procurement_flag'].sum(), 100 * top200['procurement_flag'].mean()))

    gated = df[(df['passes_specificity'] == True) & (df['procurement_flag'] != True)].head(args.topn)
    gated.to_csv(os.path.join(args.out, 'candidates_gated.tsv'), sep='\t')

    log('\n=== FINAL %d CANDIDATES (%s selector, specificity-gated) ===' % (len(gated), args.arm))
    log('  %-4s %-10s %-8s %-9s %-10s %-9s %s'
        % ('#', 'symbol', 'rank', 'vs ctrl', 'vs othCKD', 'stability', 'perturb'))
    for i, (gid, r) in enumerate(gated.iterrows(), 1):
        log('  %-4d %-10s %-8.0f %+9.2f %+10.2f %9.2f %.2f'
            % (i, r['symbol'], r['mean_rank'], r['g_DKD_vs_control'],
               r['g_DKD_vs_otherCKD'], r['mean_S'], r['mean_P']))

    # ---------------------------------------------------------------- compare with v1
    old = pd.read_csv('results/final_candidates_gated.tsv', sep='\t')
    old_syms = set(old['symbol'].dropna().astype(str))
    new_syms = set(gated['symbol'].dropna().astype(str))
    log('\n=== v1 (uncorrected) vs v2 (%s) ===' % args.arm)
    log('  shared        : %2d  %s' % (len(old_syms & new_syms),
                                       ', '.join(sorted(old_syms & new_syms))))
    log('  dropped from v1: %2d  %s' % (len(old_syms - new_syms),
                                        ', '.join(sorted(old_syms - new_syms))))
    log('  new in v2      : %2d  %s' % (len(new_syms - old_syms),
                                        ', '.join(sorted(new_syms - old_syms))))

    # where did the old top-ranked artifact genes go?
    log('\n=== rank of the immediate-early module under the two selectors ===')
    oldall = pd.read_csv('results/final_candidates.tsv', sep='\t')
    oldall = oldall.rename(columns={oldall.columns[0]: 'entrez_id'})
    oldall['entrez_id'] = oldall['entrez_id'].astype(str)
    oldall = oldall.set_index('entrez_id')
    log('  %-8s %12s %12s' % ('gene', 'v1 rank', 'v2 rank'))
    for s_ in IEG:
        gid = sym2id.get(s_)
        if gid in df.index and gid in oldall.index:
            log('  %-8s %12.0f %12.0f'
                % (s_, oldall.loc[gid, 'mean_rank'], df.loc[gid, 'mean_rank']))
    log('\nwrote %s/' % args.out)


if __name__ == '__main__':
    main()
