#!/usr/bin/env python
"""Final candidate list: RBS ranking gated by DKD specificity.

Why the gate is a pipeline STAGE and not a footnote
---------------------------------------------------
The immediate-early module (FOS/JUN/EGR1/ATF3/DUSP1/ZFP36/...) is the most cross-cohort
reproducible thing in this data, and it is an artifact: non-diabetic CKD biopsies show
145% of the DKD shift (scripts/ieg_artifact_test.py). Cases are needle biopsies and controls
are nephrectomy/living-donor tissue in EVERY cohort, so the procurement confounder is shared
by every cohort and therefore looks perfectly reproducible to any cross-cohort stability
criterion. Reproducibility across cohorts does not imply biological validity when the
confounder is shared - the only thing that separates them here is a contrast that holds
procurement constant.

ERCB's 310 non-diabetic CKD biopsies provide exactly that contrast: they were procured like
the DKD biopsies, so DKD-vs-otherCKD is confounder-matched. A gene passes only if it still
separates DKD from other CKD, in the same direction as it separates DKD from controls.

Output columns carry the full evidence chain so nothing is taken on trust.
"""
import os, sys, glob, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_cohort

ERCB = ['GSE104948', 'GSE104954']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    n1, n0 = len(a), len(b)
    if n1 < 2 or n0 < 2:
        return np.full(a.shape[1], np.nan)
    sp = np.sqrt(((n1 - 1) * a.var(0, ddof=1) + (n0 - 1) * b.var(0, ddof=1)) / (n1 + n0 - 2))
    sp[sp == 0] = np.nan
    J = 1 - 3 / (4 * (n1 + n0) - 9)
    return (a.mean(0) - b.mean(0)) / sp * J


def genome_wide_specificity():
    """Hedges g for DKD-vs-control and DKD-vs-otherCKD, all 9,900 genes, both ERCB compartments."""
    out = {}
    for c in ERCB:
        X, y, p = load_cohort(c, labelled_only=False)
        X.columns = X.columns.astype(str)
        p = p.reset_index(drop=True)
        M = X.values
        dkd = M[(p['label'] == 1).values]
        ctl = M[(p['label'] == 0).values]
        oth = M[(p['label'] == -1).values]
        out[c] = pd.DataFrame({'vs_control': hedges_g(dkd, ctl),
                               'vs_otherCKD': hedges_g(dkd, oth),
                               'otherCKD_vs_control': hedges_g(oth, ctl)},
                              index=X.columns)
        log('  %-11s DKD=%d control=%d other-CKD=%d' % (c, len(dkd), len(ctl), len(oth)))
    return out


def consensus(components_dir):
    files = sorted(glob.glob(os.path.join(components_dir, 'components_holdout_*.tsv')))
    ranks, comp = {}, {}
    for f in files:
        held = os.path.basename(f)[len('components_holdout_'):-len('.tsv')]
        d = pd.read_csv(f, sep='\t', dtype={'entrez_id': str}).set_index('entrez_id')
        ranks[held] = d['RBS'].rank(ascending=False)
        comp[held] = d
    R = pd.DataFrame(ranks)
    any_d = list(comp.values())[0]
    out = pd.DataFrame({'symbol': any_d['symbol'],
                        'mean_rank': R.mean(axis=1), 'worst_rank': R.max(axis=1),
                        'mean_RBS': pd.DataFrame({h: comp[h]['RBS'] for h in comp}).mean(axis=1),
                        'mean_S': pd.DataFrame({h: comp[h]['S'] for h in comp}).mean(axis=1),
                        'mean_R': pd.DataFrame({h: comp[h]['R'] for h in comp}).mean(axis=1),
                        'mean_P': pd.DataFrame({h: comp[h]['P'] for h in comp}).mean(axis=1)})
    return out.sort_values('mean_rank')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--components', default='results/proposed_all4')
    ap.add_argument('--min-spec', type=float, default=0.5)
    ap.add_argument('--topn', type=int, default=30)
    ap.add_argument('--out', default='results/final_candidates.tsv')
    args = ap.parse_args()

    log('=== specificity, genome-wide ===')
    spec = genome_wide_specificity()
    cons = consensus(args.components)

    vs_ctrl = pd.concat([spec[c]['vs_control'] for c in ERCB], axis=1).mean(axis=1)
    vs_oth = pd.concat([spec[c]['vs_otherCKD'] for c in ERCB], axis=1).mean(axis=1)
    oth_ctrl = pd.concat([spec[c]['otherCKD_vs_control'] for c in ERCB], axis=1).mean(axis=1)

    df = cons.join(pd.DataFrame({'g_DKD_vs_control': vs_ctrl,
                                 'g_DKD_vs_otherCKD': vs_oth,
                                 'g_otherCKD_vs_control': oth_ctrl}))
    df['passes_specificity'] = ((df['g_DKD_vs_otherCKD'].abs() >= args.min_spec) &
                                (np.sign(df['g_DKD_vs_otherCKD']) == np.sign(df['g_DKD_vs_control'])))
    # a gene whose control-contrast is driven mostly by the shared procurement effect:
    # other-CKD moves at least as far from control as DKD does, in the same direction
    df['procurement_flag'] = ((np.sign(df['g_otherCKD_vs_control']) == np.sign(df['g_DKD_vs_control'])) &
                              (df['g_otherCKD_vs_control'].abs() >= 0.8 * df['g_DKD_vs_control'].abs()))

    df = df.sort_values('mean_rank')
    df.to_csv(args.out, sep='\t')

    n_top = 200
    top = df.head(n_top)
    log('\n=== of the top %d by RBS ===' % n_top)
    log('  pass DKD-specificity gate : %d (%.0f%%)'
        % (top['passes_specificity'].sum(), 100 * top['passes_specificity'].mean()))
    log('  flagged as procurement-driven: %d (%.0f%%)'
        % (top['procurement_flag'].sum(), 100 * top['procurement_flag'].mean()))

    final = df[df['passes_specificity'] & ~df['procurement_flag']].head(args.topn)
    final.to_csv(args.out.replace('.tsv', '_gated.tsv'), sep='\t')

    log('\n=== FINAL %d CANDIDATES (RBS rank, gated on DKD specificity) ===' % len(final))
    log('  %-4s %-10s %-8s %-9s %-9s %-9s %s'
        % ('#', 'symbol', 'RBSrank', 'vs ctrl', 'vs othCKD', 'stability', 'perturb'))
    for i, (gid, r) in enumerate(final.iterrows(), 1):
        log('  %-4d %-10s %-8.0f %+9.2f %+9.2f %9.2f %.2f'
            % (i, r['symbol'], r['mean_rank'], r['g_DKD_vs_control'],
               r['g_DKD_vs_otherCKD'], r['mean_S'], r['mean_P']))

    dropped = df[df['procurement_flag']].head(12)
    log('\n=== dropped as procurement-driven (high RBS, but other-CKD moves as much or more) ===')
    for gid, r in dropped.iterrows():
        log('  %-10s RBSrank=%-6.0f DKD vs ctrl=%+.2f   otherCKD vs ctrl=%+.2f'
            % (r['symbol'], r['mean_rank'], r['g_DKD_vs_control'], r['g_otherCKD_vs_control']))

    log('\nwrote %s and %s' % (args.out, args.out.replace('.tsv', '_gated.tsv')))


if __name__ == '__main__':
    main()
