#!/usr/bin/env python
"""How much does a method-comparison verdict depend on which cohorts you happened to use?

Adding a fifth cohort reversed the sign of the RBS-versus-WGCNA reproducibility comparison
(+0.088 at four cohorts, -0.034 at five) and turned a non-significant AUROC gap into a
significant one. That is not a detail: nearly every cross-cohort benchmark in this literature
reports a point estimate from one fixed set of cohorts.

This runs every subset of the five available cohorts - 10 pairs, 10 triples, 5 quadruples and
the full set - and records, for each, which method would have been declared the winner. The
output is the distribution of verdicts a study could have reported depending on nothing but
which cohorts it had access to.

Each subset gets a full LODO evaluation of all four methods, so the only thing varying is the
cohort composition.
"""
import os, sys, itertools, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many
from bootstrap_ci import run_all_methods
from dkd_deconfound import IEG

OUT = 'results/cohort_sensitivity'
ALL5 = ['GSE30528', 'GSE96804', 'GSE104948', 'GSE142025', 'GSE294519']
METHODS = ['DEG_meta', 'WGCNA_hub', 'RBS', 'RBS_orth']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument('--K', type=int, default=50)
    ap.add_argument('--B-inner', type=int, default=30)
    ap.add_argument('--base', default='relieff')
    ap.add_argument('--power', type=int, default=6)
    ap.add_argument('--min-size', type=int, default=2)
    ap.add_argument('--seed', type=int, default=20260825)
    # 내부 선택기를 바꿔 돌릴 때 헤드라인 결과(relieff)를 덮어쓰지 않도록 출력 폴더를 받는다.
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args()
    OUT = args.out
    os.makedirs(OUT, exist_ok=True)

    data, genes = load_many(ALL5)
    genes = np.array([str(g) for g in genes])
    gpos = {g: i for i, g in enumerate(genes)}
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = {s: g for g, s in zip(gs['entrez_id'], gs['symbol']) if isinstance(s, str)}
    ieg_cols = [gpos[sym2id[s]] for s in IEG if s in sym2id and sym2id[s] in gpos]

    subsets = []
    for r in range(args.min_size, len(ALL5) + 1):
        subsets += [list(c) for c in itertools.combinations(ALL5, r)]
    log('evaluating %d cohort subsets (sizes %d-%d)\n' % (len(subsets), args.min_size, len(ALL5)))

    rows = []
    for i, sub in enumerate(subsets, 1):
        dat = {c: (data[c][0].values, data[c][1]) for c in sub}
        try:
            r = run_all_methods(dat, sub, ieg_cols, args.K, args.B_inner,
                                args.base, args.seed, args.power)
        except Exception as e:
            log('  subset %d failed: %s' % (i, str(e)[:70]))
            continue
        for m, v in r.items():
            rows.append(dict(subset='+'.join(s[3:] for s in sub), size=len(sub),
                             method=m, **v))
        if i % 5 == 0:
            log('  %d / %d subsets' % (i, len(subsets)))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'subsets.tsv'), sep='\t', index=False)

    # ---------------------------------------------------------------- verdict distribution
    log('\n=== which method would have been declared the winner? ===')
    for metric, nice in (('cross_fold_jaccard', 'reproducibility'),
                         ('external_auroc', 'external AUROC')):
        log('\n  by %s' % nice)
        log('    %-8s %-6s %s' % ('size', 'n', '  '.join('%-12s' % m for m in METHODS)))
        for size, g in df.groupby('size'):
            piv = g.pivot(index='subset', columns='method', values=metric)
            piv = piv.dropna(how='any')
            if piv.empty:
                continue
            win = piv.idxmax(axis=1).value_counts()
            log('    %-8d %-6d %s'
                % (size, len(piv),
                   '  '.join('%-12s' % ('%d (%.0f%%)' % (win.get(m, 0),
                                                         100 * win.get(m, 0) / len(piv)))
                             for m in METHODS)))

    # ---------------------------------------------------------------- paired difference
    log('\n=== RBS minus WGCNA, across subsets of each size ===')
    log('  %-6s %-8s %-30s %s' % ('size', 'n', 'reproducibility [range]', 'AUROC [range]'))
    srows = []
    for size, g in df.groupby('size'):
        piv_j = g.pivot(index='subset', columns='method', values='cross_fold_jaccard').dropna()
        piv_a = g.pivot(index='subset', columns='method', values='external_auroc').dropna()
        dj = (piv_j['RBS'] - piv_j['WGCNA_hub']) if len(piv_j) else pd.Series(dtype=float)
        da = (piv_a['RBS'] - piv_a['WGCNA_hub']) if len(piv_a) else pd.Series(dtype=float)
        log('  %-6d %-8d %+.3f  [%+.3f, %+.3f]          %+.3f  [%+.3f, %+.3f]'
            % (size, len(piv_a),
               dj.mean() if len(dj) else np.nan, dj.min() if len(dj) else np.nan,
               dj.max() if len(dj) else np.nan,
               da.mean(), da.min(), da.max()))
        srows.append(dict(size=size, n_subsets=len(piv_a),
                          jaccard_diff_mean=dj.mean() if len(dj) else np.nan,
                          jaccard_diff_min=dj.min() if len(dj) else np.nan,
                          jaccard_diff_max=dj.max() if len(dj) else np.nan,
                          auroc_diff_mean=da.mean(), auroc_diff_min=da.min(),
                          auroc_diff_max=da.max(),
                          jaccard_sign_flips=int(((dj > 0).any() and (dj < 0).any()))
                          if len(dj) else 0))
    pd.DataFrame(srows).to_csv(os.path.join(OUT, 'by_size.tsv'), sep='\t', index=False)

    log('\n=== the most and least favourable subsets for RBS ===')
    piv = df.pivot_table(index=['subset', 'size'], columns='method',
                         values='cross_fold_jaccard').dropna()
    d = (piv['RBS'] - piv['WGCNA_hub']).sort_values()
    for lbl, sel in (('WGCNA ahead by most', d.head(3)), ('RBS ahead by most', d.tail(3))):
        log('  %s:' % lbl)
        for (sub, size), v in sel.items():
            log('    %-28s (n=%d)  RBS - WGCNA = %+.3f' % (sub, size, v))
    log('\n  a benchmark that reported only one of these subsets would have reached the')
    log('  opposite conclusion from one that reported another. That is the finding.')
    log('\nwrote %s/' % OUT)


if __name__ == '__main__':
    main()
