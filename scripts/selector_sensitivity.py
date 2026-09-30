#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Is the cohort-composition instability a property of one selector pair, or general?

The headline sweep (results/cohort_sensitivity/) runs stability selection with ReliefF as
its inner selector against the connectivity comparator over all 26 cohort subsets. The
limitation that follows is real: reversals measured with one inner selector say nothing about
the others. This script repeats the same sweep with each of the other inner selectors that
Table 6 (the correction sweep) already uses, writes each to its own folder so the headline
numbers are untouched, and summarises the verdict instability per selector.

    python scripts/selector_sensitivity.py            # run missing sweeps, then summarise
    python scripts/selector_sensitivity.py --summary  # summarise what exists

Per selector, on the ten three-cohort subsets: range of the RBS minus connectivity difference
in cross-fold Jaccard, how many subsets are negative (verdict reversed), the same for external
AUROC, and the five-cohort verdict. What it does not show: anything about selectors not run,
or about the reference WGCNA implementation.
"""
import argparse
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HEADLINE = 'results/cohort_sensitivity/subsets.tsv'
OUT = 'results/cohort_sensitivity_by_selector'
# ReliefF is the headline; the rest are the inner selectors of the correction sweep.
SELECTORS = ['univariate', 'lasso', 'elastic_net', 'mrmr', 'rf', 'boruta']
NICE = {'relieff': 'ReliefF', 'univariate': 'Univariate F', 'lasso': 'LASSO',
        'elastic_net': 'Elastic net', 'mrmr': 'mRMR', 'rf': 'Random forest', 'boruta': 'Boruta'}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def diffs(path):
    d = pd.read_csv(path, sep='\t')
    p = d.pivot_table(index=['subset', 'size'], columns='method',
                      values=['cross_fold_jaccard', 'external_auroc'])
    j = (p[('cross_fold_jaccard', 'RBS')] - p[('cross_fold_jaccard', 'WGCNA_hub')]).dropna()
    a = (p[('external_auroc', 'RBS')] - p[('external_auroc', 'WGCNA_hub')]).dropna()
    return j, a


def summarise():
    rows = []
    for sel in ['relieff'] + SELECTORS:
        path = HEADLINE if sel == 'relieff' else os.path.join(OUT, sel, 'subsets.tsv')
        if not os.path.exists(path):
            log('  (no sweep yet) %s' % sel)
            continue
        j, a = diffs(path)
        j3 = j.xs(3, level='size'); a3 = a.xs(3, level='size')
        j5 = j.xs(5, level='size'); a5 = a.xs(5, level='size')
        rows.append(dict(
            selector=sel, name=NICE[sel], n3=len(j3),
            j3_lo=j3.min(), j3_hi=j3.max(), j3_neg=int((j3 < 0).sum()),
            a3_lo=a3.min(), a3_hi=a3.max(), a3_neg=int((a3 < 0).sum()),
            j5=float(j5.iloc[0]), a5=float(a5.iloc[0]),
            # verdict reversal: sign of the three-cohort majority vs the five-cohort verdict
            reversal=bool((j3 > 0).sum() > len(j3) / 2) != bool(j5.iloc[0] > 0),
        ))
    s = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    s.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)
    log('%-14s %-22s %-6s %-22s %-6s %-16s %s' % (
        'selector', '3-cohort Jaccard diff', 'neg', '3-cohort AUROC diff', 'neg', 'five-cohort J/A', 'reversal'))
    for r in rows:
        log('%-14s %+.3f to %+.3f      %d/%-3d %+.3f to %+.3f      %d/%-3d %+.3f / %+.3f   %s' % (
            r['name'], r['j3_lo'], r['j3_hi'], r['j3_neg'], r['n3'], r['a3_lo'], r['a3_hi'],
            r['a3_neg'], r['n3'], r['j5'], r['a5'], 'yes' if r['reversal'] else 'no'))
    log('\nwrote %s/summary.tsv' % OUT)
    return s


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    ap = argparse.ArgumentParser()
    ap.add_argument('--summary', action='store_true', help='summarise existing sweeps only')
    ap.add_argument('--selectors', default=','.join(SELECTORS))
    args = ap.parse_args()
    if not args.summary:
        for sel in args.selectors.split(','):
            out = os.path.join(OUT, sel)
            if os.path.exists(os.path.join(out, 'by_size.tsv')):
                log('  have %s' % sel)
                continue
            log('  sweeping %s ...' % sel)
            r = subprocess.run([sys.executable, 'scripts/cohort_sensitivity.py', '--min-size', '2',
                                '--base', sel, '--out', out])
            if r.returncode != 0:
                raise SystemExit('sweep failed for %s' % sel)
    summarise()


if __name__ == '__main__':
    main()
