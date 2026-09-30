#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""correction_sweep 이 지금까지 끝낸 selector 만 모아 sweep.tsv 를 만든다.

sweep 전체가 끝나야만 집계되면, 오래 걸리는 selector(boruta, rf) 하나 때문에 장표를
만들 수 없다. 끝난 것만으로 먼저 집계하고, 나중에 다시 돌리면 갱신되게 함.
"""
import glob
import os
import sys

import pandas as pd

OUT = 'results/correction_sweep'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    rows = []
    for f in sorted(glob.glob(os.path.join(OUT, '*', 'benchmark.tsv'))):
        sel = os.path.basename(os.path.dirname(f))
        t = pd.read_csv(f, sep='\t')
        t = t[t['K'] == 50]
        g = t.groupby('arm').agg(auroc=('external_auroc', 'mean'),
                                 ieg=('n_ieg', 'mean'),
                                 proc=('procurement_frac', 'mean'),
                                 spec=('spec_pass_frac', 'mean'),
                                 xfold=('cross_fold_jaccard', 'mean')).reset_index()
        g.insert(0, 'selector', sel)
        rows.append(g)
    if not rows:
        log('아직 끝난 selector 가 없음.')
        return 1
    df = pd.concat(rows, ignore_index=True)
    df['B'] = 40
    df.to_csv(os.path.join(OUT, 'sweep.tsv'), sep='\t', index=False)
    log('selector %d종 집계: %s' % (df['selector'].nunique(),
                                   ', '.join(sorted(df['selector'].unique()))))
    p = df[df['arm'].isin(['none', 'IEG_resid'])].pivot_table(
        index='selector', columns='arm', values='proc')
    log('')
    log('procurement 관련 feature 비율 (%)')
    log((100 * p).round(0).to_string())
    return 0


if __name__ == '__main__':
    sys.exit(main())
