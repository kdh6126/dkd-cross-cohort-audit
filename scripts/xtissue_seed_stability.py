#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""간의 작은 부호 반전(-0.008)이 시드를 바꿔도 남는가.

P1 은 간에서 "held" 로 판정됐지만 근거는 열 개 세 코호트 부분집합 가운데 하나,
그것도 -0.008 짜리 반전 하나다. 이 정도 크기는 부트스트랩 시드 하나에 흔들릴 수
있고, 흔들린다면 일반화 주장의 근거가 약해진다. 그래서 같은 부분집합들을 시드만
바꿔 여러 번 돌리고, 반전의 크기와 재현 여부를 함께 적는다.

    python scripts/xtissue_seed_stability.py --tissue liver_masld --seeds 5

부분집합마다 시드별 RBS - WGCNA 교차 폴드 Jaccard 차이를 기록한다. 판정은 하지
않는다. 반전이 몇 개 시드에서 나타나는지, 그 크기가 얼마인지만 보고한다.
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many                       # noqa: E402
from bootstrap_ci import run_all_methods             # noqa: E402
from xtissue_sweep import core_cohorts, ieg_entrez, KIDNEY5   # noqa: E402


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def summarise(out):
    """시드별·부분집합별 요약을 표 하나로 남긴다.

    계산 무작위성만 보는 점검이다. 환자 표본의 불확실성은 여기서 재지 않는다.
    """
    df = pd.read_csv(out, sep='\t')
    rows = []
    for seed, g in df.groupby('seed'):
        neg = g[g.jaccard_diff < 0]
        rows.append(dict(level='seed', key=str(seed), n=len(g),
                         n_negative=len(neg), min_diff=g.jaccard_diff.min(),
                         median_diff=g.jaccard_diff.median(),
                         max_diff=g.jaccard_diff.max(),
                         sign='mixed' if len(neg) else 'positive in all',
                         per_seed='',
                         negative_subsets=';'.join(sorted(neg.subset))))
    for sub, g in df.groupby('subset'):
        neg = g[g.jaccard_diff < 0]
        sign = ('negative in all' if len(neg) == len(g) else
                'positive in all' if len(neg) == 0 else 'mixed')
        rows.append(dict(level='subset', key=sub, n=len(g), n_negative=len(neg),
                         min_diff=g.jaccard_diff.min(), median_diff=g.jaccard_diff.median(),
                         max_diff=g.jaccard_diff.max(), sign=sign,
                         per_seed=';'.join('%+.4f' % v for v in g.sort_values('seed').jaccard_diff),
                         negative_subsets=';'.join(str(x) for x in sorted(neg.seed))))
    neg = df[df.jaccard_diff < 0]
    rows.append(dict(level='all', key='total', n=len(df), n_negative=len(neg),
                     min_diff=df.jaccard_diff.min(), median_diff=df.jaccard_diff.median(),
                     max_diff=df.jaccard_diff.max(),
                     sign='%d subset(s) ever negative' % neg.subset.nunique(),
                     per_seed='', negative_subsets=';'.join(sorted(set(neg.subset)))))
    t = pd.DataFrame(rows)
    p = out.replace('.tsv', '_summary.tsv')
    t.to_csv(p, sep='\t', index=False)
    log('wrote %s' % p)
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tissue', default='liver_masld')
    ap.add_argument('--size', type=int, default=3)
    ap.add_argument('--seeds', type=int, default=5)
    ap.add_argument('--base-seed', type=int, default=20260825)
    ap.add_argument('--K', type=int, default=50)
    ap.add_argument('--B-inner', type=int, default=30)
    ap.add_argument('--base', default='relieff')
    ap.add_argument('--power', type=int, default=6)
    ap.add_argument('--out', default=None)
    ap.add_argument('--summary-only', action='store_true',
                    help='이미 있는 표에서 요약만 다시 만든다')
    a = ap.parse_args()

    if a.tissue == 'kidney_dkd':
        cohorts, root = KIDNEY5, 'data/processed/harmonized'
    else:
        cohorts, root = core_cohorts(a.tissue), 'data/processed/xtissue/%s' % a.tissue
    out = a.out or 'results/xtissue/%s/seed_stability.tsv' % a.tissue
    os.makedirs(os.path.dirname(out), exist_ok=True)

    if a.summary_only:
        summarise(out)
        return

    data, genes = load_many(cohorts, root=root)
    genes = [str(g) for g in genes]
    pos = {g: i for i, g in enumerate(genes)}
    ieg_cols = [pos[g] for g in ieg_entrez() if g in pos]
    subsets = [list(c) for c in itertools.combinations(cohorts, a.size)]
    log('%s: %d개 부분집합 x %d개 시드' % (a.tissue, len(subsets), a.seeds))

    rows = []
    for si in range(a.seeds):
        seed = a.base_seed + 1000 * si
        for i, sub in enumerate(subsets, 1):
            dat = {c: (data[c][0].values, data[c][1]) for c in sub}
            try:
                res = run_all_methods(dat, sub, ieg_cols, a.K, a.B_inner, a.base, seed, a.power)
            except Exception as e:                       # noqa: BLE001
                log('  시드 %d 부분집합 %d 실패: %s' % (seed, i, str(e)[:80]))
                continue
            d = res['RBS']['cross_fold_jaccard'] - res['WGCNA_hub']['cross_fold_jaccard']
            da = res['RBS']['external_auroc'] - res['WGCNA_hub']['external_auroc']
            rows.append(dict(seed=seed, subset='+'.join(sub), jaccard_diff=d, auroc_diff=da))
            pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
        done = pd.DataFrame(rows)
        cur = done[done.seed == seed]
        if len(cur):
            log('  시드 %d: 음수 %d/%d, 최솟값 %+.4f'
                % (seed, int((cur.jaccard_diff < 0).sum()), len(cur), cur.jaccard_diff.min()))

    df = pd.DataFrame(rows)
    df.to_csv(out, sep='\t', index=False)
    log('')
    log('=== 시드별 요약 ===')
    for seed, g in df.groupby('seed'):
        log('  seed %d  음수 %d/%d  최솟값 %+.4f  중앙 %+.4f'
            % (seed, int((g.jaccard_diff < 0).sum()), len(g), g.jaccard_diff.min(),
               g.jaccard_diff.median()))
    log('')
    log('=== 부분집합별 (시드 간) ===')
    for sub, g in df.groupby('subset'):
        log('  %-38s 평균 %+.4f  범위 %+.4f..%+.4f  음수 %d/%d'
            % (sub, g.jaccard_diff.mean(), g.jaccard_diff.min(), g.jaccard_diff.max(),
               int((g.jaccard_diff < 0).sum()), len(g)))
    summarise(out)
    neg = df[df.jaccard_diff < 0]
    log('')
    log('전체 %d회 중 음수 %d회 (%.1f%%)' % (len(df), len(neg), 100 * len(neg) / max(1, len(df))))
    if len(neg):
        log('음수의 크기: 중앙 %+.4f, 최소 %+.4f' % (neg.jaccard_diff.median(), neg.jaccard_diff.min()))
    log('wrote %s' % out)


if __name__ == '__main__':
    main()
