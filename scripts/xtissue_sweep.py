#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""예측 P1 과 P3: 코호트 구성 스윕을 새 조직에서 원고와 같은 코드로 돌린다.

    P1  세 코호트 부분집합들에서 RBS - WGCNA_hub 의 교차 폴드 Jaccard 차이가 0 을 가로지르고,
        외부 AUROC 차이는 Jaccard 차이보다 부호가 덜 바뀐다.
    P3  보정하지 않은 RBS 상위 50 에 들어간 즉시초기 유전자 수(LODO 폴드 평균, 전체 코호트
        집합 기준)가 신장 > 간 > 대장 순이고, 대장에서는 1 이하이다.

부분집합 평가는 bootstrap_ci.run_all_methods 를 그대로 부른다. 원고의 구성 스윕
(cohort_sensitivity.py) 과 같은 함수, 같은 K, 같은 안쪽 B(30), 같은 기준 선택기(relieff),
같은 소프트 파워(6)다. 달라지는 것은 코호트와 유전자 공간뿐이다. 즉시초기 유전자 수를 세려고
return_orders 로 폴드별 순위를 함께 받는다.

신장은 --tissue kidney_dkd 로 같은 코드를 전체 다섯 코호트 한 집합에만 돌려 P3 의 기준값을 만든다.
"""
import argparse
import gzip
import itertools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many                       # noqa: E402
from bootstrap_ci import run_all_methods             # noqa: E402
from dkd_deconfound import IEG                       # noqa: E402

KIDNEY5 = ['GSE30528', 'GSE96804', 'GSE104948', 'GSE142025', 'GSE294519']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def ieg_entrez():
    sym = {}
    with gzip.open('data/raw/annotation/Homo_sapiens.gene_info.gz', 'rt', encoding='utf-8',
                   errors='replace') as f:
        hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip('\n').split('\t')
            sym[p[ix['Symbol']]] = p[ix['GeneID']]
    return [sym[s] for s in IEG if s in sym]


def core_cohorts(tissue):
    from xtissue_config import COHORTS
    tab = pd.read_csv('results/xtissue/%s/cohorts.tsv' % tissue, sep='\t')
    ok = set(tab.loc[tab['n_case'].fillna(0) > 0, 'gse'])
    return [c for c, cfg in COHORTS[tissue].items() if cfg['status'] == 'core' and c in ok]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tissue', required=True)
    ap.add_argument('--K', type=int, default=50)
    ap.add_argument('--B-inner', type=int, default=30)
    ap.add_argument('--base', default='relieff')
    ap.add_argument('--power', type=int, default=6)
    ap.add_argument('--seed', type=int, default=20260825)
    ap.add_argument('--full-only', action='store_true')
    ap.add_argument('--root', default=None, help='다른 조화 디렉터리(민감도 분석용)')
    ap.add_argument('--tag', default='sweep')
    args = ap.parse_args()

    if args.tissue == 'kidney_dkd':
        cohorts, root = KIDNEY5, 'data/processed/harmonized'
    else:
        cohorts, root = core_cohorts(args.tissue), 'data/processed/xtissue/%s' % args.tissue
    if args.root:
        root = args.root
    out = 'results/xtissue/%s/%s' % (args.tissue, args.tag)
    os.makedirs(out, exist_ok=True)

    data, genes = load_many(cohorts, root=root)
    genes = [str(g) for g in genes]
    pos = {g: i for i, g in enumerate(genes)}
    ieg_cols = [pos[g] for g in ieg_entrez() if g in pos]
    ieg_set = set(ieg_cols)
    log('%s: cohorts %s, genes %d, IEG present %d' % (args.tissue, cohorts, len(genes),
                                                     len(ieg_cols)))

    if args.full_only:
        subsets = [list(cohorts)]
    else:
        subsets = [list(c) for r in range(2, len(cohorts) + 1)
                   for c in itertools.combinations(cohorts, r)]
    rows = []
    for i, sub in enumerate(subsets, 1):
        dat = {c: (data[c][0].values, data[c][1]) for c in sub}
        try:
            res, orders = run_all_methods(dat, sub, ieg_cols, args.K, args.B_inner, args.base,
                                          args.seed, args.power, return_orders=True)
        except Exception as e:  # noqa: BLE001
            log('  subset %d failed: %s' % (i, str(e)[:100]))
            continue
        for m, v in res.items():
            n_ieg = float(np.mean([len(set(orders[m][h][:args.K]) & ieg_set) for h in sub]))
            rows.append(dict(subset='+'.join(sub), size=len(sub), method=m, n_ieg_topK=n_ieg,
                             **v))
        log('  %d/%d %s' % (i, len(subsets), '+'.join(sub)))
        pd.DataFrame(rows).to_csv(os.path.join(out, 'subsets.tsv'), sep='\t', index=False)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, 'subsets.tsv'), sep='\t', index=False)
    if args.full_only:
        return
    # P1
    rows = []
    for size, g in df.groupby('size'):
        pj = g.pivot(index='subset', columns='method', values='cross_fold_jaccard').dropna()
        pa = g.pivot(index='subset', columns='method', values='external_auroc').dropna()
        dj = pj['RBS'] - pj['WGCNA_hub']
        da = pa['RBS'] - pa['WGCNA_hub']
        rows.append(dict(size=size, n=len(pj), jaccard_min=dj.min(), jaccard_max=dj.max(),
                         jaccard_spans_zero=bool(dj.min() < 0 < dj.max()),
                         jaccard_n_negative=int((dj < 0).sum()),
                         auroc_min=da.min(), auroc_max=da.max(),
                         auroc_n_negative=int((da < 0).sum())))
    bs = pd.DataFrame(rows)
    bs.to_csv(os.path.join(out, 'by_size.tsv'), sep='\t', index=False)
    log(bs.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
