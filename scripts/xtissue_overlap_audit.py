#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""같은 연구실·같은 플랫폼 계열끼리 같은 표본을 다시 올렸는지 값 벡터로 확인한다.

신장에서 제외한 세 계열은 접근번호로는 보이지 않고 표본 수준에서만 보였다. 여기서도 같은
방법을 쓴다. 두 계열이 같은 플랫폼이면 공통 프로브에서 A 의 표본마다 B 의 모든 표본과의
상관을 구한다.

같은 칩을 다시 올린 표본은 다른 사람과의 상관 분포에서 동떨어지게 높다. 재정규화를 거치면
값이 정확히 같지는 않으므로 절대 기준 하나로 자르지 않는다. 표본마다 최대 상관이 0.995
이상이고, 그 표본의 나머지 상관 분포보다 표준편차 5배 이상 높을 때 '같은 표본' 으로 본다.
"""
import gzip
import itertools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xtissue_harmonize import matrix_path, parse_matrix   # noqa: E402

PAIRS = {
    'colon_uc': [('GSE75214', 'GSE59071', 'GPL6244'), ('GSE75214', 'GSE48958', 'GPL6244'),
                 ('GSE59071', 'GSE48958', 'GPL6244'),
                 ('GSE16879', 'GSE38713', 'GPL570'), ('GSE16879', 'GSE47908', 'GPL570'),
                 ('GSE16879', 'GSE13367', 'GPL570'), ('GSE38713', 'GSE47908', 'GPL570'),
                 ('GSE13367', 'GSE47908', 'GPL570'), ('GSE13367', 'GSE38713', 'GPL570')],
    'liver_masld': [('GSE48452', 'GSE61260', 'GPL11532'), ('GSE24807', 'GSE17470', 'GPL2895')],
}
OUT = 'results/xtissue'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def logged(df):
    v = df.values[np.isfinite(df.values)]
    if len(v) and np.percentile(v, 99) > 100:
        return np.log2(df.clip(lower=0) + 1)
    return df


def gse61260_supplement():
    """GSE61260 은 표 본문이 비어 있고 정규화 행렬이 기호 x Sample1..134 로만 올라 있다.
    감사에는 표본 라벨이 필요 없으므로 이 행렬을 그대로 쓰되, GSE48452 쪽을 같은 기호 공간으로
    바꿔 맞댄다."""
    path = 'data/raw/geo/GSE61260/GSE61260_datLiverNormalizedExpr.csv.gz'
    d = pd.read_csv(path, index_col=0)
    return [dict(title=c) for c in d.columns], d


def to_symbols(df, gpl):
    ann = pd.read_csv('data/raw/annotation/%s.annot.gz' % gpl, sep='\t', comment='!',
                      dtype=str, low_memory=False, skiprows=lambda i: False)
    return df


def gse48452_symbols():
    s, d = parse_matrix(matrix_path('GSE48452', 'GPL11532'))
    txt = gzip.open('data/raw/annotation/GPL11532.annot.gz', 'rt', encoding='utf-8',
                    errors='replace').read()
    body = txt.split('!platform_table_begin\n', 1)[-1].split('!platform_table_end', 1)[0]
    import io
    ann = pd.read_csv(io.StringIO(body), sep='\t', dtype=str, low_memory=False)
    m = dict(zip(ann['ID'], ann['Gene symbol']))
    d = d.loc[[i for i in d.index if isinstance(m.get(i), str) and '///' not in m[i]]]
    d.insert(0, 'sym', [m[i] for i in d.index])
    d['_m'] = d.drop(columns='sym').mean(axis=1)
    d = d.sort_values('_m', ascending=False).drop_duplicates('sym').drop(columns='_m')
    return s, d.set_index('sym')


def main():
    rows = []
    for domain, pairs in PAIRS.items():
        for a, b, gpl in pairs:
            try:
                if b == 'GSE61260':
                    sa, da = gse48452_symbols()
                    sb, db = gse61260_supplement()
                else:
                    sa, da = parse_matrix(matrix_path(a, gpl))
                    sb, db = parse_matrix(matrix_path(b, gpl))
            except Exception as e:  # noqa: BLE001
                log('%s vs %s: 실패 %s' % (a, b, e))
                continue
            da, db = logged(da), logged(db)
            common = da.index.intersection(db.index)
            A = da.loc[common].dropna()
            B = db.loc[A.index.intersection(db.loc[common].dropna().index)]
            A = A.loc[B.index]
            Az = (A - A.mean(axis=0)) / A.std(axis=0)
            Bz = (B - B.mean(axis=0)) / B.std(axis=0)
            R = (Az.T.values @ Bz.values) / (len(A) - 1)       # samples_A x samples_B
            # 처음에는 최대 상관이 나머지 분포보다 z 5 이상 튀는지로 판정했다. 조직 배열은
            # 서로 다른 사람끼리도 0.98 대라, 값이 완전히 같은 표본도 z 5 에 못 미쳤다. 판정을
            # 둘로 나눈다. 값이 사실상 같으면(r >= 0.9999) 같은 데이터를 다시 올린 것이다.
            # 재처리를 거쳐 값이 달라졌으면, 표본마다 최고 상관과 차선 상관의 간격이 무관한
            # 쌍보다 체계적으로 크다. 간격의 중앙값을 쌍 단위로 비교한다.
            srt = np.sort(R, axis=1)[:, ::-1]
            best, gap = srt[:, 0], srt[:, 0] - srt[:, 1]
            titles_a = {s['title'] for s in sa}
            titles_b = {s['title'] for s in sb}
            rows.append(dict(domain=domain, series_a=a, series_b=b, platform=gpl,
                             samples_a=R.shape[0], samples_b=R.shape[1], features=len(A),
                             identical_samples=int((best >= 0.9999).sum()),
                             distinct_partner_samples=int((gap >= 0.01).sum()),
                             median_best_r=float(np.median(best)),
                             median_gap=float(np.median(gap)),
                             shared_titles=len(titles_a & titles_b)))
            log('%-9s vs %-9s  %-8s identical=%d/%d  gap>=0.01=%d  median gap=%.4f  titles=%d'
                % (a, b, gpl, rows[-1]['identical_samples'], R.shape[0],
                   rows[-1]['distinct_partner_samples'], rows[-1]['median_gap'],
                   rows[-1]['shared_titles']))
    df = pd.DataFrame(rows)
    # 무관한 쌍(동일 표본이 없는 서로 다른 연구실 쌍)의 간격 중앙값을 기준으로 삼는다
    ref = df.loc[(df['identical_samples'] == 0) & (df['shared_titles'] == 0), 'median_gap']
    base = float(ref.median()) if len(ref) else float('nan')

    def verdict(r):
        if r['identical_samples'] > 0 or r['shared_titles'] > 0:
            return 'shared samples: identical values or titles'
        if r['median_gap'] >= 5 * base:
            return 'shared individuals likely: re-processed data (median gap >= 5x reference)'
        if r['distinct_partner_samples'] > 0 and r['median_best_r'] >= 0.99:
            return 'isolated high match: flag'
        return 'no evidence of shared samples'
    df['reference_median_gap'] = base
    df['verdict'] = df.apply(verdict, axis=1)
    df.to_csv(os.path.join(OUT, 'overlap_audit.tsv'), sep='\t', index=False)
    log('\nreference median gap (unrelated pairs) = %.4f' % base)
    for _, r in df.iterrows():
        log('  %-9s vs %-9s  %s' % (r['series_a'], r['series_b'], r['verdict']))


if __name__ == '__main__':
    main()
