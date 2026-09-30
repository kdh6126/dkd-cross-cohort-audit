#!/usr/bin/env python
"""보관 조건이 통제된 데이터에서도 채취 효과가 나타나는가.

논문의 가장 아픈 한계는 이것입니다 — 어떤 DKD 코호트도 허혈 시간이나 보관 조건을 기록하지
않아서, "생검 대 신장절제" 차이를 조직 처리 탓으로 돌리는 것이 추론에 머뭅니다. 보관 조건이
다른 것이 진짜 원인일 가능성을 배제하지 못합니다.

GSE163603은 그 배제를 한 단계 진행시킵니다.

    disease        Reference 9명 · 당뇨병성 신증 6명
    tissue source  Reference = Nephrectomy · DKD = Biopsy   (같은 교락이 또 있다)
    storage        66개 샘플 전부 'OCT at -80C'             (여기가 다르다)
    segment        레이저 미세절단. DKD는 간질 구획만 제출됨.

보관 조건이 전 샘플 동일하다고 명시돼 있으므로, 여기서 즉시초기 모듈 차이가 나온다면
"보관 방식이 달라서"라는 설명은 배제됩니다. 남는 후보는 채취 방식과 그에 딸린 허혈 시간
입니다. 여전히 직접 측정은 아니지만, 대안 설명 하나를 실제로 지웁니다.

비교는 간질 구획끼리만 합니다. DKD가 간질만 제출했고, Reference는 환자당 7개 구획을 냈기
때문에 전 구획을 쓰면 같은 환자가 7번 세어지는 유사반복이 됩니다.
"""
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_deconfound import IEG   # noqa: E402

ACC = 'GSE163603'
EXPR = 'data/raw/geo/GSE163603/GSE163603_processed_data.csv.gz'
META = 'data/raw/geo/GSE163603/GSE163603_series_matrix.txt.gz'
OUT = 'results/gse163603'
SEED = 0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def read_meta():
    fields = {}
    with gzip.open(META, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            p = [x.strip('"') for x in line.rstrip('\n').split('\t')]
            if line.startswith('!Sample_title'):
                fields['title'] = p[1:]
            elif line.startswith('!Sample_characteristics_ch1'):
                key = p[1].split(':')[0].strip()
                fields[key] = [x.split(': ', 1)[-1] for x in p[1:]]
    return pd.DataFrame(fields)


def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if not np.isfinite(sp) or sp == 0:
        return np.nan
    return (a.mean() - b.mean()) / sp * (1 - 3 / (4 * (na + nb) - 9))


def main():
    os.makedirs(OUT, exist_ok=True)
    m = read_meta()
    m['segment_from_title'] = m['title'].str.replace(r'^\S+\s*', '', regex=True).str.strip()
    m['patient'] = m['title'].str.extract(r'^(\S+)')[0]
    log('샘플 %d · 환자 %d' % (len(m), m['patient'].nunique()))
    log('  질환:      %s' % dict(m['disease'].value_counts()))
    log('  채취:      %s' % dict(m['tissue source'].value_counts()))
    log('  보관:      %s' % dict(m['storage'].value_counts()))
    log('')

    # 간질 구획만. DKD가 간질만 냈고, Reference 전 구획을 쓰면 유사반복이 된다.
    inter = m['segment_from_title'].str.contains('Interstitium', case=False)
    mi = m[inter].copy()
    log('간질 구획만 남김: %d 샘플, 환자 %d명' % (len(mi), mi['patient'].nunique()))
    log('  %s' % dict(mi['disease'].value_counts()))
    if mi['disease'].nunique() < 2:
        log('  두 군이 모두 있지 않습니다. 중단.')
        return 1

    X = pd.read_csv(EXPR, index_col=0)
    X.columns = [c.strip() for c in X.columns]
    cols = [c for c in mi['title'] if c in X.columns]
    missing = [c for c in mi['title'] if c not in X.columns]
    if missing:
        log('  발현표에 없는 샘플 %d개: %s' % (len(missing), ', '.join(missing[:4])))
    mi = mi[mi['title'].isin(cols)]
    Y = X[list(mi['title'])]
    log('  발현 매칭 %d 샘플, 유전자 %d' % (Y.shape[1], Y.shape[0]))

    # count 데이터이므로 CPM 로그 변환 후 유전자별 z
    cpm = np.log2(Y / Y.sum(axis=0) * 1e6 + 1)
    keep = cpm.mean(axis=1) > 1
    cpm = cpm[keep]
    Z = cpm.sub(cpm.mean(axis=1), axis=0).div(cpm.std(axis=1, ddof=1) + 1e-12, axis=0)
    log('  발현 필터 후 유전자 %d' % len(Z))

    present = [g for g in IEG if g in Z.index]
    log('')
    log('  즉시초기유전자 %d / %d 개가 이 데이터에 있습니다: %s'
        % (len(present), len(IEG), ', '.join(present)))
    if len(present) < 5:
        log('  너무 적어 모듈 점수를 만들 수 없습니다. 중단.')
        return 1

    score = Z.loc[present].mean(axis=0)
    dkd = (mi['disease'].str.contains('Diabetic')).values
    g = hedges_g(score.values[dkd], score.values[~dkd])
    t = stats.ttest_ind(score.values[dkd], score.values[~dkd], equal_var=False)

    log('')
    log('=' * 74)
    log('결과 — 즉시초기 모듈 (생검 DKD vs 신장절제 Reference)')
    log('=' * 74)
    log('  DKD(생검) %d명 평균 z      : %+.3f' % (dkd.sum(), score.values[dkd].mean()))
    log('  Reference(신절제) %d명     : %+.3f' % ((~dkd).sum(), score.values[~dkd].mean()))
    log('  Hedges g = %+.2f,  p = %.2e' % (g, t.pvalue))

    # 무작위 유전자 집합 대조 — 모듈이 특별한지 확인
    rng = np.random.default_rng(SEED)
    null = []
    for _ in range(2000):
        idx = rng.choice(len(Z), len(present), replace=False)
        s2 = Z.iloc[idx].mean(axis=0).values
        null.append(hedges_g(s2[dkd], s2[~dkd]))
    null = np.array([x for x in null if np.isfinite(x)])
    p_emp = (np.sum(np.abs(null) >= abs(g)) + 1) / (len(null) + 1)
    log('  같은 크기 무작위 유전자군의 |g| 분포에서의 위치: 경험적 p = %.4f' % p_emp)

    pd.DataFrame({'sample': mi['title'].values, 'patient': mi['patient'].values,
                  'disease': mi['disease'].values,
                  'tissue_source': mi['tissue source'].values,
                  'storage': mi['storage'].values,
                  'ieg_score': score.values}).to_csv(
        os.path.join(OUT, 'ieg_score.tsv'), sep='\t', index=False)
    pd.DataFrame({'hedges_g': [g], 'p': [t.pvalue], 'p_empirical': [p_emp],
                  'n_dkd': [int(dkd.sum())], 'n_ref': [int((~dkd).sum())],
                  'n_ieg_present': [len(present)]}).to_csv(
        os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)

    log('')
    log('=' * 74)
    log('해석')
    log('=' * 74)
    log('  이 데이터의 값어치는 크기가 아니라 메타데이터에 있습니다. 66개 샘플 전부')
    log("  'OCT at -80C'로 보관됐다고 명시돼 있습니다. 다른 어떤 코호트도 이걸 적지")
    log('  않았습니다.')
    log('')
    if p_emp < 0.05:
        log('  보관 조건이 동일한데도 모듈이 두 군을 가릅니다. "보관 방식이 달라서"라는')
        log('  대안 설명이 여기서는 성립하지 않습니다. 남는 후보는 채취 방식과 허혈 시간')
        log('  입니다.')
    else:
        log('  이 데이터에서는 모듈이 무작위 유전자군과 구분되지 않습니다. 환자 %d 대 %d로'
            % (dkd.sum(), (~dkd).sum()))
        log('  검정력이 낮으므로 부재의 증거로 읽지 않습니다.')
    log('')
    log('  한계 — 여전히 허혈 시간은 기록돼 있지 않습니다. 보관 조건 하나를 배제했을')
    log('  뿐이고, DKD 6명 대 Reference 9명으로 표본이 작습니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
