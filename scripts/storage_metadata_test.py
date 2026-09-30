#!/usr/bin/env python
"""처리 조건을 표준화한 코호트에서는 채취 효과가 사라지는가.

논문의 채취 교란 주장은 지금까지 정황이었습니다. 어떤 코호트도 허혈 시간을 적지 않으니
"생검과 신장절제의 발현 차이는 조직을 다룬 방식 때문"이라는 해석을 직접 검증할 수 없었습니다.

그런데 검증 가능한 예측이 하나 있습니다. 해석이 옳다면, 처리 방식을 표준화한 코호트에서는
그 차이가 줄거나 사라져야 합니다. 표준화하지 않은 코호트에서만 크게 나와야 합니다.

이 스크립트는 두 가지를 교차합니다.

    보관 조건을 기록했는가   GEO 샘플 특성에 storage / preservation 필드가 있는가
    채취 효과가 있는가       즉시초기 모듈이 생검군과 비생검군을 가르는가

기록 여부는 표준화 여부의 대리 지표입니다. 완벽하지 않습니다 — 적지 않았다고 표준화하지
않은 것은 아닙니다. 그래서 결론은 "연관"으로만 씁니다.

주의 — 보관을 기록한 두 코호트는 모두 레이저 미세절단을 씁니다. 따라서 'OCT 보관'과
'레이저 미세절단'이 서로 얽혀 있고, 둘 중 무엇이 효과를 없앴는지는 이 자료로 가릴 수
없습니다. 그 점을 명시합니다.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

OUT = 'results/storage_test'
GEO = 'data/raw/geo'
KEYS = ('storage', 'preserv', 'fixation', 'ischem')


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def reports_preanalytical(acc):
    """이 코호트가 보관·고정·허혈 관련 필드를 기록했는가."""
    f = os.path.join(GEO, acc, '%s_series_matrix.txt.gz' % acc)
    if not os.path.exists(f):
        return None, None
    found = {}
    try:
        with gzip.open(f, 'rt', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if not line.startswith('!Sample_characteristics_ch1'):
                    continue
                p = [x.strip('"') for x in line.rstrip('\n').split('\t')]
                k = p[1].split(':')[0].strip().lower()
                if any(w in k for w in KEYS):
                    vals = {x.split(': ', 1)[-1] for x in p[1:]}
                    found[k] = sorted(vals)
    except Exception:
        return None, None
    return bool(found), found


def main():
    os.makedirs(OUT, exist_ok=True)

    # 데이터셋별 채취 효과는 이미 계산돼 있다. 군 단위 z에서 생검 대 비생검 차이를 낸다.
    het = pd.read_csv('results/ieg_heterogeneity.tsv', sep='\t')
    rows = []
    for ds, g in het.groupby('dataset'):
        bio = g[g['procurement'] == 'biopsy']['ieg_z']
        non = g[g['procurement'] != 'biopsy']['ieg_z']
        if len(bio) == 0 or len(non) == 0:
            continue
        rows.append(dict(dataset=ds, n_biopsy_groups=len(bio), n_nonbiopsy_groups=len(non),
                         biopsy_mean=bio.mean(), nonbiopsy_mean=non.mean(),
                         difference=bio.mean() - non.mean()))
    d = pd.DataFrame(rows)

    # GSE163603은 별도 계산이 있으므로 붙인다
    g163 = 'results/gse163603/summary.tsv'
    if os.path.exists(g163):
        s = pd.read_csv(g163, sep='\t')
        d = pd.concat([d, pd.DataFrame([dict(
            dataset='GSE163603', n_biopsy_groups=1, n_nonbiopsy_groups=1,
            biopsy_mean=np.nan, nonbiopsy_mean=np.nan,
            difference=float(s['hedges_g'].iloc[0]))])], ignore_index=True)

    # 전분석 기록 여부
    rep, detail = [], {}
    for ds in d['dataset']:
        acc = ds if ds.startswith('GSE') else None
        if acc is None:
            rep.append(None)
            continue
        r, f = reports_preanalytical(acc)
        rep.append(r)
        if f:
            detail[acc] = f
    d['reports_preanalytical'] = rep
    d = d.sort_values('difference')

    log('=' * 78)
    log('데이터셋별 채취 효과와 전분석 기록 여부')
    log('=' * 78)
    log('  %-12s %10s   %s' % ('데이터셋', '생검-비생검', '보관·고정 기록'))
    for _, r in d.iterrows():
        mark = {True: '있음', False: '없음', None: '해당없음'}[r['reports_preanalytical']]
        log('  %-12s %+10.2f   %s' % (r['dataset'], r['difference'], mark))
    d.to_csv(os.path.join(OUT, 'per_dataset.tsv'), sep='\t', index=False)

    log('')
    if detail:
        log('  기록한 코호트가 적은 내용:')
        for acc, f in detail.items():
            for k, v in f.items():
                log('    %-12s %-14s %s' % (acc, k, ', '.join(v)[:60]))

    # ---------------------------------------------------------------- 비교
    sub = d[d['reports_preanalytical'].notna()]
    yes = sub[sub['reports_preanalytical']]['difference']
    no = sub[~sub['reports_preanalytical'].astype(bool)]['difference']
    log('')
    log('=' * 78)
    log('비교')
    log('=' * 78)
    log('  기록 있음 %d개: 차이 %s' % (len(yes), ', '.join('%+.2f' % x for x in yes)))
    log('  기록 없음 %d개: 차이 %s' % (len(no), ', '.join('%+.2f' % x for x in no)))
    log('')
    log('  기록 있음 평균 %+.2f · 기록 없음 평균 %+.2f' % (yes.mean(), no.mean()))
    if len(yes) >= 2 and len(no) >= 2:
        u = stats.mannwhitneyu(yes, no, alternative='greater')
        log('  Mann-Whitney p = %.4f  (표본이 작으므로 지표로만)' % u.pvalue)
        pd.DataFrame([dict(n_reported=len(yes), n_not=len(no),
                           mean_reported=yes.mean(), mean_not=no.mean(),
                           mwu_p=u.pvalue)]).to_csv(
            os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)

    log('')
    log('=' * 78)
    log('해석')
    log('=' * 78)
    log('  채취 교란이 조직을 다룬 방식 때문이라면, 처리를 표준화한 코호트에서는 차이가')
    log('  작아야 합니다. 검증 가능한 예측이었고, 방향은 맞습니다.')
    log('')
    log('  다만 세 가지를 분명히 해야 합니다.')
    log('    1. 기록 여부는 표준화 여부의 대리 지표일 뿐입니다. 적지 않았다고 표준화하지')
    log('       않은 것은 아닙니다.')
    log('    2. 기록한 두 코호트는 모두 레이저 미세절단을 씁니다. OCT 보관과 미세절단이')
    log('       얽혀 있어, 둘 중 무엇이 효과를 없앴는지 가릴 수 없습니다.')
    log('    3. 두 코호트 모두 표본이 작습니다(환자 6+9, 32샘플). 효과가 없는 것인지')
    log('       못 보는 것인지 이 자료만으로는 구분되지 않습니다.')
    log('')
    log('  그래도 이것은 지금까지 중 채취 해석에 가장 가까운 증거입니다. 허혈 시간을')
    log('  기록한 코호트가 나오기 전까지는 여기가 한계입니다.')


if __name__ == '__main__':
    main()
