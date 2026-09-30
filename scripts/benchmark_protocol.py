#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""코호트 구성 민감도를 보고 규약으로 정리한다.

기존 cohort_sensitivity 는 방법 차이의 구간과 부호 전환을 냈다. 여기서는 그 위에
**승자 빈도**를 얹는다. "어느 방법이 이겼나" 가 코호트를 어떻게 골랐느냐에 따라
얼마나 바뀌는지가, 이 자료가 다른 분야에 전할 수 있는 부분이기 때문이다.

동률 처리가 핵심이다
--------------------
처음 계산할 때 idxmax() 를 그대로 썼다. 그런데 size=2 에서는 cross-fold Jaccard 가
네 방법 모두 0.000 인 부분집합이 10개 중 2개이고, 나머지도 0.00 대 0.01 수준이다.
idxmax() 는 동률에서 첫 행을 집으므로 그 경우가 전부 첫 방법의 '승리' 로 세어졌다.
그렇게 나온 14/26 은 실제 승리가 아니다.

그래서 두 가지를 함께 낸다.

    strict   동률(최댓값이 둘 이상)인 부분집합은 승자 없음으로 빼고 센다.
    frac     동률이면 해당 방법들에 1/n 씩 나눠 준다.

그리고 크기를 섞지 않는다. 중심 분석인 size=3 을 따로 보고하고, 크기를 섞은 값은
"구성 크기에 따라 결론이 달라진다" 의 증거로만 쓴다.
"""
import os
import sys

import numpy as np
import pandas as pd

SRC = 'results/cohort_sensitivity/subsets.tsv'
OUT = 'results/cohort_sensitivity'
METRICS = ('cross_fold_jaccard', 'external_auroc')
TOL = 1e-9


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def winners(d, metric):
    """한 부분집합의 승자 목록. 최댓값을 공유하면 모두 돌려준다."""
    mx = d[metric].max()
    return list(d.loc[d[metric] >= mx - TOL, 'method'])


def tally(s, metric, sizes=None):
    sub = s if sizes is None else s[s['size'].isin(sizes)]
    strict, frac, n_tie, n_tot = {}, {}, 0, 0
    for _, d in sub.groupby('subset'):
        w = winners(d, metric)
        n_tot += 1
        if len(w) > 1:
            n_tie += 1
            for m in w:
                frac[m] = frac.get(m, 0.0) + 1.0 / len(w)
            continue
        strict[w[0]] = strict.get(w[0], 0) + 1
        frac[w[0]] = frac.get(w[0], 0.0) + 1.0
    return strict, frac, n_tie, n_tot


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    s = pd.read_csv(SRC, sep='\t')
    methods = sorted(s['method'].unique())
    rows = []

    log('=' * 74)
    log('코호트 구성 민감도 — 승자 빈도')
    log('=' * 74)
    for metric in METRICS:
        log('')
        log('[%s]' % metric)
        for sz in sorted(s['size'].unique()):
            st, fr, n_tie, n_tot = tally(s, metric, [sz])
            log('  size=%d  부분집합 %2d개 (동률 %d개)  strict: %s'
                % (sz, n_tot, n_tie,
                   ', '.join('%s %d' % (m, st[m]) for m in methods if m in st) or '없음'))
            for m in methods:
                rows.append(dict(metric=metric, size=sz, method=m,
                                 n_subsets=n_tot, n_tied=n_tie,
                                 wins_strict=st.get(m, 0), wins_frac=fr.get(m, 0.0)))
        st, fr, n_tie, n_tot = tally(s, metric)
        log('  전체    부분집합 %2d개 (동률 %d개)  strict: %s'
            % (n_tot, n_tie,
               ', '.join('%s %d' % (m, st[m]) for m in methods if m in st)))
        for m in methods:
            rows.append(dict(metric=metric, size=0, method=m, n_subsets=n_tot,
                             n_tied=n_tie, wins_strict=st.get(m, 0),
                             wins_frac=fr.get(m, 0.0)))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'winner_frequency.tsv'),
                              sep='\t', index=False)

    # ---- 두 지표가 같은 방법을 뽑는가. 동률인 부분집합은 제외한다.
    log('')
    log('=' * 74)
    log('재현성 1위와 AUROC 1위가 같은 방법인가')
    log('=' * 74)
    arows = []
    for sz in sorted(s['size'].unique()) + [0]:
        sub = s if sz == 0 else s[s['size'] == sz]
        agree = n = 0
        for _, d in sub.groupby('subset'):
            a, b = winners(d, 'cross_fold_jaccard'), winners(d, 'external_auroc')
            if len(a) > 1 or len(b) > 1:
                continue
            n += 1
            agree += int(a[0] == b[0])
        lab = '전체' if sz == 0 else 'size=%d' % sz
        log('  %-8s 판정 가능 %2d개 중 일치 %2d개 %s'
            % (lab, n, agree, '(%.0f%%)' % (100 * agree / n) if n else ''))
        arows.append(dict(size=sz, n_comparable=n, n_agree=agree,
                          frac_agree=(agree / n) if n else np.nan))
    pd.DataFrame(arows).to_csv(os.path.join(OUT, 'metric_agreement.tsv'),
                               sep='\t', index=False)

    log('')
    log('  읽는 법: size=2 는 fold 가 둘뿐이라 cross-fold Jaccard 가 거의 0 으로 붙는다.')
    log('  중심 분석은 size=3 이고, 크기를 섞은 값은 구성에 따라 결론이 달라진다는')
    log('  증거로만 쓴다.')
    log('')
    log('  %s/winner_frequency.tsv, metric_agreement.tsv 에 저장' % OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
