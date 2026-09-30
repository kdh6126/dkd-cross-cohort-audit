#!/usr/bin/env python
"""상위 관련 유전자들이 서로 얼마나 함께 움직이는가.

원고는 "the DKD programme is broadly co-expressed (mean |r| = 0.48 among top-relevance genes)"
라고 적고 있었습니다. 그런데 이 0.48 은 어느 결과 파일에도 없었습니다 — null_control.py 의
주석과 한국어 설명자료에만 있었고, 계산하는 코드가 없었습니다. 모든 수치가 코드로 재생성
된다는 것이 이 논문의 판매 논거인데, 그 논거를 무너뜨리는 숫자가 본문에 하나 있었던 셈입니다.

여기서 정의를 명시하고 계산해 저장합니다.

    상위 관련 유전자  발견 코호트 전체에서 |Hedges' g| (DKD vs 대조) 평균이 가장 큰 K개
    상관              각 코호트 안에서 그 K개 사이의 쌍별 Pearson |r|
    보고 값           코호트별 평균의 코호트 간 평균

무작위 서명이 왜 잘 작동하는지를 뒷받침하는 수치이므로, 무작위 50개 유전자 집합의 같은
값도 함께 냅니다. 상위 유전자가 무작위보다 훨씬 강하게 묶여 있어야 이 설명이 성립합니다.
"""
import os
import sys

import numpy as np
import pandas as pd

OUT = 'results/coexpression.tsv'
H = 'data/processed/harmonized'
K = 50
SEED = 0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def mean_abs_r(X):
    """열(유전자) 사이 쌍별 |Pearson r| 의 평균. 대각선은 뺀다."""
    C = np.corrcoef(X, rowvar=False)
    iu = np.triu_indices_from(C, k=1)
    v = np.abs(C[iu])
    return float(np.nanmean(v))


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    sys.path.insert(0, 'scripts')
    from cohort_sensitivity import ALL5

    g = pd.read_csv('data/processed/effect_sizes_hedges_g.tsv', sep='\t', index_col=0)
    # 코호트별 g 열의 평균 절대값으로 상위 유전자를 고른다
    gm = g.abs().mean(axis=1).sort_values(ascending=False)
    top = [str(x) for x in gm.index[:K]]

    rng = np.random.default_rng(SEED)
    rows = []
    for c in ALL5:
        f = os.path.join(H, c + '_expr.tsv')
        if not os.path.exists(f):
            continue
        E = pd.read_csv(f, sep='\t', index_col=0)
        E.index = E.index.astype(str)
        sel = [x for x in top if x in E.index]
        if len(sel) < 10:
            continue
        Xt = E.loc[sel].values.T
        pool = [x for x in E.index if x not in set(sel)]
        rnd = [mean_abs_r(E.loc[rng.choice(pool, len(sel), replace=False)].values.T)
               for _ in range(20)]
        rows.append(dict(cohort=c, n_genes=len(sel),
                         mean_abs_r_top=mean_abs_r(Xt),
                         mean_abs_r_random=float(np.mean(rnd))))

    d = pd.DataFrame(rows)
    d.loc[len(d)] = dict(cohort='MEAN', n_genes=int(d['n_genes'].mean()),
                         mean_abs_r_top=d['mean_abs_r_top'].mean(),
                         mean_abs_r_random=d['mean_abs_r_random'].mean())
    d.to_csv(OUT, sep='\t', index=False)

    log('상위 %d개 유전자의 코호트 내 평균 |r|' % K)
    log('  %-12s %10s %10s' % ('코호트', '상위', '무작위'))
    for _, r in d.iterrows():
        log('  %-12s %10.3f %10.3f' % (r['cohort'], r['mean_abs_r_top'],
                                       r['mean_abs_r_random']))
    log('')
    log('  %s 에 저장했습니다.' % OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
