#!/usr/bin/env python
"""DKD에서 본 현상이 신장병 일반의 것인가, DKD 특유의 것인가.

지금까지의 결과는 전부 DKD 하나에서 나왔습니다. "건강인을 대조로 쓰면 신호가 부풀려진다"가
DKD 특유의 문제인지, 신장병 어디서나 그런지는 확인되지 않았습니다. ST001411(신경병증)이
당뇨 합병증 쪽으로 한 번 넓혔지만 층이 달랐습니다(대사체).

ERCB는 이 질문을 같은 자료 안에서 답할 수 있게 해줍니다. 사구체와 세뇨관간질 두 구획에
DN을 포함한 진단 10종이 들어 있고, 비생검 대조군도 있습니다. 그래서 진단마다 똑같은 두
비교를 만들 수 있습니다.

    기존 방식  진단 X vs 비생검 대조군    (모든 공개 연구가 쓰는 설계)
    권고 방식  진단 X vs 나머지 생검 진단  (채취 방식이 맞춰진 설계)

DKD에서 무엇이 나왔는지는 이미 압니다. 여기서는 나머지 8~9개 진단에서도 같은 크기의 축소가
일어나는지를 봅니다.

    전부 비슷하게 축소된다면  -> 이것은 신장 생검 연구 전반의 문제입니다. 논문의 주장이
                               DKD를 넘어 넓어집니다.
    DKD만 크게 축소된다면     -> DKD 특유의 문제이고, 주장을 좁혀야 합니다.

각 진단의 표본 크기가 다르므로 통과 개수만 보지 않고 효과크기 분포로도 봅니다.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_cohort   # noqa: E402

ERCB = ['GSE104948', 'GSE104954']
OUT = 'results/ckd_generalization'
MIN_N = 12          # 이보다 작은 진단은 양쪽 비교 모두 검정력이 없다
SEED = 0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    n1, n0 = len(a), len(b)
    if n1 < 2 or n0 < 2:
        return np.full(a.shape[1], np.nan)
    sp = np.sqrt(((n1 - 1) * a.var(0, ddof=1) + (n0 - 1) * b.var(0, ddof=1)) / (n1 + n0 - 2))
    sp = np.where(sp == 0, np.nan, sp)
    return (a.mean(0) - b.mean(0)) / sp * (1 - 3 / (4 * (n1 + n0) - 9))


def norm(s):
    return ' '.join(str(s).lower().replace('/', ' ').replace('-', ' ')
                    .replace('_', ' ').split())


def main():
    os.makedirs(OUT, exist_ok=True)

    # 두 구획을 읽고, 진단명을 정규화해 공통 진단만 쓴다
    frames = []
    for c in ERCB:
        X, y, p = load_cohort(c, labelled_only=False)
        X.columns = X.columns.astype(str)
        p = p.reset_index(drop=True)
        col = 'group' if 'group' in p.columns else 'stage'
        p['dx'] = p[col].map(norm)
        frames.append((c, X, p))

    common = set(frames[0][2]['dx'].dropna()) & set(frames[1][2]['dx'].dropna())
    nonbiopsy = {'living donor', 'tumor nephrectomy'}
    diagnoses = sorted(common - nonbiopsy)
    log('두 구획 공통 진단 %d개 (비생검 대조군 제외)' % len(diagnoses))

    genes = frames[0][1].columns
    rows = []
    for dx in diagnoses:
        gs_naive, gs_matched, n_case, n_ctrl, n_other = [], [], 0, 0, 0
        for c, X, p in frames:
            M = X.values
            case = (p['dx'] == dx).values
            ctrl = p['dx'].isin(nonbiopsy).values
            other = (~case) & (~ctrl) & p['dx'].notna().values
            n_case += case.sum(); n_ctrl += ctrl.sum(); n_other += other.sum()
            if case.sum() >= 2 and ctrl.sum() >= 2:
                gs_naive.append(hedges_g(M[case], M[ctrl]))
            if case.sum() >= 2 and other.sum() >= 2:
                gs_matched.append(hedges_g(M[case], M[other]))
        if not gs_naive or not gs_matched or n_case < MIN_N:
            log('  건너뜀: %-44s (환자 %d명)' % (dx[:42], n_case))
            continue
        gn = np.nanmean(np.vstack(gs_naive), axis=0)
        gm = np.nanmean(np.vstack(gs_matched), axis=0)
        ok = np.isfinite(gn) & np.isfinite(gm)
        rows.append(dict(
            diagnosis=dx, n_case=int(n_case), n_control=int(n_ctrl), n_other=int(n_other),
            naive_large=int((np.abs(gn[ok]) >= 0.5).sum()),
            matched_large=int((np.abs(gm[ok]) >= 0.5).sum()),
            naive_median=float(np.median(np.abs(gn[ok]))),
            matched_median=float(np.median(np.abs(gm[ok]))),
            shrink=float(np.median(np.abs(gm[ok])) / np.median(np.abs(gn[ok]))),
            wilcoxon_p=float(stats.wilcoxon(np.abs(gn[ok]), np.abs(gm[ok])).pvalue)))
        pd.DataFrame({'entrez_id': genes[ok], 'g_naive': gn[ok],
                      'g_matched': gm[ok]}).to_csv(
            os.path.join(OUT, 'g_%s.tsv' % dx.replace(' ', '_')), sep='\t', index=False)

    d = pd.DataFrame(rows).sort_values('shrink')
    d.to_csv(os.path.join(OUT, 'per_diagnosis.tsv'), sep='\t', index=False)

    log('')
    log('=' * 84)
    log('진단별 — 대조군을 바꾸면 신호가 얼마나 줄어드는가')
    log('=' * 84)
    log('  %-40s %5s %8s %8s %8s' % ('진단', '환자', '|g|>0.5', '축소율', 'Wilcoxon'))
    log('  %-40s %5s %8s' % ('', '', '기존->권고'))
    for _, r in d.iterrows():
        star = ' <- DKD' if 'dn' == r['diagnosis'] or 'diabet' in r['diagnosis'] else ''
        log('  %-40s %5d %4d->%-4d %7.0f%% %9.1e%s'
            % (r['diagnosis'][:38], r['n_case'], r['naive_large'], r['matched_large'],
               100 * r['shrink'], r['wilcoxon_p'], star))

    log('')
    log('=' * 84)
    log('판정')
    log('=' * 84)
    dn = d[d['diagnosis'].str.contains('dn|diabet', case=False, regex=True)]
    others = d[~d.index.isin(dn.index)]
    log('  전체 %d개 진단의 축소율: 중앙 %.0f%% (범위 %.0f%% ~ %.0f%%)'
        % (len(d), 100 * d['shrink'].median(), 100 * d['shrink'].min(),
           100 * d['shrink'].max()))
    if len(dn):
        log('  DKD                    : %.0f%%' % (100 * dn['shrink'].iloc[0]))
        log('  나머지 %d개 진단 중앙값 : %.0f%%' % (len(others), 100 * others['shrink'].median()))
        rank = int((d['shrink'] < dn['shrink'].iloc[0]).sum()) + 1
        log('  DKD의 순위             : %d / %d (작을수록 많이 줄어듦)' % (rank, len(d)))
    log('')
    n_shrunk = int((d['shrink'] < 0.8).sum())
    log('  20%% 이상 줄어든 진단: %d / %d' % (n_shrunk, len(d)))
    if n_shrunk >= len(d) * 0.7:
        log('')
        log('  대부분의 진단에서 같은 일이 일어납니다. 이것은 DKD 특유의 문제가 아니라')
        log('  "환자를 비생검 대조군과 비교하는" 설계 전반의 문제입니다. 논문의 주장은')
        log('  DKD를 넘어 신장 생검 연구 일반으로 넓힐 수 있습니다.')
    elif len(dn) and dn['shrink'].iloc[0] < others['shrink'].median() * 0.8:
        log('')
        log('  DKD에서 유독 크게 줄어듭니다. 주장을 DKD로 좁혀야 합니다.')
    else:
        log('')
        log('  진단마다 다릅니다. 일반화도 국한도 단정할 수 없습니다.')
    log('')
    log('  주의 — 여기서 "권고 방식"의 대조군은 나머지 생검 진단 전부를 묶은 것입니다.')
    log('  ST003255에서 통합이 결과를 만든 전례가 있으므로, 이 수치는 통합 조건에서의')
    log('  값으로만 읽어야 합니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
