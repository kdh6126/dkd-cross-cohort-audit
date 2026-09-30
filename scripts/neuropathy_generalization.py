#!/usr/bin/env python
"""논문의 핵심 권고가 다른 질병·다른 층에서도 성립하는가.

논문은 "환자를 건강인과 비교하지 말고, 같은 방식으로 채취된 다른 질병과 비교하라"고
권고합니다. DKD에서는 그 권고를 따르면 신호의 상당 부분이 사라졌습니다.

ST001411은 그 권고를 더 깨끗하게 시험할 수 있는 설계입니다. 미시간대가 모은 혈장 대사체
106명분인데, 세 군이 이렇게 들어 있습니다.

    당뇨 + 신경병증   48명
    당뇨, 신경병증 없음 49명
    정상              9명

DKD 자료(ST003255)보다 나은 점이 하나 있습니다. 거기서는 대조군이 "다른 신장병"이라 서로
다른 세 질병을 묶어야 했고, 그 통합이 결과를 만들었습니다. 여기서는 대조군이 "같은 당뇨인데
합병증만 없는 사람"이라 묶을 필요가 없습니다. 기저 질환이 고정됩니다.

    비교 1  신경병증 vs 정상          — 기존 방식. 당뇨와 합병증이 섞여 있다.
    비교 2  신경병증 vs 당뇨(합병증 X) — 논문이 권고하는 방식. 당뇨가 고정된다.
    비교 3  당뇨(합병증 X) vs 정상     — 비교 1에서 빠지는 성분이 무엇인지 보여준다.

비교 1의 신호 중 얼마가 비교 2에서 살아남는지가 답입니다. 적게 살아남을수록 "건강인 대조"가
과대평가를 낳는다는 논문의 주장이 강해집니다.

주의 — 정상군이 9명뿐입니다. 비교 1과 3은 검정력이 낮으므로, 통과 개수를 직접 비교하지 않고
효과크기 분포로도 함께 봅니다.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

SRC = 'data/raw/metabolomics/ST001411/AN002361_datatable.tsv'
OUT = 'results/neuropathy'
SEED = 0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if not np.isfinite(sp) or sp == 0:
        return np.nan
    return (a.mean() - b.mean()) / sp * (1 - 3 / (4 * (na + nb) - 9))


def bh(p):
    p = np.where(np.isfinite(p), p, 1.0)
    o = np.argsort(p)
    m = len(p)
    q = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.clip(q, 0, 1)
    return out


def contrast(Z, mets, ma, mb):
    g = np.array([hedges_g(Z.loc[ma, m].values, Z.loc[mb, m].values) for m in mets])
    p = np.array([stats.ttest_ind(Z.loc[ma, m].values, Z.loc[mb, m].values,
                                  equal_var=False).pvalue for m in mets])
    return g, bh(p)


def main():
    os.makedirs(OUT, exist_ok=True)
    d = pd.read_csv(SRC, sep='\t')
    d = d.rename(columns={d.columns[0]: 'sample', d.columns[1]: 'group'})
    d['group'] = d['group'].str.replace('Group:', '', regex=False).str.strip()
    mets = [c for c in d.columns if c not in ('sample', 'group')]

    X = d[mets].apply(pd.to_numeric, errors='coerce')
    X = X.loc[:, X.notna().mean() >= 0.8]
    mets = list(X.columns)
    X = X.fillna(X.median())
    X = np.log2(X.clip(lower=X[X > 0].min().min() / 2))
    Z = (X - X.mean()) / (X.std(ddof=1) + 1e-12)

    g = d['group'].values
    neuro = g == 'Diabetic Neuropathy'
    dm_only = g == 'Diabetic non neuropathy'
    normal = g == 'Normal'
    log('군 구성: 신경병증 %d · 당뇨만 %d · 정상 %d' % (neuro.sum(), dm_only.sum(), normal.sum()))
    log('80%% 이상에서 정량된 대사물질 %d / %d' % (len(mets), len(d.columns) - 2))
    log('')

    # 파일 이름과 표의 값은 영어로 둔다. 이 결과는 공개 저장소로 나가므로, 한국어 파일명은
    # 일부 시스템에서 깨지고 표를 여는 사람이 값을 읽지 못한다. 콘솔 설명만 작업 언어로
    # 남긴다.
    specs = [('neuropathy_vs_normal', 'neuropathy vs normal (conventional)',
              '신경병증 vs 정상 (기존 방식)', neuro, normal),
             ('neuropathy_vs_diabetes_only', 'neuropathy vs diabetes only (recommended)',
              '신경병증 vs 당뇨만 (권고 방식)', neuro, dm_only),
             ('diabetes_only_vs_normal', 'diabetes only vs normal (effect of diabetes)',
              '당뇨만 vs 정상 (당뇨 자체 효과)', dm_only, normal)]
    res, gmap = [], {}
    for slug, label, ko, ma, mb in specs:
        gg, qq = contrast(Z, mets, ma, mb)
        gmap[ko] = gg
        res.append(dict(contrast=label, n_a=int(ma.sum()), n_b=int(mb.sum()),
                        n_sig=int((qq < 0.05).sum()),
                        n_large=int((np.abs(gg) >= 0.5).sum()),
                        median_abs_g=float(np.nanmedian(np.abs(gg))),
                        max_abs_g=float(np.nanmax(np.abs(gg)))))
        pd.DataFrame({'metabolite': mets, 'hedges_g': gg, 'q': qq}).to_csv(
            os.path.join(OUT, 'ST001411_%s.tsv' % slug), sep='\t', index=False)

    t = pd.DataFrame(res)
    t.to_csv(os.path.join(OUT, 'contrasts.tsv'), sep='\t', index=False)
    log('=' * 76)
    log('세 비교')
    log('=' * 76)
    log('  %-32s %5s %5s %8s %10s %10s'
        % ('비교', 'n_a', 'n_b', 'q<0.05', '|g|>=0.5', '중앙|g|'))
    for _, r in t.iterrows():
        log('  %-32s %5d %5d %8d %10d %10.3f'
            % (r['contrast'], r['n_a'], r['n_b'], r['n_sig'], r['n_large'],
               r['median_abs_g']))

    naive, matched = t.loc[0, 'n_sig'], t.loc[1, 'n_sig']
    log('')
    log('=' * 76)
    log('핵심 — 기존 방식의 신호는 권고 방식에서 얼마나 살아남는가')
    log('=' * 76)
    log('  기존 방식 유의 %d개 -> 권고 방식 유의 %d개' % (naive, matched))
    if naive:
        log('  생존율 %.0f%%' % (100 * matched / naive))
    else:
        log('  (기존 방식에서 유의한 것이 없습니다. 정상군 9명이라 검정력 부족.)')

    # 효과 크기로도 본다. 개수는 표본 크기에 좌우되므로 이쪽이 더 안정적이다.
    a = np.abs(gmap['신경병증 vs 정상 (기존 방식)'])
    b = np.abs(gmap['신경병증 vs 당뇨만 (권고 방식)'])
    ok = np.isfinite(a) & np.isfinite(b)
    log('')
    log('  효과크기 중앙값 : 기존 %.3f -> 권고 %.3f  (%.0f%%로 축소)'
        % (np.median(a[ok]), np.median(b[ok]), 100 * np.median(b[ok]) / np.median(a[ok])))
    w = stats.wilcoxon(a[ok], b[ok])
    log('  Wilcoxon p = %.2e' % w.pvalue)

    # 당뇨 자체 효과와의 상관 — 기존 방식이 무엇을 잡고 있었는지
    c = gmap['당뇨만 vs 정상 (당뇨 자체 효과)']
    ok2 = np.isfinite(gmap['신경병증 vs 정상 (기존 방식)']) & np.isfinite(c)
    r = stats.pearsonr(gmap['신경병증 vs 정상 (기존 방식)'][ok2], c[ok2])
    log('')
    log('  기존 방식 효과와 "당뇨 자체" 효과의 상관: r = %.3f (p = %.1e)' % (r[0], r[1]))
    log('  (주의 — 두 비교가 같은 정상군을 쓰므로 부풀려져 있습니다. 아래에서 보정합니다.)')

    # ------------------------------------------------ 공유 대조군 보정
    # 위 상관은 두 비교가 같은 정상군 9명을 공유하므로 부풀려져 있다. 그 9명의 표본
    # 잡음이 두 효과크기에 공통으로 들어가기 때문이다. 대조군을 겹치지 않게 쪼개서
    # 다시 잰다.
    log('')
    log('=' * 76)
    log('공유 대조군 보정 — 위 상관이 인위적으로 부풀려진 것 아닌가')
    log('=' * 76)
    idx_norm = np.where(normal)[0]
    rng = np.random.default_rng(SEED)
    rs = []
    for _ in range(200):
        perm = rng.permutation(idx_norm)
        h1, h2 = perm[:len(perm) // 2], perm[len(perm) // 2:]
        m1 = np.zeros(len(g), bool); m1[h1] = True
        m2 = np.zeros(len(g), bool); m2[h2] = True
        ga = np.array([hedges_g(Z.loc[neuro, m].values, Z.loc[m1, m].values) for m in mets])
        gb = np.array([hedges_g(Z.loc[dm_only, m].values, Z.loc[m2, m].values) for m in mets])
        ok3 = np.isfinite(ga) & np.isfinite(gb)
        if ok3.sum() > 100:
            rs.append(stats.pearsonr(ga[ok3], gb[ok3])[0])
    if rs:
        log('  정상군을 %d/%d로 겹치지 않게 나눠 200회 반복'
            % (len(idx_norm) // 2, len(idx_norm) - len(idx_norm) // 2))
        log('  겹치는 대조군 사용 시 r = %.3f' % r[0])
        log('  겹치지 않는 대조군    r = %.3f  (사분위 %.3f ~ %.3f)'
            % (np.median(rs), np.percentile(rs, 25), np.percentile(rs, 75)))
        drop = r[0] - np.median(rs)
        log('  차이 %.3f -> 공유 대조군이 상관의 %.0f%%를 만들었습니다.'
            % (drop, 100 * drop / r[0]))
        if np.median(rs) > 0.6:
            log('  보정 후에도 상관이 높습니다. "기존 방식이 잡은 것은 대부분 당뇨 자체"')
            log('  라는 결론은 유지됩니다.')
        else:
            log('  보정 후 상관이 크게 떨어집니다. 이 결론은 약화됩니다.')
        pd.DataFrame({'r_disjoint': rs}).to_csv(
            os.path.join(OUT, 'disjoint_control_r.tsv'), sep='	', index=False)

    # ---------------------------------------------------------------- 판정
    log('')
    log('=' * 76)
    log('판정')
    log('=' * 76)
    shrink = np.median(b[ok]) / np.median(a[ok])
    r_adj = float(np.median(rs)) if rs else float('nan')
    log('  DKD에서: 기존 방식 84개 -> 권고 방식 29개 (35%% 생존)')
    log('  신경병증에서: 효과크기 기준 %.0f%%로 축소' % (100 * shrink))
    log('  공유 대조군 보정 후, 기존 방식 신호와 "당뇨 자체" 효과의 공유 분산 %.0f%%'
        % (100 * r_adj ** 2))
    log('')
    if shrink < 0.8:
        log('  같은 방향입니다. "건강인을 대조로 쓰면 질병 특이 신호가 과대평가된다"는')
        log('  논문의 권고가 다른 질병·다른 합병증에서도 성립합니다.')
        log('')
        log('  그리고 이 설계가 DKD보다 깨끗합니다. 대조군이 "같은 당뇨인데 합병증만')
        log('  없는 사람"이라, ST003255에서 문제가 됐던 "여러 질병을 묶는 것"이 필요 없습니다.')
    else:
        log('  축소가 크지 않습니다. 이 질병에서는 권고의 효과가 작습니다.')
    log('')
    log('  한계 — 정상군이 9명뿐이라 기존 방식 비교의 검정력이 낮습니다. 개수보다')
    log('  효과크기 비교를 신뢰해야 합니다. 그리고 이것은 한 코호트입니다.')


if __name__ == '__main__':
    main()
