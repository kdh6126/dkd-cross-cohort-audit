#!/usr/bin/env python
"""Two ways the ST003255 result could be an artifact, tested rather than assumed.

The claim under scrutiny is that only 29 of 187 metabolites survive a DKD-versus-other-kidney-
disease contrast, so most of the naive DKD signal is not DKD-specific. Two objections deserve
an answer before that is repeated anywhere:

    반론 1  IgAN, MN, HN을 한 덩어리로 묶은 것이 과도한 단순화 아닌가.
            Pooling three diseases hides whichever one actually differs. If the survivors are
            driven by one comparator, the claim is about that disease, not about "other kidney
            disease". Tested by repeating the contrast against each disease separately and
            asking how many survivors are shared.

    반론 2  살아남은 신호가 DKD 특이성이 아니라 신기능 중증도 차이의 대리 아닌가.
            If the DKD arm simply has worse kidney function than the comparators, any marker of
            function will separate them and look DKD-specific. Tested by residualising every
            metabolite on creatinine -- a direct function marker measured in this same panel --
            and repeating the contrast. This is the same manoeuvre the transcriptomic work uses
            against the handling score, applied to a different nuisance variable.

Neither test can prove specificity. Both can refute it, which is what they are for.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

SRC = 'data/raw/metabolomics/ST003255/AN005337_datatable.tsv'
OUT = 'results/metabolomics'
SEVERITY = 'Creatinine'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    return (a.mean() - b.mean()) / (sp + 1e-12) * (1 - 3 / (4 * (na + nb) - 9))


def bh(p):
    p = np.where(np.isfinite(p), p, 1.0)
    o = np.argsort(p)
    m = len(p)
    q = p[o] * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.clip(q, 0, 1)
    return out


def load():
    d = pd.read_csv(SRC, sep='\t')
    d = d.rename(columns={d.columns[0]: 'sample', d.columns[1]: 'group'})
    d['group'] = d['group'].str.replace('Disease:', '', regex=False)
    mets = [c for c in d.columns if c not in ('sample', 'group')]
    X = d[mets].apply(pd.to_numeric, errors='coerce')
    X = X.loc[:, X.notna().mean() >= 0.8]
    X = X.fillna(X.median())
    X = np.log2(X.clip(lower=X[X > 0].min().min() / 2))
    Z = (X - X.mean()) / (X.std(ddof=1) + 1e-12)
    return Z, d['group'].values, list(X.columns)


def contrast(Z, mets, ma, mb):
    g = np.array([hedges_g(Z.loc[ma, m].values, Z.loc[mb, m].values) for m in mets])
    p = np.array([stats.ttest_ind(Z.loc[ma, m].values, Z.loc[mb, m].values,
                                  equal_var=False).pvalue for m in mets])
    return g, bh(p)


def main():
    os.makedirs(OUT, exist_ok=True)
    Z, grp, mets = load()
    dkd = grp == 'DKD'
    pooled = np.isin(grp, ['IgAN', 'MN', 'HN'])

    g0, q0 = contrast(Z, mets, dkd, pooled)
    base = {m for m, q in zip(mets, q0) if q < 0.05}
    log('기준선: DKD vs 다른 신장병(통합) 통과 %d개' % len(base))
    log('')

    # ------------------------------------------------------ 반론 1
    log('=' * 74)
    log('반론 1 — 세 질환을 묶은 것이 결과를 만든 것 아닌가')
    log('=' * 74)
    per = {}
    rows = []
    for dis in ('IgAN', 'MN', 'HN'):
        gg, qq = contrast(Z, mets, dkd, grp == dis)
        s = {m for m, q in zip(mets, qq) if q < 0.05}
        per[dis] = s
        rows.append(dict(comparator=dis, n=int((grp == dis).sum()), n_sig=len(s),
                         overlap_with_pooled=len(s & base)))
        log('  DKD vs %-5s (n=%3d) : 통과 %3d개, 통합 결과와 겹침 %3d개'
            % (dis, (grp == dis).sum(), len(s), len(s & base)))

    in_all = base & per['IgAN'] & per['MN'] & per['HN']
    in_two = {m for m in base if sum(m in per[d] for d in per) >= 2}
    in_one = {m for m in base if sum(m in per[d] for d in per) <= 1}
    log('')
    log('  통합 통과 %d개 가운데' % len(base))
    log('    세 질환 모두에서 통과 : %2d개 (%.0f%%)' % (len(in_all), 100 * len(in_all) / len(base)))
    log('    둘 이상에서 통과      : %2d개 (%.0f%%)' % (len(in_two), 100 * len(in_two) / len(base)))
    log('    한 질환 이하에서만    : %2d개 (%.0f%%)' % (len(in_one), 100 * len(in_one) / len(base)))
    log('')
    if len(in_one) / max(len(base), 1) > 0.5:
        log('  판정: 통합이 결과를 만들었습니다. "다른 신장병 전반"이라 말할 수 없습니다.')
    elif len(in_all) / max(len(base), 1) >= 0.4:
        log('  판정: 과반이 세 질환 모두에서 재현됩니다. 통합은 결과를 만들지 않았습니다.')
    else:
        log('  판정: 부분적으로만 재현됩니다. 주장을 통합 결과로만 서술해야 합니다.')
        log('  (HN은 n=24로 검정력이 낮아 단독 비교에서 불리합니다.)')

    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'robust_per_disease.tsv'),
                              sep='\t', index=False)

    # ------------------------------------------------------ 반론 2
    log('')
    log('=' * 74)
    log('반론 2 — 신기능 중증도의 대리 아닌가')
    log('=' * 74)
    if SEVERITY not in mets:
        log('  %s가 패널에 없습니다. 이 검증을 할 수 없습니다.' % SEVERITY)
        return 1

    sev = Z[SEVERITY].values
    log('  중증도 대리변수: 혈장 %s (같은 패널에서 측정됨)' % SEVERITY)
    log('  군별 평균 z: ' + ' · '.join(
        '%s %+.2f' % (d, sev[grp == d].mean()) for d in ('DKD', 'IgAN', 'MN', 'HN', 'Normal')))

    tt = stats.ttest_ind(sev[dkd], sev[pooled], equal_var=False)
    log('  DKD vs 다른 신장병의 %s 차이: g = %+.2f, p = %.1e'
        % (SEVERITY, hedges_g(sev[dkd], sev[pooled]), tt.pvalue))
    if tt.pvalue > 0.05:
        log('  -> 두 군의 신기능이 통계적으로 다르지 않습니다. 중증도 교란의 여지가 작습니다.')
    else:
        log('  -> 두 군의 신기능이 다릅니다. 아래 잔차화가 필요합니다.')

    # residualise every metabolite on creatinine, then repeat the contrast
    keep = [m for m in mets if m != SEVERITY]
    R = Z[keep].copy()
    s = (sev - sev.mean()) / (sev.std(ddof=1) + 1e-12)
    for m in keep:
        v = Z[m].values
        beta = float(np.dot(s, v) / np.dot(s, s))
        R[m] = v - beta * s

    g1, q1 = contrast(R, keep, dkd, pooled)
    after = {m for m, q in zip(keep, q1) if q < 0.05}
    before = base - {SEVERITY}
    log('')
    log('  잔차화 전 통과: %d개 (크레아티닌 제외)' % len(before))
    log('  잔차화 후 통과: %d개' % len(after))
    log('  유지율        : %.0f%%' % (100 * len(after & before) / max(len(before), 1)))
    log('')
    surv = sorted(after & before)
    log('  중증도 보정 후에도 남은 대사물질 %d개:' % len(surv))
    gmap = dict(zip(keep, g1))
    for m in sorted(surv, key=lambda x: -abs(gmap[x]))[:12]:
        log('    %-40s g = %+.2f' % (m[:38], gmap[m]))

    pd.DataFrame({'metabolite': keep, 'hedges_g_resid': g1, 'q_resid': q1,
                  'sig_before': [m in before for m in keep],
                  'sig_after': [m in after for m in keep]}).to_csv(
        os.path.join(OUT, 'robust_severity_adjusted.tsv'), sep='\t', index=False)

    log('')
    log('=' * 74)
    log('종합')
    log('=' * 74)
    frac_all = len(in_all) / max(len(base), 1)
    frac_keep = len(after & before) / max(len(before), 1)
    log('  질환별 재현 (세 질환 모두): %.0f%%' % (100 * frac_all))
    log('  중증도 보정 후 유지       : %.0f%%' % (100 * frac_keep))
    if frac_all >= 0.4 and frac_keep >= 0.5:
        log('  -> 두 반론 모두 결과를 뒤집지 못했습니다. 다만 이는 특이성의 증거이지')
        log('     전사체의 "재현성 vs 타당성" 분리를 재현한 것은 아닙니다. 그 주장은')
        log('     여러 코호트가 있어야 검증할 수 있고, 이 자료는 단일 코호트입니다.')
    else:
        log('  -> 반론이 부분적으로 성립합니다. 주장 수위를 낮춰야 합니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
