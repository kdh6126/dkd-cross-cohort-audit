#!/usr/bin/env python
"""전사체에서 찾은 후보를 대사체에서 독립 검증할 수 있는가.

가장 자연스러운 발상은 "전사체에서 찾은 마커를 대사체 데이터에서 다시 찾아본다"입니다.
분자 수준에서는 불가능합니다. 유전자와 대사물질은 다른 종류의 분자이고, 이름이 겹치는 것이
하나도 없습니다(후보 30개 vs 측정 물질 187개, 교집합 0).

가능한 것은 경로 수준입니다. 어떤 유전자가 만드는 단백질이 어떤 대사물질을 만들거나
분해한다면, 유전자 쪽 변화가 대사물질 쪽 변화로 이어져야 합니다. 그 연결이 실제로 존재하는
지점을 찾아 검증합니다.

이 두 자료 사이에 존재하는 연결은 하나입니다 — 콜라겐 대사.

    전사체 후보  COL1A2(콜라겐 자체), LUM·FMOD·THBS2·VCAN(콜라겐 결합 단백질),
                 MMP2·MMP7(콜라겐 분해 효소)
    대사체 측정  4-Hydroxyproline. 콜라겐에만 들어 있는 아미노산이라, 혈장 농도는
                 콜라겐이 얼마나 분해되고 있는지를 반영합니다.

예측은 명확합니다. 전사체 후보가 콜라겐 축적·재편을 가리킨다면, 같은 병에서 혈장
하이드록시프롤린이 올라가 있어야 합니다. 그리고 그것이 DKD 특이적이라면 다른 신장병에서는
올라가 있지 않아야 합니다.

서로 다른 환자·분자·장비이므로 정보는 독립적입니다. 다만 경로를 사후에 골랐고 코호트가
하나뿐이므로, '확증'이 아니라 '정합적 보강'으로만 읽어야 합니다.
"""
import os
import sys

import pandas as pd

OUT = 'results/metabolomics'

COLLAGEN_GENES = ['COL1A2', 'LUM', 'FMOD', 'THBS2', 'VCAN', 'MMP2', 'MMP7']
MARKER = '4-Hydroxyproline'
CONTROL = 'cis-4-Hydroxyproline'   # 콜라겐 유래가 아닌 이성질체. 음성 대조로 쓴다.

CONTRASTS = [
    ('DKD vs 정상', 'ST003255_DKD_vs_healthy.tsv'),
    ('DKD vs 다른 신장병', 'ST003255_DKD_vs_other_kidney_disease.tsv'),
    ('다른 신장병 vs 정상', 'ST003255_other_kidney_disease_vs_healthy.tsv'),
]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(f, met):
    d = pd.read_csv(os.path.join(OUT, f), sep='\t')
    r = d[d['metabolite'] == met]
    return (float(r['hedges_g'].iloc[0]), float(r['q'].iloc[0])) if len(r) else (None, None)


def per_disease(marker):
    """통합 대조군이 결과를 만들지 않았는지, 질환별로 따로 확인한다."""
    import numpy as np
    from scipy import stats
    d = pd.read_csv('data/raw/metabolomics/ST003255/AN005337_datatable.tsv', sep='\t')
    d = d.rename(columns={d.columns[0]: 's', d.columns[1]: 'g'})
    d['g'] = d['g'].str.replace('Disease:', '', regex=False)
    X = d.drop(columns=['s', 'g']).apply(pd.to_numeric, errors='coerce')
    X = X.loc[:, X.notna().mean() >= 0.8].fillna(0)
    X = X.replace(0, np.nan).fillna(X[X > 0].min().min())
    X = np.log2(X)
    Z = (X - X.mean()) / (X.std(ddof=1) + 1e-12)
    g = d['g'].values
    rows = []
    for dis in ('IgAN', 'MN', 'HN'):
        a, b = Z.loc[g == 'DKD', marker].values, Z.loc[g == dis, marker].values
        na, nb = len(a), len(b)
        sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
        gg = (a.mean() - b.mean()) / (sp + 1e-12) * (1 - 3 / (4 * (na + nb) - 9))
        rows.append(dict(comparator=dis, n=nb, hedges_g=gg,
                         p=stats.ttest_ind(a, b, equal_var=False).pvalue))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'crossomics_per_disease.tsv'), sep='\t', index=False)
    return t


def main():
    cand = pd.read_csv('results/candidates_v2/master_candidate_table.tsv', sep='\t')
    genes = set(cand[cand.columns[0]].astype(str))
    present = [g for g in COLLAGEN_GENES if g in genes]

    log('=' * 74)
    log('연결 지점: 콜라겐 대사')
    log('=' * 74)
    log('  전사체 후보 30개 중 콜라겐 관련: %d개 -> %s' % (len(present), ', '.join(present)))
    log('  대사체에서 대응하는 물질: %s (콜라겐 분해 산물)' % MARKER)
    log('  음성 대조: %s (콜라겐 유래가 아닌 이성질체)' % CONTROL)
    log('')

    log('=' * 74)
    log('검증 결과')
    log('=' * 74)
    log('  %-22s %-24s %8s %10s %s' % ('대조', '물질', 'g', 'q', '판정'))
    rows = []
    for label, f in CONTRASTS:
        for met, kind in ((MARKER, '표적'), (CONTROL, '음성대조')):
            g, q = get(f, met)
            if g is None:
                continue
            sig = q < 0.05
            rows.append(dict(contrast=label, metabolite=met, kind=kind,
                             hedges_g=g, q=q, significant=sig))
            log('  %-22s %-24s %+8.2f %10.1e %s'
                % (label, met, g, q, '유의' if sig else '-'))

    log('')
    sev = pd.read_csv(os.path.join(OUT, 'robust_severity_adjusted.tsv'), sep='\t')
    r = sev[sev['metabolite'] == MARKER]
    if len(r):
        log('  신기능(크레아티닌) 보정 후: g = %+.2f, q = %.1e, 유지 = %s'
            % (r['hedges_g_resid'].iloc[0], r['q_resid'].iloc[0], bool(r['sig_after'].iloc[0])))

    # 통합 대조군이 결과를 만들었던 전례가 있으므로, 같은 검증을 이 지표에도 적용한다.
    log('')
    log('  질환별로 따로 비교 (통합이 결과를 만들지 않았는지 확인):')
    pdz = per_disease(MARKER)
    for _, r_ in pdz.iterrows():
        log('    DKD vs %-6s (n=%3d)  g = %+.2f  p = %.1e  %s'
            % (r_['comparator'], r_['n'], r_['hedges_g'], r_['p'],
               '유의' if r_['p'] < 0.05 else '-'))
    log('    -> %d개 질환 중 %d개에서 유의' % (len(pdz), int((pdz['p'] < 0.05).sum())))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'crossomics_collagen.tsv'),
                              sep='\t', index=False)

    # ---------------------------------------------------------------- 판정
    g_dkd_h, q_dkd_h = get(CONTRASTS[0][1], MARKER)
    g_dkd_o, q_dkd_o = get(CONTRASTS[1][1], MARKER)
    g_oth_h, q_oth_h = get(CONTRASTS[2][1], MARKER)
    g_ctl_o, q_ctl_o = get(CONTRASTS[1][1], CONTROL)

    log('')
    log('=' * 74)
    log('해석')
    log('=' * 74)
    ok_dir = g_dkd_h > 0 and g_dkd_o > 0
    ok_spec = q_dkd_o < 0.05 and q_oth_h >= 0.05
    ok_ctrl = q_ctl_o >= 0.05
    log('  방향이 예측과 맞는가 (DKD에서 증가)          : %s' % ('예' if ok_dir else '아니오'))
    log('  DKD 특이적인가 (다른 신장병에서는 증가 없음)  : %s' % ('예' if ok_spec else '아니오'))
    log('  음성 대조가 조용한가                          : %s' % ('예' if ok_ctrl else '아니오'))
    log('')
    if ok_dir and ok_spec and ok_ctrl:
        log('  세 조건이 맞습니다. 전사체 후보가 가리키는 콜라겐 재편과 대사체 쪽 관찰이')
        log('  서로 어긋나지 않습니다. 다만 이것은 정합적 보강이지 독립 확증이 아닙니다.')
        log('')
        log('  수위를 낮춰야 하는 이유 셋:')
        log('    1. 경로를 사전에 정해 검정한 것이 아니라, 살아남은 물질 목록을 본 뒤')
        log('       콜라겐 연결을 찾았습니다. 사후 선택입니다.')
        log('    2. 코호트가 하나뿐이고, 대사체 쪽 독립 복제가 없습니다.')
        log('    3. 하이드록시프롤린은 신장 섬유화 지표로 이미 확립돼 있습니다')
        log('       (PubMed: hydroxyproline + kidney 700편 이상). 새 발견이 아닙니다.')
    else:
        log('  조건 일부가 어긋납니다. 이 연결은 주장으로 쓸 수 없습니다.')
    log('')
    log('  한계 — 이것은 대사물질 하나이고 경로 하나입니다. 전사체 후보 30개 중')
    log('  나머지 23개는 대사체 패널(카르니틴 40, 인지질 90, 아미노산 33)에 대응하는')
    log('  물질이 없어 검증 자체가 불가능합니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
