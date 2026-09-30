#!/usr/bin/env python
"""바이오마커 후보 대장 — 연구개발 기록용.

이 프로젝트의 목표는 검증된 바이오마커가 아니라 근거를 갖춘 후보를 제시하는 것입니다.
그래서 "신규성이 없으면 버린다"가 아니라 "각 후보에 어떤 근거가 몇 겹 쌓였는가"를 기록합니다.

근거는 서로 독립적인 축으로 셉니다. 같은 데이터를 다시 보는 것은 근거가 늘어난 것이 아니므로,
축마다 출처가 다릅니다.

    1 교차코호트 선택   5개 코호트 LODO에서 반복 선택됨            (전사체, 발굴)
    8 진단 일관성       대조 질병을 묶지 않고 하나씩 봐도 통과      (전사체, ERCB)
    2 질병 특이성       DKD vs 다른 신장병 효과가 순열 상한 초과   (전사체, ERCB)
    3 발현 패턴         DKD 우세 패턴이며 기증자 조직 탓이 아님    (전사체, 진단군별)
    4 세포 수준         KPMP 단일세포에서 검출·농축               (단일세포, 다른 환자)
    5 단백질 수준       KPMP 지역 프로테오믹스에서 확인           (단백질, 다른 환자)
    6 대사 경로         ST003255 혈장 대사체와 정합               (대사체, 다른 환자)
    7 문헌 신규성       DKD 문헌에 보고가 없거나 희소             (외부)

축 7은 근거가 아니라 "새로움"이므로 근거 점수에서 분리해 따로 표시합니다. 근거가 두꺼운데
문헌이 없으면 후속 연구 가치가 높고, 근거가 얇은데 문헌도 없으면 그냥 미확인입니다.
"""
import os
import sys

import numpy as np
import pandas as pd

SRC = 'results/candidates_v2/master_candidate_table.tsv'
OUT = 'results/candidates_v2'

# 대사체 패널에서 검증 가능한 경로. 지금은 콜라겐 회전율 하나뿐이다.
COLLAGEN_GENES = {'COL1A2', 'LUM', 'FMOD', 'THBS2', 'VCAN', 'MMP2', 'MMP7'}
COLLAGEN_NOTE = '혈장 4-하이드록시프롤린 상승과 정합 (DKD vs 다른 신장병 g=+0.72)'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    d = pd.read_csv(SRC, sep='\t')
    g = d.columns[0]

    ev = pd.DataFrame({'유전자': d[g]})
    # 축 1 — 교차코호트 선택은 이 표에 있는 것 자체가 통과를 뜻한다
    ev['교차코호트'] = True
    # 축 2 — 특이성
    ev['질병특이성'] = d['above_perm_ceiling'].astype(bool)
    # 축 3 — DKD 우세이며 기증자 조직 탓이 아님
    ev['발현패턴'] = d['pattern_class'].isin(['A', 'B']) & (~d['donor_driven'].astype(bool))
    # 축 4 — 단일세포
    ev['세포수준'] = d['sn_top_celltype'].notna()
    # 축 5 — 단백질
    ev['단백질수준'] = d['protein_glom_vs_TI_FC'].notna()
    # 축 6 — 대사경로
    ev['대사경로'] = d[g].isin(COLLAGEN_GENES)
    # 축 7 — 대조 진단을 묶지 않고 하나씩 봐도 과반에서 통과하는가.
    # 통합이 결과를 만들었을 가능성을 배제하는 축이라, 축 2와 별개로 센다.
    cons_path = 'results/specificity_sweep/candidate_consistency.tsv'
    if os.path.exists(cons_path):
        cons = pd.read_csv(cons_path, sep='	')
        half = cons['n_dx'].iloc[0] // 2 + 1
        good = set(cons[cons['n_pass'] >= half]['gene'])
        ev['진단일관성'] = d[g].isin(good)
    else:
        ev['진단일관성'] = False

    # 축 8 — ERCB 밖 독립 데이터(GSE20602)에서 DKD 쪽으로 더 크게 움직이는가.
    ext = 'results/external_validation/gse20602_nsc.tsv'
    if os.path.exists(ext):
        e = pd.read_csv(ext, sep='	', dtype={'entrez_id': str})
        e['delta'] = e['g_DKD_vs_control'].abs() - e['g_NSC_vs_control'].abs()
        thr = e.loc[~e['is_candidate'], 'delta'].quantile(0.95)   # 무작위 유전자의 95분위
        good_ids = set(e.loc[e['delta'] >= thr, 'entrez_id'])
        gs2 = pd.read_csv('data/processed/gene_space_all6.tsv', sep='	', dtype=str)
        id2sym = dict(zip(gs2['entrez_id'], gs2['symbol']))
        ev['외부독립'] = d[g].isin({id2sym.get(i) for i in good_ids})
    else:
        ev['외부독립'] = False

    axes = ['교차코호트', '질병특이성', '발현패턴', '세포수준', '단백질수준',
            '대사경로', '진단일관성', '외부독립']
    ev['근거수'] = ev[axes].sum(axis=1)

    # 신규성은 근거와 분리
    ev['문헌상태'] = d['verdict']
    ev['신규'] = d['lit_novel'].astype(bool)   # 0편인 경우만 True
    ev['DKD논문수'] = d['n_dkd']
    ev['g_vs_대조'] = d['g_DKD_vs_control'].round(2)
    ev['g_vs_다른신장병'] = d['g_DKD_vs_otherCKD'].round(2)
    ev['주요세포'] = d['sn_top_celltype'].fillna('-').astype(str).str.replace(
        r'<[^>]+>', '', regex=True).str.slice(0, 28)

    # 문헌은 0편/그 외의 이진값이 아니라 편수로 본다. 중앙값이 2.5편이므로 3편 이하를
    # 희소로 잡는다. 이진 플래그를 쓰면 논문 1편짜리(MOXD1)가 147편짜리(MMP2)와 같은
    # 칸에 들어가, 후속 가치가 가장 높은 후보가 보이지 않는다.
    SPARSE = 3

    def tier(r):
        thick = r['근거수'] >= 6
        sparse = r['DKD논문수'] <= SPARSE
        if thick and sparse:
            return 'B 근거 두꺼움 · 문헌 희소'
        if thick:
            return 'A 근거 두꺼움 · 기보고'
        if sparse:
            return 'C 근거 얇음 · 문헌 희소'
        return 'D 근거 얇음 · 기보고'

    ev['등급'] = ev.apply(tier, axis=1)
    ev['비고'] = np.where(ev['대사경로'], COLLAGEN_NOTE, '')
    ev = ev.sort_values(['근거수', '질병특이성'], ascending=False).reset_index(drop=True)
    ev.to_csv(os.path.join(OUT, 'candidate_dossier.tsv'), sep='\t', index=False)

    # ------------------------------------------------------------ 보고
    log('=' * 78)
    log('바이오마커 후보 대장 — 후보 %d개' % len(ev))
    log('=' * 78)
    log('')
    log('  근거 축별 통과 (독립적인 데이터 출처)')
    for a in axes:
        log('    %-12s %2d / %d' % (a, int(ev[a].sum()), len(ev)))
    log('')
    log('  근거 겹수 분포')
    for k, v in ev['근거수'].value_counts().sort_index(ascending=False).items():
        log('    %d겹  %2d개' % (k, v))
    log('')
    log('  등급 분포')
    for k, v in ev['등급'].value_counts().items():
        log('    %-24s %2d개' % (k, v))

    log('')
    log('=' * 78)
    log('B등급 — 근거가 두껍고 DKD 문헌이 희소한 후보 (후속 연구 우선순위)')
    log('=' * 78)
    b = ev[ev['등급'].str.startswith('B')]
    if len(b) == 0:
        log('  없음')
    for _, r in b.iterrows():
        log('  %-10s 근거 %d겹  g(대조)=%+.2f  g(다른신장병)=%+.2f  DKD논문 %d편'
            % (r['유전자'], r['근거수'], r['g_vs_대조'], r['g_vs_다른신장병'], r['DKD논문수']))
        log('             %s | %s' % (r['문헌상태'], r['주요세포']))

    log('')
    log('=' * 78)
    log('A등급 — 근거가 두껍고 이미 보고된 후보 (파이프라인 건전성의 증거)')
    log('=' * 78)
    a_ = ev[ev['등급'].str.startswith('A')]
    for _, r in a_.head(10).iterrows():
        log('  %-10s 근거 %d겹  g(다른신장병)=%+.2f  DKD논문 %d편%s'
            % (r['유전자'], r['근거수'], r['g_vs_다른신장병'], r['DKD논문수'],
               '  [대사경로 정합]' if r['대사경로'] else ''))

    log('')
    log('=' * 78)
    log('해석')
    log('=' * 78)
    log('  이 표는 "검증된 바이오마커"가 아니라 "근거가 몇 겹 쌓인 후보"의 기록입니다.')
    log('  등급 A는 파이프라인이 이미 알려진 것을 제대로 회수한다는 증거이고,')
    log('  등급 B는 같은 수준의 근거를 갖췄지만 문헌이 희소해 후속 가치가 높은 쪽입니다.')
    log('')
    log('  근거 6겹이라도 임상적으로 검증된 것은 아닙니다. 여기서 말하는 근거는 전부')
    log('  공개 데이터에 대한 계산이며, 전향적 코호트나 실험 검증은 포함돼 있지 않습니다.')
    log('')
    log('  wrote %s/candidate_dossier.tsv' % OUT)


if __name__ == '__main__':
    main()
