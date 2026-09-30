#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""공개 자료 카탈로그를 브라우저에서 본다.

    python -m streamlit run scripts/omics_catalog_app.py

이 화면이 답해야 하는 질문은 하나입니다. **같은 참여자에서 두 오믹스 층을 잰 비교를
만들 수 있는가.** 그 답이 첫 화면에 나와야 하고, 왜 그런지를 눌러서 내려갈 수 있어야
합니다. 목록을 예쁘게 보여 주는 것이 목적이 아닙니다.

판정 규칙은 하나입니다. 한 층 쌍에 대해 **사례와 대조 양쪽에 사람이 있어야** 합니다.
어느 한쪽이 0명이면 짝은 존재해도 비교가 성립하지 않습니다. 이 프로젝트가 두 번 틀린
지점이 정확히 여기이고, 그래서 사람이 세지 않고 뷰가 세게 둡니다.

숫자는 전부 db/omics_catalog.sqlite 에서 읽습니다. 화면에서 계산하지 않습니다.
"""
import os
import sqlite3
import sys

import pandas as pd
import streamlit as st

DB = 'db/omics_catalog.sqlite'
LAYERS = ['transcriptomics', 'proteomics', 'metabolomics']


@st.cache_data
def load(sql):
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    con = sqlite3.connect(os.path.join(root, DB))
    try:
        return pd.read_sql(sql, con)
    finally:
        con.close()


LOGO = 'docs/assets/etri_logo.png'


def header():
    """제목 줄. 로고 파일이 있으면 왼쪽에 놓는다.

    로고를 코드로 그리지 않습니다. 기관 마크를 손으로 흉내 내면 비슷하지만 틀린 것이
    되고, 그건 없느니만 못합니다. 파일을 그 자리에 놓으면 나타납니다.
    """
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, LOGO)
    if os.path.exists(path):
        st.image(path, width=260)
    st.title('DKD 데이터베이스')


def grid_block(sub):
    """참여자 x 오믹스 층 격자.

    표로 세면 "대조군에 단백체가 0명" 이 숫자 하나지만, 격자로 그리면 그 열이 통째로
    비어 있는 것이 보입니다. 이 프로젝트가 두 번 틀린 사실을 눈으로 확인하는 그림입니다.

    행은 참여자, 열은 층입니다. 군으로 묶어 정렬하므로 사례 구간과 대조 구간이 위아래로
    갈립니다. 층이 하나뿐인 참여자가 500명 가까이 있어 전부 그리면 격자가 아니라 띠가
    되므로, 기본은 두 층 이상 가진 사람과 그 자원의 대조군만 그립니다.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    from matplotlib.patches import FancyBboxPatch

    # 한글 라벨이 네모로 깨지지 않게 설치된 글꼴을 찾는다. 없으면 라벨을 영어로 바꾼다.
    names = {f.name for f in fm.fontManager.ttflist}
    korean = next((n for n in ('Malgun Gothic', 'NanumGothic', 'AppleGothic',
                               'Noto Sans KR', 'Gulim') if n in names), None)
    if korean:
        plt.rcParams['font.family'] = korean
    plt.rcParams['axes.unicode_minus'] = False
    lay_lab = ['Transcriptome', 'Proteome', 'Metabolome']

    st.subheader('참여자 x 오믹스 층 격자')
    res = st.selectbox('자원', sorted(sub['resource'].dropna().unique()),
                       index=sorted(sub['resource'].dropna().unique()).index('KPMP')
                       if 'KPMP' in set(sub['resource']) else 0)
    v = sub[sub['resource'] == res].copy()
    if st.checkbox('층이 하나뿐인 참여자도 포함', value=False):
        pass
    else:
        v = v[v['n_layers'] > 1]
    if v.empty:
        have = [lab for lab, col in zip(LAYERS, LAYERS)
                if int(sub[sub['resource'] == res][col].sum()) > 0]
        st.info('이 자원에는 두 층 이상 가진 참여자가 없습니다. 가진 층은 %s 하나뿐이라 '
                '같은 참여자에서 층을 잇는 비교를 만들 수 없습니다. 우리 전사체 코호트'
                '(GEO)가 전부 여기에 해당합니다 — **DKD 자료는 있지만 층이 하나입니다.**'
                % ', '.join(have))
        return

    order = {'case': 0, 'comparator': 1, 'control': 2, 'unstated': 3}
    v['o'] = v['arm'].map(order).fillna(9)
    v = v.sort_values(['o', 'disease_class', 'participant'])
    # 군만 적으면 '사례' 안에 어떤 병이 있는지 보이지 않는다. KPMP 의 사례는 CKD 와
    # AKI 이고 DKD 는 한 명도 없는데, 군만 보면 그 사실이 가려진다.
    v['block'] = v['arm'] + ' · ' + v['disease_class']
    M = v[LAYERS].to_numpy(dtype=float)
    n = len(v)

    # 칸을 띄워 그리면 하나하나가 '측정 하나' 로 읽힌다.
    # 참여자가 많으면 한 명 한 명을 줄로 그릴 수 없다. 80줄이 되면 그림이 좁고 길어져
    # 열 이름과 참여자 번호가 서로 겹치고, 무엇보다 전사체만 가진 줄이 화면을 덮어
    # 정작 봐야 할 '어느 군에서 어느 층이 비었는가' 가 묻힌다. 그럴 때는 질병 블록
    # 하나를 한 줄로 요약하고, 칸의 진하기로 그 블록에서 몇 퍼센트가 그 층을 가졌는지
    # 보여 준다. 숫자도 칸 안에 적어 어림하지 않게 한다.
    per_participant = n <= 45
    if per_participant:
        rows_lab = list(v['participant'])
        frac = M
        note = None
    else:
        gb = v.groupby('block', sort=False)
        frac = gb[LAYERS].mean().to_numpy()
        cnt = gb[LAYERS].sum().astype(int).to_numpy()
        size = gb.size().to_numpy()
        rows_lab = ['%s (%d)' % (b, s_) for b, s_ in zip(gb.groups.keys(), size)]
        note = (cnt, size)
    nr = len(rows_lab)

    fig, ax = plt.subplots(figsize=(4.9, 0.30 * nr + 1.5))
    for y in range(nr):
        for x in range(len(LAYERS)):
            f = float(frac[y, x])
            ax.add_patch(FancyBboxPatch(
                (x - 0.40, y - 0.33), 0.80, 0.66,
                boxstyle='round,pad=0,rounding_size=0.14',
                facecolor=('#f4f1ed' if f == 0 else
                           (0.757 + 0.243 * (1 - f), 0.267 + 0.733 * (1 - f),
                            0.055 + 0.945 * (1 - f))),
                edgecolor=('#e4dfd8' if f == 0 else '#c1440e'), linewidth=0.7))
            if note is not None:
                ax.text(x, y, '%d/%d' % (note[0][y, x], note[1][y]),
                        va='center', ha='center', fontsize=6.4,
                        color=('#8a8a8a' if f == 0 else
                               ('white' if f > 0.55 else '#1a1a1a')))
    ax.set_xlim(-0.62, len(LAYERS) - 0.38)
    ax.set_ylim(nr - 0.42, -0.62)
    ax.set_xticks(range(len(LAYERS)))
    ax.set_xticklabels(lay_lab, fontsize=8.5, rotation=22, ha='left',
                       rotation_mode='anchor')
    ax.xaxis.set_ticks_position('top')
    ax.set_yticks(range(nr))
    ax.set_yticklabels(rows_lab, fontsize=(6.4 if per_participant else 7.4),
                       color=('#8a8a8a' if per_participant else '#1a1a1a'))
    ax.tick_params(length=0)

    if per_participant:
        # 참여자 단위로 그릴 때만 군 경계를 긋는다. 요약 그림은 줄 자체가 블록이다.
        arms = list(v['arm'])
        blocks = list(v['block'])
        start = 0
        for i in range(1, nr + 1):
            if i == nr or blocks[i] != blocks[start]:
                ax.text(-2.15, (start + i - 1) / 2,
                        '%s (%d)' % (blocks[start], i - start),
                        rotation=90, va='center', ha='center', fontsize=7.0,
                        color=('#c1440e' if arms[start] == 'case' else '#1a1a1a'))
                if i < nr:
                    heavy = arms[i] != arms[start]
                    ax.plot([-0.62, len(LAYERS) - 0.38], [i - 0.5, i - 0.5],
                            color=('#8a8a8a' if heavy else '#e4dfd8'),
                            lw=(1.2 if heavy else 0.8), clip_on=False)
                start = i
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.tight_layout()
    st.pyplot(fig, width='content')

    # 열마다 군별로 몇 명인지. 그림 옆에 숫자가 있어야 읽는 사람이 확인할 수 있다.
    t = v.groupby(['arm', 'disease_class'])[LAYERS].sum().astype(int)
    t.columns = lay_lab
    t['참여자'] = v.groupby(['arm', 'disease_class']).size()
    st.dataframe(t, width='stretch')
    st.caption('질병별로 봅니다. 어떤 병에서 한 열이 0이면 그 병에서는 그 층을 쓰는 '
               '비교가 성립하지 않습니다. **목록에 없는 병은 그 자원에 참여자가 아예 '
               '없다는 뜻입니다** — KPMP 에는 당뇨병성 신장질환 참여자가 없습니다.')


# 근거 축. (표시 이름, 판정 함수, 근거가 되는 원본 열) 로 둔다. 판정과 표시를 한 곳에
# 모아 두면 표와 상세 화면이 어긋날 수 없다.
def _bool(v):
    """1/0/결측을 충족·미충족·평가불가로 옮긴다."""
    if v is None or (isinstance(v, float) and v != v):
        return None
    return bool(v)


EVIDENCE = [
    ('순열 상한 초과', lambda r: _bool(r['above_perm_ceiling']),
     lambda r: 'g(DKD 대 다른 CKD) = %+.2f' % r['g_DKD_vs_otherCKD']),
    ('진단 전반에서 일관', lambda r: (None if not isinstance(r['pattern_class'], str)
                                     else r['pattern_class'] in ('A', 'B')),
     lambda r: '패턴 %s' % (r['pattern_class'] or '미상')),
    ('대조군이 끌지 않음', lambda r: (None if r['donor_driven'] is None
                                     else not _bool(r['donor_driven'])),
     lambda r: '공여자 주도 여부'),
    ('단일세포 지지', lambda r: (None if r['sn_top_foldchange'] != r['sn_top_foldchange']
                                else abs(r['sn_top_foldchange']) > 1.5),
     lambda r: ('%s 에서 %.2f배' % (r['sn_top_celltype'], r['sn_top_foldchange'])
                if r['sn_top_foldchange'] == r['sn_top_foldchange'] else '조회 안 됨')),
    ('신장 부위 특이성 (질병 아님)',
     lambda r: (None if r['protein_adjP'] != r['protein_adjP']
                else r['protein_adjP'] < 0.05),
     lambda r: ('KPMP 사구체 대 세뇨관간질 단백질 q=%.3g — 부위 대비이지 '
                'DKD 대비가 아님' % r['protein_adjP']
                if r['protein_adjP'] == r['protein_adjP'] else 'KPMP 패널에 없음')),
    ('DKD 문헌에 없음', lambda r: _bool(r['lit_novel']),
     lambda r: 'DKD 논문 %d편 · 마커 제안 %d편' % (r['n_dkd'], r['n_dkd_biomarker'])),
    ('질병 대비 단백체', lambda r: (None if (r['on_somascan'] != r['on_somascan']
                                            and r['on_massspec'] != r['on_massspec'])
                                    else bool(r['sig_somascan'] or r['sig_massspec'])),
     lambda r: _protein_note(r)),
]


def _protein_note(r):
    if r['on_somascan'] != r['on_somascan'] and r['on_massspec'] != r['on_massspec']:
        return '두 코호트 모두 재지 않음'
    bits = []
    for lab, on, sig in (('SOMAscan', 'on_somascan', 'sig_somascan'),
                         ('질량분석', 'on_massspec', 'sig_massspec')):
        if r[on] == r[on] and r[on]:
            bits.append('%s %s' % (lab, '유의' if r[sig] else '유의하지 않음'))
    return ' · '.join(bits) or '두 코호트 모두 재지 않음'


MARK = {True: '\u25cf 충족', False: '\u25cb 미충족', None: '\u2014 평가불가'}


def confounding_block():
    """질병이 다른 변수와 완전히 묶여 있는가.

    사례 100명이 전부 생검이고 대조 100명이 전부 신절제이면, 어느 층에도 양쪽이 없으므로
    질병 효과와 조달 효과를 나눌 정보가 자료에 없습니다. 이때는 사후 보정도 소용이
    없고, AUC 가 높아도 무엇을 맞힌 것인지 알 수 없습니다.

    반대로 한 층 안에 양쪽이 있으면 그 층 안에서 비교할 수 있습니다. 그래서 세는 것은
    '층의 개수' 가 아니라 **사례와 대조가 함께 있는 층의 개수** 입니다.
    """
    st.subheader('질병과 다른 변수가 분리되는가')
    t = load('select * from v_confounding_summary order by axis, cohort')
    if t.empty:
        st.info('시료 수준 메타데이터가 없습니다.')
        return

    for axis, g in t.groupby('axis'):
        bad = int((g['separable_strata'] == 0).sum())
        head = '%s — 코호트 %d개 중 %d개가 완전 교란' % (axis, len(g), bad)
        with st.expander(head, expanded=(bad > 0)):
            v = g.rename(columns={'n_strata': '층 수',
                                  'separable_strata': '사례·대조가 함께 있는 층',
                                  'n_case': '사례', 'n_control': '대조'})
            v['판정'] = ['완전 교란' if x == 0 else '분리 가능'
                         for x in g['separable_strata']]
            st.dataframe(v[['cohort', '층 수', '사례·대조가 함께 있는 층',
                            '사례', '대조', '판정']],
                         width='stretch', hide_index=True)
            det = load("select * from v_confounding where axis='%s' "
                       "order by cohort, stratum" % axis)
            st.caption('층별 내역')
            st.dataframe(det.rename(columns={'stratum': '층', 'n_case': '사례',
                                             'n_control': '대조'}),
                         width='stretch', hide_index=True)

    cov = load('select * from v_kpmp_covariate')
    if not cov.empty:
        with st.expander('KPMP 안의 탐색적 공변량 — 성별과 연령대', expanded=False):
            st.caption('이 두 값은 **KPMP 에만** 있습니다. GEO 코호트는 시료 수준으로 '
                       '배포하지 않아 같은 표에 합칠 수 없습니다. 그래서 핵심 교란 보정 '
                       '변수로 올리지 않고 KPMP 안의 탐색용으로만 둡니다. 사례는 등록 '
                       '진단이 아니라 **당뇨 병력으로 고른 조작적 대리 부분집합**입니다.')
            for ax, g in cov.groupby('axis'):
                both = int(((g['n_case'] > 0) & (g['n_control'] > 0)).sum())
                st.write('**%s** — 사례·대조가 함께 있는 층 %d / %d'
                         % (ax, both, len(g)))
                st.dataframe(g.rename(columns={'stratum': '층', 'n_case': '사례',
                                               'n_control': '대조'})[['층', '사례', '대조']],
                             width='stretch', hide_index=True)
            st.caption('성별은 두 층 모두에 사례와 대조가 함께 있어, 이 분포에서는 '
                       '질병과 겹쳐 있지 않습니다 — 그래서 보정 변수로 쓰지 '
                       '않았습니다. 연령대는 다릅니다. 16개 구간 중 11개가 대조군에만 '
                       '있어 **겹침(overlap) 자체가 성립하지 않습니다**. 그 구간에는 '
                       '비교할 상대가 없으므로, "연령은 교란이 아니다" 가 아니라 '
                       '**이 자료에서는 연령을 안정적으로 분리해 보정할 수 없다**가 '
                       '맞는 서술입니다. KPMP 안에서의 관찰이고 다른 코호트로 옮길 수 '
                       '있는 결론도 아닙니다.')

    st.warning('조달(procurement)에서 여섯 코호트가 모두 완전 교란입니다. 사례는 전부 '
               '생검이고 대조는 전부 공여자나 종양신절제라, 어느 층에도 양쪽이 없습니다. '
               '이 구조에서는 높은 AUC 가 질병을 맞힌 것인지 채취 방법을 맞힌 것인지 '
               '구분되지 않습니다. 구획과 플랫폼은 같은 층 안에 양쪽이 있어 분리됩니다 '
               '— 문제는 이 자료 전체가 아니라 **조달 축 하나**입니다.')
    st.caption('나이 · 성별 · 병원 · 시퀀싱 배치도 같은 원리로 봐야 하지만, 우리가 쓴 '
               '공개 자료는 그 값을 시료 수준으로 배포하지 않아 이 표에 넣을 수 '
               '없었습니다. 없는 것을 채워 넣지 않고 없다고 적어 둡니다.')


def model_block():
    """교란을 보정하면 무엇이 달라지는가.

    좋은 보정은 두 가지를 동시에 합니다 — 교란에서 온 유전자 비율을 낮추고, 예측력은
    지킵니다. 둘 중 하나만 보면 판단이 틀립니다. 비율만 보면 아무것도 못 고르는 보정이
    최고가 되고, 예측력만 보면 교란을 그대로 둔 것이 최고가 됩니다.

    아무 보정도 하지 않은 arm 과 효과가 없는 보정도 함께 보여 줍니다. '보정을 했다' 가
    아니라 '어떤 보정이 들었다' 가 결론이기 때문입니다.
    """
    st.subheader('교란 보정에 따른 결과')
    k = st.select_slider('선택 유전자 수 K', [5, 10, 20, 50, 100], value=50)
    d = load('select * from v_correction_effect where K=%d order by procurement_frac' % k)
    if d.empty:
        st.info('이 K 에 대한 결과가 없습니다.')
        return

    base = d[d['method'] == 'none']
    b_auroc = float(base['auroc'].iloc[0]) if len(base) else float('nan')
    b_proc = float(base['procurement_frac'].iloc[0]) if len(base) else float('nan')

    v = d.copy()
    v['AUROC 변화'] = (v['auroc'] - b_auroc).round(3)
    v['교란 유전자 변화'] = (v['procurement_frac'] - b_proc).round(3)
    v['판정'] = ['보정 안 함' if m == 'none'
                 else ('효과 있음' if p <= b_proc - 0.10 else
                       ('효과 없음' if p >= b_proc - 0.02 else '부분적'))
                 for m, p in zip(v['method'], v['procurement_frac'])]
    st.dataframe(
        v.rename(columns={'method': '방법', 'auroc': 'AUROC',
                          'procurement_frac': '교란 유래 유전자 비율',
                          'reproducibility': '재현성(교차폴드 자카드)'})
        [['방법', 'AUROC', 'AUROC 변화', '교란 유래 유전자 비율',
          '교란 유전자 변화', '재현성(교차폴드 자카드)', '판정']],
        width='stretch', hide_index=True)
    st.caption('교란 유래 유전자 비율은 선택된 유전자 중 조달 신호를 따라가는 것의 '
               '비율입니다. 낮을수록 좋습니다. **AUROC 를 지키면서 이 비율을 낮춘 것만 '
               '효과가 있습니다** — ComBat 과 경험적 RUVg 는 둘 다 거의 움직이지 '
               '않았습니다.')

    st.divider()
    st.subheader('AUROC 절대값을 믿을 수 있는가')
    nl = load('select * from v_null_auroc order by K')
    st.dataframe(nl.rename(columns={'random_mean': '무작위 유전자 평균',
                                    'random_p05': '5백분위', 'random_p95': '95백분위',
                                    'permuted_mean': '라벨 섞음 평균'}),
                 width='stretch', hide_index=True)
    st.warning('질병과 무관한 **무작위 유전자 %d개로도 AUROC 가 평균 %.2f** 나옵니다'
               '(K=%d). 라벨을 섞으면 0.50 으로 떨어지므로 모형이 고장난 것은 아니고, '
               '이 자료에서는 거의 어떤 유전자 집합으로도 사례와 대조가 갈린다는 '
               '뜻입니다. 조달이 완전히 교란되어 있기 때문입니다 — 교란 구조 탭을 '
               '보세요.'
               % (int(nl['K'].iloc[-1]),
                  float(nl['무작위 유전자 평균'].iloc[-1])
                  if '무작위 유전자 평균' in nl else 0,
                  int(nl['K'].iloc[-1])))

    st.divider()
    st.subheader('선택 방법 비교')
    mm = load('select method, K, external, internal, cross_fold_jaccard, null_pct '
              'from model_method where K=%d order by external desc' % k)
    st.dataframe(mm.rename(columns={'method': '방법', 'external': '외부 AUROC',
                                    'internal': '내부 AUROC',
                                    'cross_fold_jaccard': '교차폴드 자카드',
                                    'null_pct': '귀무 분포 백분위'}),
                 width='stretch', hide_index=True)
    st.caption('귀무 분포 백분위는 같은 크기의 무작위 유전자 집합과 견준 위치입니다. '
               '50 근처면 무작위와 다르지 않다는 뜻입니다.')


def deviation_block():
    """환자별 · 유전자별 정상 대비 편차."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    st.subheader('환자별 정상 대비 편차')
    st.caption('여기 실린 여섯 코호트는 **전부 당뇨병성 신장질환**입니다 — 사례가 DN 또는 '
               'DKD 이고 대조는 신절제나 공여자입니다. 격자 탭에 DKD 가 없는 것은 그쪽이 '
               'KPMP 만 다루기 때문이고, KPMP 에는 당뇨병성 신장질환 참여자가 없습니다.')
    d = load('select * from patient_deviation')
    if d.empty:
        st.info('편차 자료가 없습니다. build_omics_catalog.py 를 먼저 돌리세요.')
        return

    dis = load('select cohort, arm, disease_group, count(*) n from sample_meta '
               "where arm='case' group by 1,2,3")
    lab = {r['cohort']: r['disease_group'] for _, r in dis.iterrows()}
    c1, c2 = st.columns([1, 1])
    coh = c1.selectbox('코호트', sorted(d['cohort'].unique()),
                       format_func=lambda x: '%s  (사례 = %s)' % (x, lab.get(x, '?')))
    show_ctl = c2.checkbox('대조군도 함께', value=True)
    v = d[d['cohort'] == coh]
    v = v[v['arm'].isin(['case', 'control'] if show_ctl else ['case'])]
    if v.empty:
        st.info('이 코호트에는 그릴 시료가 없습니다.')
        return

    M = v.pivot_table(index=['arm', 'sample'], columns='gene', values='z')
    # 열은 경로(주소)로 묶어 정렬한다. 같은 경로의 유전자가 붙어 있어야 덩어리로 읽힌다.
    addr = v.drop_duplicates('gene').set_index('gene')['pathway']
    cols = sorted(M.columns, key=lambda g: (addr.get(g, '미분류'), g))
    M = M[cols]
    M = M.sort_index(level=['arm', 'sample'], ascending=[False, True])

    hi = st.slider('색 눈금 상한 (|z|)', 1.0, 8.0, 3.0, 0.5,
                   help='몇 개의 극단값이 눈금을 끌고 가면 나머지가 전부 희게 죽습니다. '
                        '보통 z=3 이면 대조군 분포의 바깥입니다.')
    lim = hi
    fig, ax = plt.subplots(figsize=(7.2, max(2.6, 0.16 * len(M) + 1.6)))
    im = ax.imshow(M.to_numpy(), aspect='auto', cmap='RdBu_r',
                   vmin=-lim, vmax=lim, interpolation='nearest')
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, fontsize=6.2, rotation=90)
    ax.set_xlabel('유전자 (후보 %d개, 경로별로 묶음)' % len(cols), fontsize=8)
    arms = [a for a, _ in M.index]
    ax.set_yticks(range(len(M)))
    ax.set_yticklabels([s_ for _, s_ in M.index], fontsize=5.4, color='#8a8a8a')
    ax.set_ylabel('환자 한 명 = 한 줄  (GSM… 은 GEO 시료 번호)', fontsize=8)
    ax.tick_params(length=0)

    # 군이 바뀌는 자리에 선. 사례 구간과 대조 구간이 눈으로 갈려야 한다.
    start = 0
    for i in range(1, len(arms) + 1):
        if i == len(arms) or arms[i] != arms[start]:
            # 군 이름을 구간 왼쪽 바깥에 한 번씩 적는다. 줄마다 GSM 번호만 있으면
            # 어디까지가 환자이고 어디부터 대조인지 알 수 없다.
            ax.text(-0.055, 1 - (start + i) / 2 / len(arms),
                    {'case': '환자', 'control': '대조'}.get(arms[start], arms[start]),
                    transform=ax.transAxes, rotation=90, va='center', ha='center',
                    fontsize=9, color=('#c1440e' if arms[start] == 'case' else '#1a1a1a'))
            if i < len(arms):
                ax.axhline(i - 0.5, color='#1a1a1a', lw=1.2)
            start = i
    # 경로가 바뀌는 자리에도 선. 이것이 '주소' 경계다.
    prev = None
    for j, g in enumerate(cols):
        a = addr.get(g, '미분류')
        if prev is not None and a != prev:
            ax.axvline(j - 0.5, color='#1a1a1a', lw=0.8)
        prev = a
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label('z  (대조군 대비)', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    fig.tight_layout()
    st.pyplot(fig, width='content')

    # 경로별로 사례가 어느 쪽으로 쏠렸는지. 그림 옆에 숫자가 있어야 확인이 된다.
    t = (v[v['arm'] == 'case'].groupby('pathway')['z']
         .agg(['count', 'mean']).round(2)
         .rename(columns={'count': '값 개수', 'mean': '사례 평균 z'}))
    t['유전자'] = v[v['arm'] == 'case'].groupby('pathway')['gene'].nunique()
    st.dataframe(t, width='stretch')
    st.caption('경로가 "주소" 입니다. 같은 경로의 유전자를 붙여 놓았고 경계에 세로선을 '
               '그었습니다. 어느 경로에도 속하지 않는 후보는 억지로 묶지 않고 미분류로 '
               '둡니다. 대조군 행이 0 근처에서 흩어져 있어야 눈금이 제 구실을 하는 '
               '것입니다 — 대조군까지 한쪽으로 쏠리면 z 를 믿을 수 없습니다.')


def biomarker_block():
    st.subheader('바이오마커 후보 30개와 그 근거')
    st.caption('이 프로젝트가 내놓은 후보입니다. 등급은 **자료 강도**와 **문헌 신규성**의 '
               '조합이지 품질 순서가 아닙니다. 재지 않은 축은 미충족이 아니라 평가불가로 '
               '표시합니다 — 검사해서 떨어진 것과 섞으면 안 됩니다.')
    bm = load('select * from biomarker order by tier, gene')

    rows = []
    for _, r in bm.iterrows():
        d = {'유전자': r['gene'], '등급': str(r['tier'])[0]}
        for name, judge, _note in EVIDENCE:
            d[name] = MARK[judge(r)]
        d['근거 충족'] = sum(1 for name, judge, _ in EVIDENCE if judge(r) is True)
        rows.append(d)
    tab = pd.DataFrame(rows)

    c1, c2 = st.columns([1, 3])
    tiers = sorted(tab['등급'].unique())
    pick = c1.multiselect('등급', tiers, default=tiers)
    only = c2.multiselect('이 근거를 충족하는 것만', [e[0] for e in EVIDENCE])
    v = tab[tab['등급'].isin(pick)]
    for o in only:
        v = v[v[o].str.startswith('\u25cf')]
    st.write('%d개' % len(v))
    st.dataframe(v.sort_values(['등급', '근거 충족'], ascending=[True, False]),
                 width='stretch', hide_index=True)
    st.caption('**신장 부위 특이성**은 KPMP 의 사구체 대 세뇨관간질 단백질 대비입니다. '
               'KPMP 에는 당뇨병성 신장질환 참여자가 없으므로 이 축은 **질병 근거가 '
               '아닙니다** — 그 단백질이 사람 신장에서 실제로 검출되고 부위에 따라 '
               '다르게 분포한다는 것까지만 말합니다. 질병 대비 단백체 근거는 오른쪽 '
               '끝의 별도 축입니다.')

    st.divider()
    st.subheader('후보 하나를 자세히')
    g = st.selectbox('유전자', list(bm['gene']))
    r = bm[bm['gene'] == g].iloc[0]
    a, b, c = st.columns(3)
    a.metric('등급', str(r['tier'])[0])
    b.metric('g (DKD 대 다른 CKD)', '%+.2f' % r['g_DKD_vs_otherCKD'])
    c.metric('보정 전 목록에도 있었나', '있음' if r['in_uncorrected_list'] else '없음')
    det = pd.DataFrame([{'근거': name, '판정': MARK[judge(r)], '값': note(r)}
                        for name, judge, note in EVIDENCE])
    st.dataframe(det, width='stretch', hide_index=True)
    if r['rank_naive'] == r['rank_naive']:
        st.caption('대조군을 바꿨을 때 순위: 비생검 대조 %d위 -> 생검 대조 %d위. '
                   '순위가 오르는 것은 조달 비대칭이 줄었을 때 우선순위가 올라간다는 '
                   '뜻이지, 마커라는 증거가 아닙니다.'
                   % (int(r['rank_naive']), int(r['rank_matched'])))


def main():
    st.set_page_config(page_title='DKD 공개 자료 카탈로그', layout='wide')
    header()

    ds = load('select * from dataset')
    sub = load('select * from subject_layer')

    a, b, c, d = st.columns(4)
    a.metric('데이터셋', len(ds))
    b.metric('참여자·시료', len(sub))
    c.metric('오믹스 종류', ds['omics'].nunique())
    d.metric('두 층 이상 보유', int((sub['n_layers'] > 1).sum()))

    st.divider()
    (tab_c, tab_m, tab_b, tab_d, tab0, tab1, tab2, tab3,
     tab4) = st.tabs(
        ['교란 구조', '모델 결과', '바이오마커 후보', '환자별 편차', '격자',
         '참여자별 층', '질환 분류별 집계', '데이터셋', '유전자 조회'])

    with tab_d:
        st.info('격자는 "쟀는가" 만 보여 줍니다. 이 화면은 그 다음 질문에 답합니다 — '
                '**이 환자는 정상에서 어느 쪽으로 얼마나 벗어났는가.** 행이 환자, 열이 '
                '후보 유전자이고, 색은 그 코호트 대조군 대비 z 입니다. 열은 경로'
                '(주소)로 묶어 정렬했고 경계에 세로선을 그었습니다.')
        deviation_block()

    with tab_m:
        st.info('교란을 보정하면 무엇이 달라지는지를 봅니다. 좋은 보정은 **예측력은 '
                '지키면서 교란에서 온 유전자 비율을 낮춥니다**. 둘 중 하나만 보면 '
                '판단이 틀립니다 — 비율만 보면 아무것도 못 고르는 보정이 최고가 되고, '
                '예측력만 보면 교란을 그대로 둔 것이 최고가 됩니다.')
        model_block()

    with tab_c:
        st.info('모형을 돌리기 전에 보는 화면입니다. 질병이 다른 변수와 **완전히 묶여 '
                '있으면** 그 둘을 나눌 정보가 자료에 없고, 사후 보정으로도 해결되지 '
                '않습니다. 세는 것은 층의 개수가 아니라 **사례와 대조가 함께 있는 층의 '
                '개수**입니다. 0이면 완전 교란입니다.')
        confounding_block()

    with tab_b:
        st.info('이 프로젝트가 내놓은 후보 30개와 그 근거입니다. 표의 각 열이 근거 축 '
                '하나이고, 충족 · 미충족 · 평가불가를 구분합니다. **재지 않은 축을 '
                '미충족으로 칠하지 않습니다** — 검사해서 떨어진 것과 애초에 못 잰 것은 '
                '다릅니다. 아래에서 유전자 하나를 고르면 축마다 실제 값을 볼 수 '
                '있습니다.')
        biomarker_block()

    with tab0:
        st.info('참여자 한 명이 한 줄, 오믹스 층이 한 칸입니다. 채워진 칸은 그 층을 '
                '실제로 잰 것이고, 빈 칸은 재지 않은 것입니다. 군(사례·비교군·대조)으로 '
                '묶어 정렬하므로, **어떤 군에서 한 층이 통째로 비어 있는지**가 바로 '
                '보입니다. 그것이 다중 오믹스가 성립하지 않는 이유입니다.')
        grid_block(sub)

    with tab1:
        st.info('참여자 한 명이 한 줄입니다. 자원과 군으로 거르고, 특정 층을 반드시 '
                '가진 사람만 볼 수 있습니다. 다중 오믹스 후보를 직접 찾으려면 '
                '"두 층 이상만" 을 켜고 군을 대조로 바꿔 보세요 — 결과가 비면 그 자원 '
                '으로는 비교를 만들 수 없다는 뜻입니다.')
        col = st.columns(4)
        res = col[0].multiselect('자원', sorted(sub['resource'].dropna().unique()),
                                 default=sorted(sub['resource'].dropna().unique()))
        arm = col[1].multiselect('군', sorted(sub['arm'].unique()),
                                 default=sorted(sub['arm'].unique()))
        want = col[2].multiselect('반드시 가진 층', LAYERS)
        only_multi = col[3].checkbox('두 층 이상만', value=False)
        v = sub[sub['resource'].isin(res) & sub['arm'].isin(arm)]
        for w in want:
            v = v[v[w] == 1]
        if only_multi:
            v = v[v['n_layers'] > 1]
        st.write('%d명' % len(v))
        st.dataframe(v, width='stretch', hide_index=True)
        st.caption('참여자 식별자는 KPMP 공개 아틀라스가 배포하는 비식별 번호입니다. '
                   '이 화면은 로컬에서만 돕니다.')

    with tab2:
        st.info('질환 분류마다 각 층을 몇 명이 가졌는지 셉니다. **n_multilayer** 가 그 '
                '분류에서 두 층 이상 가진 사람 수이고, 0이면 그 분류에서는 다중 오믹스가 '
                '성립하지 않습니다. 군은 사례 · 비교군(다른 진단의 생검) · 대조(질병이 '
                '아닌 이유로 얻은 조직) · 미기재로 나뉩니다.')
        g = load('select * from v_layer_by_class order by resource, arm, n_subjects desc')
        st.dataframe(g, width='stretch', hide_index=True)
        st.caption('층 열의 값은 그 분류에서 해당 층을 가진 사람 수입니다. '
                   'n_multilayer 가 0이면 그 분류에서는 다중 오믹스가 성립하지 않습니다.')

    with tab3:
        st.info('받아 온 자료의 목록입니다. 접근번호 · 저장소 · 오믹스 종류 · 조직 · '
                '사례와 대조 인원 · 조달 경로를 담습니다. 조달 경로 두 열이 이 프로젝트의 '
                '핵심입니다 — 사례는 생검인데 대조가 신절제나 공여자이면 그 차이가 '
                '질병 신호처럼 보입니다.')
        o = st.multiselect('오믹스', sorted(ds['omics'].dropna().unique()),
                           default=sorted(ds['omics'].dropna().unique()))
        v = ds[ds['omics'].isin(o)]
        st.dataframe(v, width='stretch', hide_index=True)
        st.dataframe(load('select * from v_dataset_by_omics'),
                     width='stretch', hide_index=True)

    with tab4:
        st.info('9,900개 유전자 공간 안의 어떤 유전자든 코호트별 효과크기(Hedges g)를 '
                '볼 수 있습니다. 값은 각 코호트의 DKD 대 그 코호트 대조군 비교입니다. '
                '코호트마다 부호나 크기가 크게 다르면, 그 유전자는 한 코호트의 특성일 '
                '가능성이 있습니다.')
        st.subheader('유전자 조회')
        g = st.text_input('유전자 기호', 'MOXD1').strip().upper()
        if g:
            e = load("select cohort, hedges_g from gene_effect "
                     "where upper(symbol)='%s' order by cohort" % g.replace("'", ""))
            if e.empty:
                st.info('%s 는 9,900 유전자 공간에 없습니다.' % g)
            else:
                st.bar_chart(e.set_index('cohort'))
                st.dataframe(e, width='stretch', hide_index=True)
                st.caption('코호트별 Hedges g (DKD 대 그 코호트의 대조군).')


if __name__ == '__main__':
    main()
