# -*- coding: utf-8 -*-
"""알고리즘 설명과 전체 실험 결과 장표.

과제 발표에서 "무엇이 나왔나" 만큼 "무엇으로 분석했나" 를 묻습니다. 이 모듈은 그 질문에
답하는 장표를 담습니다 — 선택 알고리즘 8종과 교란 보정 6종을 각각 설명한 뒤, 같은 조건에서
돌린 비교표를 나란히 놓습니다.

    알고리즘별 비교표 (교란 보정 전)  ->  교란 보정 후 비교표

두 표를 이 순서로 두면 "알고리즘을 바꿔도 별 차이가 없고, 교란을 다루는 방식이 결과를
좌우한다" 는 이 과제의 결론이 표만 봐도 드러남.
"""
import os

from deck_style import blank, bullets, grid, head, note, panel, picture

KICK = '[다중 오믹스 DB 구축 및 DKD 바이오마커 후보 탐색]'
FIG = 'results/figures'


def s_algorithms(prs, f):
    s = blank(prs)
    head(s, KICK, '후보 선택 알고리즘 8종',
         '계열이 다른 방법을 모아 같은 조건에서 돌렸음. 마지막 하나가 제안 방법임.')
    grid(s, 0.77, 1.75, 11.80, [
        ['계열', '알고리즘', '작동 원리', '이 자료에서의 특징'],
        ['필터', 'Univariate (F-test)', '유전자 하나씩 군간 분산비로 순위',
         '가장 단순 · 부트스트랩 안정성은 최고'],
        ['임베디드', 'LASSO (L1)', '선형 모형에 L1 벌점, 계수가 0이 되며 선택',
         '희소하지만 상관된 유전자 중 하나만 남김'],
        ['임베디드', 'Elastic Net (L1+L2)', 'L1 과 L2 를 섞어 상관 그룹을 함께 남김',
         'LASSO 의 그룹 탈락 문제를 완화'],
        ['정보이론', 'mRMR', '관련성은 높이고 이미 뽑힌 것과의 중복은 낮춤',
         '중복은 줄지만 교차 코호트 재현은 낮음'],
        ['거리기반', 'ReliefF', '이웃 표본과의 차이로 특징 가중치를 갱신',
         '기존 방법 중 교차 폴드 재현성 1위'],
        ['앙상블', 'Random Forest 중요도', '트리 분할 기여도로 순위',
         '비선형을 잡지만 재현성이 가장 낮음'],
        ['앙상블', 'Boruta', '그림자(값을 섞은) 특징보다 유의하게 나은지 검정',
         '통계적 판정을 주지만 계산량이 큼'],
        ['제안', 'RBS (Robust Biomarker Score)',
         '코호트를 합치지 않고 각각 부트스트랩 선택 후 최악값으로 합산',
         '한 코호트에서만 강한 유전자를 배제'],
    ], widths=[1.3, 2.6, 4.6, 3.8], size=10.5, hl=(7,), row_h=0.40)
    note(s, 0.77, 6.15, 11.80,
         'RBS 의 핵심은 합치지 않는 것임. 코호트를 모아 한 번에 고르면 가장 큰 코호트의 '
         '특성이 목적함수를 지배함. 각 코호트에서 따로 고른 뒤 **최악의 코호트 기준**으로 '
         '합산하면 그런 유전자가 걸러지고, 효과 방향이 코호트 간에 뒤집히면 점수를 0으로 둠.')


def s_algo_compare(prs, f):
    s = blank(prs)
    head(s, KICK, '알고리즘별 비교표 — 교란 보정 전',
         '4코호트 LODO, 상위 50개 기준. 성능이 아니라 재현성이 방법을 가릅니다.')
    NAME = {'rbs': 'RBS (제안)', 'relieff': 'ReliefF', 'elastic_net': 'Elastic Net',
            'lasso': 'LASSO', 'mrmr': 'mRMR', 'boruta': 'Boruta',
            'univariate': 'Univariate', 'rf': 'Random Forest'}
    rows = [['알고리즘', '외부 AUROC', '표준편차', '내부 CV', '부트스트랩 Jaccard',
             '교차폴드 Jaccard', '무작위 대비 백분위']]
    for _, r in f['sel_table'].iterrows():
        rows.append([NAME.get(r['method'], r['method']),
                     '%.3f' % r['external'], '%.3f' % r['external_sd'],
                     '%.3f' % r['internal'], '%.3f' % r['boot_jaccard'],
                     '%.3f' % r['cross_fold_jaccard'], '%.0f%%' % (100 * r['null_pct'])])
    grid(s, 0.77, 1.80, 11.80, rows, widths=[2.2, 1.6, 1.3, 1.3, 1.9, 1.8, 2.0],
         size=11, hl=(0,), row_h=0.36)
    bullets(s, 0.77, 5.50, 11.80, 1.42, [
        (0, '외부 AUROC 로는 방법을 가릴 수 없음'),
        (1, '여덟 방법이 0.866~0.952 안에 모여 있고, 무작위 유전자도 이 범위에 들어옴'),
        (0, '교차 폴드 Jaccard 가 유일하게 벌어지는 축임'),
        (1, 'RBS 0.303 대 차순위 ReliefF 0.196. 나머지 여섯은 0.136~0.154 로 사실상 동률'),
        (1, '부트스트랩 Jaccard(같은 자료 안에서의 안정성)와는 순위가 다릅니다 — '
            'Univariate 이 0.261 로 1위지만 코호트가 바뀌면 0.140 으로 무너짐'),
    ])


def s_generative(prs, f):
    s = blank(prs)
    head(s, KICK, '생성모델·딥러닝은 이 규모에서 쓸 수 없었다',
         '증강으로 표본 부족을 메울 수 있는지 직접 시험했음.')
    g = f['gen']
    rows = [['생성 모델', '후보 서명 AUROC', '무작위 서명 AUROC', '차이', '판정']]
    for k, lab in (('vae', 'VAE (변분 오토인코더)'), ('diffusion', 'Diffusion')):
        if k in g.index:
            r = g.loc[k]
            rows.append([lab, '%.3f' % r['sig'], '%.3f' % r['rnd'],
                         '%+.3f' % r['margin'],
                         '이득 없음' if r['margin'] < 0.05 else '이득'])
    grid(s, 1.40, 1.85, 10.50, rows, widths=[3.0, 2.2, 2.2, 1.6, 1.8], size=11.5,
         hl=(0, 1), row_h=0.42)
    vae = g.loc['vae', 'margin'] if 'vae' in g.index else 0.0
    dif = g.loc['diffusion', 'sig'] if 'diffusion' in g.index else 0.0
    bullets(s, 1.40, 3.80, 10.50, 3.0, [
        (0, 'VAE 로 합성 표본을 만들어도 후보 서명이 무작위 서명을 이기지 못함'),
        (1, '차이 %+.3f — 생성기가 원자료의 상관 구조를 그대로 복제하므로 무작위 '
            '유전자도 같은 이득을 봄' % vae),
        (0, 'Diffusion 은 이 표본 수에서 학습 자체가 되지 않음'),
        (1, '후보 서명 AUROC %.3f 로 사실상 동전 던지기 — 코호트당 수십 명 규모에서 '
            '분포를 배울 수 없음' % dif),
        (0, '왜 이 과제가 깊은 모형을 주력으로 쓰지 않았는가'),
        (1, '전체 %d명, 유전자 %s개. 표본보다 변수가 38배 많은 상황에서 깊은 모형은 '
            '검증 자체가 성립하지 않음'
            % (f['n_case'] + f['n_ctrl'], format(f['genes'], ','))),
        (1, '표본 수에 맞는 방법(정규화 선형모형·앙상블·부트스트랩 안정성 선택)을 쓰고, '
            '남는 힘을 **검증 설계**에 씀'),
    ])


def s_corr_algorithms(prs, f):
    s = blank(prs)
    head(s, KICK, '교란 보정 알고리즘 6종',
         '배치 보정의 표준 도구부터 이 과제에서 만든 방법까지 같은 축으로 비교했음.')
    grid(s, 0.77, 1.75, 11.80, [
        ['알고리즘', '무엇을 가정하는가', '어떻게 제거하는가', '이 자료에서의 결과'],
        ['보정 없음', '—', '—', '기준선'],
        ['ComBat', '배치가 알려져 있고 코호트와 같다',
         '배치별 평균·분산을 경험적 베이즈로 맞춤', '교란이 코호트 안에 있어 **무효**'],
        ['RUVg (경험적 대조)', '질병과 무관한 유전자를 자료에서 찾을 수 있다',
         '가장 덜 연관된 유전자로 원치 않는 변동을 추정', '오히려 악화 — 채취 피처 27→29%'],
        ['RUVg (IEG 를 대조로)', 'IEG 모듈이 원치 않는 변동을 대표한다',
         'IEG 를 음성 대조 집합으로 지정해 인자 제거', 'IEG 는 줄지만 재현성 손실'],
        ['대리변수 보정 (SVA)', '설계를 보호한 채 숨은 인자를 찾을 수 있다',
         '잔차에서 대리변수를 추정해 함께 회귀', '채취 피처 27→6%, 특이성 최고'],
        ['핸들링 점수 잔차화', 'IEG 평균이 조직이 다루어진 정도의 대리값이다',
         '군 내부 기울기로 그 축을 각 유전자에서 제거', '채취 피처 27→5%, AUROC 손실 없음'],
    ], widths=[2.3, 3.2, 3.6, 3.2], size=10.5, hl=(4, 5), row_h=0.48)
    note(s, 0.77, 6.10, 11.80,
         '핵심은 **군 내부 기울기**임. 전체 기울기로 회귀하면 사례가 곧 생검이므로 '
         '질병 차이 자체를 지워 버림. 각 군 안에서 중심화한 뒤 기울기를 추정하면 '
         '교란된 군간 차이에 오염되지 않은 부분만 제거됨.')


def s_corr_compare(prs, f):
    s = blank(prs)
    head(s, KICK, '교란 보정 후 비교표',
         '같은 selector · 같은 fold 에서 보정만 바꿈. 4 cohorts LODO 평균, top 50.')
    t = f['corr_table']
    NAME = {'none': '보정 없음', 'ComBat': 'ComBat', 'RUVg_empirical': 'RUVg (경험적 대조)',
            'RUVg_ieg': 'RUVg (IEG 대조)', 'SVA_like': '대리변수 보정 (SVA)',
            'IEG_resid': '핸들링 점수 잔차화'}
    rows = [['보정 알고리즘', '외부 AUROC', 'IEG 개수', '채취 관련 피처',
             '특이성 통과', '배경 대비', '교차폴드']]
    hl = []
    for i, a in enumerate(('none', 'ComBat', 'RUVg_empirical', 'RUVg_ieg',
                           'SVA_like', 'IEG_resid')):
        if a not in t.index:
            continue
        r = t.loc[a]
        rows.append([NAME[a], '%.3f' % r['auroc'], '%.1f' % r['ieg'],
                     '%.0f%%' % (100 * r['proc']), '%.0f%%' % (100 * r['spec']),
                     '%.1f배' % r['enrich'], '%.3f' % r['xfold']])
        if a in ('SVA_like', 'IEG_resid'):
            hl.append(i)
    grid(s, 0.77, 1.80, 11.80, rows, widths=[2.6, 1.6, 1.3, 1.8, 1.5, 1.4, 1.4],
         size=11, hl=tuple(hl), row_h=0.38)
    bullets(s, 0.77, 4.70, 6.30, 2.2, [
        (0, 'ComBat 은 보정 없음과 숫자가 완전히 같음'),
        (1, '배치를 코호트로 잡으면 코호트 **안에** 있는 교란은 건드리지 못함'),
        (0, 'RUVg 는 대조 유전자를 무엇으로 잡느냐가 전부임'),
        (1, '경험적으로 고르면 오히려 악화, IEG 를 지정하면 개선'),
        (0, '두 방법이 들음'),
        (1, '채취 관련 피처를 27%에서 5~6%로 줄이면서 AUROC 손실이 유의하지 않음'),
    ])
    picture(s, os.path.join(FIG, 'P6_corrections.png'), 7.45, 4.80, w=5.00)


def s_conventional(prs, f):
    s = blank(prs)
    head(s, KICK, '기존 관행과의 비교',
         'DKD 문헌이 주로 쓰는 두 방식을 같은 자료·같은 관문으로 평가했음.')
    c = f['conv']
    NAME = {'DEG_meta': 'DEG 메타분석 (기존 관행)', 'WGCNA_hub': 'WGCNA 허브 (기존 관행)',
            'RBS_uncorrected': 'RBS, 보정 없음', 'RBS_corrected': 'RBS, 교란 보정'}
    rows = [['파이프라인', 'IEG 개수', '채취 관련 피처', '특이성 통과', '배경 대비']]
    hl = []
    for i, k in enumerate(('DEG_meta', 'WGCNA_hub', 'RBS_uncorrected', 'RBS_corrected')):
        if k not in c.index:
            continue
        r = c.loc[k]
        rows.append([NAME[k], '%d개' % int(r['n_ieg']),
                     '%.0f%%' % (100 * r['procurement_frac']),
                     '%.0f%%' % (100 * r['spec_pass_frac']),
                     '%.1f배' % r['spec_enrichment']])
        if k == 'RBS_corrected':
            hl.append(i)
    grid(s, 1.30, 1.85, 10.70, rows, widths=[3.4, 1.6, 2.0, 1.8, 1.6], size=11.5,
         hl=tuple(hl), row_h=0.42)
    bullets(s, 1.30, 4.25, 10.70, 2.5, [
        (0, '기존 관행도 이 자료에서는 같은 함정에 빠짐'),
        (1, 'DEG 메타분석과 WGCNA 허브 모두 상위 목록에 IEG 와 채취 관련 피처를 올림'),
        (0, '차이를 만드는 것은 알고리즘이 아니라 보정임'),
        (1, '같은 RBS 라도 보정 전후로 채취 관련 피처와 특이성 통과율이 달라짐'),
        (1, '즉 어떤 선택 알고리즘을 쓰는가보다 **교란을 어떻게 다루는가**가 결과를 '
            '더 크게 좌우함'),
    ])


def s_other_results(prs, f):
    s = blank(prs)
    head(s, KICK, '그 밖의 실험 결과 — 강건성과 다른 층에서의 확인',
         '후보를 흔들어 보고, 전사체 밖에서도 보이는지 확인했음.')
    st = f['stress']
    LAB = {'subsample': '표본 70% 추출', 'noise': '잡음 주입',
           'mask': '유전자 20% 결측', 'batch': '인공 배치 이동'}
    rows = [['섭동 방식', '후보 서명', '무작위 서명', '차이']]
    for k in ('subsample', 'noise', 'mask', 'batch'):
        if k in st.index:
            r = st.loc[k]
            rows.append([LAB[k], '%.3f' % r['sig'], '%.3f' % r['rnd'],
                         '%+.3f' % r['margin']])
    grid(s, 0.80, 1.85, 5.90, rows, widths=[2.2, 1.5, 1.5, 1.2], size=11, row_h=0.34)
    note(s, 0.80, 3.95, 5.90,
         '네 가지 섭동 모두에서 후보 서명이 무작위보다 0.09~0.11 앞섬. 섭동은 발굴이 '
         '끝난 뒤의 stress test 이고, 후보를 찾는 데는 쓰지 않았음.')
    panel(s, 7.00, 1.85, 5.70, '다른 층에서의 확인', [
        (True, '단백체 · 두 독립 코호트'),
        '  MMP7 이 두 측정 원리에서 같은 방향, CMH p = 0.030',
        (True, '단일세포 · KPMP 교란 대응 대비'),
        '  DKD 대 비당뇨 CKD  p = 0.21 (음성)',
        (True, 'GWAS · 한국인 DKD 2,532명'),
        '  후보 %d개 중 유의 %d개 (기대 1.5) — 음성' % (f['gwas_n'], f['gwas_sig']),
        (True, '대사체 · 혈장 LC-MS'),
        '  매칭 대비에서 29개 통과, 4-hydroxyproline 포함',
    ], size=10.5, row_h=0.34)
    bullets(s, 0.80, 5.45, 11.90, 1.48, [
        (0, '음성 결과를 숨기지 않았음'),
        (1, '단일세포와 GWAS 는 후보를 지지하지 않음. 그 두 검증이 신규 후보 '
            'OLFML3·PRSS23 을 걸러낸 근거이기도 함'),
        (1, 'WGCNA 허브 기준 %d가지 변형 중 IEG 가 상위 50에 든 경우는 %d건 — 교란이 '
            '허브 선택 방식 때문에 생긴 것이 아님을 확인'
            % (f['wgcna_n'], f['wgcna_ieg_top50'])),
    ])
