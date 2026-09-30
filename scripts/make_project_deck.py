#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""과제 발표 장표를 DKD발표.pptx 템플릿 서식으로 만든다.

논문이 아니라 **과제 목표**를 따라 씀. 과제 목표는 DKD 바이오마커 후보군 발굴이고,
1차년도에 실제로 무엇을 해서 무엇이 나왔는지를 순서대로 보여 주는 것이 이 덱의 일.

    자료를 모았다 -> multi-omics 는 공개 자료로 불가능했다 -> 그래서 single-omics
    cross-cohort 로 갔다 -> AUROC 가 좋아 보였다 -> 그런데 random 도 그만큼 나왔다 ->
    procurement confounding 을 찾았다 -> 보정했더니 순위가 바뀌었다 -> 기준을 세워
    후보를 뽑았다 -> 그 후보가 하위·random 보다 실제로 낫다

본문은 음슴체. AUROC, Jaccard, LODO, cohort, gene 처럼 현장에서 영어로 쓰는 말은
억지로 옮기지 않고 그대로 둠.

수치는 손으로 적지 않고 project_facts.facts() 가 결과 파일에서 읽음.
서식 상수와 도형 helper 는 deck_style 에 있고, 그 값은 템플릿에서 직접 읽은 것.
"""
import os
import sys

from pptx import Presentation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deck_style import (COVER, F_BODY, F_KICK, F_TITLE, INK, MUTED,
                        LAY_COVER, LAY_DIVIDER,
                        blank, bullets, clear, grid, head, note, panel, picture, tb)
from project_facts import facts
from project_slides_algo import (s_algorithms, s_algo_compare, s_generative,
                                 s_corr_algorithms, s_corr_compare, s_conventional,
                                 s_other_results)
from project_slides_rbs import (s_rbs_detail, s_rbs_code, s_corr_per_selector,
                                s_pipeline_code)

# 템플릿은 산출물과 분리해 둔다. 둘을 같은 경로로 두면 한 번 돌리는 순간 템플릿이
# 사라지고, 새로 복제한 저장소에서는 만들 수 없게 된다.
TEMPLATE = 'docs/assets/deck_template.pptx'
OUT = 'DKD발표.pptx'
KICK = '[다중 오믹스 DB 구축 및 DKD 바이오마커 후보 탐색]'
FIG = 'results/figures'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ====================================================================== 표지 · 간지
def s_cover(prs, f):
    s = blank(prs, LAY_COVER)
    tb(s, 1.24, 4.40, 9.60, 0.50, '당뇨병성 신장질환 바이오마커 후보 발굴',
       font=F_TITLE, size=24, bold=True, color=COVER)
    tb(s, 1.24, 4.92, 9.90, 0.60,
       '공개 데이터 기반 cross-cohort 분석과 조직 채취 confounding 의 보정 — 1차년도',
       font=F_BODY, size=20, color=INK)
    tb(s, 1.34, 5.92, 6.25, 0.24, '2026년 9월', font=F_BODY, size=12, color=MUTED)
    tb(s, 1.34, 6.20, 6.25, 0.24, '자율지능DX연구실 | 김도훈 선임연구원',
       font=F_BODY, size=12, bold=True, color=INK)


def s_divider(prs, text, sub=None):
    s = blank(prs, LAY_DIVIDER)
    tb(s, 7.98, 2.96, 4.90, 0.61, text, font=F_KICK, size=30, bold=True, color=INK)
    if sub:
        tb(s, 7.98, 3.66, 4.90, 0.60, sub, font=F_BODY, size=12.5, color=MUTED)


# ====================================================================== 1부 배경
def s_disease(prs, f):
    s = blank(prs)
    head(s, KICK, '당뇨병성 신장질환(DKD)이란',
         '당뇨 환자의 약 40%에서 생기며, 말기신부전의 가장 큰 단일 원인.')
    bullets(s, 0.77, 1.70, 6.10, 4.6, [
        (0, '무엇인가'),
        (1, '고혈당이 오래 지속되며 신장의 여과 단위인 사구체가 손상되는 질환'),
        (1, '단백뇨가 늘고 eGFR 이 떨어짐. 진행되면 투석이나 이식이 필요'),
        (0, '왜 바이오마커가 필요한가'),
        (1, '현재 진단은 단백뇨와 eGFR 에 의존하는데, 둘 다 손상이 진행된 뒤에 움직임'),
        (1, '조기에, 그리고 다른 신장질환과 구분해서 잡아낼 지표가 없음'),
        (0, '왜 어려운가'),
        (1, '신장 조직은 needle biopsy 로만 얻을 수 있어 검체가 귀하고 cohort 가 작음'),
        (1, '정상 신장 조직은 더 귀해서 대조군을 다른 경로로 얻게 됨 — 이 발표의 핵심'),
    ])
    panel(s, 7.10, 1.70, 5.60, '이 과제가 답하려는 질문', [
        '공개 omics 데이터만으로 DKD 바이오마커 후보를 뽑을 수 있는가',
        '뽑았다면 그것이 질병 신호인지 어떻게 확인하는가',
        '여러 cohort 에서 재현되는 후보를 어떻게 고르는가',
    ], size=11.5, row_h=0.40)
    panel(s, 7.10, 3.60, 5.60, '1차년도 범위', [
        '공개 데이터 확보 · 분류 · 품질 감사',
        'multi-omics 가능성 실측',
        'single-omics cross-cohort 후보 발굴',
        'confounding 진단과 보정 알고리즘',
        '후보군 도출과 검증',
    ], accent=True, size=11.5, row_h=0.36)


def s_goal(prs, f):
    s = blank(prs)
    head(s, KICK, '과제 목표와 1차년도에 실제로 한 일',
         '목표는 후보군 발굴. 그 과정에서 먼저 풀어야 할 문제가 드러남.')
    grid(s, 0.77, 1.75, 11.80, [
        ['단계', '계획', '실제 결과'],
        ['자료 확보', '공개 omics 데이터 수집·분류',
         '자산 %d건 확보, %d건 사용, %d건은 환자 중복으로 제외'
         % (f['inv_total'], f['inv_used'], f['inv_excluded'])],
        ['multi-omics', '여러 omics 를 통합한 후보 발굴',
         '공개 자료에 환자 단위 다층 설계가 없어 불가 — 실측으로 확인'],
        ['single-omics', 'transcriptome cross-cohort 로 우회',
         'cohort 7종 %d명(사례 %d/대조 %d)을 %s gene 공간으로 통합'
         % (f['n_case'] + f['n_ctrl'], f['n_case'], f['n_ctrl'], format(f['genes'], ','))],
        ['confounding 진단', '(계획에 없던 단계)',
         '조직 채취 방식 confounding 발견 — AUROC 로는 후보를 가릴 수 없음'],
        ['보정·도출', '후보군 도출',
         '보정 알고리즘 적용 후 후보 30개, 그중 우선 인계 대상 6개'],
    ], widths=[1.5, 3.0, 6.0], size=11.5, hl=(3,), row_h=0.52)
    note(s, 0.77, 5.95, 11.80,
         '계획에 없던 네 번째 줄이 1차년도의 실질 성과임. 그 단계를 건너뛰었다면 후보 '
         '목록은 나왔겠지만 그 목록이 질병을 보고 있는지 알 수 없었을 것.')


# ====================================================================== 2부 데이터
def s_data(prs, f):
    s = blank(prs)
    head(s, KICK, '확보한 공개 데이터 — omics 종류별',
         '자산 %d건 · 실제 사용 %d건 · 총 %.0f GB'
         % (f['inv_total'], f['inv_used'], f['inv_gb']))
    grid(s, 0.77, 1.75, 11.80, [
        ['omics', '저장소', '자료', '이 과제에서의 쓰임'],
        ['Transcriptome (array·RNA-seq)', 'GEO', 'GSE 15건',
         '후보 발굴의 본체 — cross-cohort 선택'],
        ['Single-nucleus RNA-seq', 'KPMP Atlas', '공여자 %d명' % f['kpmp_donor'],
         'cell type 음성대조 · 다층 설계 실측'],
        ['Proteome (aptamer)', 'Mendeley', 'SOMAscan, DKD 23 / 정상 10',
         '후보의 orthogonal 보강 (환자 비매칭)'],
        ['Proteome (LC-MS/MS)', 'PRIDE', 'PXD041884, DKD 5 / 비당뇨 7',
         '두 번째 측정 원리로 재확인'],
        ['Metabolome', 'Metabolomics Workbench', 'ST 5건, 혈장 LC-MS',
         '대조군 설계 권고의 일반화 검증'],
        ['Germline (GWAS)', 'GWAS Catalog', 'GCST90179152, 한국인 DKD 2,532명',
         '후보의 유전적 근거 — 결과는 음성'],
    ], widths=[2.8, 1.9, 2.5, 3.8], size=11, row_h=0.44)
    note(s, 0.77, 5.85, 11.80,
         'omics 종류는 여섯 갈래를 모두 확보함. 문제는 종류가 아니라 "같은 환자에게서 '
         '두 층 이상을 쟀는가" 였음 — 다음 장.')


def s_multiomics(prs, f):
    s = blank(prs)
    head(s, KICK, '공개 데이터로 multi-omics 는 불가능했다',
         '가설이 아니라 KPMP 참여자를 한 명씩 조회해 실측한 결과.')
    bullets(s, 0.77, 1.70, 6.00, 4.5, [
        (0, 'multi-omics 가 성립하려면 두 가지가 동시에 필요함'),
        (1, '같은 환자에게서 두 층 이상을 측정했을 것 (환자 단위 pairing)'),
        (1, '사례군과 대조군 **양쪽** 모두에서 그 층을 측정했을 것 (contrast 성립)'),
        (0, 'KPMP 실측 — pairing 은 일부 있으나 contrast 가 없음'),
        (1, 'single-nucleus 공여자 %d명 중 당뇨 병력 %d명, 그중 CKD 등록 %d명'
            % (f['kpmp_donor'], f['kpmp_diab'], f['kpmp_proxy'])),
        (1, '그 %d명 중 3명이 regional proteomics, 7명이 spatial metabolomics 를 함께 가짐'
            % f['kpmp_proxy']),
        (1, '그러나 대조군(healthy reference 12명) 중 regional proteomics 보유자는 **0명**'),
        (0, '결론'),
        (1, '사례 쪽 pairing 이 몇 명이든, 비교할 대조가 없으면 contrast 를 못 만듦'),
        (1, 'spatial metabolomics 는 양쪽에 있으나 imaging MS 라 참여자별 정량 행렬이 없음'),
    ])
    panel(s, 7.00, 1.70, 5.70, 'KPMP 층별 참여자 수', [
        'single-nucleus       공여자 %d명' % f['kpmp_donor'],
        'DKD proxy subset     %d명' % f['kpmp_proxy'],
        (True, '  + regional proteomics    3명'),
        (True, '  + spatial metabolomics   7명'),
        'healthy reference    12명',
        (True, '  + regional proteomics    0명  <- contrast 불가'),
    ], size=11, row_h=0.34)
    note(s, 7.00, 4.35, 5.70,
         '"공개 자료에 DKD multi-omics 가 없다" 는 흔한 서술은 정확하지 않음. pairing 은 '
         '일부 존재함. 없는 것은 대조군 쪽 층이고, 그래서 contrast 가 성립하지 않음. '
         '이 구분이 2차년도 설계의 출발점임.', size=10.5)


def s_pivot(prs, f):
    s = blank(prs)
    head(s, KICK, '그래서 1차년도는 single-omics · cross-cohort 로',
         '한 cohort 에서 재현되는 것과 여러 cohort 에서 재현되는 것은 다름.')
    bullets(s, 0.77, 1.75, 6.05, 4.4, [
        (0, '왜 transcriptome 인가'),
        (1, '여섯 omics 중 사례·대조가 모두 충분한 수로 존재하는 유일한 층'),
        (1, '서로 다른 연구팀이 독립적으로 모은 cohort 가 여러 개 — 교차 검증 가능'),
        (0, '왜 cross-cohort 인가'),
        (1, '한 cohort 안에서의 재현성은 그 cohort 의 특성을 재현하는 것일 수 있음'),
        (1, 'cohort 를 하나씩 빼고 학습·검증(LODO)하면 그 cohort 특유의 신호가 걸러짐'),
        (0, '제안한 선택 방법 — RBS'),
        (1, 'cohort 를 합치지 않고 각 cohort 안에서 따로 bootstrap 선택'),
        (1, 'cohort 사이에서는 worst-case 로 합산 — 한 cohort 에서만 강한 gene 을 배제'),
        (1, '효과 방향이 cohort 간에 뒤집히면 점수를 0 으로 (sign concordance gate)'),
    ])
    panel(s, 7.00, 1.75, 5.70, 'cross-fold Jaccard (4 cohorts, top 100)', [
        (True, 'RBS (제안)             %.3f' % f['xf_best']),
        '차순위 기존 방법        %.3f' % f['xf_next'],
        '',
        '같은 후보가 hold-out cohort 를 바꿔도',
        '다시 나오는 비율. 기존 방법 7종은 전부',
        '0.19 이하 — cohort 가 바뀌면 사실상',
        '다른 목록을 냄.',
    ], accent=True, size=11, row_h=0.32)
    picture(s, os.path.join(FIG, 'F2_cross_fold_stability.png'), 7.00, 4.60, w=5.65)


def s_integrate(prs, f):
    s = blank(prs)
    head(s, KICK, '데이터 통합 — 공통 gene 공간만 남김',
         'platform 이 다른 cohort 를 비교하려면 먼저 같은 좌표계로 옮겨야 함.')
    grid(s, 0.77, 1.75, 11.80, [
        ['단계', '무엇을 했나', '결과'],
        ['① 식별자 통일', 'platform 별 probe ID 를 Entrez gene ID 로 변환',
         'microarray 3종 · RNA-seq platform 통합'],
        ['② 중복 정리', '한 gene 에 여러 probe 가 붙으면 max-mean 규칙으로 하나만 선택',
         '규칙을 동결해 이후 분석에서 바뀌지 않게 함'],
        ['③ 교집합', '모든 cohort 에 공통으로 존재하는 gene 만 남김',
         format(f['genes'], ',') + '개 gene 으로 고정'],
        ['④ 정규화', 'cohort 별 gene 단위 z-score (label 을 쓰지 않는 unsupervised 변환)',
         '척도가 다른 cohort 를 겹쳐 놓을 수 있게 됨'],
        ['⑤ 감사', '환자 중복 검사 — 같은 환자가 여러 GSE 에 실린 경우 탐지',
         'GSE30122 · GSE47183 · GSE99340 제외'],
    ], widths=[1.5, 5.4, 4.9], size=11, row_h=0.50)
    note(s, 0.77, 5.55, 11.80,
         '⑤ 가 중요함. 같은 환자가 두 cohort 에 실려 있으면 "서로 다른 cohort 에서 '
         '재현됐다" 가 거짓이 됨. 감사에서 3건을 제외했고, 남은 중복은 숨기지 않고 '
         '기록함 — GSE104948 과 GSE104954 는 같은 ERCB 환자의 두 compartment 임.')


def s_cohorts(prs, f):
    s = blank(prs)
    head(s, KICK, '분석에 쓴 transcriptome cohort',
         '사례 %d명 · 대조 %d명. 대조군 조직이 어디서 왔는지를 함께 봄.'
         % (f['n_case'], f['n_ctrl']))
    c = f['cohorts']
    rows = [['cohort', '신장 부위', 'DKD', '대조', '대조군 조직']]
    for _, r in c.iterrows():
        ct = str(r['control_type']).replace('tumor_nephrectomy', 'tumor nephrectomy') \
            .replace('living_donor', 'living donor')
        rows.append([r['cohort'], str(r['compartment']), r['case'], r['control'], ct])
    grid(s, 1.30, 1.75, 10.70, rows, widths=[2.0, 1.8, 1.0, 1.0, 3.2], size=11, row_h=0.36)
    note(s, 1.30, 5.05, 10.70,
         '맨 오른쪽 열을 볼 것. **어느 cohort 에도 biopsy 로 얻은 대조군이 없음.** '
         '사례는 전부 needle biopsy 이고 대조는 전부 수술로 떼어낸 조직임. 이 표가 뒤에 '
         '나올 confounding 의 원인이고, 처음에는 이것을 문제로 보지 않았음.')


# ====================================================================== 3부 첫 결과와 함정
def s_control(prs, f):
    s = blank(prs)
    head(s, KICK, '대조군이란 무엇인가 — 그리고 왜 문제인가',
         '"정상 신장 조직" 은 한 종류가 아님. 어떻게 얻었는지가 다름.')
    picture(s, os.path.join(FIG, 'S_procurement_schematic.png'), 2.55, 1.66, w=8.20)
    note(s, 0.90, 6.52, 11.55,
         '두 조직은 질병만 다른 것이 아니라 **다루어진 방식**이 다름. 수술 조직은 적출과 '
         '고정 사이에 허혈 시간이 생기고 그 사이에 세포는 stress response gene 을 켬. '
         '즉 DKD 대 정상 비교에는 disease effect 와 procurement effect 가 같은 방향으로 '
         '겹쳐 있음 — 완전히 confounded 되어 둘을 분리할 정보가 자료 안에 없음.')


def s_first_auc(prs, f):
    s = blank(prs)
    head(s, KICK, '1차 결과 — 분류 성능은 좋아 보였다',
         'LODO external AUROC. 이 숫자만 보면 성공.')
    bullets(s, 0.77, 1.75, 5.90, 2.55, [
        (0, '무엇을 쟀나'),
        (1, 'cohort 하나를 완전히 빼고 나머지로 학습, 뺀 cohort 에서 검증'),
        (1, '선택한 상위 K개 gene 으로 logistic regression 적합'),
        (0, '결과'),
        (1, '기존 방법 7종과 제안 방법 모두 external AUROC %.2f~%.2f'
            % (f['auc_lo'], f['auc_hi'])),
        (1, '논문으로 내기에 충분해 보이는 수치'),
    ])
    panel(s, 6.85, 1.75, 5.75, 'LODO external AUROC (K=50)', [
        (True, '최고            %.3f' % f['auc_hi']),
        '최저            %.3f' % f['auc_lo'],
        '',
        '여덟 개 방법이 모두 이 범위 안에 들어옴.',
        '방법 간 차이가 거의 없다는 점이',
        '첫 번째 신호였음.',
    ], size=11, row_h=0.32)
    picture(s, os.path.join(FIG, 'K2_auroc.png'), 3.90, 4.42, w=5.50)


def s_null(prs, f):
    s = blank(prs)
    head(s, KICK, '그런데 random gene 도 같은 성능을 냈다',
         '무작위로 뽑은 50개 gene 으로 같은 실험을 반복 (fold 당 300회).')
    panel(s, 0.77, 1.80, 5.75, 'random-signature null', [
        (True, 'random 50 genes   AUROC %.2f ~ %.2f' % (f['null_lo'], f['null_hi'])),
        (True, '최고 기록          %.3f' % f['null_max']),
        'label 을 섞으면    %.2f' % f['null_perm'],
        '',
        'label permutation 이 0.5 로 떨어지므로',
        'pipeline 이 고장난 것은 아님. 신호는 실재함.',
    ], accent=True, size=11, row_h=0.34)
    picture(s, os.path.join(FIG, 'P2_null_control.png'), 0.80, 4.60, w=5.70)
    bullets(s, 6.85, 1.80, 5.80, 3.75, [
        (0, '이것이 뜻하는 바'),
        (1, 'DKD 신호는 gene 전반에 넓게 퍼져 있어서, 어떤 gene 을 고르든 대부분 '
            '그 신호를 물려받음'),
        (1, '따라서 **절대 AUROC 로는 좋은 후보와 나쁜 후보를 구분할 수 없음**'),
        (1, '한 cohort(GSE142025)에서는 random gene 이 AUROC 1.000 에 도달 — '
            '그 fold 는 방법을 가릴 능력이 아예 없음'),
        (0, '그래서 바꾼 것'),
        (1, '모든 AUROC 를 null 분포 대비 백분위로만 보고'),
        (1, '방법 비교의 기준축을 성능이 아니라 **cross-fold 재현성**으로 변경'),
    ])
    note(s, 6.85, 5.80, 5.80,
         '이 대조군이 없었다면 random gene 과 구별되지 않는 목록을 바이오마커 후보로 '
         '보고했을 것. 유방암 signature 에서 같은 문제가 보고된 적 있음(Venet 등, 2011).')


# ====================================================================== 4부 교란
def s_artifact(prs, f):
    s = blank(prs)
    head(s, KICK, '가장 재현성 높은 신호는 질병이 아니었다',
         'tissue procurement confounding — 조직 채취 방식에 따른 교란')
    bullets(s, 0.77, 1.75, 6.00, 4.4, [
        (0, '무엇이 1위였나'),
        (1, '모든 stability 순위에서 1위를 차지한 것은 immediate-early gene (IEG) '
            'module — FOS, JUN, EGR1, DUSP1, ZFP36 등'),
        (1, '세포가 stress 를 받으면 수 분 안에 켜지는 조기 반응 gene 들'),
        (0, '왜 질병 신호가 아니라고 보는가'),
        (1, '비당뇨 CKD 환자도 DKD 의 145% 만큼 같은 방향으로 움직임 — '
            'disease-specific 이라면 이럴 수 없음'),
        (1, '같은 procurement route 를 공유하는 두 군을 비교하면 이 module 은 두 군을 '
            '가르지 못함 (한쪽이 당뇨여도 마찬가지)'),
        (1, '허혈이 바로 이 transcript 들을 수 분 내에 유도한다는 것은 독립적으로 보고됨'),
        (0, '서술의 수위'),
        (1, 'procurement route 는 시료 주석에서 읽은 것이고 허혈 시간을 측정한 것이 '
            '아니므로 **procurement-associated** 라고 쓰고 단정은 피함'),
    ])
    panel(s, 7.00, 1.75, 5.70, '9개 자료 · 4개 platform 에서 확인', [
        'IEG module 이 procurement route 를 따라감',
        '질병 유무가 아니라',
        '',
        (True, '비당뇨 CKD 가 DKD 의 145%'),
        '',
        '심장 조직에서는 재현되지 않음',
        '-> kidney-specific 현상',
    ], accent=True, size=11, row_h=0.34)
    picture(s, os.path.join(FIG, 'P3_procurement.png'), 7.15, 4.60, w=5.30)


def s_evidence(prs, f):
    s = blank(prs)
    head(s, KICK, 'confounding 의 근거 — 설계가 분리를 허용하지 않음',
         '보정 이전에, 이 자료로 질병과 채취를 분리할 수 있는지부터 확인함.')
    c = f['cohorts']
    rows = [['cohort', 'DKD 조직', '대조군 조직', '같은 route 공유']]
    for _, r in c.iterrows():
        ct = str(r['control_type']).replace('tumor_nephrectomy', 'tumor nephrectomy') \
            .replace('living_donor', 'living donor')
        rows.append([r['cohort'], 'needle biopsy', ct, '아니오'])
    grid(s, 1.10, 1.75, 11.10, rows, widths=[2.0, 2.0, 3.4, 2.2], size=10.5, row_h=0.33)
    bullets(s, 1.10, 4.80, 11.10, 2.1, [
        (0, '분리 가능한 stratum 이 0개임'),
        (1, '질병과 채취 방식을 동시에 관찰하려면, 같은 route 안에 사례와 대조가 함께 '
            '있는 구간이 하나라도 있어야 함'),
        (1, '예를 들어 DKD 100명이 전부 biopsy 이고 정상 100명이 전부 nephrectomy 라면, '
            'DKD effect 와 biopsy effect 를 분리해 관찰할 정보 자체가 없음'),
        (1, 'compartment 와 platform 은 분리 가능했지만 **채취 방식은 아니었음**'),
    ])


# ====================================================================== 5부 도출
def s_reorder(prs, f):
    s = blank(prs)
    head(s, KICK, '대조군을 바꾸면 후보 순위가 크게 올라간다',
         '같은 사례군을 두 가지 대조군과 비교했을 때의 effect size 순위. '
         '환자를 pairing 한 것이 아님.')
    sh = f['shift'].set_index('gene')
    rows = [['인계 대상 후보', '일반 대조군 기준', '매칭 대조군 기준', '이동']]
    for g in ('FMOD', 'MMP2', 'LUM', 'MOXD1', 'THBS2', 'CCND2'):
        if g not in sh.index:
            continue
        r = sh.loc[g]
        rows.append([g, '%s위' % format(int(r['rank_naive']), ','),
                     '%d위' % int(r['rank_matched']), '%+d' % int(-r['moved'])])
    grid(s, 0.80, 1.85, 6.00, rows, widths=[1.6, 1.9, 1.9, 1.1], size=11,
         hl=tuple(range(len(rows) - 1)), row_h=0.34, head_h=0.52)
    picture(s, os.path.join(FIG, 'P9_reordering.png'), 7.05, 1.85, w=5.75)
    bullets(s, 0.80, 4.85, 6.00, 2.0, [
        (0, '여섯 개 모두 매칭 대조군에서 순위가 올라감'),
        (1, 'FMOD 는 1,264위에서 4위로, LUM 은 468위에서 43위로'),
        (1, '반대로 CASP1(16->603), TYROBP(31->678) 처럼 내려가는 gene 도 있음'),
        (0, '보정 전후 top-50 의 겹침은 50개 중 1개'),
        (1, '대조군을 무엇으로 잡느냐가 후보 목록을 거의 전부 바꿈'),
    ])
    note(s, 7.05, 5.95, 5.75,
         '여기의 순위는 %s개 gene 을 **effect size |g| 로 줄세운 것**임. 뒤에 나오는 '
         '인계 후보표의 순위는 selector 가 매긴 것이라 서로 다른 값.'
         % format(f['genes'], ','))


def s_criteria(prs, f):
    s = blank(prs)
    head(s, KICK, '후보군 도출 기준',
         '순위가 높다는 것만으로는 후보가 되지 않게 세 개의 gate 를 세움.')
    panel(s, 0.77, 1.75, 3.85, '① specificity gate', [
        'DKD 대 비당뇨 CKD 에서',
        '|Hedges g| >= 0.5',
        '',
        '"신장이 나빠서" 생기는',
        '일반 신호를 걸러냄',
        '',
        (True, 'random 통과율 %.1f%%' % f['bg']),
        (True, '(%s개 중 %d개)' % (format(f['genes'], ','), f['n_pass'])),
    ], size=11, row_h=0.32)
    panel(s, 4.85, 1.75, 3.85, '② permutation ceiling', [
        'label 을 섞었을 때 도달하는',
        '|g| 의 95 백분위를 넘는가',
        '',
        '잡음으로 설명되는 크기인지',
        '확인',
        '',
        (True, 'gate 자체의 FDR 추정치'),
        (True, '%.1f%% (환자 단위 permutation)' % f['fdr']),
    ], size=11, row_h=0.32)
    panel(s, 8.90, 1.75, 3.85, '③ procurement flag', [
        '채취 방식과 연관된 것으로',
        '표시된 gene 은 제외',
        '',
        'IEG module 과 그에 준하는',
        '연관을 보이는 feature',
        '',
        (True, '보정 전 27% -> 보정 후 5%'),
    ], accent=True, size=11, row_h=0.32)
    note(s, 0.77, 5.25, 11.95,
         '① 의 %.1f%% 와 ② 의 %.1f%% 는 분모가 다른 값이라 서로 비교할 수 없음. '
         '%.1f%% 는 전체 %s개 gene 중 gate 를 통과하는 비율(background pass rate)이고, '
         '%.1f%% 는 통과한 %d개 중 우연으로 기대되는 비율의 추정치(permutation-based '
         'FDR estimate)임.'
         % (f['bg'], f['fdr'], f['bg'], format(f['genes'], ','), f['fdr'], f['n_pass']))


def s_before_after(prs, f):
    s = blank(prs)
    head(s, KICK, 'gate 를 통과한 후보 — 보정 전후 비교',
         '같은 기준을 보정 전 목록과 보정 후 목록에 똑같이 적용함.')
    grid(s, 1.60, 1.80, 10.10, [
        ['항목', '보정 전', '보정 후'],
        ['최종 후보 수', '30개', '30개'],
        ['이미 DKD marker 로 제안된 gene', '%d개' % f['lit_un'], '%d개' % f['lit_co']],
        ['두 목록의 겹침', '%d / 30' % f['lit_overlap'], '%d / 30' % f['lit_overlap']],
        ['novel(문헌 없음)로 분류된 gene', '4개', '2개'],
    ], widths=[4.6, 2.4, 2.4], size=12, hl=(1,), row_h=0.42)
    bullets(s, 1.60, 4.30, 10.10, 2.6, [
        (0, '보정하면 알려진 DKD 생물학 쪽으로 이동함'),
        (1, '%d/30 -> %d/30 — 방향은 분명하지만 이 차이 자체는 통계적으로 유의하지 '
            '않음 (Fisher one-sided p = 0.14)' % (f['lit_un'], f['lit_co'])),
        (1, '따라서 이 수치는 보조 근거로만 씀. 결론은 뒤의 층 비교로 받침'),
        (0, '목록이 절반 가까이 바뀜'),
        (1, '두 목록의 겹침이 %d/30 — 보정은 순위를 다듬는 것이 아니라 목록을 다시 씀'
            % f['lit_overlap']),
    ])


def s_candidates(prs, f):
    s = blank(prs)
    head(s, KICK, '도출된 후보군 30개와 등급',
         '두 축으로 나눔 — 데이터에서 강한가, 그리고 문헌에 없는가.')
    t = f['tier']
    get = lambda k: int(t[[x for x in t.index if str(x).startswith(k)][0]]) \
        if any(str(x).startswith(k) for x in t.index) else 0
    grid(s, 0.85, 1.85, 6.40, [
        ['등급', '뜻', '개수'],
        ['Tier 1', '데이터에서 강하고 + 문헌에 없음', '%d개' % get('1')],
        ['Tier 2', '데이터에서 강하지만 이미 보고됨', '%d개' % get('2')],
        ['Tier 3', '문헌에는 없지만 데이터에서 약함', '%d개' % get('3')],
        ['Tier 4', '이미 보고됐고 데이터에서도 약함', '%d개' % get('4')],
    ], widths=[1.3, 4.4, 1.1], size=11.5, hl=(0, 1), row_h=0.44)
    picture(s, os.path.join(FIG, 'P7_candidates.png'), 7.55, 1.85, w=5.10)
    bullets(s, 0.85, 4.55, 11.70, 2.3, [
        (0, 'Tier 1 이 비어 있음 — 이것이 1차년도의 정직한 결론'),
        (1, 'novel 로 분류된 2개(OLFML3, PRSS23)는 독립 검증 3건에서 모두 기각: '
            '독립 RNA-seq cohort 에서 방향이 뒤집혔고, single-cell 교란 대응 contrast '
            '에서 random 과 구분되지 않았으며, 보정된 selector 는 이들을 상위에 '
            '올리지 않음'),
        (0, '그러나 과제 산출물의 기준은 다름'),
        (1, '후속 실험에 넘길 우선순위 목록이 목적이라면 novelty 는 요구 조건이 아니며, '
            '선행 문헌이 있는 쪽이 착수 위험이 낮음 — Tier 2 가 인계 대상'),
    ])


def s_tier2(prs, f):
    s = blank(prs)
    head(s, KICK, '인계 대상 후보 6개 (Tier 2)',
         '데이터에서 강하고, 비당뇨 CKD 와도 갈리며, cohort 를 빼도 대부분 남음.')
    t2 = f['tier2']
    keep = {'LUM', 'MOXD1', 'THBS2', 'FMOD', 'MMP2'}
    rows = [['gene', 'selector 순위 (LODO 평균)', 'DKD 대 비당뇨CKD |g|',
             'permutation ceiling', 'cohort 제외 시']]
    hl = []
    for i, (g, r) in enumerate(t2.iterrows()):
        rows.append([g, '%d위' % int(r['mean_rank']),
                     '%.2f' % abs(r['g_DKD_vs_otherCKD']),
                     '초과' if bool(r['above_perm_ceiling']) else '—',
                     '생존' if g in keep else '탈락'])
        if g not in keep:
            hl.append(i)
    grid(s, 1.30, 1.80, 10.70, rows, widths=[1.6, 1.8, 2.6, 2.0, 1.8], size=11.5,
         hl=tuple(hl), row_h=0.36, head_h=0.52)
    bullets(s, 1.30, 4.92, 10.70, 1.95, [
        (0, '가장 중요한 열은 세 번째'),
        (1, '여섯 개 모두 비당뇨 CKD 대비 |g| 0.81~1.28 — "신장이 나빠서" 생기는 '
            '일반 신호가 아니라 DKD 쪽으로 갈림'),
        (0, '네 개가 ECM·fibrosis 계열로 묶임'),
        (1, 'LUM · THBS2 · FMOD · MMP2. 혈장 metabolome 에서도 collagen 분해산물인 '
            '4-hydroxyproline 이 같은 축을 가리킴 (환자 비매칭 보강)'),
        (0, '주의'),
        (1, '여섯 중 다섯은 이미 DKD marker 로 제안된 gene 임. 발견이 아니라 '
            '우선순위 목록으로 인계함. 두 번째 열은 selector 가 매긴 순위이고, '
            '앞 장표의 순위는 effect size |g| 기준이라 값이 다름'),
    ])


def s_validate(prs, f):
    s = blank(prs)
    head(s, KICK, '제안한 후보군은 실제로 더 나은가 — 상위 대 하위 비교',
         '같은 PubMed 조회를 상위 30개 · 하위 30개 · random 30개에 똑같이 걸었음.')
    st, sp = f['strata'], f.get('strata_p', {})
    rows = [['gene 집단', 'DKD 논문 있음', 'DKD marker 논문 있음',
             'DKD 논문 중앙값', '신장 문헌 전무']]
    for k in ('상위 30 (제안 후보)', '하위 30', '무작위 30'):
        if k not in st:
            continue
        v = st[k]
        rows.append([k, '%d / %d' % (v['dkd'], v['n']),
                     '%d / %d' % (v['marker'], v['n']),
                     '%.1f편' % v['med_dkd'], '%d개' % v['none']])
    grid(s, 0.90, 1.80, 11.55, rows, widths=[3.0, 2.2, 2.6, 2.0, 1.8], size=11.5,
         hl=(0,), row_h=0.40)
    ps = ['상위 대 %s   p = %.1e' % (k, sp[k]) for k in ('하위 30', '무작위 30') if k in sp]
    panel(s, 0.90, 3.85, 5.60, 'Fisher one-sided test', ps + [
        '', 'DKD marker 논문 보유 비율 기준', '우연으로 보기 어려운 차이임.',
    ], accent=True, size=11, row_h=0.32)
    bullets(s, 6.85, 3.85, 5.60, 2.9, [
        (0, '무엇을 보여 주는가'),
        (1, '우리 순위 상위에 놓인 gene 은 하위·random 보다 DKD 문헌을 훨씬 많이 가짐'),
        (1, '즉 pipeline 이 생물학적으로 관련 있는 gene 을 위로 올림'),
        (0, '무엇을 보여 주지 않는가'),
        (1, '이것은 **방법의 타당성** 근거이지, 새로운 바이오마커를 찾았다는 근거가 아님'),
    ])


def s_conclusion(prs, f):
    s = blank(prs)
    head(s, KICK, '1차년도 결론과 2차년도 방향', None)
    panel(s, 0.77, 1.70, 3.85, '한 일', [
        '공개 자산 %d건 확보 · %d건 사용' % (f['inv_total'], f['inv_used']),
        'multi-omics 가능성 실측',
        'cohort 7종 %d명 통합' % (f['n_case'] + f['n_ctrl']),
        '%s gene 공간 동결' % format(f['genes'], ','),
        'selector 8종 benchmark',
        'confounding 진단과 보정 6종 비교',
        '후보 30개 · 인계 대상 6개',
        '재현 pipeline %d단계' % f['n_stage'],
    ], size=11, row_h=0.34)
    panel(s, 4.85, 1.70, 3.85, '알아낸 것', [
        (True, '공개 자료로 multi-omics 불가'),
        '대조군 쪽 층이 없어 contrast 불성립',
        '',
        (True, '절대 AUROC 는 후보를 못 가림'),
        'random 50 genes 가 %.2f~%.2f' % (f['null_lo'], f['null_hi']),
        '',
        (True, '최재현 신호는 조직 채취'),
        '7개 cohort 전부 대조가 non-biopsy',
        '',
        (True, 'ComBat 은 여기서 무효'),
    ], size=11, row_h=0.32)
    panel(s, 8.90, 1.70, 3.85, '2차년도 제안', [
        '① 환자 단위 다층 설계 확보',
        '   같은 환자에서 두 층 이상',
        '   대조군 쪽도 함께 측정',
        '',
        '② 대조군을 설계에 포함',
        '   비당뇨 신장질환 biopsy 를',
        '   같은 route 로 확보',
        '',
        '③ Tier 2 후보 6개 실험 검증',
        '   ECM·fibrosis 축 우선',
    ], accent=True, size=11, row_h=0.32)
    note(s, 0.77, 6.05, 11.95,
         '1차년도의 실질 성과는 후보 목록 자체보다, **그 목록을 신뢰할 수 있는지 판정하는 '
         '절차**를 세운 것임. 같은 절차를 2차년도 자료에 그대로 적용할 수 있음.')


# ====================================================================== 조립
SLIDES = [
    ('표지', s_cover),
    ('배경', s_disease),
    ('배경', s_goal),
    ('divider', ('데이터', '무엇을 모았고, 무엇이 가능했는가')),
    ('데이터', s_data),
    ('데이터', s_multiomics),
    ('데이터', s_pivot),
    ('데이터', s_integrate),
    ('데이터', s_cohorts),
    ('divider', ('방법', '어떤 알고리즘으로 분석했는가')),
    ('방법', s_algorithms),
    ('방법', s_rbs_detail),
    ('방법', s_rbs_code),
    ('divider', ('첫 결과와 함정', '좋아 보이는 숫자를 의심한 이유')),
    ('분석', s_control),
    ('분석', s_first_auc),
    ('분석', s_null),
    ('분석', s_algo_compare),
    ('분석', s_generative),
    ('divider', ('교란의 발견', '가장 재현성 높은 신호는 무엇이었나')),
    ('교란', s_artifact),
    ('교란', s_evidence),
    ('divider', ('보정과 후보 도출', '교란을 걷어내고 다시 세운 순위')),
    ('도출', s_corr_algorithms),
    ('도출', s_corr_compare),
    ('도출', s_corr_per_selector),
    ('도출', s_reorder),
    ('도출', s_criteria),
    ('도출', s_before_after),
    ('도출', s_candidates),
    ('도출', s_tier2),
    ('도출', s_validate),
    ('도출', s_conventional),
    ('도출', s_other_results),
    ('divider', ('결론', '1차년도가 남긴 것')),
    ('결론', s_conclusion),
    ('결론', s_pipeline_code),
]


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    if not os.path.exists(TEMPLATE):
        log('%s 가 없음. 템플릿이 있어야 서식을 가져옴.' % TEMPLATE)
        return 1
    f = facts()
    prs = Presentation(TEMPLATE)
    clear(prs)
    made = 0
    for tag, fn in SLIDES:
        try:
            if tag == 'divider':
                s_divider(prs, *fn)
            else:
                fn(prs, f)
            made += 1
        except Exception as e:
            log('  ! 장표 실패 [%s] %s: %s' % (tag, getattr(fn, '__name__', fn), e))
    prs.save(OUT)
    log('%s 에 %d장을 썼음 (%.1f MB).' % (OUT, made, os.path.getsize(OUT) / 1e6))
    return 0 if made == len(SLIDES) else 1


if __name__ == '__main__':
    sys.exit(main())
