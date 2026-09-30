#!/usr/bin/env python
"""Deck slides for the three extensions: single-cell, confounder-aware selection, germline."""
from make_ppt import (add_slide, header, txt, rule, table, picture,
                      INK, MUTED, ACCENT, BLUE)


def s_singlecell(prs, d):
    s = add_slide(prs)
    header(s, 'EXTENSION 1', '단일세포 — 아티팩트를 자체 음성 대조군과 함께 확증',
           'KPMP h5ad 원본 사용. API는 dmr 슬라이스에 통계를 아예 주지 않습니다. '
           '도너 단위 pseudobulk(817개 프로파일), 세포 단위 검정 아님.')
    picture(s, 'F10_singlecell_control.png', 0.6, 2.1, w=11.6)
    rows = [['대비', '비교 대상', 'IEG 평균 g', 'IEG vs 랜덤'],
            ['DKD(27) vs reference(12)', '생검 vs 비생검', '-0.843', 'p = 6.0e-76'],
            ['DKD(27) vs 비당뇨CKD(10)', '생검 vs 생검', '+0.082', 'p = 0.21 (n.s.)'],
            ['reference(12) vs nephrectomy(6)', '비생검 vs 비생검', '-0.133', 'p = 0.89 (n.s.)']]
    table(s, 0.6, 5.75, 7.4, rows, col_w=[3.0, 2.0, 1.3, 1.5], size=10)
    txt(s, 8.3, 5.85, 4.5, 1.4, [
        'KPMP에서 DKD는 생검, 두 대조군은 모두 비생검입니다.',
        '모듈이 생검-vs-비생검만 가르고 그 외에는 아무것도 가르지 않습니다.',
        '',
        '후보는 정직하게 나쁩니다: 교란 통제 대비에서 +0.145(5%)로 랜덤(+0.127, 4%)과 차이 없음.'],
        size=10, color=MUTED, line_spacing=1.22)


def s_deconfound(prs, d):
    s = add_slide(prs)
    header(s, 'EXTENSION 2', '교란 보정을 선택 목적함수 안으로',
           '기존에는 사후 게이트였고 비당뇨 CKD 생검이 필요해 11개 중 2개 데이터셋에서만 가능했습니다. '
           'IEG 점수는 발현 행렬만으로 계산됩니다.')
    picture(s, 'F11_deconfound.png', 0.6, 2.15, w=12.0)
    rows = [['arm', 'K', '채취 플래그', '특이성 통과', 'AUROC', '교차폴드 J'],
            ['보정 없음', '10', '55%', '18%', '0.932', '0.228'],
            ['잔차화(orth)', '10', '5%', '42%', '0.860', '0.204'],
            ['보정 없음', '50', '28%', '32%', '0.952', '0.303'],
            ['잔차화(orth)', '50', '6%', '42%', '0.936', '0.248']]
    table(s, 0.6, 5.5, 6.8, rows, col_w=[1.8, 0.7, 1.3, 1.3, 1.1, 1.3], size=10,
          hl_rows=(2, 4))
    txt(s, 7.7, 5.6, 5.1, 1.6, [
        'K=10에서 채취 플래그 55%→5%, 특이성 통과 18%→42%.',
        '',
        '비용도 명시: AUROC 0.932→0.860(K=10), 교차폴드 0.303→0.248(K=50). '
        'IEG 점수가 case status와 상관되므로 실제 질환 신호도 일부 제거됩니다 — 공짜가 아닙니다.',
        '',
        '소프트 버전(penalty)은 더 나쁩니다: 채취 제거는 덜 되고 AUROC는 0.75까지 떨어집니다.'],
        size=9.5, color=MUTED, line_spacing=1.2)


def s_gwas(prs, d):
    s = add_slide(prs)
    header(s, 'EXTENSION 3', 'Germline 레이어 — 명확한 음성',
           'GCST90179152, 한국인 DKD (2,532 cases / 31,347 controls). 변이 4,976,559개. '
           'p값·좌표만 있어 위치 기반 유전자 수준 검정만 가능합니다.')
    rows = [['항목', '값'],
            ['empirical p < 0.05 후보', '1 / 30  (우연 기대 1.5)'],
            ['genome-wide 유의(5e-8) 도달', '0'],
            ['0.05 꼬리 이항검정', 'p = 0.785'],
            ['균등분포 대비 KS 검정', 'D = 0.120,  p = 0.736']]
    table(s, 0.6, 2.2, 6.4, rows, col_w=[3.0, 2.6], size=11.5)
    txt(s, 0.6, 4.0, 6.4, 1.6, [
        '두 가지 함정을 처리했습니다:',
        '· min-p는 창 안 SNP 수에 편향 → 배경 유전자를 SNP 개수로 매칭',
        '· LD로 변이가 독립이 아님 → 매칭된 경험적 널이 1차적으로 흡수',
        '  (MAGMA가 아니며 그렇게 제시하지 않습니다)'],
        size=10.5, line_spacing=1.3)
    txt(s, 7.3, 2.25, 5.5, 0.4, '순위 상위와 Tier 1', size=13, bold=True, color=BLUE)
    rows2 = [['유전자', 'Tier', 'min p', 'empirical p'],
             ['TNNT2', '4', '8.1e-04', '0.049'],
             ['C1orf21', '4', '7.4e-04', '0.050'],
             ['PRSS23', '1', '7.3e-04', '0.094'],
             ['OLFML3', '1', '1.7e-01', '0.926']]
    table(s, 7.3, 2.7, 5.5, rows2, col_w=[1.6, 0.8, 1.3, 1.4], size=10.5, hl_rows=(3, 4))
    txt(s, 7.3, 4.4, 5.5, 1.8, [
        '해석은 부풀리지 않습니다. 2,532 cases면 검정력이 낮고, 발현 기반 후보가 common '
        'germline risk를 가질 의무도 없습니다 — 신장 형질의 GWAS 신호는 대부분 조절성이고 '
        '작용 유전자와 멀리 떨어져 있습니다.',
        '',
        '이 레이어는 지지도 반박도 추가하지 않습니다.'],
        size=10, color=MUTED, line_spacing=1.25)
    rule(s, 0.6, 6.3, 12.2)
    txt(s, 0.6, 6.45, 12.2, 0.7, [
        '세 확장의 순효과: 방법론(아티팩트 3자원 3플랫폼 복제 + 보정을 목적함수로 일반화)은 강해졌고, '
        '후보(교란 통제 대비에서 랜덤과 동등, germline 지지 없음, GSE175759에서 부호 반전)는 약해졌습니다.',
        '→ 논문은 방법론을 앞세우고, 후보는 파이프라인이 자기 출력을 승격시키지 않는 예시로 씁니다.'],
        size=11, line_spacing=1.3)


EXTRA_SLIDES2 = [s_singlecell, s_deconfound, s_gwas]
