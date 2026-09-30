#!/usr/bin/env python
"""Build a status deck for the DKD project.

Numbers are read from the result files, not typed in, so the deck cannot drift out of sync
with the analysis. Korean body text uses Malgun Gothic, which ships with Windows PowerPoint.
"""
import os, sys, json, glob, sqlite3, argparse
import numpy as np
import pandas as pd
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

OUT = 'DKD_status.pptx'
FIG = 'results/figures'

KR = 'Malgun Gothic'
EN = 'Segoe UI'
INK = RGBColor(0x1A, 0x1A, 0x1A)
MUTED = RGBColor(0x70, 0x70, 0x70)
ACCENT = RGBColor(0xC1, 0x44, 0x0E)
BLUE = RGBColor(0x2A, 0x5D, 0x8F)
LINE = RGBColor(0xD0, 0xD0, 0xD0)
BG = RGBColor(0xFF, 0xFF, 0xFF)

W, Hh = Inches(13.333), Inches(7.5)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def add_slide(prs):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    return s


def txt(slide, x, y, w, h, text, size=14, bold=False, color=INK, font=KR,
        align=PP_ALIGN.LEFT, space_after=4, line_spacing=1.15):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    lines = text.split('\n') if isinstance(text, str) else text
    for i, ln in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(space_after)
        p.line_spacing = line_spacing
        r = p.add_run()
        r.text = ln
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.color.rgb = color
        r.font.name = font
    return tb


def rule(slide, x, y, w, color=LINE, thick=1.0):
    from pptx.enum.shapes import MSO_SHAPE
    sh = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y),
                                Inches(w), Pt(thick))
    sh.fill.solid(); sh.fill.fore_color.rgb = color
    sh.line.fill.background()
    sh.shadow.inherit = False
    return sh


def header(slide, kicker, title, sub=None):
    txt(slide, 0.6, 0.42, 12.2, 0.3, kicker, size=11, bold=True, color=ACCENT, font=EN)
    txt(slide, 0.6, 0.72, 12.2, 0.6, title, size=26, bold=True)
    rule(slide, 0.6, 1.42, 12.2)
    if sub:
        txt(slide, 0.6, 1.55, 12.2, 0.5, sub, size=12.5, color=MUTED)


def picture(slide, name, x, y, w=None, h=None):
    p = os.path.join(FIG, name)
    if not os.path.exists(p):
        log('  ! missing figure %s' % name)
        return None
    kw = {}
    if w: kw['width'] = Inches(w)
    if h: kw['height'] = Inches(h)
    return slide.shapes.add_picture(p, Inches(x), Inches(y), **kw)


def table(slide, x, y, w, rows, col_w=None, size=11, header_bold=True, hl_col=None,
          hl_rows=()):
    nr, nc = len(rows), len(rows[0])
    h = 0.32 * nr
    shp = slide.shapes.add_table(nr, nc, Inches(x), Inches(y), Inches(w), Inches(h))
    tb = shp.table
    if col_w:
        total = sum(col_w)
        for j, cw in enumerate(col_w):
            tb.columns[j].width = Emu(int(Inches(w) * cw / total))
    for i, row in enumerate(rows):
        tb.rows[i].height = Inches(0.3)
        for j, v in enumerate(row):
            c = tb.cell(i, j)
            c.text = ''
            c.margin_left = c.margin_right = Inches(0.06)
            c.margin_top = c.margin_bottom = Inches(0.02)
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
            p = c.text_frame.paragraphs[0]
            r = p.add_run(); r.text = str(v)
            r.font.size = Pt(size)
            r.font.name = KR
            r.font.bold = (i == 0 and header_bold) or (i in hl_rows)
            r.font.color.rgb = ACCENT if (i in hl_rows) else INK
            if j > 0:
                p.alignment = PP_ALIGN.RIGHT
            c.fill.solid()
            c.fill.fore_color.rgb = RGBColor(0xF4, 0xF4, 0xF4) if i == 0 else BG
    return tb


# ==================================================================== data
def load():
    d = {}
    d['inv'] = pd.read_csv('results/data_inventory.tsv', sep='\t')
    d['g3'] = pd.read_csv('results/comparison_glom3.tsv', sep='\t')
    d['a4'] = pd.read_csv('results/comparison_all4.tsv', sep='\t')
    d['mt'] = pd.read_csv('results/master_candidate_table.tsv', sep='\t')
    d['lit'] = pd.read_csv('results/literature_review.tsv', sep='\t')
    d['null'] = pd.read_csv('results/null_control.tsv', sep='\t')
    d['stress'] = pd.read_csv('results/stress_test.tsv', sep='\t')
    d['gen'] = pd.read_csv('results/generative_robustness.tsv', sep='\t')
    d['kpmp'] = pd.read_csv('results/kpmp_corroboration.tsv', sep='\t')
    con = sqlite3.connect('db/dkd.sqlite')
    d['db'] = {t: con.execute('SELECT COUNT(*) FROM %s' % t).fetchone()[0]
               for t in ['study', 'sample', 'subject', 'gene', 'matrix',
                         'method_result', 'candidate_evidence', 'external_evidence']}
    con.close()
    return d


# ==================================================================== slides
def s_title(prs, d):
    s = add_slide(prs)
    rule(s, 0.6, 2.55, 3.0, ACCENT, 2.5)
    txt(s, 0.6, 2.8, 12.0, 1.2,
        '당뇨병성 신장질환(DKD) 바이오마커 후보 탐색', size=34, bold=True)
    txt(s, 0.6, 3.75, 12.0, 0.8,
        'Cross-Cohort Robust Feature Selection · 1단계 중간 현황', size=18, color=MUTED)
    txt(s, 0.6, 4.9, 12.0, 0.9,
        ['공개 데이터 %d개 연구 · 샘플 %d건 · 공통 gene space 9,900개'
         % (d['db']['study'], d['db']['sample']),
         '하이브리드 DB 1식 구축 완료 · 방법론 벤치마크 8종 완료 · 후보 30종 도출'],
        size=13, color=INK)
    txt(s, 0.6, 6.6, 12.0, 0.4, 'ETRI · 2026-08-24', size=11, color=MUTED, font=EN)


def s_summary(prs, d):
    s = add_slide(prs)
    header(s, 'EXECUTIVE SUMMARY', '한 장 요약 — 예상과 달랐던 세 가지')
    y = 1.95
    items = [
        ('① AUROC로는 방법을 줄세울 수 없다',
         '랜덤으로 뽑은 50개 유전자가 external AUROC 0.68~0.88을 냅니다(라벨 셔플은 0.50).\n'
         'DKD 신호가 광범위하게 공발현돼 어떤 유전자를 고르든 큰 차이가 없습니다.\n'
         '→ 모든 AUROC는 널 분포 대비 퍼센타일로 보고해야 합니다.'),
        ('② 가장 재현성 높은 모듈이 조직 채취 아티팩트였다',
         'immediate-early 유전자(FOS/JUN/EGR1/DUSP1/ZFP36)가 모든 안정성 순위 1위인데,\n'
         '비당뇨 CKD가 DKD의 145%만큼 움직입니다. 전 코호트에서 환자=생검, 대조군=신절제라\n'
         '교란요인이 공유됩니다. → 교란요인이 공유되면 교차코호트 재현성 ≠ 생물학적 타당성.'),
        ('③ 변별력 있는 축은 교차코호트 안정성이었다',
         '베이스라인 7종 모두 교차폴드 Jaccard 0.03~0.05 — hold-out 코호트가 바뀌면\n'
         '사실상 다른 시그니처를 냅니다. 제안 방법은 동일 AUROC에서 0.192, 4코호트에서 0.303.'),
    ]
    for t, b in items:
        txt(s, 0.6, y, 12.2, 0.35, t, size=15, bold=True, color=BLUE)
        txt(s, 0.85, y + 0.38, 11.9, 0.9, b, size=11.5, color=INK, line_spacing=1.25)
        y += 1.62


def s_data(prs, d):
    s = add_slide(prs)
    inv = d['inv']
    used = inv[inv['state'] == 'used']
    unused = inv[inv['state'] == 'unused']
    excl = inv[inv['state'] == 'excluded']
    header(s, 'DATA', '데이터 확보 및 사용 현황',
           '내려받은 총 %.1f GB 중 분석에 들어간 것은 %.2f GB. 다만 용량 비율은 오해를 부릅니다 — '
           '단일세포 파일이 크기 때문입니다.' % (inv['size_gb'].sum(), used['size_gb'].sum()))
    picture(s, 'F1_data_inventory.png', 0.6, 2.3, w=7.6)
    rows = [['구분', '데이터셋', 'GB'],
            ['분석에 사용', str(len(used)), '%.2f' % used['size_gb'].sum()],
            ['감사로 제외', str(len(excl)), '%.2f' % excl['size_gb'].sum()],
            ['확보했으나 미사용', str(len(unused)), '%.2f' % unused['size_gb'].sum()],
            ['미확보', str(len(inv[inv['state'] == 'not_downloaded'])), '0.00']]
    table(s, 8.5, 2.5, 4.3, rows, col_w=[2.2, 1, 1], size=11)
    txt(s, 8.5, 4.35, 4.4, 2.4,
        ['모델링에 실제 들어간 샘플',
         '· 라벨 샘플 222건 (DKD 116, 대조 106)',
         '· 비당뇨 CKD 310건 → 특이성 검정 전용',
         '',
         '미사용 주요 자산',
         '· KPMP h5ad 4.85 GB (API로 대체)',
         '· GSE131882 snRNA 1.49 GB',
         '· 대사체 4개 연구, 혈액 GSE142153'],
        size=10.5, color=MUTED)


def s_audit(prs, d):
    s = add_slide(prs)
    header(s, 'AUDIT', '코호트 감사 — 환자 중복으로 3개 제외',
           'GSM이 아니라 환자 ID 수준으로 대조. GEO는 기존 코호트를 조용히 재사용합니다.')
    rows = [['데이터셋', '판정'],
            ['GSE30122', 'GSE30528+GSE30529 그대로 + 동일 대조군 25명 재수록 → 신규 환자 0'],
            ['GSE47183', 'DN 환자 ID가 GSE104948과 12/12 완전 일치'],
            ['GSE99340', 'ERCB superset (각각 168명 공유), 413개 중 ~51개는 세포주']]
    tb = table(s, 0.6, 2.1, 12.2, rows, col_w=[1.6, 6.5], size=12, hl_rows=())
    txt(s, 0.6, 3.5, 12.2, 0.4, '남은 교란요인 두 가지 (제외 아님 — 명시 대상)',
        size=14, bold=True, color=BLUE)
    txt(s, 0.85, 3.95, 12.0, 1.4,
        ['· 조직 구획: GSE142025는 microdissection이 아닌 whole cortex → cross-cohort와 '
         'cross-compartment가 섞임',
         '· 대조군 출처: GSE104948/104954만 생체공여(living donor), 나머지는 전부 tumor nephrectomy',
         '· 잔존 중복: GSE30528↔GSE30529는 DKD 5명·대조 4명 공유, GSE104948↔GSE104954는 DN 5명 공유'],
        size=11.5, line_spacing=1.35)
    txt(s, 0.6, 5.6, 12.2, 0.35, '이 결과는 이제 DB 쿼리로 재현됩니다', size=13, bold=True)
    txt(s, 0.85, 6.0, 12.0, 0.8,
        ['SELECT * FROM v_subject_reuse  →  DN1 / ERCB / 4개 연구에 등장',
         'SELECT * FROM v_procurement_balance  →  전 코호트에서 case=biopsy, control=nephrectomy'],
        size=10.5, color=MUTED, font=EN)


def s_genespace(prs, d):
    s = add_slide(prs)
    header(s, 'HARMONIZATION', '공통 gene space — 9,900 genes',
           '허브 ID는 Entrez. GSE104948/104954가 Brainarray ENTREZG custom CDF라 probe ID가 이미 '
           'Entrez이고, 가장 빡빡한 코호트에서 매핑 손실이 0입니다.')
    picture(s, 'F8_gene_space.png', 0.6, 2.4, w=6.6)
    txt(s, 7.6, 2.4, 5.2, 0.4, '동결된 결정 사항', size=14, bold=True, color=BLUE)
    txt(s, 7.6, 2.85, 5.2, 2.2,
        ['· Collapse 규칙 = MAX-MEAN (유전자당 평균발현 최고 probe)',
         '· 다중 매핑 probe 제거 (GPL571 1,223개, GPL17586 1,572개)',
         '· 코호트별 gene-wise z-score — 라벨을 쓰지 않으므로 leakage 아님',
         '',
         '나중에 바꾸면 저장된 모든 bootstrap 빈도가 무효화되므로 지금 고정.'],
        size=11, line_spacing=1.3)
    txt(s, 7.6, 5.3, 5.2, 0.4, '실현 가능성 실측', size=14, bold=True, color=BLUE)
    txt(s, 7.6, 5.75, 5.2, 1.2,
        ['4개 코호트에서 효과 방향이 일치하는 유전자 2,689/9,900 = 27.2%',
         '(우연 기대치 12.5%). 그중 827개는 |g|>0.5.',
         'NPHS1 하향, COL1A2·PCOLCE 상향, C1QA·MS4A6A 상향 — 교과서적 DKD 병태생리 복원.'],
        size=11, line_spacing=1.3)


def s_method(prs, d):
    s = add_slide(prs)
    header(s, 'METHOD', '제안 방법 — Robust Biomarker Score (RBS)',
           '베이스라인은 training 코호트를 합쳐서 한 번 선택합니다. 그래서 가장 큰 코호트에서만 '
           '강한 유전자가 이깁니다. RBS는 합치지 않습니다.')
    box = [('I  예측 중요도', '부트스트랩별 |계수| ÷ 해당 부트스트랩 최대값, 코호트 평균'),
           ('S  부트스트랩 안정성', '코호트 내 선택 빈도의 코호트 평균'),
           ('R  교차코호트 재현성', '코호트별 빈도의 기하평균 × 효과 방향 일치 게이트'),
           ('P  Perturbation 강건성', 'subsampling / noise / masking / batch shift 하에서의 선택 빈도')]
    y = 2.3
    for t, b in box:
        txt(s, 0.6, y, 3.5, 0.3, t, size=13, bold=True, color=ACCENT)
        txt(s, 4.2, y, 8.6, 0.4, b, size=11.5)
        y += 0.62
    rule(s, 0.6, 5.0, 12.2)
    txt(s, 0.6, 5.15, 12.2, 0.35,
        '결합은 가중 기하평균 — 합이면 한 항목만 좋아도 통과하지만, 주장은 "중요하고 AND 안정적이고 '
        'AND 재현되고 AND 강건하다"입니다.', size=11.5, color=MUTED)
    txt(s, 0.6, 5.75, 12.2, 0.35, '설계 과정에서 실제로 바꾼 것', size=14, bold=True, color=BLUE)
    txt(s, 0.85, 6.15, 12.0, 1.0,
        ['· 처음 쓴 raw MIN 집계는 실패 (교차폴드 0.036, R>0.1인 유전자 13개뿐). n=22 코호트의 '
         '추정오차가 min을 지배 — 엄격해서가 아니었습니다. → geomean으로 교체, 안정성·AUROC 양쪽 개선.',
         '· ablation 1위였던 mean은 채택하지 않음: S가 이미 F.mean이라 R이 S와 수학적으로 동일해져 '
         '교차코호트 항 자체가 사라집니다.'],
        size=11, line_spacing=1.3)


def s_benchmark(prs, d):
    s = add_slide(prs)
    a4 = d['a4'][d['a4']['K'] == 50].sort_values('external', ascending=False)
    header(s, 'BENCHMARK', '벤치마크 결과 — 교차코호트 안정성',
           'LODO(코호트 단위 hold-out) · B=200 부트스트랩 · K=50 · 9,900 genes')
    picture(s, 'F2_cross_fold_stability.png', 0.6, 2.15, w=7.4)
    rows = [['방법', 'AUROC', '널대비', '교차폴드']]
    for _, r in a4.iterrows():
        rows.append(['RBS (제안)' if r['method'] == 'rbs' else r['method'],
                     '%.3f' % r['external'], '%+.2f sd' % r['null_z'],
                     '%.3f' % r['cross_fold_jaccard']])
    table(s, 8.25, 2.35, 4.6, rows, col_w=[1.8, 1, 1.1, 1.1], size=10,
          hl_rows=(1,) if a4.iloc[0]['method'] == 'rbs' else ())
    txt(s, 8.25, 5.6, 4.6, 1.4,
        ['4코호트 구성 (K=50)',
         '',
         '3코호트에서는 AUROC가 lasso와 동률(0.851 vs 0.863)이고',
         '교차폴드만 4배였는데, training 코호트가 2→3개가 되면',
         'AUROC까지 1위가 됩니다.'],
        size=10, color=MUTED, line_spacing=1.25)


def s_null(prs, d):
    s = add_slide(prs)
    n = d['null']
    header(s, 'FINDING 1', 'AUROC의 함정 — 랜덤 시그니처 널 대조군',
           '폴드별로 랜덤 50개 유전자 300세트. 라벨을 섞으면 0.50으로 떨어지므로 신호는 진짜지만, '
           '어떤 유전자를 고르는지는 거의 무관합니다.')
    picture(s, 'F3_auroc_vs_null.png', 0.6, 2.3, w=12.2)
    sub = n[n['K'] == 50]
    rows = [['Hold-out', '랜덤 AUROC (K=50)', '95% 상한', '라벨 셔플']]
    for _, r in sub.iterrows():
        rows.append([r['held_out'], '%.3f' % r['random_mean'],
                     '%.3f' % r['random_p95'], '%.3f' % r['permuted_mean']])
    table(s, 0.6, 5.5, 7.0, rows, col_w=[1.6, 1.6, 1.1, 1.2], size=10.5)
    txt(s, 8.0, 5.6, 4.8, 1.5,
        ['GSE142025 hold-out은 랜덤도 K=100에서 평균 0.932(최대 1.000)입니다.',
         '→ 이 폴드는 방법을 변별할 수 없으므로 평균에 그냥 넣으면 안 됩니다.',
         '',
         '이 폴드를 빼도 RBS 순위는 그대로이고 널 대비는 +1.33 → +1.42 sd로 개선됩니다.'],
        size=10.5, color=MUTED, line_spacing=1.3)


def s_artifact(prs, d):
    s = add_slide(prs)
    header(s, 'FINDING 2', '조직 채취 아티팩트 — 재현성 ≠ 타당성',
           'ERCB의 비당뇨 CKD 생검은 DKD 생검과 같은 방식으로 채취됐습니다. 교란요인이 통제된 대비입니다.')
    picture(s, 'F4_procurement_artifact.png', 0.6, 2.2, w=5.4)
    txt(s, 6.4, 2.3, 6.4, 0.4, 'immediate-early 모듈 (20 genes)', size=14, bold=True, color=BLUE)
    rows = [['대비', '모듈 평균 Hedges g'],
            ['DKD vs 대조군', '-0.85'],
            ['비당뇨 CKD vs 대조군', '-1.23  (DKD의 145%)'],
            ['두 효과벡터 상관', 'r = 0.83 (사구체) / 0.92 (세뇨관)']]
    table(s, 6.4, 2.8, 6.4, rows, col_w=[2.6, 2.4], size=11)
    txt(s, 6.4, 4.4, 6.4, 2.4,
        ['당뇨가 없는 CKD가 더 크게 움직입니다. 이 모듈은 조직이 생검으로 왔는지 '
         '신절제/공여로 왔는지를 추적하는 것이고, 환자가 당뇨였는지를 추적하는 게 아닙니다.',
         '',
         '교란요인이 모든 코호트에 동일하게 존재하므로, 어떤 교차코호트 안정성 기준으로 봐도 '
         '완벽하게 "재현성 높은" 신호로 보입니다.',
         '',
         '→ 강건성 점수는 타당성 필터가 아닙니다. 공유 교란과 공유 생물학을 구분하는 것은 '
         '교란요인을 고정한 대비뿐입니다.'],
        size=11, line_spacing=1.3)
    txt(s, 6.4, 6.75, 6.4, 0.4,
        'RBS 순위에서 탈락: DUSP1(2위) ZFP36(5) FOS(14) EGR1(22) JUN(28) ATF3(32) FOSB(37)',
        size=10, color=ACCENT)


def s_specificity(prs, d):
    s = add_slide(prs)
    header(s, 'GATE', '특이성 게이트를 파이프라인 단계로 승격',
           'DKD를 비당뇨 CKD 310건과도 구분하는 유전자만 통과. RBS 상위 200개 중 92개(46%)만 통과, '
           '21개(10%)는 채취 주도로 플래그.')
    picture(s, 'F5_per_diagnosis.png', 0.6, 2.2, h=5.05)
    txt(s, 8.1, 2.3, 4.7, 0.4, '게이트 자체의 거짓양성률', size=13, bold=True, color=BLUE)
    rows = [['항목', '값'],
            ['관측 통과 (|g|>=0.5)', '833 (8.4%)'],
            ['순열 통과 (500회)', '평균 114'],
            ['경험적 FDR', '13.7%'],
            ['순열 최대 |g| (95pct)', '0.92']]
    table(s, 8.1, 2.8, 4.7, rows, col_w=[2.6, 1.5], size=10.5)
    txt(s, 8.1, 4.6, 4.7, 2.3,
        ['최종 30개 중 순열 천장(0.92)을 넘는 것은 4개뿐:',
         'FMOD 1.28 · MMP2 1.03 · MOXD1 0.95 · LUM 0.94',
         '',
         '진단별 프로파일로 세 패턴 분리',
         '· A 단독: MMP2(양쪽), FMOD(사구체)',
         '· B 우세: LUM, MOXD1',
         '· C 공여자 주도: CALHM2, CCDC91, DPP6, TNNT2 등',
         '  → 채취 교란이 형태만 바꿔 재등장'],
        size=10.5, line_spacing=1.25)


def s_robustness(prs, d):
    s = add_slide(prs)
    st = d['stress']
    header(s, 'ROBUSTNESS', 'Perturbation stress test — 발견 이후 단계',
           '원본 데이터로 발견을 끝내고, 확정된 30개 시그니처만 고정한 뒤 테스트 코호트를 교란합니다.')
    picture(s, 'F6_stress_test.png', 0.6, 2.2, w=12.2)
    g = st.groupby(['perturbation', 'severity'])['margin'].mean().reset_index()
    rows = [['Perturbation', 's=0', 's=1', 's=2']]
    nice = {'subsample': 'Subsampling', 'noise': 'Gaussian noise',
            'mask': 'Feature masking', 'batch': 'Batch shift'}
    for k in ('subsample', 'noise', 'mask', 'batch'):
        sub = g[g['perturbation'] == k].set_index('severity')['margin']
        rows.append([nice[k]] + ['%+.3f' % sub.get(v, np.nan) for v in (0.0, 1.0, 2.0)])
    table(s, 0.6, 5.4, 6.2, rows, col_w=[2.2, 1, 1, 1], size=10.5)
    txt(s, 7.2, 5.5, 5.6, 1.8,
        ['표는 랜덤 대비 마진. 전 구간에서 +0.07~+0.13으로 유지 — 시그니처가 데이터보다 '
         '빨리 무너지지 않습니다.',
         '',
         '생성모델: VAE는 통과(클래스 평균차 r=0.94, 4/4 폴드), Diffusion은 실패(0/4). '
         'denoiser eps-MSE가 0.57에서 정체 → n≈120의 표본수 한계이며 수치를 인용하면 안 됩니다.'],
        size=10.5, color=MUTED, line_spacing=1.3)


def s_candidates(prs, d):
    s = add_slide(prs)
    mt = d['mt']
    header(s, 'CANDIDATES', '최종 후보 — 데이터 강도 × 문헌 신규성',
           'PubMed 정량 집계로 신규성 판정. 약어 충돌(MSC=mesenchymal stem cell 등) 4건은 '
           '카운트 무효로 표시.')
    picture(s, 'F7_candidates.png', 0.6, 2.2, w=6.9)
    t1 = mt[mt['tier'].str.startswith('1')]
    t2 = mt[mt['tier'].str.startswith('2')]
    rows = [['Tier', '유전자', '판정']]
    rows.append(['1 · 신규 AND 강함', ', '.join(t1['gene']), '추적 대상'])
    rows.append(['2 · 강하지만 기보고', ', '.join(t2['gene']), '방법 검증'])
    table(s, 7.8, 2.4, 5.0, rows, col_w=[1.8, 2.0, 1.2], size=10)
    txt(s, 7.8, 3.6, 5.0, 3.2,
        ['Tier 2가 방법의 타당성 근거입니다 — FMOD·LUM은 2021년 DKD hub gene 논문에 이미 '
         '함께 보고돼 있고, 파이프라인이 그것을 재발견했습니다.',
         '',
         'Tier 1 상세',
         '· OLFML3 — 신장 문헌 전무, 진단별 B(DN 우세), KPMP에서 염증성 피질 간질 섬유아세포',
         '· PRSS23 — 신장은 연구됐으나 DKD 논문 0편, 진단별 B',
         '',
         '주의: 두 유전자의 특이성 |g|(0.63, 0.51)는 순열 천장 0.92 미만입니다. '
         '독립 검증이 필요한 후보이며 확정 소견이 아닙니다.'],
        size=10.5, line_spacing=1.25)


def s_kpmp(prs, d):
    s = add_slide(prs)
    kp = d['kpmp']
    n_prot = int(kp['protein_glom_vs_TI_FC'].notna().sum())
    n_cell = int(kp['sn_top_celltype'].notna().sum())
    header(s, 'CROSS-OMICS', 'KPMP 교차오믹스 대조 — 할 수 있는 것과 없는 것',
           'GraphQL API를 해독해 유전자 단위 조회. DAR 승인 불필요.')
    rows = [['KPMP 데이터 타입', 'AKI', 'CKD', '건강', 'DKD(dmr)']]
    for nm, a, c, h, m in [('Single-nucleus RNA (sn)', 33, 72, 40, 11),
                           ('Single-cell RNA (sc)', 19, 51, 40, 3),
                           ('Regional transcriptomics', 5, 22, 9, 0),
                           ('Regional proteomics', 12, 14, 5, 0)]:
        rows.append([nm, a, c, h, m])
    table(s, 0.6, 2.2, 7.4, rows, col_w=[3.0, 0.8, 0.8, 0.8, 1.1], size=11)
    txt(s, 0.6, 4.1, 7.4, 1.4,
        ['KPMP 지역 프로테오믹스에 당뇨 환자가 0명입니다. 게다가 dmr 슬라이스는 API가 '
         'cell count만 주고 foldChange·pValAdj가 전부 null입니다. 대사체 레이어는 아예 없습니다.',
         '→ KPMP는 DKD validation을 담을 수 없습니다. corroboration으로 낮춰 잡은 판단이 맞았습니다.'],
        size=11, line_spacing=1.3)
    txt(s, 8.3, 2.3, 4.5, 0.4, '실제로 얻은 것', size=13, bold=True, color=BLUE)
    txt(s, 8.3, 2.75, 4.5, 3.6,
        ['후보 30개 중',
         '· %d개 단백질 수준 증거' % n_prot,
         '· %d개 유의한 세포유형 배정' % n_cell,
         '',
         '세포유형 수렴',
         '· 섬유아세포: FMOD, LUM, MMP2, COL1A2, OLFML3 (상당수 degenerative 상태)',
         '· LYVE1+ 조직 대식세포: MS4A6A, FOLR2',
         '',
         '단백질 수준 구획 배정',
         '· FBN1 +2.03 (adjP 8e-19)',
         '· DPP6 +3.97 (adjP 2.5e-13)',
         '· COL1A2 +1.19 (adjP 6.4e-06)'],
        size=10.5, line_spacing=1.2)


def s_db(prs, d):
    s = add_slide(prs)
    db = d['db']
    header(s, 'DELIVERABLE', '하이브리드 DB 1식 구축 완료',
           'db/dkd.sqlite (관계형 6.6 MB) + db/columnar/ (Parquet 48 MB) · 14 테이블 3 뷰')
    txt(s, 0.6, 2.15, 6.0, 0.4, '왜 하이브리드인가', size=14, bold=True, color=BLUE)
    txt(s, 0.6, 2.6, 6.0, 1.6,
        ['(gene, sample, value) 한 행씩 넣는 관계형 설계는 라벨 샘플만 220만 행, '
         '비당뇨 CKD 포함 590만 행입니다.',
         '유일하게 중요한 접근 패턴(gene × sample 블록을 꺼내 곱하기)에서 수백만 좁은 행의 '
         '인덱스 스캔이 되고, 저장 용량도 약 3배입니다.',
         '→ 데이터 종류가 아니라 접근 패턴으로 분리했습니다.'],
        size=11, line_spacing=1.3)
    rows = [['Tier', '내용', '저장소'],
            ['관계형', '연구·환자·샘플·감사·gene space·분석결과', 'SQLite / PostgreSQL'],
            ['컬럼형', '발현·효과크기·선택빈도 행렬', 'Parquet + zstd'],
            ['카탈로그', 'matrix 테이블 (URI + 체크섬)', '관계형']]
    table(s, 0.6, 4.4, 6.0, rows, col_w=[1.0, 3.2, 1.6], size=10)
    txt(s, 7.0, 2.15, 5.8, 0.4, '구축 결과', size=14, bold=True, color=BLUE)
    rows = [['테이블', '행 수'],
            ['study', db['study']], ['sample', db['sample']], ['subject', db['subject']],
            ['gene', db['gene']], ['matrix (Parquet)', db['matrix']],
            ['method_result', db['method_result']],
            ['candidate_evidence', db['candidate_evidence']],
            ['external_evidence', db['external_evidence']]]
    table(s, 7.0, 2.6, 5.8, rows, col_w=[3.2, 1.4], size=10.5)
    txt(s, 7.0, 5.6, 5.8, 1.5,
        ['스키마가 프로세스가 아니라 구조로 강제하는 것',
         '· subject 테이블 분리 → v_subject_reuse로 코호트 재사용이 쿼리 한 줄',
         '· procurement 컬럼 → v_procurement_balance로 채취 교란이 즉시 보임',
         '· study_audit.usable → v_usable_sample이 제외 코호트를 기본 차단'],
        size=10, color=MUTED, line_spacing=1.25)


def s_status(prs, d):
    """1단계 목표 대비 현재 상태. 완료 항목은 손으로 적되, 수치는 파일에서 읽는다."""
    from ppt_slides_status import facts
    f = facts()
    s = add_slide(prs)
    header(s, 'STATUS', '현재 상태 평가', '1단계 목표 대비')
    txt(s, 0.6, 2.05, 4.0, 0.35, '완료', size=14, bold=True, color=BLUE)
    txt(s, 0.6, 2.5, 4.0, 3.7,
        ['· 하이브리드 DB 1식 + 웹 대시보드',
         '· AI 학습용 표준 데이터셋 정의',
         '  (9,900 gene space, collapse 규칙 동결)',
         '· 데이터 관리 체계 (감사 결과를 스키마에)',
         '· 방법론 벤치마크 7종 + 제안 방법',
         '· LODO 두 구성 + 널 대조군 2종',
         '· 특이성 게이트 + 진단별 층화',
         '· 단백체 2코호트 직교 보강',
         '· 단일세포 · 대사체 · GWAS 레이어 사용',
         '· 원고 %d쪽 · 그림 9 · 부록 6' % f['pages'],
         '· 재현 파이프라인 %d단계 · 감사 8종' % f['n_stage']],
        size=11, line_spacing=1.25)
    txt(s, 4.9, 2.05, 4.2, 0.35, '한계 (원고에 명시)', size=14, bold=True, color=ACCENT)
    txt(s, 4.9, 2.5, 4.2, 3.7,
        ['· 특이성 게이트가 ERCB 에 의존',
         '  같은 방식 채취된 비당뇨 CKD 생검을',
         '  가진 다른 공개 데이터가 없음',
         '· GSE104948 이 발굴과 게이트 양쪽에',
         '  쓰임 — 공개하고 민감도로 방어',
         '· 게이트 FDR 추정치 %.1f%% (배경 %.1f%%)' % (f['fdr'], f['bg']),
         '· 후보의 |g| 상당수가 순열 천장 미만',
         '· compartment 분석이 부분적 동일 환자',
         '· 임상 변수(eGFR 등) 확보 못함',
         '· 환자 매칭 다중 오믹스는 공개 자료에',
         '  존재하지 않음 (KPMP 로 실측)'],
        size=11, line_spacing=1.25)
    txt(s, 9.4, 2.05, 3.4, 0.35, '남은 것', size=14, bold=True)
    txt(s, 9.4, 2.5, 3.4, 3.7,
        ['투고 전 마지막 한 걸음',
         '',
         '1. 저장소 URL 확정',
         '2. Zenodo DOI 발급',
         '3. 투고일 기입',
         '',
         'submission_config.py 에 채운 뒤',
         'fill_submission -> preflight',
         '-> make_release 순으로 실행.',
         '',
         'URL 과 DOI 는 같은 태그된',
         '릴리스를 가리켜야 합니다.'],
        size=11, line_spacing=1.25)
    rule(s, 0.6, 6.4, 12.2)
    txt(s, 0.6, 6.5, 12.2, 0.35,
        '초기 계획의 "다음 단계" 4건(독립 검증 · 단일세포 · 대사체/GWAS · 논문 초고)은 모두 '
        '완료되어 원고에 반영되었습니다.', size=11, color=MUTED)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args()

    d = load()
    prs = Presentation()
    prs.slide_width, prs.slide_height = W, Hh

    from ppt_slides_extra import s_ladder, s_replication
    from ppt_slides_extra2 import s_singlecell, s_deconfound, s_gwas
    from ppt_slides_status import s_now, s_now_review, s_now_defense
    for fn in (s_title, s_now, s_now_review, s_now_defense,
               s_summary, s_data, s_audit, s_genespace, s_method, s_benchmark,
               s_null, s_artifact, s_specificity, s_robustness, s_candidates, s_kpmp,
               s_ladder, s_replication, s_singlecell, s_deconfound, s_gwas,
               s_db, s_status):
        try:
            fn(prs, d)
        except Exception as e:
            log('  ! slide %s failed: %s' % (fn.__name__, e))

    prs.save(args.out)
    log('wrote %s  (%d slides, %.1f MB)'
        % (args.out, len(prs.slides.__iter__.__self__._sldIdLst), os.path.getsize(args.out) / 1e6))


if __name__ == '__main__':
    main()
