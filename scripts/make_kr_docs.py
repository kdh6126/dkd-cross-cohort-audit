#!/usr/bin/env python
"""Build the Korean explainer deck and Word document from scripts/kr_content.py.

One content module feeds both formats so they cannot drift. Numbers in kr_content.py are
traceable to files under results/; SOURCES there names them.

    python scripts/make_kr_docs.py            # both
    python scripts/make_kr_docs.py --pptx     # deck only
    python scripts/make_kr_docs.py --docx     # document only
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr_content as C
import kr_detail as X   # noqa: E402

FONT = '맑은 고딕'
INK = (0x1a, 0x1a, 0x1a)
MUTED = (0x5a, 0x5a, 0x5a)
ACCENT = (0xa8, 0x32, 0x2a)
RULE = (0xd8, 0xd8, 0xd8)
BAND = (0xf2, 0xf0, 0xec)
FIGDIR = 'results/figures'

SW, SH = 13.333, 7.5           # slide size, inches
ML, MR = 0.85, 0.85            # margins
CW = SW - ML - MR              # content width


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ===================================================================== PPTX

def build_pptx(path='DKD_설명자료.pptx'):
    from pptx import Presentation
    from pptx.util import Inches, Pt
    from pptx.dml.color import RGBColor
    from pptx.enum.shapes import MSO_SHAPE

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(SW), Inches(SH)
    counter = {'n': 0}

    def blank():
        return prs.slides.add_slide(prs.slide_layouts[6])

    def tb(s, l, t, w, h):
        box = s.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
        box.text_frame.word_wrap = True
        return box.text_frame

    def run(para, text, size, bold=False, colour=INK):
        r = para.add_run()
        r.text = text
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.name = FONT
        r.font.color.rgb = RGBColor(*colour)
        return r

    def rect(s, l, t, w, h, colour):
        sh = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(l), Inches(t), Inches(w), Inches(h))
        sh.fill.solid()
        sh.fill.fore_color.rgb = RGBColor(*colour)
        sh.line.fill.background()
        sh.shadow.inherit = False
        return sh

    def header(s, kicker, title):
        """Kicker line plus title plus a hairline, the same on every content slide."""
        if kicker:
            f = tb(s, ML, 0.42, CW, 0.35)
            run(f.paragraphs[0], kicker, 12, bold=True, colour=ACCENT)
        f = tb(s, ML, 0.72, CW, 0.85)
        run(f.paragraphs[0], title, 26, bold=True)
        rect(s, ML, 1.52, CW, 0.02, RULE)

    def footer(s):
        counter['n'] += 1
        f = tb(s, SW - 1.35, SH - 0.6, 0.9, 0.35)
        p = f.paragraphs[0]
        p.alignment = 2
        run(p, str(counter['n']), 10, colour=MUTED)

    def callout(s, text, top):
        """Muted band for the 'so what' line that closes most slides."""
        h = 0.42 + 0.24 * (len(text) // 78)
        rect(s, ML, top, CW, h, BAND)
        f = tb(s, ML + 0.22, top + 0.08, CW - 0.44, h - 0.16)
        run(f.paragraphs[0], text, 13, bold=True)

    # ---------------------------------------------------------- slide types

    def slide_title():
        s = blank()
        rect(s, 0, 0, SW, 0.28, ACCENT)
        f = tb(s, ML, 2.25, CW, 1.3)
        run(f.paragraphs[0], C.TITLE, 40, bold=True)
        f2 = tb(s, ML, 3.55, CW, 0.6)
        run(f2.paragraphs[0], C.SUBTITLE, 17, colour=MUTED)
        f3 = tb(s, ML, 6.35, CW, 0.4)
        run(f3.paragraphs[0], C.FOOTER, 11, colour=MUTED)

    def slide_divider(num, text):
        s = blank()
        rect(s, 0, 0, 0.22, SH, ACCENT)
        f = tb(s, 1.4, 3.0, CW - 1.0, 0.7)
        run(f.paragraphs[0], num, 15, bold=True, colour=ACCENT)
        f2 = tb(s, 1.4, 3.4, CW - 1.0, 1.0)
        run(f2.paragraphs[0], text, 32, bold=True)
        footer(s)

    def slide_bullets(kicker, title, lead, bullets, note=None):
        s = blank()
        header(s, kicker, title)
        top = 1.75
        if lead:
            f = tb(s, ML, top, CW, 0.5)
            run(f.paragraphs[0], lead, 15, colour=MUTED)
            top += 0.55
        f = tb(s, ML, top, CW, SH - top - 1.4)
        for i, b in enumerate(bullets):
            p = f.paragraphs[0] if i == 0 else f.add_paragraph()
            p.space_after = Pt(13)
            run(p, '· ' + b, 15)
        if note:
            callout(s, note, SH - 1.25)
        footer(s)
        return s

    def slide_pairs(kicker, title, lead, pairs, note=None, label_w=2.6):
        """Label on the left, explanation on the right. Used for the beats of a finding."""
        s = blank()
        header(s, kicker, title)
        top = 1.75
        if lead:
            f = tb(s, ML, top, CW, 0.45)
            run(f.paragraphs[0], lead, 15, colour=MUTED)
            top += 0.5
        avail = SH - top - (1.35 if note else 0.75)
        row_h = min(0.92, avail / max(len(pairs), 1))
        size = 14 if len(pairs) <= 5 else 12.5
        for k, (lab, txt) in enumerate(pairs):
            y = top + k * row_h
            lf = tb(s, ML, y, label_w, row_h)
            run(lf.paragraphs[0], lab, size, bold=True, colour=ACCENT)
            rf = tb(s, ML + label_w, y, CW - label_w, row_h)
            run(rf.paragraphs[0], txt, size)
        if note:
            callout(s, note, SH - 1.25)
        footer(s)
        return s

    def slide_table(kicker, title, headers, rows, widths=None, note=None):
        from pptx.util import Inches as In
        s = blank()
        header(s, kicker, title)
        nr, nc = len(rows) + 1, len(headers)
        # 행 높이를 0.42 로 고정하면 행이 많은 표가 슬라이드를 넘어간다. 실제로 데이터셋
        # 표(15행)가 1인치 넘쳤다. 남은 공간에 맞춰 줄이고, 좁아지면 글자도 줄인다.
        avail = SH - 1.78 - (1.35 if note else 0.62)
        row_h = min(0.42, avail / max(nr, 1))
        hdr_pt = 12 if row_h >= 0.36 else (11 if row_h >= 0.30 else 10)
        body_pt = 11 if row_h >= 0.36 else (10 if row_h >= 0.30 else 9)
        gt = s.shapes.add_table(nr, nc, In(ML), In(1.78), In(CW), In(row_h * nr)).table
        for rw in gt.rows:
            rw.height = In(row_h)
        if widths:
            tot = sum(widths)
            for j, wd in enumerate(widths):
                gt.columns[j].width = int(In(CW) * wd / tot)
        for j, h in enumerate(headers):
            c = gt.cell(0, j)
            c.text = h
            for p in c.text_frame.paragraphs:
                for r in p.runs:
                    r.font.size = Pt(hdr_pt)
                    r.font.bold = True
                    r.font.name = FONT
        for i, row in enumerate(rows, 1):
            for j, v in enumerate(row):
                c = gt.cell(i, j)
                c.text = str(v)
                c.text_frame.word_wrap = True
                for p in c.text_frame.paragraphs:
                    for r in p.runs:
                        r.font.size = Pt(body_pt)
                        r.font.name = FONT
        if note:
            callout(s, note, min(1.78 + row_h * nr + 0.18, SH - 1.25))
        footer(s)
        return s

    def slide_figure(kicker, title, png, caption):
        from PIL import Image
        s = blank()
        header(s, kicker, title)
        fp = os.path.join(FIGDIR, png)
        if os.path.exists(fp):
            box_w, box_h = CW - 1.2, 3.9
            with Image.open(fp) as im:
                ar = im.width / im.height
            w = min(box_w, box_h * ar)
            h = w / ar
            s.shapes.add_picture(fp, Inches(ML + (CW - w) / 2), Inches(1.8 + (box_h - h) / 2),
                                 width=Inches(w), height=Inches(h))
        else:
            log('  ! figure missing: %s' % fp)
        f = tb(s, ML, 5.85, CW, 1.0)
        run(f.paragraphs[0], caption, 12.5, colour=MUTED)
        footer(s)
        return s

    # ---------------------------------------------------------- deck

    slide_title()

    slide_divider('요약', '한 장으로 보는 결론')
    slide_pairs('요약', '무엇을 했고 무엇이 나왔나', C.SUMMARY_LEAD, C.SUMMARY, label_w=2.4)

    # 현재 상황. 설명자료는 "무엇을 알아냈나" 를 다루지만, 받아 보는 쪽은 대개
    # "지금 어디까지 왔나" 를 먼저 묻습니다. 그 답을 앞쪽에 둡니다.
    for head, items in C.status_now():
        slide_bullets('현재 상황', head,
                      C.STATUS_NOW_LEAD if head == '끝난 것' else None, items)
    slide_pairs('현재 상황', '나올 만한 반론과 답', None,
                C.status_defense(), label_w=3.4)

    slide_divider('1', '먼저 알아야 할 용어')
    half = (len(C.GLOSSARY) + 1) // 2
    for idx, chunk in enumerate((C.GLOSSARY[:half], C.GLOSSARY[half:])):
        slide_pairs('용어', '용어 정리 (%d/2)' % (idx + 1),
                    C.GLOSSARY_INTRO if idx == 0 else None, chunk, label_w=3.1)

    slide_divider('2', '왜 이 연구를 했는가')
    for head, lead, bullets, note in C.WHY:
        slide_bullets('왜 했는가', head, lead, bullets, note)

    slide_divider('3', '어떤 데이터를 썼는가')
    slide_bullets('데이터', '사용한 데이터의 범위', None, C.DATA_INTRO)
    slide_pairs('데이터', '저장소와 접근번호 — GSE, ST 는 무엇인가', None,
                X.REPOSITORIES, label_w=3.2)
    slide_table('데이터', 'GEO 전사체 15건 — 공개된 원래 제목과 출처',
                ['접근번호', '용도', '공개된 제목', '측정 장비', '표본', '출처 논문'],
                X.dataset_rows_geo(), widths=[1.1, 1.4, 3.0, 1.9, 0.6, 2.0])
    slide_table('데이터', 'Metabolomics Workbench 대사체 5건',
                ['접근번호', '용도', '공개된 제목', '측정 방식', '', '출처'],
                X.dataset_rows_mw(), widths=[1.1, 1.4, 4.6, 1.4, 0.2, 1.3])
    slide_table('데이터', '전체 수량', ['항목', '수', '설명'],
                [list(r) for r in X.totals()], widths=[2.4, 1.1, 6.5])
    slide_pairs('데이터', '신장 부위 — 사구체와 세뇨관간질은 무엇이 다른가', None,
                X.COMPARTMENT, label_w=3.4)
    slide_pairs('데이터', '대조군이란 무엇인가', None, X.CONTROL_MEANING, label_w=3.2)
    slide_pairs('데이터', 'ERCB — 이 자료가 왜 특별한가', None, X.ERCB, label_w=2.8)
    slide_table('데이터', '데이터셋 지도 — 이름이 비슷해 헷갈리는 것 정리',
                ['역할', '데이터셋', '무엇에 쓰나'],
                [list(r) for r in C.DATASET_MAP], widths=[2.1, 3.4, 4.5],
                note=C.DATASET_NAMING)
    slide_table('데이터', '전사체 코호트 7종, 258명분',
                ['코호트', '신장 부위', 'DKD', '대조군', '용도'],
                [list(r) for r in C.DISCOVERY_COHORTS],
                widths=[2, 2, 1, 1, 3.4], note=C.COHORT_NOTE)
    slide_table('데이터', '같은 환자가 중복되어 제외한 3건',
                ['접근번호', '제외 사유'], [list(r) for r in C.EXCLUDED], widths=[1.6, 8.4])
    slide_bullets('데이터', '왜 이 정리가 중요했는가', None, C.AUDIT_STORY)
    slide_pairs('데이터', '샘플 이름에서 환자 번호를 어떻게 복원했나', None,
                X.SUBJECT_RECOVERY, label_w=3.2)
    slide_pairs('데이터', '데이터베이스 뷰는 어떻게 보나', None, X.DB_VIEWS, label_w=3.0)
    slide_bullets('데이터', '유전자 이름을 통일한 방법', None, C.GENESPACE)
    slide_pairs('데이터', 'Entrez 번호로 통일한 이유와 방법', None, X.ENTREZ_WHY, label_w=3.6)
    slide_pairs('데이터', '9,900개는 정한 숫자인가', None, X.GENESPACE_9900, label_w=3.4)
    slide_pairs('데이터', '무엇으로 측정한 값인가', None, X.MEASUREMENT, label_w=2.8)

    slide_divider('4', '어떤 방법을 비교했는가')
    slide_pairs('방법', '네 방법은 무엇을 하려는 것인가', None, X.METHODS_PURPOSE,
                label_w=3.4)
    slide_figure('방법', '그림 — 한 번 돌릴 때의 절차 (LODO)', 'K1_lodo_flow.png',
                 '코호트 하나를 통째로 빼고 나머지로 유전자를 고른 뒤, 빼둔 코호트에서 '
                 '채점합니다. 빼는 코호트를 바꿔가며 반복합니다.')
    slide_pairs('방법', '발굴 코호트를 어떻게 쓰는가', None, X.LODO_FLOW, label_w=3.6)
    slide_table('방법', '비교한 네 가지 방법',
                ['방법', '무엇인가', '어디서 쓰나', '강점과 약점'],
                [list(r) for r in C.METHODS_DETAIL], widths=[1.4, 3.5, 2.5, 2.6])
    slide_bullets('방법', 'RBS 안에서 쓴 기본 선택자', None, [C.BASE_SELECTORS])
    slide_table('방법', '네 방법의 원 논문', ['방법', '원 정의', '출처'],
                [list(r) for r in X.METHOD_REFS], widths=[1.8, 4.0, 4.2])
    slide_table('방법', '기본 선택자의 원 논문', ['선택자', '무엇을 하는가', '출처'],
                [list(r) for r in X.BASE_SELECTOR_REFS], widths=[1.6, 3.8, 4.6])
    slide_table('방법', '비교한 교란 보정 6가지', ['보정', '무엇인가', '보통 어디에 쓰나',
                                                     '여기서 무엇이 나왔나'],
                [list(r) for r in C.CORRECTIONS_DETAIL], widths=[1.7, 3.0, 2.7, 2.6])
    slide_table('방법', '교란 보정은 실제로 무엇을 하는가',
                ['보정', '어떻게 작동하는가', '원 논문'],
                [list(r) for r in X.CORRECTION_HOWTO], widths=[1.8, 4.6, 3.6])
    slide_pairs('방법', 'AUROC란 무엇인가', None, X.AUROC_EXPLAIN, label_w=3.2)
    slide_figure('방법', '그림 — AUROC를 세는 법', 'K2_auroc.png',
                 '환자와 대조군을 짝지어 환자 점수가 높은 짝을 셉니다. 그 비율이 ROC 곡선 '
                 '아래 넓이와 정확히 같습니다.')
    slide_pairs('방법', '우리는 AUROC를 어떻게 쟀나', None, X.AUROC_HOWWE, label_w=3.6)
    slide_pairs('방법', '새 알고리즘 개발이 의미가 있는가', None, C.METHOD_VS_DESIGN,
                label_w=3.6)

    slide_divider('5', '무엇을 알아냈는가')
    figmap = {
        '발견 1': ('P2_null_control.png',
                  '무작위로 고른 50개 유전자가 조합에 따라 AUROC 0.68~0.88에 도달합니다. '
                  '라벨을 뒤섞으면 0.50입니다.'),
        '발견 2': ('P1_cohort_sensitivity.png',
                  '(a) 부분집합 크기별 승률 (b) 개별 조합의 재현성 차이 — 부호가 바뀝니다 '
                  '(c) 같은 조합을 외부 AUROC로 보면 방향이 안정적입니다.'),
        '발견 3': ('P3_procurement.png',
                  '(a) 신장 데이터셋 9개의 즉시초기 모듈 수준 (b) 같은 모듈을 네 가지 비교로 '
                  '본 것 — 심장에서는 재현되지 않습니다.'),
        '발견 4': ('P5_block_size.png',
                  '덩어리 크기 비율을 바꾸면 허브 기준의 우위가 부호까지 뒤집힙니다. '
                  '강건성이 아니라 크기 효과입니다.'),
        '발견 5': ('P6_corrections.png',
                  '(a) 채취 관련 유전자 비율 (b) 특이성 통과율을 배경 대비 배수로 '
                  '(c) 보정의 AUROC 비용.'),
    }
    extramap = {'발견 1': ('AUROC로 알 수 없다는 것이 표본이 적어서인가',
                           X.AUROC_WHY_UNINFORMATIVE)}
    readmap = {
        '발견 1': C.FIG_NULL, '발견 2': C.FIG_SENS, '발견 3': C.FIG_PROC,
        '발견 4': C.FIG_BLOCK, '발견 5': C.FIG_CORR,
    }
    for kicker, title, lead, pairs, note in C.FINDINGS:
        slide_pairs(kicker, title, lead, pairs, note)
        if kicker in figmap:
            png, cap = figmap[kicker]
            slide_figure(kicker, '그림 — ' + title, png, cap)
            # 캡션 상자는 서너 줄이 한계라 축과 볼 곳을 설명할 수 없다. 읽는 법을
            # 따로 한 장 붙인다.
            if kicker in readmap:
                slide_pairs(kicker, '그림 읽는 법 — ' + title, None,
                            readmap[kicker], label_w=3.2)
            if kicker in extramap:
                t2, block = extramap[kicker]
                slide_pairs(kicker, t2, None, block, label_w=3.6)
        if kicker == '발견 3':
            slide_figure(kicker, '그림 — 채취 방식을 맞추면 신호가 사라진다',
                         'P4_singlecell_control.png',
                         'KPMP 단일세포. DKD vs 비생검 대조군은 p = 6×10⁻⁷⁶, '
                         'DKD vs 비당뇨 신장병 생검은 p = 0.21로 전혀 갈리지 않습니다.')
            slide_pairs(kicker, '그림 읽는 법 — 채취 방식을 맞추면 신호가 사라진다', None,
                        C.FIG_SC, label_w=3.4)

    slide_divider('6', '바이오마커 후보 대장')
    slide_table('산출물', '이미 알려진 DKD 바이오마커는 무엇인가',
                ['연구', '무엇을 보고했나', '출처'],
                [list(r) for r in X.KNOWN_BIOMARKERS], widths=[2.4, 5.0, 2.6])
    slide_pairs('산출물', '근거 축 여덟은 우리가 만든 것인가', None, X.EVIDENCE_ORIGIN,
                label_w=3.4)
    slide_table('산출물', '근거를 어떤 축으로 세었는가',
                ['근거 축', '무엇을 확인하는가', '통과'],
                [list(r) for r in C.EVIDENCE_AXES], widths=[2.2, 6.4, 1.4],
                note=C.EVIDENCE_NOTE)
    slide_pairs('산출물', '등급 번호는 품질 순서가 아니다', None, X.TIER_MEANING,
                label_w=3.8)
    slide_table('산출물', '후보 30개의 등급 분포', ['등급', '개수', '내용'],
                [list(r) for r in C.TIERS], widths=[2.6, 0.9, 6.5])
    slide_pairs('산출물', 'B등급 후보 — MOXD1', C.DELIVERABLE_LEAD, C.MOXD1, label_w=3.4)
    slide_pairs('산출물', 'MOXD1 — 어디서 나왔고 무엇으로 검증했나', None,
                C.MOXD1_PROVENANCE, label_w=3.8)
    slide_table('산출물', 'MOXD1 실측값', ['항목', '값', '뜻'],
                [list(r) for r in X.MOXD1_VALUES], widths=[2.2, 2.4, 5.4])
    slide_pairs('산출물', '후보를 교차 오믹스로 검증할 수 있는가', None,
                X.MOXD1_CROSSCHECK, label_w=3.6)
    slide_pairs('산출물', '떨어진 후보와 그 교훈', None, C.DROPPED, label_w=3.0)
    slide_pairs('산출물', '우리 게이트도 같은 문제가 있는지 확인했다', None,
                C.SELFCHECK, label_w=3.8)

    slide_pairs('산출물', 'GSE175759의 역할 — 두 용도를 섞지 말 것', None,
                C.GSE175759_ROLE, label_w=3.6)
    slide_pairs('스트레스 테스트', '검증을 돌리기 전에 정한 기준', C.TIERING_RULE,
                C.TIERING, label_w=3.4)
    slide_pairs('CKD', 'CKD는 다중 오믹스가 되는가', None, C.KPMP_MODALITY, label_w=3.2)
    for k in (0, 4):
        slide_pairs('CKD', 'DKD만의 문제인가 (%d/2)' % (k // 4 + 1), None,
                    C.CKD_GENERAL[k:k + 4], label_w=3.6)
    for k in (0, 5):
        slide_pairs('채취 교란', '처리를 표준화하면 효과가 사라지는가 (%d/2)' % (k // 5 + 1),
                    None, C.STORAGE_TEST[k:k + 5], label_w=4.0)
    for k in (0, 5, 8):
        chunk = C.STRESS[k:k + (5 if k == 0 else 3)]
        slide_pairs('스트레스 테스트', '후보가 살아남는지 시험한다 (%d/3)' % (k // 4 + 1),
                    C.STRESS_LEAD if k == 0 else None, chunk, label_w=3.6)
    for k in (0, 5):
        slide_pairs('일반화', '권고가 다른 합병증에서도 성립하는가 (%d/2)' % (k // 5 + 1),
                    None, C.NEUROPATHY[k:k + 5], label_w=3.4)
    slide_pairs('보조자료', '표본이 작아 보조로 남긴 것', None, C.SUPPLEMENTARY, label_w=3.0)
    s_v = blank()
    header(s_v, '스트레스 테스트', '결과 — 셋 중 둘 통과')
    f_v = tb(s_v, ML, 2.3, CW, 2.2)
    run(f_v.paragraphs[0], C.STRESS_VERDICT, 16, bold=True)
    footer(s_v)

    slide_divider('7', '교차 오믹스 — 다른 층에서의 점검')
    slide_pairs('교차 오믹스', '먼저 용어를 구분한다', C.MO_LEAD, C.CROSSOMICS_DEF, label_w=4.2)
    slide_pairs('교차 오믹스', '층위별로 실제 무엇을 했는가', None, C.CROSSOMICS_DONE,
                label_w=3.2)
    slide_table('교차 오믹스', '어느 데이터와 어느 데이터를, 무엇을 근거로 이었는가',
                ['층', '데이터 쌍', '이을 수 있는 근거와 한계'],
                [list(r) for r in C.CROSSOMICS_PAIRS], widths=[1.5, 3.3, 5.2])
    slide_pairs('교차 오믹스', '층 이름을 정확히 — 무엇이 진짜 다른 층인가', None,
                X.LAYER_NAMING, label_w=3.4)
    slide_pairs('교차 오믹스', 'CKD 에서는 가능한가', None, X.CKD_MULTIOMICS, label_w=3.6)
    slide_pairs('교차 오믹스', '단백질 층 교차검증 (1) 두 코호트가 무엇인가', None,
                X.PROTEOME_RESULT[:6], label_w=3.8)
    slide_pairs('교차 오믹스', '단백질 층 교차검증 (2) 무엇이 나왔나', None,
                X.PROTEOME_RESULT[6:11], label_w=3.8)
    slide_pairs('교차 오믹스', '단백질 층 교차검증 (3) 검정을 어떻게 방어했나', None,
                X.PROTEOME_RESULT[11:], label_w=3.8)
    slide_figure('교차 오믹스', '그림 — 단백질 층에서의 후보 성적',
                 'K3_proteome.png',
                 '왼쪽은 후보 12개가 단백질에서 어느 쪽으로 얼마나 움직였는가입니다. 채워진 표시가 q<0.05 이고, 동그라미는 SOMAscan 네모는 질량분석입니다. MMP7 은 양쪽에 다 있습니다. 오른쪽은 그 비율을 배경 단백질과 견준 것입니다.')
    slide_pairs('교차 오믹스', '정말 그 유전자의 단백질을 쟀는가', None,
                X.IDENTITY_CAVEAT, label_w=3.4)
    slide_pairs('교차 오믹스', '그래서 교차 검증이 가능한가', None, C.CROSSOMICS_VERDICT,
                label_w=3.0)
    slide_pairs('교차 오믹스', '이 접근의 한계', None, C.CROSSOMICS_LIMIT, label_w=3.0)
    slide_pairs('검증', '자가검증에서 무엇이 나왔나', None, X.SELF_AUDIT, label_w=3.8)
    slide_pairs('교차 오믹스', 'DKD에서 2차 층이 어려웠던 이유', None, C.MO_WHY_HARD,
                label_w=3.4)
    slide_bullets('다중 오믹스', '그 판단은 틀렸습니다', None, C.MO_WRONG,
                  '검색어를 질병명으로 잡은 것이 원인이었습니다. 데이터가 없던 것이 아닙니다.')
    slide_pairs('다중 오믹스', 'ST003255 — 찾아낸 자료', None, C.MO_DATASET, label_w=2.4)
    slide_pairs('다중 오믹스', '실제로 분석한 결과', None, C.MO_RESULT, label_w=3.4)
    for lab, ttl in (('MO_LESSON', '제대로 말하면 무엇인가'),):
        slide_pairs('다중 오믹스', ttl, None, C.MO_LESSON, label_w=3.6)
    slide_pairs('다중 오믹스', '전사체 후보를 대사체에서 확인할 수 있는가', None,
                C.MO_COLLAGEN, label_w=3.6)
    slide_pairs('다중 오믹스', '그래서 새 바이오마커인가 — 아닙니다', None,
                C.COLLAGEN_CAVEAT, label_w=3.4)
    slide_table('다중 오믹스', '질병별로 무엇이 가능한가',
                ['대상', '판정', '근거'], [list(r) for r in C.MO_OPTIONS],
                widths=[2, 1.7, 6.3])
    slide_pairs('다중 오믹스', '다음 단계 권고', None, C.MO_RECOMMEND, label_w=3.2)

    slide_divider('8', '우리의 기술적 기여')
    s_ = blank()
    header(s_, '기여', '새 알고리즘은 없다. 그럼 무엇이 기여인가')
    f_ = tb(s_, ML, 1.8, CW, 4.6)
    run(f_.paragraphs[0], C.CONTRIB_LEAD, 15, colour=MUTED)
    footer(s_)
    for lab, what, honest in C.CONTRIB:
        s2 = blank()
        header(s2, '기여', lab)
        f1 = tb(s2, ML, 1.8, CW, 2.6)
        run(f1.paragraphs[0], what, 15)
        rect(s2, ML, 4.7, CW, 0.02, RULE)
        f2 = tb(s2, ML, 4.9, CW, 1.4)
        run(f2.paragraphs[0], '정직하게 — ' + honest, 13, colour=MUTED)
        footer(s2)
    s3 = blank()
    header(s3, '기여', '무엇이 기여인지 한 줄로')
    f3 = tb(s3, ML, 2.4, CW, 2.4)
    run(f3.paragraphs[0], C.CONTRIB_VERDICT, 17, bold=True)
    footer(s3)

    slide_divider('9', '정리')
    slide_pairs('정리', '기억할 다섯 가지', None, C.TAKEAWAYS, label_w=3.6)
    slide_table('정리', '현재 상태', ['항목', '상태'], [list(r) for r in C.STATUS],
                widths=[1.7, 8.3])
    slide_bullets('정리', '남은 일', None, C.NEXT)

    slide_divider('부록', '참고문헌')
    for gname, items in C.REFERENCES:
        slide_bullets('참고문헌', gname, None, items)
    slide_bullets('참고문헌', '데이터 출처', None, [
        'GEO — https://www.ncbi.nlm.nih.gov/geo/',
        'KPMP Atlas — https://atlas.kpmp.org/ (인용 문구 준수 필요)',
        'Metabolomics Workbench ST003255 — https://www.metabolomicsworkbench.org/',
        'NHGRI-EBI GWAS Catalog — https://www.ebi.ac.uk/gwas/'])

    prs.save(path)
    log('  wrote %s (%d slides)' % (path, len(prs.slides._sldIdLst)))
    return path


# ===================================================================== DOCX

def build_docx(path='DKD_설명자료.docx'):
    from docx import Document
    from docx.shared import Pt, RGBColor, Inches
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml.ns import qn

    doc = Document()
    st = doc.styles['Normal']
    st.font.name = FONT
    st.font.size = Pt(10.5)
    st.element.rPr.rFonts.set(qn('w:eastAsia'), FONT)

    def H(text, level):
        h = doc.add_heading(text, level=level)
        for r in h.runs:
            r.font.name = FONT
            r._element.rPr.rFonts.set(qn('w:eastAsia'), FONT)
            r.font.color.rgb = RGBColor(*INK)
        return h

    def P(text, bullet=False, muted=False, bold=False):
        p = doc.add_paragraph(style='List Bullet' if bullet else None)
        r = p.add_run(text)
        r.font.name = FONT
        r._element.rPr.rFonts.set(qn('w:eastAsia'), FONT)
        r.font.bold = bold
        if muted:
            r.font.size = Pt(9)
            r.font.color.rgb = RGBColor(*MUTED)
        return p

    def PAIR(label, text):
        p = doc.add_paragraph()
        a = p.add_run(label + ' — ')
        a.font.bold = True
        a.font.name = FONT
        a._element.rPr.rFonts.set(qn('w:eastAsia'), FONT)
        a.font.color.rgb = RGBColor(*ACCENT)
        b = p.add_run(text)
        b.font.name = FONT
        b._element.rPr.rFonts.set(qn('w:eastAsia'), FONT)
        return p

    def NOTE(text):
        p = doc.add_paragraph()
        r = p.add_run('▶ ' + text)
        r.font.bold = True
        r.font.name = FONT
        r._element.rPr.rFonts.set(qn('w:eastAsia'), FONT)
        return p

    def TBL(headers, rows):
        t = doc.add_table(rows=1, cols=len(headers))
        t.style = 'Light Grid Accent 1'
        for j, h in enumerate(headers):
            c = t.rows[0].cells[j]
            c.text = h
            for p in c.paragraphs:
                for r in p.runs:
                    r.font.bold = True
                    r.font.name = FONT
                    r._element.rPr.rFonts.set(qn('w:eastAsia'), FONT)
                    r.font.size = Pt(9.5)
        for row in rows:
            cells = t.add_row().cells
            for j, v in enumerate(row):
                cells[j].text = str(v)
                for p in cells[j].paragraphs:
                    for r in p.runs:
                        r.font.name = FONT
                        r._element.rPr.rFonts.set(qn('w:eastAsia'), FONT)
                        r.font.size = Pt(9.5)
        doc.add_paragraph()
        return t

    def FIG(png, caption):
        fp = os.path.join(FIGDIR, png)
        if not os.path.exists(fp):
            log('  ! figure missing: %s' % fp)
            return
        doc.add_picture(fp, width=Inches(6.1))
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
        P(caption, muted=True)

    H(C.TITLE, 0)
    P(C.SUBTITLE, bold=True)
    P(C.FOOTER, muted=True)

    H('요약 — 한 장으로 보는 결론', 1)
    P(C.SUMMARY_LEAD)
    for lab, txt in C.SUMMARY:
        PAIR(lab, txt)

    H('현재 상황', 1)
    P(C.STATUS_NOW_LEAD, muted=True)
    for head, items in C.status_now():
        H(head, 2)
        for it in items:
            P(it, bullet=True)
    H('나올 만한 반론과 답', 2)
    for q, a in C.status_defense():
        PAIR(q, a)

    H('1. 먼저 알아야 할 용어', 1)
    P(C.GLOSSARY_INTRO)
    for term, expl in C.GLOSSARY:
        PAIR(term, expl)

    H('2. 왜 이 연구를 했는가', 1)
    for head, lead, bullets, note in C.WHY:
        H(head, 2)
        if lead:
            P(lead)
        for b in bullets:
            P(b, bullet=True)
        if note:
            NOTE(note)

    H('3. 어떤 데이터를 썼는가', 1)
    for b in C.DATA_INTRO:
        P(b)
    H('저장소와 접근번호 — GSE, ST 는 무엇인가', 2)
    for lab, txt in X.REPOSITORIES:
        PAIR(lab, txt)
    H('GEO 전사체 15건 — 공개된 원래 제목과 출처', 2)
    TBL(['접근번호', '용도', '공개된 제목', '측정 장비', '표본', '출처 논문'],
        X.dataset_rows_geo())
    H('Metabolomics Workbench 대사체 5건', 2)
    TBL(['접근번호', '용도', '공개된 제목', '측정 방식', '', '출처'], X.dataset_rows_mw())
    H('전체 수량', 2)
    TBL(['항목', '수', '설명'], [list(r) for r in X.totals()])
    H('신장 부위 — 사구체와 세뇨관간질은 무엇이 다른가', 2)
    for lab, txt in X.COMPARTMENT:
        PAIR(lab, txt)
    H('대조군이란 무엇인가', 2)
    for lab, txt in X.CONTROL_MEANING:
        PAIR(lab, txt)
    H('ERCB — 이 자료가 왜 특별한가', 2)
    for lab, txt in X.ERCB:
        PAIR(lab, txt)
    H('데이터셋 지도 — 이름이 비슷해 헷갈리는 것 정리', 2)
    P(C.DATASET_MAP_LEAD)
    TBL(['역할', '데이터셋', '무엇에 쓰나'], [list(r) for r in C.DATASET_MAP])
    P(C.DATASET_NAMING, muted=True)
    H('전사체 코호트 7종, 258명분', 2)
    TBL(['코호트', '신장 부위', 'DKD', '대조군', '용도'], [list(r) for r in C.DISCOVERY_COHORTS])
    P(C.COHORT_NOTE, muted=True)
    H('같은 환자가 중복되어 제외한 3건', 2)
    TBL(['접근번호', '제외 사유'], [list(r) for r in C.EXCLUDED])
    H('왜 이 정리가 중요했는가', 2)
    for b in C.AUDIT_STORY:
        P(b, bullet=True)
    H('샘플 이름에서 환자 번호를 어떻게 복원했나', 2)
    for lab, txt in X.SUBJECT_RECOVERY:
        PAIR(lab, txt)
    H('데이터베이스 뷰는 어떻게 보나', 2)
    for lab, txt in X.DB_VIEWS:
        PAIR(lab, txt)
    H('유전자 이름을 통일한 방법', 2)
    for b in C.GENESPACE:
        P(b, bullet=True)
    H('Entrez 번호로 통일한 이유와 방법', 2)
    for lab, txt in X.ENTREZ_WHY:
        PAIR(lab, txt)
    H('9,900개는 정한 숫자인가', 2)
    for lab, txt in X.GENESPACE_9900:
        PAIR(lab, txt)
    H('무엇으로 측정한 값인가', 2)
    for lab, txt in X.MEASUREMENT:
        PAIR(lab, txt)

    H('4. 어떤 방법을 비교했는가', 1)
    P(C.METHODS_LEAD)
    H('네 방법은 무엇을 하려는 것인가', 2)
    for lab, txt in X.METHODS_PURPOSE:
        PAIR(lab, txt)
    FIG('K1_lodo_flow.png', '한 번 돌릴 때의 절차 (LODO)')
    H('발굴 코호트를 어떻게 쓰는가', 2)
    for lab, txt in X.LODO_FLOW:
        PAIR(lab, txt)
    TBL(['방법', '무엇인가', '어디서 쓰나', '강점과 약점'],
        [list(r) for r in C.METHODS_DETAIL])
    H('네 방법의 원 논문', 2)
    TBL(['방법', '원 정의', '출처'], [list(r) for r in X.METHOD_REFS])
    H('기본 선택자의 원 논문', 2)
    TBL(['선택자', '무엇을 하는가', '출처'], [list(r) for r in X.BASE_SELECTOR_REFS])
    H('RBS 안에서 쓴 기본 선택자', 2)
    P(C.BASE_SELECTORS)
    H('비교한 교란 보정 6가지', 2)
    TBL(['보정', '무엇인가', '보통 어디에 쓰나', '여기서 무엇이 나왔나'],
        [list(r) for r in C.CORRECTIONS_DETAIL])
    H('교란 보정은 실제로 무엇을 하는가', 2)
    TBL(['보정', '어떻게 작동하는가', '원 논문'], [list(r) for r in X.CORRECTION_HOWTO])
    H('AUROC란 무엇인가', 2)
    for lab, txt in X.AUROC_EXPLAIN:
        PAIR(lab, txt)
    FIG('K2_auroc.png', 'AUROC를 세는 법')
    H('우리는 AUROC를 어떻게 쟀나', 2)
    for lab, txt in X.AUROC_HOWWE:
        PAIR(lab, txt)
    H('AUROC로 알 수 없다는 것이 표본이 적어서인가', 2)
    for lab, txt in X.AUROC_WHY_UNINFORMATIVE:
        PAIR(lab, txt)
    H('새 알고리즘 개발이 의미가 있는가', 2)
    for lab, txt in C.METHOD_VS_DESIGN:
        PAIR(lab, txt)

    H('5. 무엇을 알아냈는가', 1)
    figmap = {
        '발견 1': [('P2_null_control.png',
                   '무작위로 고른 50개 유전자가 조합에 따라 AUROC 0.68~0.88에 도달합니다.')],
        '발견 2': [('P1_cohort_sensitivity.png',
                   '재현성 차이는 부호가 바뀌고, 외부 AUROC는 방향이 안정적입니다.')],
        '발견 3': [('P3_procurement.png',
                   '즉시초기 모듈이 채취 방식을 갈라냅니다. 심장에서는 재현되지 않습니다.'),
                  ('P4_singlecell_control.png',
                   'KPMP 단일세포. 채취 방식을 맞추면 신호가 사라집니다(p = 0.21).')],
        '발견 4': [('P5_block_size.png',
                   '덩어리 크기 비율에 따라 허브 기준의 우위가 뒤집힙니다.')],
        '발견 5': [('P6_corrections.png',
                   '보정별로 남는 교란, 특이성, AUROC 비용을 함께 봅니다.')],
    }
    for kicker, title, lead, pairs, note in C.FINDINGS:
        H('%s. %s' % (kicker, title), 2)
        if lead:
            P(lead)
        for lab, txt in pairs:
            PAIR(lab, txt)
        if note:
            NOTE(note)
        for png, cap in figmap.get(kicker, []):
            FIG(png, cap)
        # 그림을 넣었으면 읽는 법도 같이 넣는다. 슬라이드와 같은 내용이다.
        for rk, block in (('발견 1', C.FIG_NULL), ('발견 2', C.FIG_SENS),
                          ('발견 3', C.FIG_PROC), ('발견 4', C.FIG_BLOCK),
                          ('발견 5', C.FIG_CORR)):
            if kicker == rk:
                H('그림 읽는 법 — ' + title, 3)
                for lab, txt in block:
                    PAIR(lab, txt)
        if kicker == '발견 3':
            H('그림 읽는 법 — 채취 방식을 맞추면 신호가 사라진다', 3)
            for lab, txt in C.FIG_SC:
                PAIR(lab, txt)

    H('6. 바이오마커 후보 대장', 1)
    P(C.DELIVERABLE_LEAD)
    H('이미 알려진 DKD 바이오마커는 무엇인가', 2)
    TBL(['연구', '무엇을 보고했나', '출처'], [list(r) for r in X.KNOWN_BIOMARKERS])
    H('근거 축 여덟은 우리가 만든 것인가', 2)
    for lab, txt in X.EVIDENCE_ORIGIN:
        PAIR(lab, txt)
    H('근거를 어떤 축으로 세었는가', 2)
    TBL(['근거 축', '무엇을 확인하는가', '통과'], [list(r) for r in C.EVIDENCE_AXES])
    P(C.EVIDENCE_NOTE, muted=True)
    H('등급 번호는 품질 순서가 아니다', 2)
    for lab, txt in X.TIER_MEANING:
        PAIR(lab, txt)
    H('후보 30개의 등급 분포', 2)
    TBL(['등급', '개수', '내용'], [list(r) for r in C.TIERS])
    H('B등급 후보 — MOXD1', 2)
    for lab, txt in C.MOXD1:
        PAIR(lab, txt)
    H('MOXD1 — 어디서 나왔고 무엇으로 검증했나', 2)
    for lab, txt in C.MOXD1_PROVENANCE:
        PAIR(lab, txt)
    H('MOXD1 실측값', 2)
    TBL(['항목', '값', '뜻'], [list(r) for r in X.MOXD1_VALUES])
    H('후보를 교차 오믹스로 검증할 수 있는가', 2)
    for lab, txt in X.MOXD1_CROSSCHECK:
        PAIR(lab, txt)
    H('떨어진 후보와 그 교훈', 2)
    for lab, txt in C.DROPPED:
        PAIR(lab, txt)
    H('우리 게이트도 같은 문제가 있는지 확인했다', 2)
    for lab, txt in C.SELFCHECK:
        PAIR(lab, txt)

    H('GSE175759의 역할 — 두 용도를 섞지 말 것', 2)
    for lab, txt in C.GSE175759_ROLE:
        PAIR(lab, txt)
    H('검증을 돌리기 전에 정한 기준', 2)
    P(C.TIERING_RULE)
    for lab, txt in C.TIERING:
        PAIR(lab, txt)
    H('CKD는 다중 오믹스가 되는가', 2)
    for lab, txt in C.KPMP_MODALITY:
        PAIR(lab, txt)
    H('DKD만의 문제인가 — 신장병 8종에서 확인', 2)
    for lab, txt in C.CKD_GENERAL:
        PAIR(lab, txt)
    H('처리를 표준화하면 채취 효과가 사라지는가', 2)
    for lab, txt in C.STORAGE_TEST:
        PAIR(lab, txt)
    H('스트레스 테스트 — 후보가 살아남는지 시험한다', 2)
    P(C.STRESS_LEAD)
    for lab, txt in C.STRESS:
        PAIR(lab, txt)
    P(C.STRESS_VERDICT, bold=True)
    H('권고가 다른 합병증에서도 성립하는가 — 당뇨병성 신경병증', 2)
    for lab, txt in C.NEUROPATHY:
        PAIR(lab, txt)
    H('표본이 작아 보조로 남긴 것', 2)
    for lab, txt in C.SUPPLEMENTARY:
        PAIR(lab, txt)

    H('7. 교차 오믹스 — 다른 층에서의 점검', 1)
    P(C.MO_LEAD)
    H('먼저 용어를 구분한다', 2)
    for lab, txt in C.CROSSOMICS_DEF:
        PAIR(lab, txt)
    H('층위별로 실제 무엇을 했는가', 2)
    for lab, txt in C.CROSSOMICS_DONE:
        PAIR(lab, txt)
    H('이 접근의 한계', 2)
    H('어느 데이터와 어느 데이터를, 무엇을 근거로 이었는가', 2)
    TBL(['층', '데이터 쌍', '이을 수 있는 근거와 한계'],
        [list(r) for r in C.CROSSOMICS_PAIRS])
    H('층 이름을 정확히 — 무엇이 진짜 다른 층인가', 2)
    for lab, txt in X.LAYER_NAMING:
        PAIR(lab, txt)
    H('CKD 에서는 가능한가', 2)
    for lab, txt in X.CKD_MULTIOMICS:
        PAIR(lab, txt)
    H('단백질 층에서 실제로 재현됐다', 2)
    for lab, txt in X.PROTEOME_RESULT:
        PAIR(lab, txt)
    FIG('K3_proteome.png', '단백질 층에서의 후보 성적')
    H('정말 그 유전자의 단백질을 쟀는가', 2)
    for lab, txt in X.IDENTITY_CAVEAT:
        PAIR(lab, txt)
    H('자가검증에서 무엇이 나왔나', 2)
    for lab, txt in X.SELF_AUDIT:
        PAIR(lab, txt)
    H('그래서 교차 검증이 가능한가', 2)
    for lab, txt in C.CROSSOMICS_VERDICT:
        PAIR(lab, txt)
    H('이 접근의 한계', 2)
    for lab, txt in C.CROSSOMICS_LIMIT:
        PAIR(lab, txt)
    H('DKD에서 2차 층이 어려웠던 이유', 2)
    for lab, txt in C.MO_WHY_HARD:
        PAIR(lab, txt)
    H('그 판단은 틀렸습니다', 2)
    for b in C.MO_WRONG:
        P(b, bullet=True)
    H('ST003255 — 찾아낸 자료', 2)
    for lab, txt in C.MO_DATASET:
        PAIR(lab, txt)
    H('실제로 분석한 결과', 2)
    for lab, txt in C.MO_RESULT:
        PAIR(lab, txt)
    H('제대로 말하면 무엇인가', 2)
    for lab, txt in C.MO_LESSON:
        PAIR(lab, txt)
    H('전사체 후보를 대사체에서 확인할 수 있는가', 2)
    for lab, txt in C.MO_COLLAGEN:
        PAIR(lab, txt)
    H('그래서 새 바이오마커인가 — 아닙니다', 2)
    for lab, txt in C.COLLAGEN_CAVEAT:
        PAIR(lab, txt)
    H('질병별로 무엇이 가능한가', 2)
    TBL(['대상', '판정', '근거'], [list(r) for r in C.MO_OPTIONS])
    H('다음 단계 권고', 2)
    for lab, txt in C.MO_RECOMMEND:
        PAIR(lab, txt)

    H('8. 우리의 기술적 기여', 1)
    P(C.CONTRIB_LEAD)
    for lab, what, honest in C.CONTRIB:
        H(lab, 2)
        P(what)
        P('정직하게 — ' + honest, muted=True)
    P(C.CONTRIB_VERDICT, bold=True)

    H('9. 정리', 1)
    for lab, txt in C.TAKEAWAYS:
        PAIR(lab, txt)
    H('현재 상태', 2)
    TBL(['항목', '상태'], [list(r) for r in C.STATUS])
    H('남은 일', 2)
    for b in C.NEXT:
        P(b, bullet=True)

    H('부록 A — 참고문헌', 1)
    for gname, items in C.REFERENCES:
        H(gname, 2)
        for it in items:
            P(it, bullet=True)

    H('부록 A-2 — 데이터 출처', 1)
    TBL(['자원', 'URL'], [
        ['GEO', 'https://www.ncbi.nlm.nih.gov/geo/'],
        ['KPMP Atlas', 'https://atlas.kpmp.org/'],
        ['Metabolomics Workbench (ST003255)', 'https://www.metabolomicsworkbench.org/'],
        ['NHGRI-EBI GWAS Catalog', 'https://www.ebi.ac.uk/gwas/']])

    H('부록 B — 이 자료의 수치가 나온 파일', 1)
    P('모든 수치는 아래 파일에서 직접 읽었습니다. 파이프라인을 다시 돌리면 같은 값이 나옵니다.')
    TBL(['항목', '파일'], [[k, v] for k, v in sorted(C.SOURCES.items())])

    doc.save(path)
    log('  wrote %s' % path)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pptx', action='store_true')
    ap.add_argument('--docx', action='store_true')
    a = ap.parse_args()
    both = not (a.pptx or a.docx)
    if a.pptx or both:
        build_pptx()
    if a.docx or both:
        build_docx()


if __name__ == '__main__':
    main()
