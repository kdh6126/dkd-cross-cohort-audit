#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Briefings in Bioinformatics 판을 만든다.

BMC 판(submission/dkd-manuscript.tex)이 감사를 거친 원본이다. BiB 는 문제 해결 프로토콜을
2,000-5,000 단어로 받으므로 본문은 따로 줄여 썼고(submission_bib/bib-manuscript.tex), 줄이며
본문에서 빠진 방법과 결과는 BMC 판의 문장을 그대로 떼어 Supplementary Notes 로 묶는다.
손으로 옮겨 적으면 두 판이 어긋나므로, 보충 자료는 여기서 매번 BMC 판에서 다시 만든다.

    1. 그림을 BiB 번호로 복사한다
    2. 참고문헌 양식: Vancouver 에서 저자 넷 이상이면 셋 뒤 et al. (BiB 규정)
    3. Supplementary Notes S1-S10 을 BMC 판에서 잘라 조립한다
    4. 둘 다 컴파일한다
    5. 검사: BiB 본문의 수치가 전부 BMC 판에 있는가, 단어 수, 미해결 참조, overfull
    6. Supplementary File 1-7 을 BiB 표기로 다시 만들고 BMC 잔재가 없는지 읽어 본다
"""
import os
import re
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
DST = 'submission_bib'
BMC = 'submission/dkd-manuscript.tex'
B = chr(92)
E = re.escape(B)
NL = chr(10)

# 모두 BiB 본문 폭(526pt)으로 그린 판이다. BMC 폭으로 그린 것을 넣으면 1.41 배로
# 늘어나 그림 글자가 본문보다 커진다.
FIGS = [('Fig1', 'P1_cohort_sensitivity_bib'), ('Fig2', 'P3_procurement_bib'),
        ('Fig3', 'P6_corrections_bib'), ('Fig4', 'P9_reordering_bib'),
        ('Fig5', 'P10_cross_tissue_summary')]
# 보충자료 그림, 등장 순서. 원본 이름은 BMC 판에서 잘라 올 때 되살아나므로 여기서 다시 바꾼다.
SUPP_FIGS = [('FigS1', 'P2_null_control'), ('FigS2', 'P4_singlecell_control'),
             ('FigS3', 'P5_block_size'), ('FigS4', 'P7_candidates'),
             ('FigS5', 'P10_cross_tissue')]
FIGDIR = 'figures'

# BiB 본문에 들어간 그림과 표. 보충 자료에서는 떼어 내고 본문 번호로 가리킨다.
MAIN_FIG = {'fig:sensitivity': 1, 'fig:procurement': 2, 'fig:corrections': 3,
            'fig:reordering': 4, 'fig:xtissue': 5}
MAIN_TAB = {'tab:overlap': 1, 'tab:cohorts': 2, 'tab:comparator': 3, 'tab:checklist': 4,
            'tab:selectors': 5, 'tab:corrections': 6, 'tab:sweep': 7, 'tab:ckdgen': 8,
            'tab:xtissue': 9}
# 보충 자료에 제목 줄이 들어오지 않는 절. 가리키는 곳을 말로 적는다.
SEC_TEXT = {
    'sec:cohorts': 'the main-text section ' + B + 'textit{Cohort assembly and audit}',
    'sec:methods': 'the main-text ' + B + 'textit{Materials and methods}',
    'sec:results': 'the main-text ' + B + 'textit{Results}',
    'sec:corrections': 'the main-text section ' + B + 'textit{Confounder score, corrections and control definition}',
    'sec:sensitivity': 'the main-text section ' + B + 'textit{Reproducibility verdicts depend on cohort composition}',
    'sec:procurement': 'Supplementary Note~S3',
    'sec:xtissue': 'the main-text section ' + B + 'textit{Both mechanisms outside kidney}',
    'sec:discussion': 'the main-text ' + B + 'textit{Discussion}',
    'sec:limitations': 'Supplementary Note~S8',
    'sec:conclusions': 'the main-text ' + B + 'textit{Conclusion}',
    'sec:background': 'the main-text ' + B + 'textit{Introduction}',
}
# BMC 판 문장 가운데 보충 자료로 옮기면 사실이 아니게 되는 것
REWORD = [
    ('and we report it in the main text rather' + NL + 'than as supplementary material because it bounds the claim.',
     'and the main text summarises it because it bounds the claim.'),
]

NOTES = [
    ('S1', 'Harmonisation, selection strategies and evaluation',
     [(E + r'subsection\{Harmonisation\}', E + r'subsection\{Null controls and uncertainty\}')]),
    ('S2', 'Null controls, uncertainty and the random-signature null',
     [(E + r'subsection\{Null controls and uncertainty\}', E + r'subsection\{Confounder score'),
      (E + r'subsection\{Absolute AUROC is uninformative', E + r'subsection\{Reproducibility verdicts')]),
    ('S3', 'Procurement evidence: per-dataset values, single-nucleus control and heart',
     [(E + r'textbf\{A note on what is and is not measured\.\}', E + r'subsection\{Why one method escaped')]),
    ('S4', 'Why the connectivity criterion escaped the confounder',
     [(r'The simulation uses a generative model', E + r'subsection\{Predictions fixed before the liver'),
      (E + r'subsection\{Why one method escaped', E + r'subsection\{Which corrections work\}')]),
    ('S5', 'Correction benchmark details and the seven inner selectors',
     [(E + r'subsection\{Which corrections work\}', E + r'subsection\{The control-group problem')]),
    ('S6', 'Control-definition analysis by diagnosis and external checks',
     [(E + r'subsection\{The control-group problem', E + r'subsection\{Both mechanisms outside kidney')]),
    ('S7', 'Candidates and protein-level corroboration',
     [(E + r'subsection\{What the pipeline returned', E + r'section\{Discussion\}')]),
    ('S8', 'Cohort overlap and further limitations',
     [(r'Residual overlap that could not be removed is disclosed', r'Five cohorts entered the discovery rotation'),
      (r'The confounder-matched evidence rests on ERCB and KPMP', E + r'section\{Conclusions\}')]),
    ('S9', 'Implementation and regeneration',
     [(E + r'subsection\{Implementation\}', E + r'section\{Results\}'),
      (E + r'textbf\{Regeneration\.\}', E + r'bmhead\{Competing interests\}')]),
    ('S10', 'Cross-tissue test with hash-fixed predictions: full methods and liver control substitution',
     [(E + r'subsection\{Predictions fixed before the liver and colon data were read\}', E + r'subsection\{Implementation\}'),
      # 끝 표지를 begin{figure} 로 두면 stage_figures 가 그림을 옮긴 뒤에는 후보 절까지
      # 딸려 들어온다. 실제로 S7 의 후보 절이 S10 에 중복되고 sec:candidates 라벨이
      # 두 번 정의됐다. 다음 절 제목으로 끊는다.
      (E + r'textbf\{Control substitution in liver \(post hoc\)\.\}',
       E + r'subsection\{What the pipeline returned')]),
]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def figures():
    # BiB 폭 전용 판을 먼저 그린다.
    for cmd in ([sys.executable, 'scripts/paper_figures.py', '--bib'],
                [sys.executable, 'scripts/reordering_figure.py', '--bib']):
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace')
        if r.returncode != 0:
            raise SystemExit('BiB 폭 그림 생성 실패: %s%s%s' % (' '.join(cmd), NL, r.stderr[-800:]))
    fdir = os.path.join(DST, FIGDIR)
    os.makedirs(fdir, exist_ok=True)
    for new, src in FIGS + SUPP_FIGS:
        shutil.copy2(os.path.join('results/figures', src + '.pdf'), os.path.join(fdir, new + '.pdf'))
    # 예전에는 루트에 두었다. 남아 있으면 preflight 가 낯선 파일로 잡는다.
    for f in os.listdir(DST):
        if re.match(r'Fig\d+\.pdf$', f):
            os.remove(os.path.join(DST, f))


def bst():
    s = open('submission/sn-vancouver-num.bst', encoding='utf-8').read()
    a = 'nameptr #6' + NL + '          #1 + ='
    b = 'numnames #6' + NL + '          > and'
    if s.count(a) != 1 or s.count(b) != 1:
        raise SystemExit('Vancouver 양식의 저자 수 규칙을 찾지 못했습니다')
    s = s.replace(a, 'nameptr #3' + NL + '          #1 + =').replace(b, 'numnames #3' + NL + '          > and')
    open(os.path.join(DST, 'bib-vancouver.bst'), 'w', encoding='utf-8').write(s)
    shutil.copy2('submission/dkd-references.bib', os.path.join(DST, 'dkd-references.bib'))


def cut(src, start, end):
    m = re.search(start, src)
    if not m:
        raise SystemExit('시작 표지를 찾지 못했습니다: ' + start)
    n = re.search(end, src[m.start() + 1:])
    if not n:
        raise SystemExit('끝 표지를 찾지 못했습니다: ' + end)
    return src[m.start():m.start() + 1 + n.start()]


def drop_main_floats(t):
    def keep(m):
        lab = re.search(E + r'label\{([^}]+)\}', m.group(0))
        if lab and (lab.group(1) in MAIN_FIG or lab.group(1) in MAIN_TAB):
            return ''
        return m.group(0)
    return re.sub('(?s)' + E + r'begin\{(figure|table)\}.*?' + E + r'end\{\1\}', keep, t)


def supplement():
    src = open(BMC, encoding='utf-8').read()
    body = []
    for num, title, parts in NOTES:
        body.append(B + 'suppnote{' + title + '}' + B + 'label{note:' + num + '}')
        for a, b in parts:
            body.append(drop_main_floats(cut(src, a, b)))
    # 본문 Fig. 6 은 요약판이다. 접근번호 라벨이 붙은 코호트별 판은 BMC 판의 그림 환경을
    # 그대로 가져와 S10 끝에 둔다.
    full = re.search('(?s)' + E + r'begin\{figure\}(?:(?!' + E + r'end\{figure\}).)*?'
                     + E + r'label\{fig:xtissue\}.*?' + E + r'end\{figure\}', src).group(0)
    full = full.replace(B + 'label{fig:xtissue}', B + 'label{fig:xtissue-full}')
    full = re.sub(r'% source: [A-Za-z0-9_]+' + NL, '', full)
    full = re.sub(E + r'includegraphics\[[^]]*\]\{[^}]+\}',
                  lambda m: B + 'includegraphics[width=' + B + 'textwidth]{P10_cross_tissue}', full)
    body.append(full)
    t = NL.join(body)
    t = re.sub(r'(?m)^%%=+%%$', '', t)
    # 그림 파일: BMC 판은 FigN 로 번호가 매겨져 있으므로 원본 이름으로 되돌린다
    t = re.sub(r'% source: ([A-Za-z0-9_]+)' + NL + r'((?:.*' + NL + r')*?)(.*?' + E + r'includegraphics\[[^]]*\]\{)Fig\d+\}',
               lambda m: m.group(2) + m.group(3) + m.group(1) + '}', t)
    for new, src in SUPP_FIGS:
        o = B + 'includegraphics[width=' + B + 'textwidth]{' + src + '}'
        if t.count(o) != 1:
            raise SystemExit('보충자료 그림 %s 가 %d 번 나옵니다' % (src, t.count(o)))
        t = t.replace(o, B + 'includegraphics[width=' + B + 'textwidth]{' + new + '}')
    t = re.sub(r'Additional' + r'\s+' + r'file(~|\s)', lambda m: 'Supplementary File' + m.group(1), t)
    t = t.replace('Additional files', 'Supplementary Files')
    for o, n in REWORD:
        if o not in t:
            raise SystemExit('고칠 문장을 찾지 못했습니다: ' + o[:40])
        t = t.replace(o, n)
    labels = set(re.findall(E + r'label\{([^}]+)\}', t))
    unresolved = []

    def fix_ref(m):
        word, lab = m.group(1), m.group(2)
        if lab in labels:
            return m.group(0)
        if lab in MAIN_FIG:
            return 'Figure~%d of the main text' % MAIN_FIG[lab]
        if lab in MAIN_TAB:
            return 'Table~%d of the main text' % MAIN_TAB[lab]
        if lab in SEC_TEXT:
            return SEC_TEXT[lab]
        unresolved.append(lab)
        return 'the main text'
    t = re.sub(r'Eq\.~\(' + E + r'ref\{([^}]+)\}\)', lambda m: m.group(0) if m.group(1) in labels else 'Eq.~(1) of the main text', t)
    t = re.sub(r'(Sections?~|Table~|Figure~|Fig\.~)' + E + r'ref\{([^}]+)\}', fix_ref, t)
    t = re.sub(r'(?:and|,)~' + E + r'ref\{([^}]+)\}',
               lambda m: m.group(0) if m.group(1) in labels else 'and ' + SEC_TEXT.get(m.group(1), 'the main text'), t)
    left = [x for x in re.findall(E + r'ref\{([^}]+)\}', t) if x not in labels]
    if unresolved or left:
        log('  보충 자료에서 말로 바꾼 참조: %s' % ', '.join(sorted(set(unresolved + left))))

    pre = [
        B + 'documentclass[11pt]{article}',
        B + 'usepackage[a4paper,margin=22mm]{geometry}',
        B + 'usepackage[T1]{fontenc}', B + 'usepackage{lmodern}',
        B + 'usepackage{graphicx,booktabs,array,amsmath,amssymb,url}',
        B + 'usepackage[numbers,square,sort&compress]{natbib}',
        B + 'usepackage[hidelinks]{hyperref}',
        B + 'graphicspath{{' + FIGDIR + '/}}',
        B + 'providecommand{' + B + 'bibcommenthead}{}',
        B + 'providecommand{' + B + 'botrule}{' + B + 'bottomrule}',
        B + 'renewcommand{' + B + 'footnotetext}[1]{' + B + 'par' + B + 'smallskip{' + B + 'footnotesize #1}}',
        B + 'renewcommand' + B + 'thesection{S' + B + 'arabic{section}}',
        B + 'renewcommand' + B + 'thefigure{S' + B + 'arabic{figure}}',
        B + 'renewcommand' + B + 'thetable{S' + B + 'arabic{table}}',
        B + 'renewcommand' + B + 'theequation{S' + B + 'arabic{equation}}',
        B + 'setlength{' + B + 'parskip}{4pt}',
        # DEG_meta 문장처럼 끊을 수 없는 낱말이 줄 폭을 7pt 넘겼다. 표준 처방.
        B + 'emergencystretch=3em',
        B + 'begin{document}',
        B + 'begin{center}{' + B + 'Large Supplementary Data}' + B + B + '[4pt]',
        'Auditing cross-cohort feature-selection benchmarks: cohort composition and '
        'tissue-context-dependent confounding across kidney, liver and colon transcriptomes' + B + B + '[4pt]',
        B + 'textit{Briefings in Bioinformatics}' + B + 'end{center}',
        'Supplementary Notes S1--S10 carry the methods and results condensed out of the main text. '
        'Their wording is taken from the full-length version of the analysis, so section, figure '
        'and table numbers prefixed with S refer to this document, and others to the main text. '
        'Supplementary Files 1--7 are supplied separately.',
        # \section 을 재정의하면 thebibliography 의 \section*{References} 까지 함께 망가져
        # "Supplementary Note S11. *" 라는 가짜 절이 생겼다. 전용 명령을 따로 만든다.
        B + 'newcommand' + B + 'suppnote[1]{' + B + 'refstepcounter{section}' + B + 'par' + B + 'bigskip'
        + B + 'noindent{' + B + 'large' + B + 'bfseries Supplementary Note ' + B + 'thesection. #1}'
        + B + 'par' + B + 'medskip}',
    ]
    doc = NL.join(pre) + NL + t + NL + B + 'bibliographystyle{bib-vancouver}' + NL \
        + B + 'bibliography{dkd-references}' + NL + B + 'end{document}' + NL
    # section 재정의는 begin{document} 뒤에 둬야 서문의 제목이 영향을 받지 않는다
    open(os.path.join(DST, 'bib-supplement.tex'), 'w', encoding='utf-8').write(doc)


def compile_tex(name):
    run = lambda *c: subprocess.run(c, cwd=DST, capture_output=True)
    run('pdflatex', '-interaction=nonstopmode', name + '.tex')
    run('bibtex', name)
    run('pdflatex', '-interaction=nonstopmode', name + '.tex')
    p = run('pdflatex', '-interaction=nonstopmode', name + '.tex')
    logtxt = open(os.path.join(DST, name + '.log'), encoding='utf-8', errors='replace').read()
    errs = [l for l in logtxt.splitlines() if l.startswith('!')]
    undef = re.findall(r"Reference `([^']+)' .*undefined", logtxt)
    cites = re.findall(r"Citation `([^']+)' .*undefined", logtxt)
    over = re.findall(r'Overfull \\hbox \(([\d.]+)pt', logtxt)
    pages = re.search(r'Output written on .*?\((\d+) pages', logtxt)
    # OUP 템플릿은 표제면 장식 때문에 11.38pt overfull 을 원본 샘플에서도 낸다. 우리 내용이
    # 아니므로 그 값은 기대값으로 본다.
    log('  %-18s pages %s, errors %d, undefined refs %d, undefined cites %d, overfull %d%s'
        % (name, pages.group(1) if pages else '?', len(errs), len(set(undef)), len(set(cites)),
           len(over), (' (max %.0fpt)' % max(map(float, over))) if over else ''))
    for l in errs[:5]:
        log('    ' + l)
    return not errs and not undef and not cites


def words(tex):
    s = open(tex, encoding='utf-8').read()
    body = s[s.index(B + 'section{Introduction}'):s.index(B + 'section{Supplementary data}')]
    body = re.sub(r'(?m)%.*$', '', body)
    body = re.sub('(?s)' + E + r'begin\{(figure|table)\*?\}.*?' + E + r'end\{\1\*?\}', ' ', body)
    body = re.sub('(?s)' + E + r'begin\{equation\}.*?' + E + r'end\{equation\}', ' ', body)
    body = re.sub(E + r'(cite|ref|label)\{[^}]*\}', ' X ', body)
    body = re.sub(E + r'[a-zA-Z]+\*?', ' ', body)
    n_body = len(re.sub(r'[{}$~]', ' ', body).split())
    ab = re.search('(?s)' + E + r'abstract\{(.*?)\}' + NL + NL, s).group(1)
    n_ab = len(re.sub(r'[{}$~]', ' ', re.sub(E + r'[a-zA-Z]+', ' ', ab)).split())
    kp = len(re.findall(E + r'item ', s[s.index('Key Points'):s.index(B + 'maketitle')]))
    kw = len(re.search(E + r'keywords\{([^}]*)\}', s).group(1).split(';'))
    return n_body, n_ab, kp, kw


def numbers_traceable(tex):
    """BiB 본문의 수치가 모두 BMC 판에 있는가. 줄여 쓰다가 숫자를 잘못 옮기는 것을 막는다."""
    s = open(tex, encoding='utf-8').read()
    s = s[s.index(B + 'abstract{'):s.index(B + 'section{Conflicts of interest}')]
    s = re.sub(r'(?m)%.*$', '', s)
    s = re.sub(E + r'includegraphics\[[^]]*\]', '', s)
    s = re.sub(E + r'tabcolsep=[\d.]+pt', '', s)
    s = re.sub(r'p\{\d+mm\}', '', s)
    # 열 지정 문법은 수치가 아니다. \multicolumn{2}{c}{..} 의 2 를 본문 수치로 읽고
    # BMC 판에 없다고 막은 적이 있다.
    s = re.sub(E + r'multicolumn\{\d+\}\{[^}]*\}', '', s)
    s = re.sub(E + r'cmidrule(\(lr\))?\{[\d-]+\}', '', s)
    s = re.sub(E + r'(cite|ref|label)\{[^}]*\}', '', s)
    ref = re.sub(r'\s+', ' ', open(BMC, encoding='utf-8').read())
    ref_nums = set(re.findall(r'\d[\d,{}]*(?:\.\d+)?', ref))
    ref_nums = {x.replace('{,}', ',') for x in ref_nums}
    missing = []
    for tok in re.findall(r'(?<![A-Za-z\d])\d[\d,{}]*(?:\.\d+)?', s):
        tok = tok.replace('{,}', ',').rstrip(',')
        if tok not in ref_nums and tok.replace(',', '') not in ref_nums:
            missing.append(tok)
    return sorted(set(missing))


def leftovers():
    """BiB 로 나가는 파일에 BMC 판 표기가 남았는가. PDF 는 글자층을, XLSX 는 모든 셀을 읽는다."""
    import pymupdf
    import openpyxl
    # 참고문헌의 학술지명 BMC Bioinformatics 는 잔재가 아니므로 BMC 를 따로 찾지 않는다
    pat = re.compile(r'Additional\s+file|Additional_file|Sections?\s?\d')
    hits = []
    for f in sorted(os.listdir(DST)):
        p = os.path.join(DST, f)
        texts = []
        if f.startswith('Supplementary_File') and f.endswith('.pdf'):
            texts = [pg.get_text() for pg in pymupdf.open(p)]
        elif f.startswith('Supplementary_File') and f.endswith('.xlsx'):
            wb = openpyxl.load_workbook(p, read_only=True)
            texts = [str(c) for ws in wb.worksheets for row in ws.iter_rows(values_only=True)
                     for c in row if isinstance(c, str)]
        elif f in ('bib-manuscript.pdf', 'bib-supplement.pdf'):
            texts = [pg.get_text() for pg in pymupdf.open(p)]
        for t in texts:
            for m in pat.finditer(t):
                hits.append('%s: %s' % (f, t[max(0, m.start() - 30):m.end() + 20].replace(NL, ' ')))
    return sorted(set(hits))[:20]


# 핵심 주장은 숫자가 어딘가 있는지가 아니라, 측정량·비교군·값·방향이 함께 맞아야 한다.
# 아래 각 항목은 (설명, BiB 본문에서 찾을 정규식, BMC 본문에서 찾을 정규식) 이다. 숫자만
# 보는 검사는 부호나 분모가 바뀌어도 통과했다.
CLAIMS = [
    ('3코호트 Jaccard 구간(안정성-연결성)',
     r'cross-fold agreement[^.]*?ranged from \$-0\.008\$ to \$\+0\.174\$',
     r'cross-fold agreement between\s+stability selection and the connectivity comparator ranged from \$-0\.008\$ to \$\+0\.174\$'),
    ('신장 3코호트 Jaccard 구간',
     r'RBS \$-\$ connectivity in cross-fold agreement spanned\s+\$-0\.106\$ to \$\+0\.170\$',
     r'RBS \$-\$ connectivity in cross-fold agreement spans \$-0\.106\$ to \$\+0\.170\$'),
    ('외부 AUROC 이점은 평균이 양수',
     r'mean advantage[^.]*?positive at every subset size|mean external AUROC advantage',
     r'mean external AUROC advantage|advantage is positive at every size'),
    ('P3 상위 50 즉시초기 유전자 수',
     r'8\.0 immediate-early\s+genes in kidney, 0\.6 in liver and 0\.0 in colon',
     r'8\.0 immediate-early\s+genes in kidney, 0\.6 in liver and 0\.0 in colon'),
    ('비대칭 코호트에서 모듈이 사례에서 더 낮음',
     r'lower in cases in 9 of 10 cohorts across kidney and liver',
     r'lower in cases in 9 of 10 cohorts across\s+the two tissues'),
    ('대장 대칭 코호트에서 모듈이 사례에서 더 높음',
     r'higher in inflamed\s+cases in 9 of 10 cohorts \(\$\+0\.51\$ to \$\+3\.43\$\)',
     r'higher in inflamed cases in 9 of 10 cohorts \(\$\+0\.51\$ to \$\+3\.43\$\)'),
    ('효과 일치도(대장>신장>간)',
     r'0\.66\s+in colon against 0\.22 in kidney and 0\.15 in liver',
     r'was 0\.66 in colon against 0\.22 in kidney and 0\.15 in liver'),
    ('대장 축소 재표집에서 반전 2/5, 최대 0.013',
     r'two of five replicates, at most 0\.013|crossed zero in two of them by at most 0\.013',
     r'two of five\s*replicates, at most 0\.013|crossed zero in two of them, by at most 0\.013'),
    ('보정 전후 절차 유래 유전자 비율',
     r'procurement-flagged genes from 27\\% to 5--6\\%|from 27\\% to 5\\%',
     r'from 27\\% to 5--6\\%'),
    ('간의 반전이 다섯 시드에서 같은 부분집합에 재현',
     r'five bootstrap\s+seeds reproduced it in the same subset[^.]*?\$-0\.017\$ to \$-0\.008\$',
     r'five bootstrap seeds reproduced\s+it in the same subset[^.]*?\$-0\.017\$ to \$-0\.008\$'),
    ('대조군 치환 후 남는 효과 크기',
     r'retained only 14--56\\% of the naive',
     r'retained only 14--56\\% of the naive'),
]


def claims_match():
    """측정량과 방향까지 함께 대조한다. 한쪽에만 있으면 실패로 본다."""
    bib = re.sub(r'\s+', ' ', open(os.path.join(DST, 'bib-manuscript.tex'), encoding='utf-8').read())
    bmc = re.sub(r'\s+', ' ', open(BMC, encoding='utf-8').read())
    bad = []
    for label, p_bib, p_bmc in CLAIMS:
        ok_b = re.search(re.sub(r'\\s\+', ' ', p_bib), bib) is not None
        ok_m = re.search(re.sub(r'\\s\+', ' ', p_bmc), bmc) is not None
        log('  %-38s BiB %s  BMC %s' % (label, 'OK' if ok_b else '***', 'OK' if ok_m else '***'))
        if not (ok_b and ok_m):
            bad.append(label)
    return bad


def main():
    os.makedirs(DST, exist_ok=True)
    figures()
    bst()
    supplement()
    ok = compile_tex('bib-manuscript')
    ok = compile_tex('bib-supplement') and ok
    # 본문이 번호로 가리키는 보충 그림이 실제로 그 번호인가
    aux = open(os.path.join(DST, 'bib-supplement.aux'), encoding='utf-8', errors='replace').read()
    m = re.search(E + r'newlabel\{fig:xtissue-full\}\{\{([^}]+)\}', aux)
    said = re.findall(r'Supplementary Figure~(S\d+)', open(os.path.join(DST, 'bib-manuscript.tex'), encoding='utf-8').read())
    got = m.group(1) if m else '?'
    log('  본문이 부른 보충 그림 %s, 실제 코호트별 그림 %s %s' % (sorted(set(said)), got,
        'OK' if set(said) <= {got} else '*** 불일치'))
    ok = ok and set(said) <= {got}
    n_body, n_ab, kp, kw = words(os.path.join(DST, 'bib-manuscript.tex'))
    log('  본문 %d단어 (BiB 문제 해결 프로토콜 2,000-5,000) %s' % (n_body, 'OK' if n_body <= 5000 else '*** 초과'))
    log('  초록 %d단어, Key Points %d개 (3-5) %s, 키워드 %d개 (최대 6) %s'
        % (n_ab, kp, 'OK' if 3 <= kp <= 5 else '***', kw, 'OK' if kw <= 6 else '***'))
    miss = numbers_traceable(os.path.join(DST, 'bib-manuscript.tex'))
    log('  BMC 판에 없는 수치: %s' % (', '.join(miss) if miss else '없음'))
    log('  핵심 주장 대조 (측정량·비교군·값·방향):')
    bad_claims = claims_match()
    ok = ok and not bad_claims
    # 부록은 BMC 판을 복사하지 않고 BiB 표기로 다시 만든다. 복사본은 제목이 Additional file 이고
    # 본문 절 번호가 BMC 판 것이었다.
    r = subprocess.run([sys.executable, 'scripts/additional_files.py', '--target', 'bib'],
                       capture_output=True, env=dict(os.environ, PYTHONIOENCODING='utf-8'))
    if r.returncode != 0:
        log(r.stderr.decode('utf-8', 'replace')[-800:])
        ok = False
    left = leftovers()
    log('  보충 파일의 BMC 잔재: %s' % ('; '.join(left) if left else '없음'))
    ok = ok and not left
    stale = os.path.join(DST, 'P7_candidates_full.pdf')
    if os.path.exists(stale):
        os.remove(stale)
    for f in os.listdir(DST):
        if f.endswith(('.aux', '.blg', '.out', '.log')):
            os.remove(os.path.join(DST, f))
    # 제출 위생도 여기서 본다. submission/ 만 검사하고 BiB 폴더를 검사했다고 착각한 적이 있다.
    r2 = subprocess.run([sys.executable, 'scripts/preflight.py', '--target', 'bib'],
                        capture_output=True, env=dict(os.environ, PYTHONIOENCODING='utf-8'))
    tail = r2.stdout.decode('utf-8', 'replace') + r2.stderr.decode('utf-8', 'replace')
    ready = [l for l in tail.splitlines() if 'READY' in l]
    log('  제출 위생(BiB 폴더): %s' % (ready[0].strip() if ready else '판정 없음'))

    return 0 if ok and n_body <= 5000 and not miss else 1


if __name__ == '__main__':
    sys.exit(main())
