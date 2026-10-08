#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""GigaScience 판을 만든다 (OUP Transfer Desk, 2026-10-06, BiB 데스크 리젝 뒤).

BMC 판(submission/dkd-manuscript.tex)이 감사를 거친 전체 본문이다. GigaScience Research 논문의 절 순서:

    Abstract(Background/Results/Conclusions, 250 단어) - Keywords - Background - Data Description -
    Analyses - Discussion - Potential Implications - Methods - Availability of Source Code and
    Requirements - Data Availability - Declarations - References - Additional Files

BiB 판처럼 본문을 다시 쓰지 않는다. BMC 본문을 절 단위로 떼어 이 순서로 다시 붙이고, 초록만 250 단어로
줄여 쓴다(새 수치 없음; 아래 검사가 BMC 본문에 없는 수치를 잡는다).

조판은 저널의 LaTeX 템플릿(oup-contemporary.cls, \\journal{gigascience}; scripts/templates/gigascience/)
으로 한다 (2026-10-08, 사용자 요청). 그 클래스는 절 번호를 찍지 않으므로 Section~\\ref 는 절 제목으로
바꾼다(BiB 판과 같은 방식); 두 단 조판이라 그림·표는 전부 figure*/table* 로 둔다. 참고문헌 양식은
클래스가 정한다(vancouver-authoryear.bst, num-refs).

    python scripts/giga_build.py
"""
import os
import re
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
sys.path.insert(0, ROOT)
import submission_config as cfg  # noqa: E402

SRC = 'submission'
BMC = os.path.join(SRC, 'dkd-manuscript.tex')
DST = 'submission_giga'
NAME = 'giga-manuscript'
TEMPLATE = os.path.join('scripts', 'templates', 'gigascience')
TEMPLATE_FILES = ['oup-contemporary.cls', 'vancouver-authoryear.bst', 'gigascience-logo.pdf', 'oup.pdf']
B = chr(92)
E = re.escape(B)
NL = chr(10)
N_FIG = 10

TITLE = cfg.PROJECT['release_title']

# 250 단어 한도. BMC 초록(350 단어)의 주장·한정·수치만 남겼다. 새 수치는 넣지 않는다.
ABSTRACT = {
    'Background': (
        "Cross-cohort feature selection is widely used to derive reproducible transcriptomic "
        "biomarker signatures from public cohorts. This assumes that reproducibility reflects "
        "disease biology and that the benchmark measures the methods rather than the cohorts. We "
        "tested both in diabetic kidney disease, then in liver and colon against predictions fixed "
        "in advance."),
    'Results': (
        "Over all 26 subsets of five discovery cohorts, reproducibility verdicts were unstable: a "
        "stability-selection strategy won cross-fold agreement in 9 of 10 three-cohort subsets yet "
        "lost on all five, while its external AUROC advantage stayed positive. The most "
        "reproducible signal tracked procurement, not disease: an immediate-early gene module "
        "separated biopsy from surgical tissue in 8 of 9 kidney datasets, and a biopsy control "
        "retained only 14--56" + B + "% of the naive signal across seven diagnoses. Residualising "
        "on a handling score cut procurement-driven genes from 27" + B + "% to 5" + B + "% with no "
        "AUROC loss detected on a wide interval; batch correction by cohort was inert. In liver "
        "and colon the module was lower in biopsy cases than surgical controls in 9 of 10 such "
        "cohorts, verdict instability recurred in liver but not colon, and the module rose with "
        "inflammation, failing one prediction."),
    'Conclusions': (
        "Without a cohort-composition sensitivity analysis, a cross-cohort reproducibility "
        "comparison describes the cohorts rather than the methods, and should be reported "
        "separately from predictive performance. A shared pre-analytical confounder becomes "
        "reproducible wherever the design is shared; removing it rewrites the candidate list, 14 "
        "of 30 genes surviving. We release the harmonised cohorts, the schema-enforced audit "
        "database, the candidate-evidence ledger and the regenerating pipeline."),
}

# GigaScience 가 요구하는 목록. 값은 저장소·requirements.txt 와 같아야 한다 (아래에서 대조).
SOURCE_CODE = [
    ('Project name', 'dkd-cross-cohort-audit'),
    ('Project home page', B + 'url{%s}' % cfg.REPOSITORY['url']),
    ('Operating system(s)', 'Platform independent (developed and run on Windows 11)'),
    ('Programming language', 'Python 3.14; R 4.6.1 for the reference WGCNA comparator only'),
    ('Other requirements', 'Python packages pinned in requirements.txt (numpy 2.5.2, scipy 1.18.1, '
                           'pandas 3.0.5, scikit-learn 1.9.0, matplotlib 3.11.1, pyarrow 25.0.1, '
                           'h5py 3.16.0); R package WGCNA 1.74; a LaTeX distribution to rebuild '
                           'the manuscript'),
    ('License', 'MIT'),
    ('Restrictions to use by non-academics', 'None'),
    ('Archived version', B + 'url{%s} (version %s)' % (cfg.REPOSITORY['doi'], cfg.PROJECT['version'])),
]

AI_SECTION_TITLE = 'Disclosure of use of AI-assisted tools including generative AI'

# 빼는 그림과, 그 그림을 가리키던 문장의 대체문 (2026-10-08). 개요 도식은 글로, 후보 매트릭스는
# 전체 30개를 담은 Additional file 5 로 돌린다. 수치는 바뀌지 않는다.
DROP_FIGS = {
    'fig:overview': [
        ('Figure~' + B + 'ref{fig:overview} sets out the study workflow: audit the benchmark input,',
         'The study proceeds in four steps: audit the benchmark input,'),
    ],
    'fig:candidates': [
        ('Figure~' + B + 'ref{fig:candidates} shows what stands behind each of the 30, and why none',
         'Additional file~5 shows what stands behind each of the 30, and why none'),
        ('Candidate evidence over all 30 specificity-gated candidates: Figure~' + B + 'ref{fig:candidates} extended to include the 18 tier-4 candidates omitted from the main text for legibility.',
         'Candidate evidence over all 30 specificity-gated candidates: for each candidate the selection and '
         'control tests, the human-kidney context, the literature status, and the classification by current '
         'evidence strength and prior DKD literature.'),
    ],
}
# 두 단 조판의 전체 폭(488pt)에 BMC 폭(372pt)으로 그린 그림을 넣으면 1.31 배로 커진다. 0.9 폭이면 1.18 배,
# 그림 최소 글자 6.0pt 가 7.1pt 로 저널 하한(7pt)을 넘는다 (0.85 는 6.7pt 로 미달).
FIG_WIDTH = '0.9'

# 저널 요구: 범례 바로 아래 'Alt text:'. 수치는 범례·본문에 있는 것만 쓴다 (빌드가 대조한다).
ALT = {
    'fig:overview': (
        "Flow diagram of the study. Top row: audit the benchmark input for patient overlap and control "
        "procurement; then two sensitivity analyses, cohort composition and control definition; then "
        "diagnosis and correction of a procurement-associated immediate-early programme; then an "
        "auditable candidate-evidence ledger. A lower band separates supporting layers that corroborate "
        "or falsify candidates from participant-matched integrative multi-omics, which is marked as not "
        "claimed."),
    'fig:null': (
        "Dot plot of external AUROC by held-out cohort for three groups: random 50-gene signatures, "
        "permuted-label controls and the selected signatures. Random signatures form a band from about "
        "0.68 to 0.88 that differs by fold; permuted labels sit near 0.50; the selected signatures are "
        "shown against that band for each fold."),
    'fig:sensitivity': (
        "Three panels. Left, stacked bars by benchmark size (three, four, five cohorts) of the "
        "percentage of subsets in which each selection method wins on cross-fold reproducibility. "
        "Middle, one point per subset of stability selection minus the connectivity comparator in "
        "cross-fold agreement, spanning minus 0.106 to plus 0.170 with points on both sides of zero at "
        "every size except five. Right, the mean external AUROC difference with its range, positive at "
        "every size."),
    'fig:procurement': (
        "Two panels. Left, a dumbbell chart of nine kidney datasets in which biopsy-procured groups "
        "(orange) sit below nephrectomy- or donor-procured groups (blue) on the immediate-early module "
        "score in eight of nine; GSE162830 shows no gap. Right, horizontal bars of Hedges g for four "
        "contrasts: in KPMP, diabetic kidney disease against non-biopsy reference tissue minus 0.84 and "
        "against other biopsies plus 0.08; in heart failure, disease against donor minus 0.23 and "
        "ischaemic against idiopathic plus 0.33."),
    'fig:singlecell': (
        "Box plots of the immediate-early module score in KPMP single-nucleus pseudobulk by donor. "
        "Diabetic kidney disease against non-diabetic chronic kidney disease, both needle biopsies, "
        "shows overlapping distributions (p equals 0.21); diabetic kidney disease against non-biopsy "
        "reference tissue shows a clear separation."),
    'fig:blocksize': (
        "Two line panels over the ratio of disease-block size to confounder-block size in a simulation. "
        "Left, the percentage of the top 50 features drawn from the confounder block for each criterion "
        "as the ratio changes. Right, the difference between the connectivity and stability criteria, "
        "crossing zero where the ratio equals one. A marker shows the real-data regime, 19 "
        "immediate-early genes against 705 top-connectivity genes."),
    'fig:corrections': (
        "Three panels of horizontal bars for six correction strategies. Left, the percentage of the top "
        "50 selected features flagged as procurement-driven: 27 percent with no correction or ComBat, "
        "5 to 6 percent with handling-score residualisation or surrogate-variable analysis. Middle, "
        "specificity pass rate as fold enrichment over the 8.3 percent background. Right, external "
        "AUROC, 0.926 for no correction and residualisation and 0.915 for surrogate-variable analysis."),
    'fig:reordering': (
        "Two panels. Left, a slope chart of 46 genes drawn between their rank against non-biopsy "
        "controls and their rank against other-diagnosis biopsies in diabetic nephropathy, on a "
        "logarithmic rank axis; blue lines rise under the matched control and orange lines fall, 23 "
        "each. Right, horizontal bars of the number of genes shared by the two top-50 lists for seven "
        "ERCB biopsy diagnoses."),
    'fig:xtissue': (
        "Three panels. Left, a dot plot of Hedges g for the handling score, case minus control, with 95 "
        "percent bootstrap intervals, one row per cohort, split into cohorts with biopsy cases against "
        "surgical or donor controls (above the dashed line) and cohorts with both arms obtained by the "
        "same procedure (below); hollow points show the inflammation-adjusted value. Middle, bars of the "
        "mean number of immediate-early genes in the stability-selection top 50 per tissue, with and "
        "without correction. Right, per-subset reproducibility differences for kidney, liver and colon "
        "with reversed verdicts highlighted."),
    'fig:candidates': (
        "Two panels. Left, a matrix of 12 candidate genes (rows) against evidence criteria (columns) "
        "grouped as selection and control tests, human-kidney context and literature status; orange "
        "cells are criteria met, pale cells criteria not met, dashes unavailable assessments. Right, all "
        "30 candidates placed on a grid of current evidence strength against prior diabetic kidney "
        "disease literature, showing that no candidate is both strong and novel."),
}

# 저널 상한 15 단어. 제목만 바꾸고 범례는 그대로 둔다.
TITLE_FIX = {
    'The reproducibility verdict depends on which cohorts are in the benchmark; the mean predictive advantage stays positive.':
        'The reproducibility verdict depends on which cohorts are benchmarked; the predictive advantage stays positive.',
    'Outside kidney, the immediate-early module follows procurement where it is asymmetric and inflammation where it is not, and reaches the most reproducible genes mainly in kidney.':
        'Outside kidney, the module follows procurement where asymmetric and inflammation where matched.',
}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def between(s, start_pat, end_pat, name):
    a = re.search(start_pat, s)
    if not a:
        raise SystemExit('%s: 시작 표지를 못 찾음 %r' % (name, start_pat))
    b = re.search(end_pat, s[a.start():])
    if not b:
        raise SystemExit('%s: 끝 표지를 못 찾음 %r' % (name, end_pat))
    return s[a.start():a.start() + b.start()]


def front_matter():
    """템플릿 매크로: \\author[affil\\authfn{n}]{}, \\affil[]{}, \\authnote{}, \\papercat{}."""
    lines = [B + 'title{' + TITLE + '}', '']
    corr = [a for a in cfg.AUTHORS if a.get('corresponding')]
    for a in cfg.AUTHORS:
        mark = B + 'authfn{1}' if a.get('corresponding') else ''
        lines.append(B + 'author[%d%s]{%s %s}' % (a['affiliation'], mark, a['given'], a['family']))
    lines.append('')
    for num, af in sorted(cfg.AFFILIATIONS.items()):
        lines.append(B + 'affil[%d]{%s, %s, %s, %s}' % (num, af['division'], af['organisation'], af['city'], af['country']))
    lines.append('')
    lines.append(B + 'authnote{' + B + 'authfn{1}Corresponding author. E-mail: ' + '; '.join(a['email'] for a in corr) + '}')
    emails = '; '.join('%s %s, %s' % (a['given'], a['family'], a['email']) for a in cfg.AUTHORS if a.get('email'))
    lines.append(B + 'authnote{Author e-mail addresses: ' + emails + '.}')
    lines.append('')
    lines.append(B + 'papercat{Research}')
    lines.append(B + 'runningauthor{%s et al.}' % cfg.AUTHORS[0]['family'])
    lines.append(B + 'jvolume{00}' + B + 'jnumber{0}' + B + 'jyear{%s}' % cfg.PROJECT['year'])
    return NL.join(lines)


def abstract_words():
    plain = lambda t: re.sub(r'\\[a-zA-Z]+|[${}~]', '', t).replace('--', ' ').split()
    return sum(len(plain(v)) for v in ABSTRACT.values())


def numbers_in(text):
    return set(re.findall(r'(?<![\w.])\d+(?:[.,]\d+)?(?![\w])', text))


def drop_figures(body):
    """DROP_FIGS 의 그림 환경을 지우고 참조 문장을 대체문으로; 남은 그림은 등장 순서로 Fig1.. 로 다시 이름 붙인다.
    돌려주는 사전은 {새 이름: BMC 이름} 이다 (파일 복사용)."""
    for lab, pairs in DROP_FIGS.items():
        # 다른 그림 환경을 넘어가지 않도록 \end{figure} 를 건너뛰지 못하게 한다 (처음엔 2-10번 그림을 통째로 지웠다)
        pat = re.compile(r'(?s)' + E + r'begin\{figure\}(?:(?!' + E + r'end\{figure\}).)*?' + E + r'label\{' + re.escape(lab) + r'\}' + NL + E + r'end\{figure\}' + NL + r'?')
        if len(pat.findall(body)) != 1:
            raise SystemExit('그림 %s 환경을 %d 번 찾음' % (lab, len(pat.findall(body))))
        body = pat.sub('', body, count=1)
        for old, new in pairs:
            if body.count(old) != 1:
                raise SystemExit('그림 %s 참조 문장을 못 찾음: %s' % (lab, old[:60]))
            body = body.replace(old, new)
        if re.search(E + r'ref\{' + re.escape(lab) + r'\}', body):
            raise SystemExit('그림 %s 를 가리키는 참조가 남아 있음' % lab)
    names = {}
    def rename(m):
        new = 'Fig%d' % (len(names) + 1)
        names[new] = m.group(2)
        return m.group(1) + new + '}'
    body = re.sub(r'(' + E + r'includegraphics\[[^]]*\]\{)(Fig\d+)\}', rename, body)
    body = body.replace(B + 'includegraphics[width=' + B + 'textwidth]{Fig',
                        B + 'includegraphics[width=' + FIG_WIDTH + B + 'textwidth]{Fig')
    return body, names


def figure_text(body):
    """범례 아래 Alt text 삽입, 긴 제목 축약, 제목·범례 길이 검사 (15 / 300 단어)."""
    bad = []
    for lab, alt in ALT.items():
        if lab in DROP_FIGS:
            continue
        pat = B + 'label{' + lab + '}' + NL
        if body.count(pat) != 1:
            raise SystemExit('그림 라벨 %s 를 %d 번 찾음' % (lab, body.count(pat)))
        body = body.replace(pat, pat + B + 'smallskip' + B + 'noindent{' + B + 'footnotesize Alt text: ' + alt + '}' + NL)
    for old, new in TITLE_FIX.items():
        # 원문의 제목은 줄바꿈으로 감겨 있다
        pat = re.compile(re.escape(B + 'textbf{') + r'\s+'.join(map(re.escape, old.split())) + re.escape('}'))
        if len(pat.findall(body)) != 1:
            raise SystemExit('그림 제목을 못 찾음: %s' % old[:50])
        body = pat.sub(lambda _: B + 'textbf{' + new + '}', body, count=1)
    # figure 환경 안에서만 찾는다; 표의 \caption 에서 시작해 다음 그림 라벨까지 이어 붙던 것을 막는다
    for blk in re.finditer(r'(?s)' + E + r'begin\{figure\}(.*?)' + E + r'end\{figure\}', body):
        m = re.search(r'(?s)' + E + r'caption\{(.*)\}' + E + r'label\{(fig:[^}]+)\}', blk.group(1))
        if not m:
            continue
        cap, lab = m.group(1), m.group(2)
        title = re.search(E + r'textbf\{([^}]*)\}', cap)
        nt = len(title.group(1).split()) if title else 0
        nl = len(re.sub(r'\\[a-zA-Z]+|[{}$]', '', cap).split())
        if nt > 15 or nl > 300:
            bad.append('%s 제목 %d 단어, 범례 %d 단어' % (lab, nt, nl))
    if bad:
        raise SystemExit('그림 제목/범례 길이 초과: ' + '; '.join(bad))
    return body


def section_refs(body):
    """클래스가 절 번호를 찍지 않으므로 Section~\\ref{} 를 절 제목으로 바꾼다 (BiB 판과 같은 방식)."""
    titles = {}
    for m in re.finditer(E + r'(?:sub)?section\*?\{([^}]*)\}' + E + r'label\{(sec:[^}]+)\}', body):
        titles[m.group(2)] = ' '.join(m.group(1).split())
    missing = sorted(set(re.findall(E + r'ref\{(sec:[^}]+)\}', body)) - set(titles))
    if missing:
        raise SystemExit('제목을 모르는 절 라벨: %s' % ', '.join(missing))
    body = re.sub(r'Sections~' + E + r'ref\{(sec:[^}]+)\}', lambda m: 'the sections ' + B + 'textit{' + titles[m.group(1)] + '}', body)
    body = re.sub(r'Section~' + E + r'ref\{(sec:[^}]+)\}', lambda m: 'the section ' + B + 'textit{' + titles[m.group(1)] + '}', body)
    body = re.sub(r'(and|,)~' + E + r'ref\{(sec:[^}]+)\}', lambda m: m.group(1) + ' ' + B + 'textit{' + titles[m.group(2)] + '}', body)
    body = re.sub(E + r'ref\{(sec:[^}]+)\}', lambda m: B + 'textit{' + titles[m.group(1)] + '}', body)
    # 문장 첫머리의 the -> The. 줄 첫머리는 안 된다: 원문이 줄바꿈으로 감겨 있어 "of\nSection" 이
    # 문장 중간에서도 줄 첫머리에 온다 (그렇게 해서 "of The section" 이 일곱 군데 생겼다).
    body = re.sub(r'(\A|(?<=\n\n)|(?<=[.!?:]\s)|(?<=[.!?:]\n))the (section|sections) ' + E + 'textit',
                  lambda m: m.group(1) + 'The ' + m.group(2) + ' ' + B + 'textit', body)
    return body


def wide_floats(body):
    """두 단 조판: 그림·표는 전부 단 너비를 넘으므로 figure*/table* 로."""
    # 치환 문자열을 그대로 주면 re 가 '\b' 를 백스페이스로 읽는다; 람다로 넘긴다
    body = re.sub(E + r'begin\{figure\}\[[^]]*\]', lambda _: B + 'begin{figure*}[!t]', body)
    body = body.replace(B + 'end{figure}', B + 'end{figure*}')
    body = re.sub(E + r'begin\{table\}\[[^]]*\]', lambda _: B + 'begin{table*}[!t]', body)
    body = body.replace(B + 'end{table}', B + 'end{table*}')
    return body


def addfile_formats(text):
    """저널 요구: 파일 이름, 형식(확장자), 제목, 설명. 실제 파일의 확장자를 읽어 넣는다."""
    def fmt(m):
        n = int(m.group(1))
        ext = [f for f in os.listdir(DST) if f.startswith('Additional_file_%d.' % n)]
        if len(ext) != 1:
            raise SystemExit('Additional file %d 의 파일을 못 찾음' % n)
        return B + 'textbf{Additional file %d (%s).}' % (n, ext[0].rsplit('.', 1)[1].upper())
    out, k = re.subn(E + r'textbf\{Additional file (\d+)\.\}', fmt, text)
    if k != 7:
        raise SystemExit('Additional file 항목 %d 개 (7 이어야 함)' % k)
    return out


def build_tex():
    s = open(BMC, encoding='utf-8').read()
    body, fig_names = drop_figures(s[s.index(B + 'maketitle'):])
    body = figure_text(body)

    background = between(body, E + r'section\{Background\}', E + r'section\{Methods\}', 'Background')
    methods = between(body, E + r'section\{Methods\}', E + r'section\{Results\}', 'Methods')
    results = between(body, E + r'section\{Results\}', E + r'section\{Discussion\}', 'Results')
    discussion = between(body, E + r'section\{Discussion\}', E + r'section\{Conclusions\}', 'Discussion')
    conclusions = between(body, E + r'section\{Conclusions\}', E + r'section\*\{List of abbreviations\}', 'Conclusions')
    abbrev = between(body, E + r'section\*\{List of abbreviations\}', E + r'backmatter', 'Abbreviations')
    addfiles = between(body, E + r'bmhead\{Supplementary information\}', E + r'section\*\{Declarations\}', 'Additional files')
    decl = between(body, E + r'section\*\{Declarations\}', E + r'bibliography\{', 'Declarations')

    # Methods 를 셋으로: 머리 문단, 자료 절(코호트·조정 -> Data Description), 나머지
    m_head = between(methods, E + r'section\{Methods\}', E + r'subsection\{Cohort assembly and audit\}', 'Methods head')
    m_data = between(methods, E + r'subsection\{Cohort assembly and audit\}', E + r'subsection\{Selection methods\}', 'Methods data')
    m_rest = methods[methods.index(B + 'subsection{Selection methods}'):]

    # 머리 문단의 "the sections below" 는 자료 절이 앞으로 옮겨가면 틀린 말이 된다
    old = ('Section~' + B + 'ref{sec:cohorts} and' + NL + 'Section~' + B + 'ref{sec:harmonisation} cover assembly and harmonisation, Section~' + B + 'ref{sec:selectors}')
    if m_head.count(old) != 1:
        raise SystemExit('Methods 머리 문단의 참조 문장이 바뀌었습니다; giga_build.py 를 맞추십시오')
    # 번호 없는 절 제목으로는 "Section X the four strategies" 식 생략문이 읽히지 않으므로 문장을 다시 쓴다
    end = 'how each reported number is regenerated from source data.'
    a, b = m_head.index(old), m_head.index(end) + len(end)
    m_head = (m_head[:a]
              + 'Cohort assembly and harmonisation are described under Data description'
              + ' (Sections~' + B + 'ref{sec:cohorts} and~' + B + 'ref{sec:harmonisation}). The subsections'
              + ' below cover the four strategies (Section~' + B + 'ref{sec:selectors}), the nulls and the'
              + ' sensitivity sweeps (Section~' + B + 'ref{sec:nulls}), the confounder score and the corrections'
              + ' benchmarked against it (Section~' + B + 'ref{sec:corrections}), the test of both mechanisms in'
              + ' liver and colon (Section~' + B + 'ref{sec:xtmethods}), and how each reported number is'
              + ' regenerated from source data (Section~' + B + 'ref{sec:implementation}).'
              + m_head[b:])
    old2 = 'The study has two layers and the sections below follow them.'
    if m_head.count(old2) != 1:
        raise SystemExit('Methods 첫 문장이 바뀌었습니다')
    m_head = m_head.replace(old2, 'The study has two layers.')

    # 선언부에서 자료 가용성 블록을 떼어 GigaScience 의 Data Availability 절로
    avail = between(decl, E + r'bmhead\{Availability of data and materials\}', E + r'bmhead\{Competing interests\}', 'Availability')
    decl = decl.replace(avail, '')
    avail_body = avail[avail.index(NL):]  # 제목 줄 제거

    # AI 문단은 Methods 의 Implementation 절에 있다. 끝에 다시 둔다.
    ai = re.search(r'(?s)' + E + r'textbf\{Use of AI tools\.\}(.*?)(?=' + NL + NL + r')', m_rest)
    if not ai:
        raise SystemExit('Use of AI tools 문단을 못 찾음')
    ai_text = ai.group(1).strip()
    # 저널은 이 문단을 끝의 Disclosure 절에 요구한다; Methods 에 같은 문단을 두 번 두지 않는다
    m_rest = m_rest.replace(ai.group(0), '', 1)

    results = results.replace(B + 'section{Results}', B + 'section{Analyses}', 1)
    conclusions = conclusions.replace(B + 'section{Conclusions}', B + 'section{Potential implications}', 1)

    src_list = NL.join('  ' + B + 'item ' + B + 'textbf{%s:} %s' % kv for kv in SOURCE_CODE)

    pre = [
        '%% GigaScience Research article, OUP contemporary template (oup-contemporary.cls,',
        '%% \\journal{gigascience}). Generated by scripts/giga_build.py from',
        '%% submission/dkd-manuscript.tex; edit that file and rebuild, do not edit this one.',
        B + 'PassOptionsToPackage{hyphens}{url}',   # 저장소 주소가 하이픈에서 끊기게
        B + 'documentclass[a4paper,num-refs]{oup-contemporary}',
        B + 'journal{gigascience}',
        B + 'usepackage{graphicx,amsmath,array,multirow}',
        B + 'graphicspath{{figures/}}',
        B + 'providecommand{' + B + 'bibcommenthead}{}',
        B + 'providecommand{' + B + 'botrule}{' + B + 'bottomrule}',
        B + 'providecommand{' + B + 'bmhead}[1]{' + B + 'subsection{#1}}',
        B + 'renewcommand{' + B + 'footnotetext}[1]{' + B + 'par' + B + 'smallskip{' + B + 'footnotesize #1}}',
        B + 'emergencystretch=3em',
        B + 'hyphenation{OpenAI ChatGPT Anthropic}',   # 'Ope-nAI' 로 끊긴 적이 있다
        '',
        front_matter(),
        '',
        B + 'begin{document}',
        B + 'begin{frontmatter}',
        B + 'maketitle',
        B + 'begin{abstract}',
    ]
    for k, v in ABSTRACT.items():
        pre.append(B + 'textbf{%s:} %s' % (k, v) + B + 'par')
    pre += [
        B + 'end{abstract}',
        B + 'begin{keywords}',
        '; '.join(cfg.PROJECT['keywords']),
        B + 'end{keywords}',
        B + 'end{frontmatter}',
        '',
    ]

    parts = pre + [
        background,
        B + 'section{Data description}' + B + 'label{sec:data}' + NL,
        m_data,
        results,
        discussion,
        conclusions,
        m_head, m_rest,
        B + 'section{Availability of source code and requirements}' + NL + B + 'begin{itemize}' + NL
        + src_list + NL + B + 'end{itemize}' + NL,
        B + 'section{Data availability}' + NL + NL
        + 'The data sets supporting the results of this article are available in the Zenodo repository, '
        + B + 'url{%s}' % cfg.REPOSITORY['doi'] + ' ' + B + 'cite{kim2026zenodo}, with the primary public records listed below.'
        + avail_body,
        abbrev.replace(B + 'section*{List of abbreviations}', B + 'section{Declarations}' + NL + B + 'bmhead{List of abbreviations}'),
        decl.replace(B + 'section*{Declarations}' + NL, ''),
        B + 'bmhead{' + AI_SECTION_TITLE + '}' + NL + NL + ai_text
        + ' The tools are cited as ' + B + 'cite{claudecode2026,chatgpt2026}.' + NL,
        B + 'section{Additional files}' + NL + NL + addfile_formats(addfiles.replace(B + 'bmhead{Supplementary information}' + NL, '')),
        B + 'bibliography{dkd-references,giga-extra}',
        B + 'end{document}',
    ]
    out = NL.join(parts)
    out = out.replace(B + 'backmatter', '')
    out = out.replace('Fig.~' + B + 'ref', 'Figure~' + B + 'ref')   # BMC 판의 'Fig. 2' 와 'Figure 2' 혼용을 통일
    out = wide_floats(section_refs(narrow_equations(out)))
    # 점검표(Table 6)의 셋째 열은 BMC 판에서 "Section 3.2" 가 들어가던 28mm 폭이다. 절 제목이 들어가면
    # 네 줄로 감기므로, 두 단 전체 폭(172mm)에 맞춰 넓힌다.
    spec_old = B + 'begin{tabular}{@{}p{46mm}p{46mm}p{28mm}@{}}'
    if out.count(spec_old) != 1:
        raise SystemExit('점검표의 열 규격을 %d 번 찾음' % out.count(spec_old))
    out = out.replace(spec_old, B + 'begin{tabular}{@{}p{50mm}p{56mm}p{60mm}@{}}')
    return out, s, fig_names


def narrow_equations(body):
    """한 줄에 두 식을 \qquad 로 나란히 둔 수식은 한 단 폭을 넘는다. gathered 로 두 줄에 둔다."""
    def fix(m):
        inner = m.group(1)
        if B + 'qquad' not in inner:
            return m.group(0)
        lab = re.search(E + r'label\{[^}]+\}', inner)
        content = inner.replace(lab.group(0), '') if lab else inner
        # 치환 문자열의 '\\' 는 re 가 '\' 하나로 줄인다; 람다로 넘긴다
        content = re.sub(r',\s*' + E + r'qquad\s*', lambda _: ', ' + B + B + NL, content).strip()
        return (B + 'begin{equation}' + B + 'begin{gathered}' + NL + content + NL + B + 'end{gathered}'
                + (lab.group(0) if lab else '') + B + 'end{equation}')
    return re.sub(r'(?s)' + E + r'begin\{equation\}(.*?)' + E + r'end\{equation\}', fix, body)


def source_zip(fig_names):
    """저널에 올릴 LaTeX 소스 묶음: 컴파일에 필요한 것만 (tex, bbl, bib 둘, 템플릿 네 파일, 그림)."""
    import zipfile
    out = os.path.join(DST, NAME + '-source.zip')
    names = [NAME + '.tex', NAME + '.bbl', 'dkd-references.bib', 'giga-extra.bib'] + TEMPLATE_FILES
    names += ['figures/' + new + '.png' for new in sorted(fig_names, key=lambda n: int(n[3:]))]
    # 저널 시스템(Editorial Manager)은 하위 폴더가 있는 묶음을 거부한다 (2026-10-08). 그림을 루트에
    # 평평하게 넣고, 묶음 안의 tex 만 그림 경로를 루트로 돌린다. 폴더의 tex 와의 차이는 그 한 줄뿐이다.
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as z:
        for n in names:
            path = os.path.join(DST, n)
            if not os.path.exists(path):
                raise SystemExit('소스 묶음에 넣을 파일이 없습니다: %s' % path)
            if n == NAME + '.tex':
                tex = open(path, encoding='utf-8').read()
                gp = B + 'graphicspath{{figures/}}'
                if tex.count(gp) != 1:
                    raise SystemExit('graphicspath 줄을 못 찾음')
                z.writestr(n, tex.replace(gp, B + 'graphicspath{{./}}'))
            else:
                z.write(path, os.path.basename(n))
    return out, len(names)


def compile_tex():
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    ok = True
    for cmd in (['pdflatex', '-interaction=nonstopmode', NAME], ['bibtex', NAME],
                ['pdflatex', '-interaction=nonstopmode', NAME], ['pdflatex', '-interaction=nonstopmode', NAME]):
        r = subprocess.run(cmd, cwd=DST, capture_output=True, env=env)
        if r.returncode != 0 and cmd[0] == 'bibtex':
            log(r.stdout.decode('utf-8', 'replace')[-600:])
            ok = False
    logtxt = open(os.path.join(DST, NAME + '.log'), encoding='utf-8', errors='replace').read()
    undef = sorted(set(re.findall(r"(?:Reference|Citation) `([^']+)' on page", logtxt)))
    errs = re.findall(r'(?m)^! (.*)$', logtxt)
    over = [float(x) for x in re.findall(r'Overfull ' + E + r'hbox \(([\d.]+)pt too wide', logtxt)]
    # 템플릿이 초록 위아래에 그리는 장식선은 \hsize+4pc+6.64pt 로 일부러 54.64pt 넘친다
    over = [o for o in over if abs(o - 54.64) > 0.05]
    shutil.copy2(os.path.join(DST, NAME + '.log'), os.path.join('reports', 'giga-manuscript.log'))
    log('  컴파일: 오류 %d, 미해결 참조 %d %s, 넘침 >10pt %d (최대 %.1fpt)' % (
        len(errs), len(undef), undef[:6], sum(1 for o in over if o > 10), max(over or [0])))
    if errs:
        log('   ' + errs[0][:200])
    return ok and not errs and not undef and not [o for o in over if o > 10]


def main():
    os.makedirs(DST, exist_ok=True)
    os.makedirs(os.path.join(DST, 'figures'), exist_ok=True)
    tex, bmc, fig_names = build_tex()
    open(os.path.join(DST, NAME + '.tex'), 'w', encoding='utf-8').write(tex)

    # 초록 검사: 길이와 수치 출처
    n = abstract_words()
    nums = numbers_in(' '.join(ABSTRACT.values())) - numbers_in(bmc)
    log('  초록 %d 단어 (한도 250) %s; BMC 본문에 없는 수치: %s' % (n, 'OK' if n <= 250 else '*** 초과', sorted(nums) or '없음'))
    ok = n <= 250 and not nums
    alt_nums = numbers_in(' '.join(v for k, v in ALT.items() if k not in DROP_FIGS)) - numbers_in(bmc)
    log('  Alt text %d개; BMC 본문에 없는 수치: %s' % (len(ALT) - len(DROP_FIGS), sorted(alt_nums) or '없음'))
    ok = ok and not alt_nums
    # 코드 요구사항 목록의 판 번호는 requirements.txt 와 같아야 한다
    req = open('requirements.txt', encoding='utf-8').read()
    for pkg, ver in re.findall(r'([a-z-]+) (\d+\.\d+\.\d+)', SOURCE_CODE[4][1]):
        if '%s==%s' % (pkg, ver) not in req:
            log('  *** requirements.txt 와 다름: %s %s' % (pkg, ver))
            ok = False

    for f in os.listdir(os.path.join(DST, 'figures')):
        os.remove(os.path.join(DST, 'figures', f))
    # 그림은 PNG(600 dpi)로 올린다 (사용자 요청 2026-10-08; BiB 때 ScholarOne 이 PDF/EPS 변환에 실패했다).
    # tex 는 확장자 없이 부르므로 PNG 만 있으면 그것을 쓴다. 글자 크기 검사는 벡터 원본을 봐야 하므로
    # 새 이름 -> 원본 PDF 대응을 figures/sources.json 에 남긴다.
    import json
    for new, old in fig_names.items():
        src_pdf = os.path.join(SRC, old + '.pdf')
        r = subprocess.run(['pdftoppm', '-r', '600', '-png', '-singlefile', src_pdf,
                            os.path.join(DST, 'figures', new)], capture_output=True)
        if r.returncode != 0 or not os.path.exists(os.path.join(DST, 'figures', new + '.png')):
            raise SystemExit('PNG 변환 실패: %s' % src_pdf)
    json.dump({new: os.path.join(SRC, old + '.pdf').replace(os.sep, '/') for new, old in fig_names.items()},
              open(os.path.join(DST, 'figures', 'sources.json'), 'w', encoding='utf-8'), indent=1)
    log('  그림 %d개 PNG 600 dpi (BMC %s 제외), 폭 %s' % (len(fig_names), ', '.join(sorted(set('Fig%d' % i for i in range(1, N_FIG + 1)) - set(fig_names.values()))), FIG_WIDTH + B + 'textwidth'))
    shutil.copy2(os.path.join(SRC, 'dkd-references.bib'), os.path.join(DST, 'dkd-references.bib'))
    for f in TEMPLATE_FILES:
        shutil.copy2(os.path.join(TEMPLATE, f), os.path.join(DST, f))
    stale = os.path.join(DST, 'bib-vancouver.bst')   # 이전 article 판이 쓰던 양식
    if os.path.exists(stale):
        os.remove(stale)
    ok = compile_tex() and ok

    r = subprocess.run([sys.executable, 'scripts/additional_files.py', '--target', 'giga'],
                       capture_output=True, env=dict(os.environ, PYTHONIOENCODING='utf-8'))
    if r.returncode != 0:
        log(r.stderr.decode('utf-8', 'replace')[-800:])
        ok = False
    # additional_files.py 가 부록 5번에 넣는 그림 원본을 여기 두고 간다
    stale = os.path.join(DST, 'P7_candidates_full.pdf')
    if os.path.exists(stale):
        os.remove(stale)
    for f in os.listdir(DST):
        if f.endswith(('.aux', '.blg', '.out', '.log')):
            os.remove(os.path.join(DST, f))
    zp, nz = source_zip(fig_names)
    log('  LaTeX 소스 묶음 %s (%d 파일)' % (zp, nz))
    pages = re.search(r'Pages:\s+(\d+)', subprocess.run(['pdfinfo', os.path.join(DST, NAME + '.pdf')],
                                                      capture_output=True, text=True).stdout or '')
    log('  %s/%s.pdf %s쪽' % (DST, NAME, pages.group(1) if pages else '?'))
    log('GigaScience 판 %s' % ('OK' if ok else '*** 실패'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
