#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""GigaScience 판을 만든다 (OUP Transfer Desk, 2026-10-06, BiB 데스크 리젝 뒤).

BMC 판(submission/dkd-manuscript.tex)이 감사를 거친 전체 본문이다. GigaScience 는 초기 투고 형식이
자유롭지만 Research 논문의 절 순서를 정해 두었다:

    Abstract(Background/Results/Conclusions, 250 단어) - Keywords - Background - Data Description -
    Analyses - Discussion - Potential Implications - Methods - Availability of Source Code and
    Requirements - Data Availability - Declarations - References - Additional Files

BiB 판처럼 본문을 다시 쓰지 않는다. BMC 본문을 절 단위로 떼어 이 순서로 다시 붙이고, 초록만 250 단어로
줄여 쓴다(새 수치 없음; 아래 검사가 BMC 본문에 없는 수치를 잡는다). 절 번호는 그대로 둔다(상호참조와
부록이 번호를 쓴다). 더블 스페이스, 쪽 번호, 그림은 본문 자리에 넣고 파일로도 따로 둔다.

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

# AI 사용: 본문 Methods 의 문단을 그대로 끝에 다시 두고, 저널 규정대로 도구를 참고문헌으로 인용한다.
# 모델·사용 기간·사용자는 저자가 확인한다 (giga-extra.bib 의 대괄호 자리).
AI_SECTION_TITLE = 'Disclosure of use of AI-assisted tools including generative AI'


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


def author_block():
    marks = {}
    lines = []
    names = []
    for a in cfg.AUTHORS:
        star = B + 'textsuperscript{*}' if a.get('corresponding') else ''
        names.append('%s %s%stextsuperscript{%d}%s' % (a['given'], a['family'], B, a['affiliation'], star))
    lines.append(', '.join(names[:-1]) + ' and ' + names[-1] + B + B + '[6pt]')
    for num, af in sorted(cfg.AFFILIATIONS.items()):
        lines.append('{' + B + 'small ' + B + 'textsuperscript{%d}%s, %s, %s, %s}' % (
            num, af['division'], af['organisation'], af['city'], af['country']) + B + B)
    corr = next(a for a in cfg.AUTHORS if a.get('corresponding'))
    lines.append('{' + B + 'small ' + B + 'textsuperscript{*}Corresponding author: %s %s, ' % (corr['given'], corr['family'])
                 + B + 'href{mailto:%s}{%s}}' % (corr['email'], corr['email']) + B + B + '[2pt]')
    emails = '; '.join('%s %s, %s' % (a['given'], a['family'], a['email']) for a in cfg.AUTHORS if a.get('email'))
    lines.append('{' + B + 'small Author e-mail addresses: ' + emails + '}')
    return NL.join(lines)


def abstract_words():
    plain = lambda t: re.sub(r'\\[a-zA-Z]+|[${}~]', '', t).replace('--', ' ').split()
    return sum(len(plain(v)) for v in ABSTRACT.values())


def numbers_in(text):
    return set(re.findall(r'(?<![\w.])\d+(?:[.,]\d+)?(?![\w])', text))


def build_tex():
    s = open(BMC, encoding='utf-8').read()
    body = s[s.index(B + 'maketitle'):]

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
    m_head = m_head.replace(old, 'Cohort assembly and harmonisation are described under Data description'
                            ' (Sections~' + B + 'ref{sec:cohorts} and~' + B + 'ref{sec:harmonisation}); Section~' + B + 'ref{sec:selectors}')
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

    results = results.replace(B + 'section{Results}', B + 'section{Analyses}', 1)
    conclusions = conclusions.replace(B + 'section{Conclusions}', B + 'section{Potential implications}', 1)

    src_list = NL.join('  ' + B + 'item ' + B + 'textbf{%s:} %s' % kv for kv in SOURCE_CODE)

    pre = [
        '%% GigaScience Research article. Generated by scripts/giga_build.py from',
        '%% submission/dkd-manuscript.tex; edit that file and rebuild, do not edit this one.',
        B + 'documentclass[11pt,a4paper]{article}',
        B + 'PassOptionsToPackage{hyphens}{url}',   # 저장소 주소가 하이픈에서 끊기게
        B + 'usepackage[a4paper,margin=25mm]{geometry}',
        B + 'usepackage[T1]{fontenc}', B + 'usepackage{lmodern}', B + 'usepackage{cmap}',
        B + 'usepackage{graphicx,booktabs,array,multirow,amsmath,amssymb,amsfonts,textcomp,url,xcolor}',
        B + 'usepackage{setspace}',
        B + 'usepackage[numbers,square,sort&compress]{natbib}',
        B + 'usepackage[hidelinks]{hyperref}',
        B + 'graphicspath{{figures/}}',
        B + 'providecommand{' + B + 'bibcommenthead}{}',
        B + 'providecommand{' + B + 'botrule}{' + B + 'bottomrule}',
        B + 'providecommand{' + B + 'bmhead}[1]{' + B + 'subsection*{#1}}',
        B + 'renewcommand{' + B + 'footnotetext}[1]{' + B + 'par' + B + 'smallskip{' + B + 'footnotesize #1}}',
        B + 'setlength{' + B + 'parskip}{4pt}',
        B + 'emergencystretch=3em',
        B + 'renewcommand{' + B + 'topfraction}{0.9}' + B + 'renewcommand{' + B + 'bottomfraction}{0.7}',
        B + 'renewcommand{' + B + 'textfraction}{0.08}' + B + 'renewcommand{' + B + 'floatpagefraction}{0.75}',
        B + 'raggedbottom',
        B + 'begin{document}',
        B + 'doublespacing',
        B + 'begin{center}',
        '{' + B + 'Large ' + B + 'bfseries ' + TITLE + '}' + B + B + '[10pt]',
        author_block(),
        B + 'end{center}',
        B + 'bigskip',
        B + 'section*{Abstract}',
    ]
    for k, v in ABSTRACT.items():
        pre.append(B + 'textbf{%s:} %s' % (k, v) + NL)
    pre.append(B + 'noindent' + B + 'textbf{Keywords:} ' + '; '.join(cfg.PROJECT['keywords']) + NL)

    parts = pre + [
        background,
        B + 'section{Data description}' + B + 'label{sec:data}' + NL,
        m_data,
        results,
        discussion,
        conclusions,
        m_head, m_rest,
        B + 'section*{Availability of source code and requirements}' + NL + B + 'begin{itemize}' + NL
        + src_list + NL + B + 'end{itemize}' + NL,
        B + 'section*{Data availability}' + avail_body,
        abbrev.replace(B + 'section*{List of abbreviations}', B + 'section*{Declarations}' + NL + B + 'bmhead{List of abbreviations}'),
        decl.replace(B + 'section*{Declarations}' + NL, ''),
        B + 'bmhead{' + AI_SECTION_TITLE + '}' + NL + NL + ai_text
        + ' The tools are cited as ' + B + 'cite{claudecode2026,chatgpt2026}.' + NL,
        B + 'section*{Additional files}' + NL + NL + addfiles.replace(B + 'bmhead{Supplementary information}' + NL, ''),
        B + 'bibliographystyle{bib-vancouver}',
        B + 'bibliography{dkd-references,giga-extra}',
        B + 'end{document}',
    ]
    out = NL.join(parts)
    out = out.replace(B + 'backmatter', '')
    return out, s


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
    log('  컴파일: 오류 %d, 미해결 참조 %d %s, 넘침 >10pt %d (최대 %.1fpt)' % (
        len(errs), len(undef), undef[:6], sum(1 for o in over if o > 10), max(over or [0])))
    if errs:
        log('   ' + errs[0][:200])
    return ok and not errs and not undef and not [o for o in over if o > 10]


def main():
    os.makedirs(DST, exist_ok=True)
    os.makedirs(os.path.join(DST, 'figures'), exist_ok=True)
    tex, bmc = build_tex()
    open(os.path.join(DST, NAME + '.tex'), 'w', encoding='utf-8').write(tex)

    # 초록 검사: 길이와 수치 출처
    n = abstract_words()
    nums = numbers_in(' '.join(ABSTRACT.values())) - numbers_in(bmc)
    log('  초록 %d 단어 (한도 250) %s; BMC 본문에 없는 수치: %s' % (n, 'OK' if n <= 250 else '*** 초과', sorted(nums) or '없음'))
    ok = n <= 250 and not nums
    # 코드 요구사항 목록의 판 번호는 requirements.txt 와 같아야 한다
    req = open('requirements.txt', encoding='utf-8').read()
    for pkg, ver in re.findall(r'([a-z-]+) (\d+\.\d+\.\d+)', SOURCE_CODE[4][1]):
        if '%s==%s' % (pkg, ver) not in req:
            log('  *** requirements.txt 와 다름: %s %s' % (pkg, ver))
            ok = False

    for i in range(1, N_FIG + 1):
        shutil.copy2(os.path.join(SRC, 'Fig%d.pdf' % i), os.path.join(DST, 'figures', 'Fig%d.pdf' % i))
    shutil.copy2(os.path.join(SRC, 'dkd-references.bib'), os.path.join(DST, 'dkd-references.bib'))
    shutil.copy2(os.path.join('submission_bib', 'bib-vancouver.bst'), os.path.join(DST, 'bib-vancouver.bst'))
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
    pages = re.search(r'Pages:\s+(\d+)', subprocess.run(['pdfinfo', os.path.join(DST, NAME + '.pdf')],
                                                      capture_output=True, text=True).stdout or '')
    log('  %s/%s.pdf %s쪽' % (DST, NAME, pages.group(1) if pages else '?'))
    log('GigaScience 판 %s' % ('OK' if ok else '*** 실패'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
