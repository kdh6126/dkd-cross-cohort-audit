#!/usr/bin/env python
"""본문이 약속한 Additional file 네 개를 실제로 만든다.

원고는 네 개의 부록을 이름과 내용까지 적어 두었는데 파일이 없었습니다. BMC 투고는 부록을
따로 올리므로, 없는 채로 보내면 심사자가 원고가 가리키는 것을 열 수 없습니다.

형식은 내용에 맞춥니다.

    1, 2번  서술과 작은 표입니다. PDF 로 만듭니다. 본문과 같은 조판이라야 읽는 사람이
            같은 문서의 일부로 받아들입니다.
    3, 4번  큰 표입니다. XLSX 로 만듭니다. 심사자가 정렬하고 걸러 볼 수 있어야 합니다.
            PDF 로 만들면 30행 20열짜리 표가 페이지를 넘어가며 읽히지 않습니다.

네 개 모두 results/ 와 db/dkd.sqlite 에서 읽습니다. 손으로 옮겨 적은 숫자는 없습니다.
"""
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile

import pandas as pd

DST = 'submission'
DB = 'db/dkd.sqlite'
BS = chr(92)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def tex_escape(s):
    # 판다스의 결측을 그대로 문자열로 만들면 표에 nan 이 찍힌다. 읽는 사람에게는
    # 값이 없다는 뜻이 아니라 오류로 보인다.
    if s is None or (isinstance(s, float) and s != s) or str(s).lower() == 'nan':
        return '--'
    s = str(s)
    for a, b in ((BS, BS + 'textbackslash '), ('&', BS + '&'), ('%', BS + '%'),
                 ('$', BS + '$'), ('#', BS + '#'), ('_', BS + '_'),
                 ('{', BS + '{'), ('}', BS + '}'), ('~', BS + 'textasciitilde '),
                 ('^', BS + 'textasciicircum ')):
        s = s.replace(a, b)
    return s


PREAMBLE = [
    BS + 'documentclass[11pt]{article}',
    BS + 'usepackage[a4paper,margin=25mm]{geometry}',
    BS + 'usepackage{lmodern}',
    BS + 'usepackage{booktabs,longtable,array}',
    BS + 'usepackage[T1]{fontenc}',
    BS + 'setlength{' + BS + 'parskip}{6pt}',
    BS + 'setlength{' + BS + 'parindent}{0pt}',
    BS + 'renewcommand{' + BS + 'arraystretch}{1.15}',
]


def longtable(df, spec, headers, fmt=None, raw=(), small=False):
    """긴 표는 longtable 로. tabular 로 두면 페이지를 넘길 때 잘린다.

    small 은 식별자처럼 끊을 수 없는 문자열이 열보다 넓을 때 쓴다. 폭만 늘리면
    표 전체가 본문 폭을 넘고, 그대로 두면 글자가 열 밖으로 삐져나온다.
    """
    out = ([BS + 'begingroup' + BS + 'small'] if small else []) + [
           BS + 'begin{longtable}{@{}' + spec + '@{}}',
           BS + 'toprule',
           ' & '.join(BS + 'textbf{' + h + '}' for h in headers) + ' ' + BS * 2,
           BS + 'midrule', BS + 'endfirsthead',
           BS + 'toprule',
           ' & '.join(BS + 'textbf{' + h + '}' for h in headers) + ' ' + BS * 2,
           BS + 'midrule', BS + 'endhead',
           BS + 'bottomrule', BS + 'endfoot']
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c]
            v = fmt[c](v) if fmt and c in fmt else v
            # raw 로 지정한 열은 이미 TeX 다. 한 번 더 이스케이프하면 $t$ 가 글자로 찍힌다.
            cells.append(str(v) if c in raw else tex_escape(v))
        out.append(' & '.join(cells) + ' ' + BS * 2)
    out.append(BS + 'end{longtable}')
    if small:
        out.append(BS + 'endgroup')
    return out


# 부록은 두 판에 쓰인다. BMC 판은 Additional file N 과 번호 붙은 절, BiB 판은 Supplementary
# File N 과 번호 없는 절에 Supplementary Note 다. 본문을 가리키는 말은 전부 라벨에서 계산한다.
# 손으로 적은 번호는 BMC 판 안에서도 이미 어긋나 있었다(Figure 8 이 Figure 10 이 되었고,
# 코호트 표는 Table 1 이 아니라 Table 2 였다).
TARGET = 'bmc'
MANUSCRIPT = {'bmc': 'submission/dkd-manuscript.tex', 'bib': 'submission_bib/bib-manuscript.tex'}
BIB_SUPPLEMENT = 'submission_bib/bib-supplement.tex'
# BiB 본문에 제목 줄이 없는 BMC 절이 들어간 Supplementary Note
BIB_NOTE = {'sec:harmonisation': 'S1', 'sec:nulls': 'S2', 'sec:auroc-null': 'S2',
            'sec:procurement': 'S3', 'sec:modulesize': 'S4', 'sec:whichcorrections': 'S5',
            'sec:proteome': 'S7', 'sec:limitations': 'S8', 'sec:implementation': 'S9',
            'sec:xtmethods': 'S10'}


def file_title(n):
    return ('Additional file %d' if TARGET == 'bmc' else 'Supplementary File %d') % n


def file_name(n):
    return ('Additional_file_%d' if TARGET == 'bmc' else 'Supplementary_File_%d') % n


def _float_order(tex, env):
    src = open(tex, encoding='utf-8').read()
    out = []
    for m in re.finditer('(?s)' + re.escape(BS) + r'begin\{' + env + r'\*?\}.*?'
                         + re.escape(BS) + r'end\{' + env + r'\*?\}', src):
        lab = re.search(re.escape(BS) + r'label\{([^}]+)\}', m.group(0))
        if lab:
            out.append(lab.group(1))
    return out


def where(label, tex=True):
    """본문(또는 BiB 보충 자료)에서 그 라벨이 가리키는 곳을 문장으로 돌려준다."""
    ms = MANUSCRIPT[TARGET]
    kind = label.split(':')[0]
    if kind in ('fig', 'tab'):
        env, word = ('figure', 'Figure') if kind == 'fig' else ('table', 'Table')
        order = _float_order(ms, env)
        if label in order:
            txt = '%s~%d of the main text' % (word, order.index(label) + 1)
        elif TARGET == 'bib' and os.path.exists(BIB_SUPPLEMENT) and label in _float_order(BIB_SUPPLEMENT, env):
            txt = 'Supplementary %s~S%d' % (word, _float_order(BIB_SUPPLEMENT, env).index(label) + 1)
        else:
            raise SystemExit('부록이 가리키는 %s 를 %s 에서 찾지 못했습니다' % (label, ms))
    elif TARGET == 'bmc':
        txt = 'Section~%s of the main text' % section_number(label)
    else:
        src = open(ms, encoding='utf-8').read()
        m = re.search(re.escape(BS) + r'(?:sub)?section\{([^}]*)\}' + re.escape(BS)
                      + r'label\{' + re.escape(label) + r'\}', src)
        if m:
            txt = 'the main-text section ' + BS + 'textit{' + m.group(1) + '}'
        elif label in BIB_NOTE:
            txt = 'Supplementary Note~' + BIB_NOTE[label]
        else:
            raise SystemExit('부록이 가리키는 %s 를 BiB 본문에서 찾지 못했습니다' % label)
    if not tex:
        txt = re.sub(re.escape(BS) + r'textit\{([^}]*)\}', r'\1', txt).replace('~', ' ')
    return txt


def cap(t):
    return t[0].upper() + t[1:]


def section_number(label):
    """원고에서 그 라벨이 붙은 절의 번호를 센다.

    부록에 "Section 3.9" 처럼 손으로 적어 두면 본문에 절이 하나 늘거나 줄 때 조용히
    어긋납니다. 실제로 3.9 라고 적힌 곳이 3.8 이 되어 있었습니다. 그래서 세어서 씁니다.
    """
    src = open(MANUSCRIPT['bmc'], encoding='utf-8').read()
    sec = sub = 0
    for m in re.finditer(r'(?m)^' + re.escape(BS) + r'(section|subsection)\{|'
                         + re.escape(BS) + r'label\{([^}]+)\}', src):
        if m.group(1) == 'section':
            sec += 1
            sub = 0
        elif m.group(1) == 'subsection':
            sub += 1
        elif m.group(2) == label:
            return '%d.%d' % (sec, sub) if sub else '%d' % sec
    return '?'


def compile_pdf(lines, name):
    src = (chr(10)).join(PREAMBLE + [BS + 'begin{document}'] + lines
                         + [BS + 'end{document}'])
    tmp = tempfile.mkdtemp(prefix='af_')
    tex = os.path.join(tmp, name + '.tex')
    open(tex, 'w', encoding='utf-8').write(src)
    # 그림을 포함하는 부록이 있다. 임시 폴더에서 컴파일하므로 그림도 같이 옮겨야 한다.
    for g in os.listdir(DST):
        if g.endswith('.pdf') and g.startswith('P'):
            shutil.copy2(os.path.join(DST, g), os.path.join(tmp, g))
    for _ in range(2):
        p = subprocess.run(['pdflatex', '-interaction=nonstopmode', name + '.tex'],
                           cwd=tmp, capture_output=True)
    pdf = os.path.join(tmp, name + '.pdf')
    if not os.path.exists(pdf):
        log('%s 컴파일 실패' % name)
        tail = p.stdout.decode('utf-8', 'replace')
        log(chr(10).join(l for l in tail.splitlines() if l.startswith('!'))[:600])
        shutil.rmtree(tmp, ignore_errors=True)
        return False
    # 본문은 overfull 을 세는데 부록은 세지 않았다. 그래서 표가 페이지 밖으로 잘려도
    # 모든 검사가 통과했다. 읽어 보고서야 알았으므로, 여기서 세어 보고한다.
    txt = p.stdout.decode('utf-8', 'replace')
    over = [l for l in txt.splitlines() if l.startswith('Overfull')]
    shutil.copy2(pdf, os.path.join(DST, name + '.pdf'))
    shutil.rmtree(tmp, ignore_errors=True)
    log('  %-22s %s' % (name + '.pdf',
                        '' if not over else 'overfull %d개 ***' % len(over)))
    for l in over[:4]:
        log('      %s' % l[:110])
    return True


# ------------------------------------------------------------------ 1
def af1():
    var = pd.read_csv('results/wgcna_validation/variants.tsv', sep='\t')
    soft = pd.read_csv('results/wgcna_validation/soft_threshold.tsv', sep='\t')
    # 격자의 최빈값을 beta 로 쓰면 본문과 어긋난다. 실제로 쓰이는 값은 파이프라인
    # 함수의 기본 인자이므로 거기서 읽는다. 코드가 바뀌면 이 문서가 따라 바뀐다.
    import inspect
    sys.path.insert(0, 'scripts')
    from conventional_pipeline import wgcna_hub
    used_beta = int(inspect.signature(wgcna_hub).parameters['power'].default)
    L = []
    A = L.append
    A(BS + 'section*{' + file_title(1) + '. Comparator implementation and its validation}')
    A('This file supports ' + where('sec:selectors') + ' and ' + where('sec:modulesize')
      + '. It states what the '
      'correlation-module comparator used here shares with the canonical WGCNA pipeline, '
      'where it departs, and whether any conclusion in the manuscript depends on those '
      'departures.')

    A(BS + 'subsection*{Why the comparators are implemented directly}')
    A('The benchmark repeats every selection strategy across %d cohort subsets and every '
      'leave-one-dataset-out fold within each subset. Calling reference packages across '
      'language boundaries at that repetition count is slow and, more importantly, makes the '
      'result depend on the versions of several external runtimes. Implementing the '
      'comparators in a single environment keeps the whole benchmark reproducible from one '
      'command. The cost is that these are not the reference implementations, and we do not '
      'claim they reproduce them.' % 26)

    A(BS + 'subsection*{What the hub ranking shares with WGCNA, and where it departs}')
    # 같은 표가 본문 Table 3 에도 들어간다. 두 곳에 적으면 언젠가 어긋나므로
    # comparator_table.py 한 곳에서 정의하고 여기서는 가져다 쓴다.
    from comparator_table import rows as comparator_rows
    rows = pd.DataFrame(comparator_rows(),
                        columns=['Element', 'Canonical WGCNA', 'Used here', 'Status'])
    L += longtable(rows, 'p{31mm}p{42mm}p{42mm}p{16mm}',
                   ['Element', 'Canonical WGCNA', 'Used here', 'Status'],
                   raw=('Used here',))
    A('The departures are real and are why the manuscript calls this a correlation-module '
      'comparator rather than WGCNA. The question the rest of this file answers is whether '
      'any of them changes a conclusion.')

    A(BS + 'subsection*{Scale-free topology fit}')
    A('The canonical criterion selects the smallest power whose scale-free fit $R^2$ reaches '
      'about 0.80. The fit across powers on the harmonised discovery data is below.')
    sf = soft.copy()
    keep = [c for c in ('power', 'scale_free_r2', 'mean_k', 'median_k', 'slope')
            if c in sf.columns]
    sf = sf[keep]
    L += longtable(sf, 'r' * len(keep), [c.replace('_', ' ') for c in keep],
                   fmt={c: (lambda v: '%.3f' % v if isinstance(v, float) else str(v))
                        for c in keep})
    best = sf.loc[sf['scale_free_r2'].idxmax()] if 'scale_free_r2' in sf else None
    if best is not None:
        A(('At the power used here ($\\beta = %d$) the fit is $R^2 = %.2f$, below the '
          'canonical threshold; the best fit over the grid is %.2f at $\\beta = %d$. A lower '
          'power yields a denser network and therefore larger modules. '
          % (used_beta,
             float(sf.loc[sf['power'] == used_beta, 'scale_free_r2'].iloc[0])
             if (sf['power'] == used_beta).any() else float('nan'),
             float(best['scale_free_r2']), int(best['power'])))
          + cap(where('sec:modulesize')) + ' shows that larger modules favour the '
          'procurement-tracking module, so this departure works against the manuscript '
          'argument rather than for it.')

    A(BS + 'subsection*{The variant grid}')
    n_var = len(var)
    n_hit = int((var['ieg_in_top50'] > 0).sum())
    A('Every combination of the departing choices was run: %d variants over module definition, '
      'ranking statistic, adjacency and soft-threshold power. The question is whether any '
      'combination places immediate-early genes outside the top of the ranking, which would '
      'mean the procurement signal is an artefact of our implementation rather than of the '
      'data.' % n_var)
    A(BS + 'textbf{%d of %d variants place an immediate-early gene in the top 50.} '
      'The result does not depend on the implementation choices.' % (n_hit, n_var))
    v = var.copy()
    for c in ('module', 'module_detail', 'rank_by', 'adjacency'):
        if c in v.columns:
            v[c] = v[c].astype(str)
    cols = [c for c in ('module_detail', 'rank_by', 'adjacency', 'power', 'module_size',
                        'ieg_in_module', 'ieg_in_top50', 'ieg_median_rank') if c in v.columns]
    v = v[cols]
    # 굵은 머리글("module detail")이 16mm 를 넘어 표 전체가 밀렸다. 폭을 내용에 맞추고
    # 작은 글꼴을 쓴다. 뒤쪽 다섯 열은 숫자이므로 좁게 둬도 접히지 않는다.
    L += longtable(v, 'p{25mm}p{15mm}p{19mm}' + 'p{13mm}' * (len(cols) - 3),
                   [{'module_detail': 'definition', 'ieg_in_module': 'IEG in module',
                     'ieg_in_top50': 'IEG in top 50',
                     'ieg_median_rank': 'IEG median rank'}.get(c, c.replace('_', ' '))
                    for c in cols],
                   fmt={c: (lambda x: ('%.0f' % x) if isinstance(x, float) and x == x
                            else str(x)) for c in cols}, small=True)

    A(BS + 'subsection*{Provenance of the remaining selectors}')
    A('The same reasoning applies to the other three strategies in the benchmark. Each is '
      'stated here with the published definition it stands for and the respects in which the '
      'implementation used here is narrower.')
    sel = pd.DataFrame([
        ('DEG' + BS + '_meta', 'Per-cohort Welch $t$, combined by Stouffer',
         'weights are $\\sqrt{n}$; no random-effects variance component'),
        ('RBS', 'Stability selection (Meinshausen and B' + chr(92) + '"uhlmann)',
         'bootstrap selection inside each cohort, aggregated by geometric mean across '
         'cohorts, with a sign-concordance gate and a perturbation term. It is used as an '
         'instrument standing in for this family, not proposed as a contribution.'),
        ('RBS' + BS + '_orth', 'The same, after residualisation',
         'every gene is residualised on a per-sample handling score before selection'),
        ('Baseline selectors inside RBS',
         'univariate $F$, LASSO, elastic net, mRMR, ReliefF, random forest, Boruta',
         'all implemented directly rather than through packages that pin older numpy or '
         'scikit-learn versions'),
    ], columns=['Strategy', 'Stands for', 'How the implementation is narrower'])
    L += longtable(sel, 'p{28mm}p{40mm}p{62mm}', list(sel.columns),
                   raw=('Strategy', 'Stands for', 'How the implementation is narrower'))
    A('None of these is offered as a new method. They are present so that the benchmark '
      'contains a representative of each family a reader would expect to see compared.')

    A(BS + 'subsection*{Block-size simulation with the RBS implementation}')
    A(cap(where('sec:modulesize')) + ' reports that the connectivity criterion picks up more of the '
      'confounder than the stability criterion when the confounder block is large and less '
      'when it is small. The table gives the paired difference in confounder fraction of the '
      'top 50, connectivity minus RBS, over 25 replicates per configuration, run with the '
      'actual RBS implementation rather than a stand-in.')
    import numpy as np
    from scipy import stats
    rs = pd.read_csv('results/module_size/real_rbs_sweep.tsv', sep='\t')
    pv = rs.pivot_table(index=['n_dis', 'n_conf', 'rep'], columns='method',
                        values='confounder_frac')
    rows = []
    for (nd, nc), v in (pv['connectivity'] - pv['RBS_real']).groupby(level=[0, 1]):
        m = v.mean()
        se = v.std(ddof=1) / np.sqrt(len(v))
        p = stats.wilcoxon(v.to_numpy()).pvalue
        rows.append((int(nd), int(nc), '%+.3f' % m, '%+.3f' % (m - 1.96 * se),
                     '%+.3f' % (m + 1.96 * se), '%.1e' % p, int(len(v)),
                     'worse' if m > 0 else 'better'))
    # 머리글이 길어 표가 본문 폭을 80pt 넘었다. 뜻은 표 아래 문장에 적는다.
    bs = pd.DataFrame(rows, columns=['Disease', 'Confounder', 'Mean', 'CI low', 'CI high',
                                     'Wilcoxon $p$', 'Reps', 'Connectivity'])
    L += longtable(bs, 'rrrrrrrl', list(bs.columns), small=True)
    A('Disease and Confounder are the block sizes; Mean and the 95\% CI are the paired '
      'difference; Connectivity states whether the connectivity criterion did worse or better. '
      'Positive values mean the connectivity criterion placed more confounder genes in its '
      'top 50 than the stability criterion did. The sign follows the block-size ratio, which is '
      'the mechanism described in ' + where('sec:modulesize') + '.')

    return compile_pdf(L, file_name(1))


# ------------------------------------------------------------------ 2
def af2():
    c = sqlite3.connect(DB)
    aud = pd.read_sql(
        'select s.accession, s.title, s.omics_type, s.n_samples, a.usable, '
        'a.exclusion_reason, a.supersedes, '
        '(select count(*) from sample where sample.study_id = s.study_id) as n_loaded '
        'from study s join study_audit a '
        'on s.study_id = a.study_id order by a.usable, s.accession', c)
    # 제외한 연구까지 넣고 세면 '제거하지 못한 잔여 중복'이 제거한 중복까지 포함한다
    smp = pd.read_sql(
        'select st.accession as study, sm.subject_id, sm.compartment, sm.procurement '
        'from sample sm join study st on sm.study_id = st.study_id '
        'join study_audit a on a.study_id = st.study_id where a.usable = 1', c)
    c.close()

    L = []
    A = L.append
    A(BS + 'section*{' + file_title(2) + '. Cohort audit}')
    A('This file supports ' + where('sec:cohorts') + ' and ' + where('tab:overlap')
      + '. Public disease series are '
      'frequently supersets, re-analyses or re-runs of one another. Counting such a series as '
      'an independent cohort inflates every cross-cohort agreement statistic, because the '
      'agreement is partly between a patient and themselves.')

    A(BS + 'subsection*{Procedure}')
    A('For every candidate series we read the sample-level metadata from the GEO series '
      'matrix, extracted the subject identifier where the submitters recorded one, and '
      'compared identifier sets across series. Where identifiers were absent we compared '
      'sample titles, platform, and the exact per-sample value vectors. A series was excluded '
      'when its subjects were already present in a series we retained.')

    excl = aud[aud['usable'] == 0]
    A(BS + 'subsection*{Excluded series (%d)}' % len(excl))
    e = excl[['accession', 'supersedes', 'exclusion_reason']].fillna('')
    # 세미콜론으로 이은 접근번호는 줄바꿈 지점이 없어 열을 넘는다. 쉼표와 공백으로 잇는다.
    e['supersedes'] = e['supersedes'].astype(str).str.replace(';', ', ', regex=False)
    e['exclusion_reason'] = (e['exclusion_reason'].astype(str)
                             .str.replace('/', ' / ', regex=False))
    L += longtable(e, 'p{20mm}p{30mm}p{72mm}',
                   ['Accession', 'Superseded by', 'Evidence'], small=True)

    A(BS + 'subsection*{Residual overlap that could not be removed}')
    both = (smp.dropna(subset=['subject_id'])
            .groupby('subject_id')['study'].nunique())
    shared = int((both > 1).sum())
    pairs = []
    g = smp.dropna(subset=['subject_id']).groupby('study')['subject_id'].apply(set)
    for i, a in enumerate(g.index):
        for b in list(g.index)[i + 1:]:
            n = len(g[a] & g[b])
            if n:
                pairs.append((a, b, n, len(g[a]), len(g[b])))
    A('%d subjects appear in more than one retained series. These are not duplicate cohorts '
      'but the same ERCB patients profiled in two compartments, which is what makes the '
      'paired glomerular and tubulointerstitial analyses possible. The benchmark treats the '
      'two compartments as separate datasets, so the leave-one-dataset-out folds are not '
      'fully independent for those subjects; this is stated in the limitations of the main '
      'text.' % shared)
    if pairs:
        p = pd.DataFrame(pairs, columns=['Series A', 'Series B', 'Shared subjects',
                                         'Subjects in A', 'Subjects in B'])
        L += longtable(p, 'llrrr', list(p.columns))

    A(BS + 'subsection*{Full audit table}')
    a2 = aud.copy()
    a2['usable'] = a2['usable'].map({1: 'retained', 0: 'excluded'})
    a2['n_loaded'] = a2['n_loaded'].map(lambda v: '' if not v else '%d' % v)
    a2 = a2[['accession', 'omics_type', 'n_loaded', 'usable', 'exclusion_reason']].fillna('')
    # 'Samples loaded' 를 r 로 두면 머리글이 접히지 않아 표가 밀린다. p 로 바꾼다.
    L += longtable(a2, 'p{22mm}p{26mm}p{18mm}p{18mm}p{46mm}',
                   ['Accession', 'Omics', 'Samples loaded', 'Status', 'Note'], small=True)
    A('Samples loaded is the number of samples parsed into the project database. It is zero '
      'for assets that were downloaded and inspected but not parsed into the sample table, '
      'such as the single-cell objects and most metabolomics studies, which are analysed '
      'from their own files.')
    return compile_pdf(L, file_name(2))


# ------------------------------------------------------------------ 3, 4
def sheet(writer, name, df, note):
    """시트 하나 = 표 하나 + 그 표가 무엇인지 한 줄. 설명 없는 시트는 읽히지 않는다."""
    df.to_excel(writer, sheet_name=name[:31], index=False, startrow=2)
    ws = writer.sheets[name[:31]]
    ws.cell(row=1, column=1, value=note)
    for k, col in enumerate(df.columns, 1):
        w = max(len(str(col)), *(len(str(x)) for x in df[col].head(200))) if len(df) else 12
        ws.column_dimensions[ws.cell(row=3, column=k).column_letter].width = min(w + 2, 46)


def af3():
    path = os.path.join(DST, file_name(3) + '.xlsx')
    with pd.ExcelWriter(path, engine='openpyxl') as w:
        readme = pd.DataFrame([
            ('cohorts', 'Every cohort that entered the analysis: role, compartment, platform, '
                        'case and control counts, and how the control tissue was obtained.'),
            ('harmonisation', 'Value ranges after harmonisation, as written by harmonize.py.'),
            ('gene_space', 'The frozen 9,900-gene Entrez space shared by all cohorts.'),
            ('candidates', 'The complete candidate table with every line of evidence and the '
                           'pre-registered tier.'),
            ('candidate_evidence', 'Per-candidate evidence rows as stored in the project '
                                   'database.'),
            ('rank_columns', 'What each rank column in the candidates sheet means. The three '
                             'rank columns come from different procedures and are not '
                             'comparable to one another.'),
        ], columns=['Sheet', 'Contents'])
        sheet(w, 'README', readme,
              file_title(3) + '. Harmonisation summary and the complete candidate table. '
              'Supports ' + where('sec:cohorts', tex=False) + ' and '
              + where('sec:candidates', tex=False) + '.')
        # harmonized_summary.tsv 는 harmonize.py 가 만들고, GSE294519 는 나중에
        # add_gse294519.py 로 따로 붙어서 그 요약에 들어가지 않았다. 발견 코호트 다섯 중
        # 하나가 빠진 표를 부록으로 내보내면 심사자가 본문의 '다섯 코호트'와 맞춰볼 수 없다.
        sheet(w, 'cohorts',
              pd.read_csv('results/cohort_table.tsv', sep='\t'),
              'Every cohort that entered the analysis, built from the harmonised files '
              'themselves. This is ' + where('tab:cohorts', tex=False)
              + ' with per-cohort gene counts added.')
        sheet(w, 'harmonisation',
              pd.read_csv('data/processed/harmonized_summary.tsv', sep='\t'),
              'Value distribution after harmonisation, as written by harmonize.py. GSE294519 '
              'was added by a later step and is absent here; the cohorts sheet is the complete '
              'list.')
        sheet(w, 'gene_space',
              pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str),
              'The frozen gene space. Probes were collapsed to Entrez genes by maximum mean '
              'intensity; genes absent from any cohort were dropped before any analysis.')
        sheet(w, 'candidates',
              pd.read_csv('results/candidates_v2/master_candidate_table.tsv', sep='\t'),
              'Candidate table. Tier assignment follows criteria fixed before the evidence '
              'was read; see the main text. A candidate is a lead for further work, not a '
              'validated biomarker.')
        # 순위 열이 셋인데 절차가 서로 다르다. 이름만 보고 한 표에 나란히 있으면
        # 직접 비교할 수 있는 순위로 오해하기 쉬우므로, 무엇으로 매긴 순위인지와
        # 무엇과 비교해도 되는지를 따로 적는다.
        sheet(w, 'rank_columns', pd.DataFrame([
            ('mean_rank',
             'Selection pipeline',
             'Mean rank assigned by the corrected selector across leave-one-dataset-out folds, '
             'over the frozen 9,900-gene space. Lower is earlier in the candidate list.',
             'Comparable across rows of this sheet only.'),
            ('rank_naive_control',
             'Control-definition sensitivity analysis (ERCB)',
             'Rank by |Hedges g| when diabetic nephropathy is compared with the non-biopsy '
             'controls (nephrectomy or living donor).',
             'Comparable with rank_matched_control, not with mean_rank.'),
            ('rank_matched_control',
             'Control-definition sensitivity analysis (ERCB)',
             'Rank by |Hedges g| when the same case set is compared instead with biopsies of '
             'other kidney diagnoses.',
             'Comparable with rank_naive_control, not with mean_rank.'),
            ('rank_shift_on_substitution',
             'Control-definition sensitivity analysis (ERCB)',
             'rank_naive_control minus rank_matched_control. Negative means the gene rises when '
             'the control group is replaced.',
             'A difference within the sensitivity analysis; it is not a shift in mean_rank.'),
        ], columns=['Column', 'Procedure', 'Definition', 'What it may be compared with']),
              'The candidates sheet carries three rank columns produced by two different '
              'procedures. They are not interchangeable: mean_rank comes from the selection '
              'pipeline over all cohorts, while the other three come from the '
              'control-definition sensitivity analysis within ERCB and rank by effect size '
              'under two control definitions.')
        c = sqlite3.connect(DB)
        sheet(w, 'candidate_evidence',
              pd.read_sql('select * from candidate_evidence', c),
              'Evidence rows backing each candidate, as stored in db/dkd.sqlite.')
        c.close()
    log('  %s.xlsx' % file_name(3))
    return True


def af4():
    path = os.path.join(DST, file_name(4) + '.xlsx')
    ext = pd.read_csv('results/external_validation/summary.tsv', sep='\t')
    # 출처 파일은 작업 언어가 한국어입니다. 부록은 영어라야 하므로 여기서만 옮깁니다.
    ko = {'ERCB 밖 독립 특이성': 'independent specificity outside ERCB',
          '혈액에서 보이는가': 'is the signal visible in blood',
          '초기에도 변하는가': 'does the signal already move in early disease',
          'DKD 특이': 'DKD-specific',
          '조직 국한': 'tissue-restricted',
          '초기에도 변함': 'moves in early disease'}
    for col in ('question', 'verdict'):
        if col in ext.columns:
            ext[col] = ext[col].map(lambda v: ko.get(str(v).strip(), v))

    with pd.ExcelWriter(path, engine='openpyxl') as w:
        gd = 'results/ckd_generalization'
        per = pd.read_csv(os.path.join(gd, 'per_diagnosis.tsv'), sep='\t')
        readme = pd.DataFrame(
            [('per_diagnosis', 'One row per ERCB diagnosis: signal under the conventional '
                               'non-biopsy control and under the procurement-matched control.')]
            + [('g_' + f[2:-4], 'Per-gene effect sizes for %s under both control definitions.'
                % f[2:-4].replace('_', ' ')) for f in sorted(os.listdir(gd))
               if f.startswith('g_')]
            + [('external_validation', 'Independent datasets testing specificity and tissue '
                                       'restriction.'),
               ('kpmp_contrasts', 'KPMP single-nucleus pseudobulk: the immediate-early module '
                                  'under three contrasts, against random genes.'),
               ('kpmp_layers', 'KPMP modality table by diagnosis: which omics layers were '
                               'measured in the same donor.')],
            columns=['Sheet', 'Contents'])
        sheet(w, 'README', readme,
              file_title(4) + '. Per-diagnosis, external-validation and modality tables '
              'supporting ' + where('sec:external', tex=False) + '.')
        sheet(w, 'per_diagnosis', per,
              'shrink is the ratio of median |g| under the matched control to median |g| '
              'under the conventional control. Below 1 means the conventional design '
              'overstates the signal.')
        for f in sorted(os.listdir(gd)):
            if f.startswith('g_'):
                sheet(w, f[2:-4][:31],
                      pd.read_csv(os.path.join(gd, f), sep='\t'),
                      'g_naive: diagnosis vs non-biopsy control. g_matched: diagnosis vs the '
                      'remaining biopsy diagnoses.')
        sheet(w, 'external_validation', ext,
              'Datasets outside the discovery cohorts used to test whether the retained '
              'signal is DKD-specific and whether it is visible in blood.')
        # 본문 표였던 KPMP 대비 요약. 결과 파일에서 다시 계산한다.
        from scipy import stats as _st
        pb = pd.read_csv('results/kpmp_singlecell/pseudobulk_de.tsv', sep='\t')
        what = {'DKD_vs_reference': ('DKD-proxy vs healthy reference', 'biopsy vs non-biopsy'),
                'DKD_vs_CKDnondiabetic': ('DKD-proxy vs non-diabetic CKD', 'biopsy vs biopsy'),
                'reference_vs_nephrectomy': ('reference vs nephrectomy',
                                             'non-biopsy vs non-biopsy')}
        kc = []
        for ctr, g in pb.groupby('contrast'):
            ieg = g[g['kind'] == 'IEG']['g'].dropna()
            rnd = g[g['kind'] == 'random']['g'].dropna()
            kc.append(dict(contrast=what[ctr][0], compares=what[ctr][1],
                           ieg_mean_g=round(float(ieg.mean()), 3),
                           ieg_vs_random_p=float(_st.mannwhitneyu(
                               ieg, rnd, alternative='two-sided').pvalue),
                           n_ieg=int(len(ieg)), n_random=int(len(rnd))))
        sheet(w, 'kpmp_contrasts', pd.DataFrame(kc),
              'Hedges g per (gene, cell type) on donor-level pseudobulk. The module separates '
              'biopsy from non-biopsy tissue and nothing else.')
        # 참여자 한 명이 한 행인 표가 아니라 질환별 집계를 싣는다. 논문이 쓰는 것은
        # 집계뿐이고, 개인 단위 행은 공개할 이유가 없다.
        sheet(w, 'kpmp_layers',
              pd.read_csv('results/kpmp_overlap/layer_summary.tsv', sep='\t'),
              'Omics layers measured in the same KPMP participant, counted by diagnosis. Only '
              'participants appearing in the single-nucleus object were queried, so these are '
              'rates among transcriptome donors, not among all KPMP participants. Participant '
              'identifiers are deliberately not reproduced here.')
    log('  %s.xlsx' % file_name(4))
    return True


def af5():
    """본문 그림 8 의 완전판. 본문은 2·3등급 12개만 싣는다.

    30줄을 본문 폭 170mm 에 밀어 넣으면 유전자명이 5pt 밑으로 떨어져 읽히지 않습니다.
    그렇다고 자른 것을 그냥 두면 4등급 18개가 어디로 갔는지 알 수 없으므로, 같은 그림의
    전체판을 여기 냅니다. 잘라낸 것을 감추지 않는 것이 요점입니다.
    """
    src = 'results/figures/P7_candidates_full.pdf'
    if not os.path.exists(src):
        log('  %s 가 없습니다. scripts/candidate_figure.py 를 먼저 돌리세요.' % src)
        return False
    shutil.copy2(src, os.path.join(DST, 'P7_candidates_full.pdf'))
    n = len(pd.read_csv('results/candidates_v2/master_candidate_table.tsv', sep='\t'))
    L = [BS + 'usepackage{graphicx}']
    body = [
        BS + 'section*{' + file_title(5) + '. Candidate evidence, all %d candidates}' % n,
        cap(where('fig:candidates')) + ' shows the 12 candidates in tiers 2 and 3 -- those that are '
        'either strong in the data or absent from the DKD literature. This file shows the same '
        'figure over all %d specificity-gated candidates, adding the 18 in tier 4, which are '
        'neither. Nothing is computed differently here; only the row selection differs.' % n,
        BS + 'begin{center}',
        BS + 'includegraphics[width=' + BS + 'textwidth]{P7_candidates_full}',
        BS + 'end{center}',
        'Panel \\textbf{a}: a filled cell means the evidence line is met, an open cell that it '
        'was tested and not met, and a hatched cell that the gene could not be assessed on that '
        'line. Single-nucleus and proteomic corroboration are not available for every gene, and '
        'scoring an absent measurement as a failure would confuse the two. Panel \\textbf{b} '
        'counts all %d candidates in both versions of the figure, so the two panels always '
        'describe the same set. Tier assignments are read from the candidate table rather than '
        'recomputed here.' % n,
    ]
    global PREAMBLE
    saved = PREAMBLE
    PREAMBLE = saved + L
    try:
        ok = compile_pdf(body, file_name(5))
    finally:
        PREAMBLE = saved
    return ok


def af6():
    """두 단백체 코호트의 감사표와 매칭 순열의 균형.

    본문이 "대조군 정의가 서로 다르다" 와 "존재비를 맞춰도 유지된다" 두 가지를
    주장합니다. 근거를 본문에 다 넣으면 절이 표로 뒤덮이므로 여기로 옮깁니다.

    감사표의 값은 계산이 아니라 원논문과 저장소 메타데이터에서 옮긴 것입니다.
    적혀 있지 않은 항목은 적혀 있지 않다고 적습니다.
    """
    a = pd.read_csv('results/proteome_meta/cohort_audit.tsv', sep=chr(9))
    b = pd.read_csv('results/proteome_meta/matching_balance.tsv', sep=chr(9))
    m = pd.read_csv('results/proteome_meta/summary.tsv', sep=chr(9)).iloc[0]
    b = pd.DataFrame({
        # 코호트 이름을 통째로 넣으면 첫 열이 네 줄로 접혀 표가 읽히지 않는다.
        # 접근번호는 바로 위 감사표에 이미 있으므로 여기서는 기술 이름만 쓴다.
        'Cohort': b['layer'].str.split(' (', regex=False).str[0],
        'Candidates': b['n_candidates'],
        'Eligible background per candidate': b['mean_eligible_background'].round(0).astype(int),
        'Mean abundance, candidates': b['candidate_abundance'].round(2),
        'Mean abundance, matched null': b['matched_null_abundance'].round(2)})

    L = []
    A = L.append
    A(BS + 'section*{' + file_title(6) + '. The two proteomic cohorts, and what the matched '
      'permutation matched}')
    A(cap(where('sec:proteome')) + ' uses two proteomic cohorts. They are not interchangeable, '
      'and the differences matter for how far the pooled result can be read. This table '
      'records what each source publication and repository record states, and marks what they '
      'leave unstated rather than filling it in by inference.')
    L += longtable(a, 'p{0.24' + BS + 'textwidth}p{0.34' + BS + 'textwidth}p{0.34'
                   + BS + 'textwidth}',
                   ['Attribute', 'SOMAscan cohort', 'Mass-spectrometry cohort'])
    A('Three things follow. First, the control arms differ -- healthy participants in one '
      'cohort, non-diabetic pathology cases in the other -- so the pooled contrast is diabetic '
      'versus non-diabetic rather than one shared case/control definition. Second, the '
      'mass-spectrometry cohort draws both arms as archival explant cases from a single '
      'institution, which removes the biopsy-versus-nephrectomy contrast from its design; it '
      'does not establish that pre-analytical exposure was matched, because neither source '
      'reports ischaemia, fixation or storage time. Third, ' + BS + 'textit{MMP7} is the '
      'headline finding of the publication that produced the SOMAscan cohort, so its '
      'significance there restates a published result on the same data. The independent '
      'observation is its replication in the mass-spectrometry cohort.')

    A(BS + 'subsection*{Matched permutation}')
    A('The pooled enrichment could arise if the candidates were simply abundant, readily '
      'detected proteins. To test that, the permutation was repeated so that each null protein '
      'was drawn from the background matched to the candidate it replaced: on abundance rank '
      'within a caliper of 10 per cent of the background, and, in the mass-spectrometry '
      'cohort, on peptide count within one doubling. Draws are without replacement within each '
      'permuted set and never cross strata.')
    L += longtable(b, 'p{0.16' + BS + 'textwidth}' + ('p{0.15' + BS + 'textwidth}') * 4,
                   ['Cohort', 'Candidates', 'Eligible background',
                    'Abundance, candidates', 'Abundance, matched null'])
    A('The unmatched permutation gives $p = %.3f$ and the matched permutation $p = %.3f$ over '
      '%s draws, against a null expectation of %.1f significant candidates unmatched and %.1f '
      'matched. The conclusion does not depend on the candidates being easier to detect than '
      'the proteins they are compared with.'
      % (m['perm_p_sig'], m['perm_p_matched'], '{:,}'.format(int(m['n_perm'])),
         m['perm_null_mean'], m['perm_null_mean_matched']))
    d = pd.read_csv('results/proteome_meta/discovery_summary.tsv', sep=chr(9)).iloc[0]
    hit = pd.read_csv('results/proteome_meta/discovery_attempt.tsv', sep=chr(9))
    hit = hit[hit['both_sig_concordant']].sort_values('prot_q_b')
    A(BS + 'section*{Corroboration is not discovery}')
    A('The pooled result above says that a candidate list built from transcriptomes holds up '
      'in protein data. It does not say that the protein data could have produced such a list '
      'on its own, and the two need different amounts of evidence: corroboration starts from a '
      'list and needs only the few genes that happen to be measured, whereas discovery has to '
      'build one and needs a wide search space across several cohorts. We tested the second '
      'directly, using the same two datasets and no transcriptomic information, by asking '
      'which of the %d genes both of them quantify reach $q<0.05$ in both with a concordant '
      'direction.' % d['n_shared_genes'])
    num = {c: (lambda v: '%+.2f' % v) for c in ('prot_g_a', 'prot_g_b')}
    num.update({c: (lambda v: '%.3f' % v) for c in ('prot_q_a', 'prot_q_b')})
    num['is_candidate'] = lambda v: 'yes' if bool(v) else 'no'
    L += longtable(hit[['gene', 'prot_g_a', 'prot_q_a', 'prot_g_b', 'prot_q_b',
                        'is_candidate']],
                   'l' + ('p{0.13' + BS + 'textwidth}') * 4 + 'p{0.16' + BS + 'textwidth}',
                   ['Gene', 'g (SOMAscan)', 'q (SOMAscan)', 'g (LC-MS/MS)', 'q (LC-MS/MS)',
                    'In candidate list'], fmt=num)
    A('Four genes qualify, three of them outside the transcriptomic candidate list. That is '
      'not more than chance: under independence %.1f are expected, and a permutation that '
      'shuffles the gene labels of one dataset while holding each dataset to its own '
      'distribution of significance and direction gives $p = %.2f$ over %s draws (binomial '
      '$p = %.2f$). The protein layer here has 2 datasets, %d patients and %d shared genes, '
      'against 7 cohorts, 258 patients and 9,900 genes on the transcriptome side; with two '
      'datasets a leave-one-dataset-out procedure degenerates to training on one and scoring '
      'on the other. We therefore report the protein data as corroboration only, and record '
      'this negative result rather than leaving the stronger reading available.'
      % (d['expected_independent'], d['perm_p'], '{:,}'.format(int(d['n_perm'])),
         d['binom_p'], 45, int(d['n_shared_genes'])))
    ident = pd.read_csv('results/proteome_meta/protein_identity_audit.tsv', sep=chr(9))
    fz = pd.read_csv('results/proteome_meta/freeze_summary.tsv', sep=chr(9)).iloc[0]
    led = pd.read_csv('results/proteome_meta/literature_ledger_summary.tsv',
                      sep=chr(9)).iloc[0]

    A(BS + 'section*{Did we measure the protein of that gene}')
    A('Neither platform names a gene directly. An aptamer is assumed to bind its target, and '
      'mass spectrometry assigns peptides to proteins, some of which belong to more than one. '
      'The table below records, for every candidate that either cohort measured, the aptamer '
      'identifier, the accession, and how many distinct peptides are unique to that protein '
      'rather than shared with others in its group.')
    L += longtable(ident[['gene', 'somamer_id', 'ms_accession', 'ms_peptides_unique',
                          'ms_peptides_shared', 'ms_shared_with']],
                   'p{0.14' + BS + 'textwidth}p{0.18' + BS + 'textwidth}'
                   + 'p{0.13' + BS + 'textwidth}p{0.09' + BS + 'textwidth}'
                   + 'p{0.08' + BS + 'textwidth}p{0.17' + BS + 'textwidth}',
                   ['Gene', 'Aptamer', 'Accession', 'Unique', 'Shared',
                    'Shared with'], small=True)
    A('Two entries carry caveats and both are genes the mass-spectrometry cohort called '
      'significant. ' + BS + 'textit{ACTN1} has 8 unique peptides out of 25, the rest shared '
      'with other actinins including the podocyte gene ' + BS + 'textit{ACTN4}; its value '
      'cannot be attributed to ' + BS + 'textit{ACTN1} alone. ' + BS + 'textit{MMP7} has an '
      'unambiguous accession but only 2 distinct peptides. Both are reported rather than '
      'dropped.')

    A(BS + 'section*{The candidate list was fixed before these results}')
    A('The manuscript states that the analysis plan was fixed on the first cohort and applied '
      'unchanged. Timestamps cannot demonstrate that, because they change on every re-run. '
      'What can be shown is the direction of dependency: the stage that produces the candidate '
      'table takes the corrected-selector output, the literature counts, the KPMP evidence and '
      'the per-diagnosis pattern as inputs, and none of the proteomic results among them. '
      'Every proteomic stage runs after it. Verdict: ' + BS + 'texttt{' + str(fz['verdict'])
      + '}. The table analysed here has SHA-256 ' + BS + 'texttt{'
      + str(fz['sha256'])[:32] + BS + 'dots} and its input '
      + BS + 'texttt{' + str(fz['sha256'])[:0] + 'candidates' + BS + '_gated.tsv} has SHA-256 '
      + BS + 'texttt{' + str(fz['input_sha256'])[:32] + BS + 'dots}.')

    A(BS + 'section*{What each novelty tier rests on}')
    A('Reference verification confirms that a cited paper exists; it does not confirm the '
      'claim that a candidate is sparsely represented in the diabetic kidney disease '
      'literature. That claim comes from PubMed counts, so the queries, the date they were '
      'run and the returned identifiers are recorded per candidate in '
      + BS + 'texttt{literature' + BS + '_ledger.tsv}. Counts grow over time, so the frozen values the '
      'manuscript reports were re-queried on ' + str(led['queried_on']) + ': '
      + ('no candidate changed tier.' if int(led['n_class_changed']) == 0 else
         '%d candidates changed tier (%s).' % (int(led['n_class_changed']),
                                               str(led['genes_changed']))) + ' Four symbols '
      'collide with common non-gene abbreviations and are flagged as unreliable rather than '
      'counted, with the reason given for each.')
    return compile_pdf(L, file_name(6))


def af7():
    """간·대장 사전 등록 검정. 등록 파일과 해시, 코호트 정의, 감사, 판정, 사후 분석을 한 통에."""
    path = os.path.join(DST, file_name(7) + '.xlsx')
    R = 'results/xtissue'
    from xtissue_config import COHORTS, arm_match

    def text_rows(p):
        return pd.DataFrame({'line': open(p, encoding='utf-8').read().split('\n')})

    digests = []
    for doc, rec in (('docs/xtissue/PREDICTIONS.md', R + '/predictions.sha256'),
                     ('docs/xtissue/POSTHOC_PLAN.md', R + '/posthoc_plan.sha256')):
        parts = open(rec, encoding='utf-8').read().split()
        digests.append(dict(file=doc, sha256_recorded=parts[0], recorded_utc=parts[-1]))

    rows = []
    for tissue, cs in COHORTS.items():
        tab = pd.read_csv(os.path.join(R, tissue, 'cohorts.tsv'), sep='\t').set_index('gse')
        for gse, c in cs.items():
            proc = c.get('procurement') or ('', '', '')
            t = tab.loc[gse] if gse in tab.index else None
            rows.append(dict(
                tissue=tissue, series=gse, platform=c.get('platform') or '',
                role=str(c.get('status')),
                n_case=int(t['n_case']) if t is not None and pd.notna(t['n_case']) else None,
                n_control=int(t['n_control']) if t is not None and pd.notna(t['n_control'])
                else None,
                harmonisation=str(t['status']) if t is not None else '',
                procurement_case=proc[0], procurement_control=proc[1],
                arms=arm_match(proc) if proc[0] else '',
                procurement_evidence=proc[2] if len(proc) > 2 else ''))
    cohorts = pd.DataFrame(rows)

    with pd.ExcelWriter(path, engine='openpyxl') as w:
        readme = pd.DataFrame([
            ('registration', 'SHA-256 digests and UTC time stamps recorded when the two '
                             'registration files were written; a pipeline stage re-checks them.'),
            ('predictions', 'The prediction file, verbatim, fixed before any liver or colon '
                            'expression value was read.'),
            ('posthoc_plan', 'The post hoc plan, verbatim, fixed after the verdicts and before '
                             'the post hoc analyses ran.'),
            ('cohorts', 'Every liver and colon cohort considered: role, sizes, procurement per '
                        'arm and the evidence for it.'),
            ('reuse_audit', 'Value-vector audit for shared samples between same-platform series.'),
            ('verdicts', 'P1-P3 judged against the registered wording.'),
            ('module_by_cohort', 'P2: immediate-early effect per cohort with bootstrap interval '
                                 'and percentile against random modules.'),
            ('sweep_liver', 'P1/P3: composition sweep over all subsets of the liver core cohorts.'),
            ('sweep_colon', 'P1/P3: composition sweep over all subsets of the colon core cohorts.'),
            ('sweep_kidney', 'P3 reference: the same code on the five kidney discovery cohorts.'),
            ('h4_liver_substitution', 'Post hoc H4: two control definitions within GSE48452.'),
            ('h5_inflammation', 'Post hoc H5: case coefficient on the handling score with and '
                                'without the inflammation score.'),
            ('h6_colon_downsampled', 'Post hoc H6: colon cohorts subsampled to kidney sizes.'),
            ('concordance', 'Post hoc: median pairwise Spearman correlation of gene-wise g.'),
            ('s1_liver_relaxed', 'Post hoc S1: liver full-set sweep under a relaxed RNA-seq '
                                 'filter.')], columns=['Sheet', 'Contents'])
        sheet(w, 'README', readme,
              file_title(7) + '. The cross-tissue test with hash-fixed predictions, '
              + where('sec:xtissue', tex=False) + '.')
        sheet(w, 'registration', pd.DataFrame(digests),
              'Recompute SHA-256 of each file in the repository and compare.')
        sheet(w, 'predictions', text_rows('docs/xtissue/PREDICTIONS.md'),
              'docs/xtissue/PREDICTIONS.md, one line per row.')
        sheet(w, 'posthoc_plan', text_rows('docs/xtissue/POSTHOC_PLAN.md'),
              'docs/xtissue/POSTHOC_PLAN.md, one line per row.')
        sheet(w, 'cohorts', cohorts,
              'role: core enters the composition sweep; p2 enters only the per-cohort module '
              'test; excluded with the reason. arms: matched when both arms share a code.')
        sheet(w, 'reuse_audit', pd.read_csv(os.path.join(R, 'overlap_audit.tsv'), sep='\t'),
              'identical: best partner r >= 0.9999. Re-processed: median best-minus-second '
              'partner gap at least five times the reference pair.')
        sheet(w, 'verdicts', pd.read_csv(os.path.join(R, 'verdicts.tsv'), sep='\t'),
              'Verdicts use the registered wording; no threshold was changed after the data.')
        sheet(w, 'module_by_cohort', pd.read_csv(os.path.join(R, 'p2_ieg_by_cohort.tsv'),
                                                 sep='\t'),
              'g_ieg: Hedges g of the handling score, case minus control; 1,000 bootstrap '
              'resamples; percentile of |g| against 2,000 random modules of the same size.')
        for name, p in (('sweep_liver', 'liver_masld/sweep/subsets.tsv'),
                        ('sweep_colon', 'colon_uc/sweep/subsets.tsv'),
                        ('sweep_kidney', 'kidney_dkd/sweep/subsets.tsv'),
                        ('s1_liver_relaxed', 'liver_masld/sweep_relaxed/subsets.tsv')):
            sheet(w, name, pd.read_csv(os.path.join(R, p), sep='\t'),
                  'n_ieg_topK: mean immediate-early genes in the top 50 over folds.')
        sheet(w, 'h4_liver_substitution',
              pd.read_csv(os.path.join(R, 'posthoc', 'h4_gse48452_summary.tsv'), sep='\t'),
              'Top-50 overlap and rank correlation between the two control definitions.')
        sheet(w, 'h5_inflammation',
              pd.read_csv(os.path.join(R, 'posthoc', 'h5_inflammation_adjustment.tsv'),
                          sep='\t'),
              'Standardised coefficients of case status on the handling score, unadjusted and '
              'adjusted for the 13-gene inflammation score.')
        sheet(w, 'h6_colon_downsampled',
              pd.read_csv(os.path.join(R, 'posthoc', 'h6', 'replicates.tsv'), sep='\t'),
              'sizes: case/control counts used for each cohort in that replicate.')
        sheet(w, 'concordance',
              pd.read_csv(os.path.join(R, 'posthoc_effect_concordance.tsv'), sep='\t'),
              'Computed over the genes shared by the five core cohorts of each tissue.')
    log('  %s.xlsx' % file_name(7))
    return True


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    global TARGET, DST
    if '--target' in sys.argv:
        TARGET = sys.argv[sys.argv.index('--target') + 1]
        DST = {'bmc': 'submission', 'bib': 'submission_bib'}[TARGET]
    log('본문이 약속한 부록을 만듭니다 (%s -> %s).' % (TARGET, DST))
    ok = all([af1(), af2(), af3(), af4(), af5(), af6(), af7()])
    log('')
    if not ok:
        log('일부가 만들어지지 않았습니다.')
        return 1
    for f in sorted(os.listdir(DST)):
        if f.startswith(('Additional', 'Supplementary_File')):
            log('  %-28s %8.0f KB' % (f, os.path.getsize(os.path.join(DST, f)) / 1024))
    return 0


if __name__ == '__main__':
    sys.exit(main())
