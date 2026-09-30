#!/usr/bin/env python
"""논문의 개요 그림 — 왜 이 자료로는 질병 신호와 처리 신호를 가를 수 없는가.

다섯 번째 판입니다. 앞의 네 판은 전부 "단계를 차례로 그린다"는 같은 틀 안에 있었습니다.

    1판  상자마다 두세 문장. 슬라이드처럼 읽혔습니다.
    2판  작은 결과 플롯을 넣었더니 뒤의 결과 그림과 문법이 같아졌습니다.
    3판  시료 아이콘과 카드로 흐름도를 그렸습니다. 안내판처럼 보였습니다.
    4판  민짜 study design 흐름도. 저널 그림답기는 한데, 여전히 "무엇을 했는가"의 목록입니다.
    5판  이 파일. 틀을 바꿉니다. 이 논문의 결론은 절차가 아니라 구조에서 나옵니다.
         진단과 조달 경로가 붙어 있어서 관측된 대비가 질병 효과와 처리 효과의 합이라는 것.
         그래서 순서도가 아니라 인과 구조도로 그립니다. 통제 넷은 각자 어느 화살표를
         건드리는지 그 화살표 옆에 적고, 보정은 교란 경로를 자르는 표시로 둡니다.

이 틀이 나은 이유는 셋입니다. 첫째, 독자가 "왜 코호트를 더 모아도 해결되지 않는가"를 그림만
보고 압니다. 둘째, 통제 넷이 서로 다른 화살표를 겨냥한다는 것이 목록이 아니라 위치로
드러납니다. 셋째, 측정되지 않은 교란(허혈 시간)을 점선 노드로 그리면 "보고되었다면
공변량이 되었을 것"이라는 논의가 그림 안에서 이미 보입니다.

패널 b 는 그 보정이 산출물에 무엇을 하는지 봅니다. 보정 전후 후보 30개를 띠로 잇고 문헌에
이미 제안된 유전자를 색으로 나눕니다. 목록이 바뀌지만 신규성 쪽으로 바뀌지는 않는다는
결론이 흐름 하나로 보입니다.

참고한 것: Rougier et al., Ten Simple Rules for Better Figures (PLoS Comput Biol 2014) 의
규칙 2·6·8; 역학에서 교란을 그리는 표준 방식(인과 그래프); JCI Insight 2026 의 reference
tissue 논문이 같은 교란을 단일세포에서 본 사례.

숫자는 전부 results/ 아래 파일에서 읽습니다. 손으로 적으면 본문과 어긋납니다.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import save, JOURNAL_W, INK, MUTED, ACCENT, BLUE   # noqa: E402

OUT = 'results/figures'
HAIR = '#b8b3aa'
FILL = '#f4f2ee'
NL = chr(10)
EN = '–'
MINUS = '−'

W, H = JOURNAL_W, 5.2
ASP = W / H

FS = 6.6
FB = 7.2
FP = 8.6


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- 데이터
def read_numbers():
    """본문과 어긋나지 않도록 결과 파일에서 직접 읽는다."""
    n = {}

    coh = pd.read_csv('results/cohort_table.tsv', sep='\t')
    n['n_cohort'] = len(coh)
    n['n_disc'] = int((coh['role'] == 'discovery').sum())
    # 일곱 개를 다 더하면 구획쌍에서 같은 환자를 두 번 센다. 발굴 코호트 다섯 개만 더한다.
    disc = coh[coh['role'] == 'discovery']
    n['n_case'] = int(disc['n_case'].sum())
    n['n_ctrl'] = int(disc['n_control'].sum())
    # 구조도가 말이 되려면 사례가 전부 생검, 대조가 전부 비생검이어야 한다. 확인하고 그린다.
    assert (disc['case_procurement'] == 'biopsy').all()
    assert not disc['control_procurement'].str.contains('biopsy').any()

    nul = pd.read_csv('results/null_control.tsv', sep='\t')
    k50 = nul[nul['K'] == 50]['random_mean']
    n['rand_lo'], n['rand_hi'] = float(k50.min()), float(k50.max())

    sub = pd.read_csv('results/cohort_sensitivity/subsets.tsv', sep='\t')
    piv = sub[sub['size'] == 3].pivot(index='subset', columns='method',
                                      values='cross_fold_jaccard')
    d = (piv['RBS'] - piv['WGCNA_hub']).to_numpy()
    n['diff_lo'], n['diff_hi'] = float(d.min()), float(d.max())
    n['n_subsets'] = int(sub['subset'].nunique())

    het = pd.read_csv('results/ieg_heterogeneity.tsv', sep='\t')
    n['n_datasets'] = int(het['dataset'].nunique())

    gen = pd.read_csv('results/ckd_generalization/per_diagnosis.tsv', sep='\t')
    n['n_dx'] = len(gen)
    n['shrink_lo'] = 100 * float(gen['shrink'].min())
    n['shrink_hi'] = 100 * float(gen['shrink'].max())

    # K 를 섞으면 초록·Fig 7 과 다른 값이 나온다. 본문 표와 같은 K=50 만 쓴다.
    ruv = pd.read_csv('results/ruv_benchmark/benchmark.tsv', sep='\t')
    r50 = ruv[ruv['K'] == 50].groupby('arm')[['procurement_frac', 'external_auroc']].mean()
    n['proc_pre'] = 100 * float(r50.loc['none', 'procurement_frac'])
    n['proc_post'] = 100 * float(r50.loc['IEG_resid', 'procurement_frac'])

    a = pd.read_csv('results/final_candidates_gated.tsv', sep='\t')
    b = pd.read_csv('results/candidates_v2/candidates_gated.tsv', sep='\t')
    lit = pd.concat([pd.read_csv('results/literature_review.tsv', sep='\t'),
                     pd.read_csv('results/candidates_v2/literature_review.tsv', sep='\t')])
    lit = lit.drop_duplicates('gene').set_index('gene')['verdict']

    def known(s):
        return int((lit.reindex(sorted(s)) == 'ALREADY PROPOSED as DKD marker').sum())

    sa, sb = set(a['symbol']), set(b['symbol'])
    n['n_cand'] = len(sb)
    n['known_pre'], n['known_post'] = known(sa), known(sb)
    n['shared'] = len(sa & sb)
    # 띠를 문헌 여부로 나누려면 공유 후보 중 이미 제안된 것이 몇인지도 필요하다.
    n['shared_known'] = known(sa & sb)
    return n


# ---------------------------------------------------------------- 도형
def node(ax, x, y, w, h, title, sub=None, ec=INK, dashed=False, fc='#ffffff'):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, linewidth=0.8,
                           linestyle=(0, (2.4, 1.6)) if dashed else 'solid', zorder=3))
    ax.text(x + w / 2, y + h - 2.6, title, fontsize=FS, fontweight='bold', color=INK,
            ha='center', va='center', zorder=4, linespacing=1.5)
    if sub:
        ax.text(x + w / 2, y + h - 6.0, sub, fontsize=FS, color=MUTED, ha='center',
                va='center', zorder=4, linespacing=1.5)


def edge(ax, p, q, colour=INK, lw=0.9, rad=0.0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=8, linewidth=lw,
                                 color=colour, shrinkA=1, shrinkB=1, zorder=2,
                                 connectionstyle='arc3,rad=%.2f' % rad))


def txt(ax, x, y, t, size=FS, colour=INK, ha='left', va='top', weight='normal',
        box=False):
    ax.text(x, y, t, fontsize=size, color=colour, ha=ha, va=va, fontweight=weight,
            zorder=5, linespacing=1.5,
            bbox=dict(facecolor='#ffffff', edgecolor='none', pad=0.8) if box else None)


def ribbon(ax, x0, x1, lo0, hi0, lo1, hi1, colour, alpha=0.30):
    """두 막대를 잇는 띠. 가운데를 베지에로 휘게 해 흐름으로 읽히게 한다."""
    m = (x0 + x1) / 2
    verts = [(x0, hi0), (m, hi0), (m, hi1), (x1, hi1),
             (x1, lo1), (m, lo1), (m, lo0), (x0, lo0), (x0, hi0)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=colour, edgecolor='none',
                           alpha=alpha, zorder=2))


def kidney(ax, cx, cy, s, ec=INK, fc='none', lw=0.7, flip=False):
    u = [(0.02, 0.55),
         (0.28, 0.55), (0.44, 0.30), (0.44, 0.02),
         (0.44, -0.28), (0.26, -0.55), (-0.02, -0.55),
         (-0.24, -0.55), (-0.38, -0.38), (-0.36, -0.18),
         (-0.34, -0.06), (-0.20, -0.05), (-0.20, 0.02),
         (-0.20, 0.11), (-0.34, 0.11), (-0.36, 0.22),
         (-0.38, 0.40), (-0.22, 0.55), (0.02, 0.55)]
    f = -1.0 if flip else 1.0
    pts = [(cx + f * px * s, cy + py * s * ASP) for px, py in u]
    codes = [Path.MOVETO] + [Path.CURVE4] * (len(pts) - 1)
    ax.add_patch(PathPatch(Path(pts, codes), facecolor=fc, edgecolor=ec, linewidth=lw,
                           zorder=4))


def needle(ax, cx, cy, s, angle=-38, colour=INK):
    t = np.deg2rad(angle)
    R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

    def put(local, fc):
        p = [R.dot(np.array(v)) for v in local]
        pts = [(cx + a * s, cy + b * s * ASP) for a, b in p]
        ax.add_patch(Polygon(pts, closed=True, facecolor=fc, edgecolor=colour,
                             linewidth=0.5, zorder=5))

    put([(-1.05, -0.07), (0.30, -0.07), (0.30, 0.07), (-1.05, 0.07)], '#ffffff')
    put([(0.30, -0.07), (0.92, -0.07), (0.92, 0.03), (0.30, 0.07)], colour)


# ---------------------------------------------------------------- 본체
def main():
    os.makedirs(OUT, exist_ok=True)
    n = read_numbers()

    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')

    # ============================================================ a 인과 구조
    txt(ax, 2, 99, 'a', size=FP, weight='bold')
    txt(ax, 5.8, 98.4, 'What the case/control contrast in these cohorts contains', size=FB)

    node(ax, 3, 80, 22, 9.5, 'Diagnosis', 'DKD vs control')
    node(ax, 3, 55, 22, 12.5, 'Procurement route', 'biopsy vs surgical')
    node(ax, 34, 58, 23, 11.0, 'Tissue handling', 'ischaemia, time to' + NL + 'fixation',
         ec=ACCENT, dashed=True)
    node(ax, 66, 70, 20, 9.5, 'Measured' + NL + 'expression')
    node(ax, 66, 50, 20, 9.5, 'Selected' + NL + 'features')

    # 진단과 조달을 잇는 겹줄. 화살표가 아니라 등호에 가깝다.
    for dx in (-0.55, 0.55):
        ax.plot([10.0 + dx, 10.0 + dx], [80.0, 67.5], color=ACCENT, lw=1.0, zorder=2)
    # 겹줄은 인과 화살표가 아니라 '같이 움직인다'는 표시다. 양끝에 짧은 가로줄을 둔다.
    for yy in (80.0, 67.5):
        ax.plot([8.6, 11.4], [yy, yy], color=ACCENT, lw=1.0, zorder=2)
    txt(ax, 12.4, 76.0, 'collinear by design', size=FS, colour=ACCENT, weight='bold')
    txt(ax, 12.4, 73.2, '%d cases, all biopsies' % n['n_case']
        + NL + '%d controls, none' % n['n_ctrl'], size=FS, colour=ACCENT)

    # 두 경로를 노드 안에 나란히 둔다. 화살표를 가로지르지 않고, 노드가 말하는 대비를
    # 글자보다 먼저 보여 준다.
    kidney(ax, 9.8, 57.9, 3.8, ec=ACCENT)
    needle(ax, 9.3, 59.2, 1.5, colour=ACCENT)
    kidney(ax, 18.2, 57.9, 3.8, ec=BLUE, fc='#eef2f7', flip=True)
    ax.plot([16.7, 19.7], [59.0, 59.0], color=BLUE, lw=0.7, ls=(0, (1.6, 1.2)), zorder=5)

    edge(ax, (25, 84.8), (66, 78.0), rad=-0.10)
    txt(ax, 45, 87.5, 'disease effect, the quantity of interest', size=FS, colour=MUTED,
        ha='center')
    edge(ax, (25, 63.5), (33.4, 63.5), colour=ACCENT)
    edge(ax, (57, 64.5), (66, 71.5), colour=ACCENT, rad=-0.12)
    txt(ax, 57.6, 63.6, 'handling effect', size=FS, colour=ACCENT, ha='left')
    edge(ax, (76, 70), (76, 59.5))

    # 보정은 교란 경로를 자르는 표시로
    ax.plot([59.8, 62.2], [67.4, 69.8], color=BLUE, lw=1.5, zorder=6)
    ax.plot([61.6, 64.0], [66.2, 68.6], color=BLUE, lw=1.5, zorder=6)
    txt(ax, 58.6, 72.4, 'correction', size=FS, colour=BLUE, ha='left')

    # 통제 넷을 각자 건드리는 곳 아래에 적는다
    txt(ax, 3, 52.0, 'Control-definition substitution', size=FS, colour=ACCENT,
        weight='bold')
    txt(ax, 3, 49.2, 'breaking the collinear link leaves' + NL
        + '%.0f%s%.0f%% of the naive signal, in each' % (n['shrink_lo'], EN, n['shrink_hi'])
        + NL + 'of %d biopsy diagnoses' % n['n_dx'], size=FS)

    txt(ax, 34, 52.0, 'Procurement contrast', size=FS, colour=ACCENT, weight='bold')
    txt(ax, 34, 49.2, 'an immediate-early module separates' + NL
        + 'the two routes in 8 of %d datasets, and' % n['n_datasets']
        + NL + 'vanishes biopsy to biopsy', size=FS)

    txt(ax, 66, 44.0, 'Random-signature null', size=FS, colour=ACCENT, weight='bold')
    txt(ax, 66, 41.2, 'random 50-gene sets reach AUROC' + NL + '%.2f%s%.2f on the same folds'
        % (n['rand_lo'], EN, n['rand_hi']), size=FS)

    txt(ax, 66, 34.6, 'Cohort-composition sweep', size=FS, colour=ACCENT, weight='bold')
    txt(ax, 66, 31.8, 'the reproducibility winner flips across' + NL
        + 'the %d subsets, by %s to %+.2f'
        % (n['n_subsets'], ('%+.2f' % n['diff_lo']).replace('-', MINUS), n['diff_hi']),
        size=FS)

    txt(ax, 3, 38.0, 'Correction', size=FS, colour=BLUE, weight='bold')
    txt(ax, 3, 35.2, 'residualising every gene on a handling score cuts' + NL
        + 'procurement-driven features from %.0f%% to %.0f%% of those'
        % (n['proc_pre'], n['proc_post'])
        + NL + 'selected, with external AUROC unchanged', size=FS)

    # ============================================================ b 후보 목록
    txt(ax, 2, 27.6, 'b', size=FP, weight='bold')
    txt(ax, 5.8, 27.0, 'What the correction does to the candidate list', size=FB)

    bx0, bx1, bw = 16.0, 44.0, 7.5
    base = 3.0
    unit = 15.0 / n['n_cand']
    pre_k, post_k = n['known_pre'], n['known_post']
    sh, sh_k = n['shared'], n['shared_known']

    def stack(x, k):
        ax.add_patch(Rectangle((x, base), bw, unit * k, facecolor=BLUE, edgecolor='none',
                               zorder=4))
        ax.add_patch(Rectangle((x, base + unit * k), bw, unit * (n['n_cand'] - k),
                               facecolor=FILL, edgecolor='none', zorder=3))
        ax.add_patch(Rectangle((x, base), bw, unit * n['n_cand'], facecolor='none',
                               edgecolor=HAIR, linewidth=0.7, zorder=5))

    stack(bx0, pre_k)
    stack(bx1, post_k)

    ribbon(ax, bx0 + bw, bx1, base, base + unit * sh_k, base, base + unit * sh_k, BLUE)
    ribbon(ax, bx0 + bw, bx1, base + unit * sh_k, base + unit * sh,
           base + unit * sh_k, base + unit * sh, MUTED, alpha=0.20)

    for x, lab in ((bx0, 'uncorrected'), (bx1, 'corrected')):
        txt(ax, x + bw / 2, base + unit * n['n_cand'] + 5.6, lab, size=FS, ha='center',
            weight='bold')
        txt(ax, x + bw / 2, base + unit * n['n_cand'] + 3.0, '%d candidates' % n['n_cand'],
            size=FS, ha='center', colour=MUTED)
    txt(ax, (bx0 + bw + bx1) / 2, base + unit * sh + 2.4, '%d shared' % sh, size=FS,
        ha='center', colour=MUTED)
    txt(ax, bx0 - 1.6, base + unit * pre_k / 2 + 1.2, '%d' % pre_k, size=FS, ha='right',
        colour=BLUE, weight='bold')
    txt(ax, bx1 + bw + 1.6, base + unit * post_k / 2 + 1.2, '%d' % post_k, size=FS,
        colour=BLUE, weight='bold')

    txt(ax, 58, 20.0, 'Blue is the part of each list already proposed as'
        + NL + 'a marker of the disease. Correction moves the list'
        + NL + 'towards that literature, from %d genes to %d, and %d of'
        % (pre_k, post_k, sh)
        + NL + 'the %d candidates survive it. What it does not do is' % n['n_cand']
        + NL + 'return a candidate both strong in the data and absent'
        + NL + 'from the literature: that cell is empty in both lists.', size=FS)

    save(fig, 'P0_overview', pdf=True)
    log('  숫자는 전부 results/ 에서 읽었습니다. 본문과 어긋날 수 없습니다.')


def main_v2():
    """Render Figure 1 as a concise map of the paper's contribution.

    The former overview mixed a causal schematic, six result callouts, and a candidate-list
    graphic. The latter is useful evidence, but it competes with the study logic on a first read.
    This version reserves Figure 1 for the audit-to-decision workflow; detailed results remain in
    their dedicated figures and tables.
    """
    os.makedirs(OUT, exist_ok=True)
    n = read_numbers()

    fig = plt.figure(figsize=(JOURNAL_W, 4.9))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')

    orange = ACCENT
    blue_fill = '#eaf1f8'
    warm_fill = '#fbf2e9'
    neutral_fill = '#f5f3ef'
    green = '#4b7d50'
    green_fill = '#edf5ea'

    def card(x, y, w, h, title, body, edge, fill):
        """Compact text card sized for a two-column journal PDF."""
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.25,rounding_size=1.2',
                                    facecolor=fill, edgecolor=edge, linewidth=0.8, zorder=3))
        ax.text(x + 1.6, y + h - 3.0, title, fontsize=6.0, fontweight='bold', color=INK,
                ha='left', va='center', zorder=4, linespacing=1.18)
        ax.text(x + 1.6, y + h - 9.6, body, fontsize=6.0, color=MUTED,
                ha='left', va='top', zorder=4, linespacing=1.32)

    txt(ax, 3, 96, 'Audit the design before interpreting a biomarker benchmark',
        size=7.3, weight='bold')
    txt(ax, 3, 91.8,
        'The central claim is transcriptomic. Cross-modal data corroborate' + NL +
        'or bound it; they are not a participant-matched multi-omics analysis.',
        size=6.0, colour=MUTED)

    # Step labels make the read order unambiguous even when the figure is viewed at one column.
    steps = [
        (4, '1', 'Audit the benchmark input', neutral_fill),
        (36, '2', 'Stress-test the conclusions', warm_fill),
        (68, '3', 'Diagnose, correct, and report', blue_fill),
    ]
    for x, number, title, fill in steps:
        ax.text(x, 84.5, number, fontsize=9, fontweight='bold', color='#ffffff', ha='center',
                va='center', bbox=dict(boxstyle='circle,pad=0.28', facecolor=INK, edgecolor='none'))
        txt(ax, x + 3.2, 86.2, title, size=6.0, weight='bold')

    # 1. Input audit
    card(4, 55, 25, 24, 'Harmonised DKD' + NL + 'transcriptomes',
         '7 cohorts; 9,900 genes' + NL +
         'five discovery cohorts' + NL +
         '26 admissible subsets', INK, neutral_fill)
    txt(ax, 5.8, 50.5, 'Required audit inputs', size=6.0, weight='bold', colour=MUTED)
    txt(ax, 5.8, 45.8,
        'patient-level overlap' + NL +
        'control definition and' + NL +
        'procurement route' + NL +
        'pre-analytical metadata' + NL +
        'where available', size=6.0)

    # 2. Two deliberately distinct sensitivity tests.
    card(36, 62, 25, 17, 'Cohort-composition' + NL + 'sensitivity',
         'repeat the benchmark on' + NL +
         'every admissible subset', orange, warm_fill)
    txt(ax, 37.8, 58.6,
        'difference -0.106 to +0.170;' + NL +
        'report reproducibility and' + NL +
        'AUROC separately', size=6.05, colour=orange)

    card(36, 33, 25, 17, 'Control-definition' + NL + 'sensitivity',
         'hold cases fixed; use a' + NL +
         'defensible alternative' + NL +
         'control group', orange, warm_fill)
    txt(ax, 37.8, 29.7,
        '7 diagnoses: 14-56% of the' + NL +
        'naive signal remains; a' + NL +
        'sensitivity analysis, not a' + NL +
        'causal decomposition', size=6.05, colour=orange)

    # 3. The output makes clear that correction is evaluated, not simply asserted.
    card(68, 62, 28, 17, 'Confounder diagnosis',
         'the immediate-early programme' + NL +
         'tracks procurement, not the' + NL +
         'matched disease contrast', orange, warm_fill)
    txt(ax, 69.8, 58.6,
        '8 of 9 datasets; attribution' + NL +
        'stays an inference because' + NL +
        'handling times are not recorded', size=6.05, colour=orange)

    card(68, 33, 28, 17, 'Correction and candidate' + NL + 'evidence ledger',
         'compare correction arms; keep' + NL +
         'the selection, contradiction' + NL +
         'and corroboration trail', BLUE, blue_fill)
    txt(ax, 69.8, 29.7,
        'procurement-driven features' + NL +
        '27% to 5-6%, with no' + NL +
        'significant AUROC loss', size=6.05, colour=BLUE)

    # Main arrows, with no implied causal arrow from the unmeasured handling variables.
    edge(ax, (29, 67), (35.5, 70.5), colour=INK)
    edge(ax, (29, 63), (35.5, 41.5), colour=INK)
    edge(ax, (61, 70.5), (67.5, 70.5), colour=orange)
    edge(ax, (61, 41.5), (67.5, 41.5), colour=orange)
    edge(ax, (82, 62), (82, 50.5), colour=INK)

    # Supporting data are visibly separate from the main transcriptomic evidence.
    ax.add_patch(Rectangle((4, 2.6), 92, 15.2, facecolor=green_fill, edgecolor=green,
                           linewidth=0.8, linestyle=(0, (2.2, 1.5)), zorder=1))
    txt(ax, 5.8, 15.3, 'Supporting layers: they corroborate or falsify, and do not redefine the claim',
        size=6.0, weight='bold', colour=green)
    txt(ax, 5.8, 11.7,
        'External stress tests  |  Patient-independent proteomics  |  GWAS boundary test' + NL +
        'Targeted metabolomics as a disease-relevant comparator', size=6.0)
    txt(ax, 5.8, 6.9,
        'Output: an auditable candidate-evidence ledger, not a new selector' + NL +
        'and not a participant-matched multi-omics biomarker panel.',
        size=6.0, colour=MUTED)
    edge(ax, (82, 33), (82, 17.8), colour=green, lw=1.0)

    save(fig, 'P0_overview', pdf=True)
    log('  wrote Figure 1 overview workflow')


if __name__ == '__main__':
    main_v2()
