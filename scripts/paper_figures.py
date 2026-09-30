#!/usr/bin/env python
"""The manuscript figures, rebuilt around what the analysis actually supports.

P0  framework overview - main pipeline, control analyses and exploratory extensions
P1  cohort-composition sensitivity - the verdict depends on which cohorts you have
P2  random-signature null - absolute AUROC carries little information here
P3  the procurement axis across nine kidney datasets, plus the heart counter-example
P4  KPMP single-nucleus with its biopsy-vs-biopsy negative control
P5  module-size bias in simulation, with the real selector
P6  corrections, reported as enrichment over the gate's 8.3% background rate
"""
import os
import shutil
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import JOURNAL_W, BIB_W, save, INK, MUTED, ACCENT, BLUE, GREY

OUT = 'results/figures'
# 그림을 그릴 폭. BMC 는 372pt, BiB 는 526pt 다. 같은 그림을 두 폭으로 그려야
# 어느 쪽에서도 확대·축소 없이 들어가고 글자 크기가 본문과 맞는다.
W = JOURNAL_W
SUF = ''
PANEL_BELOW = False
GREEN = '#4a7c3f'
METHODS = ['DEG_meta', 'WGCNA_hub', 'RBS', 'RBS_orth']
NICE = {'DEG_meta': 'DEG meta', 'WGCNA_hub': 'connectivity', 'RBS': 'stability (RBS)',
        'RBS_orth': 'stability + correction'}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def panels_below(fig, axes, pad=0.085, letters='abcdefg', letter_y=None):
    """패널 제목을 떼고 (a)(b)(c) 를 그림 아래에 둔다.

    좁은 패널에서는 제목이 옆 패널을 덮거나 잘린다. 설명은 캡션이 맡고
    그림에는 글자만 남긴다. tight_layout 뒤의 실제 위치를 읽어 가운데 맞춘다.
    """
    if not PANEL_BELOW:
        return
    for ax in axes:
        # set_title('') 는 가운데 제목만 지운다. loc='left' 로 쓴 제목은 딴 객체라 남는다.
        for loc in ('left', 'center', 'right'):
            ax.set_title('', loc=loc)
    fig.tight_layout(rect=(0, pad, 1, 1))
    y = pad * 0.18 if letter_y is None else letter_y
    for ax, ch in zip(axes, letters):
        bb = ax.get_position()
        fig.text((bb.x0 + bb.x1) / 2, y, ch, ha='center', va='bottom',
                 fontsize=8, fontweight='bold', color=INK)


# ---------------------------------------------------------------- P0
def p0_overview():
    fig, ax = plt.subplots(figsize=(13, 7.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    def box(x, y, w, h, title, lines, fc, ec=INK, lw=1.0):
        patch = FancyBboxPatch((x, y), w, h,
                               boxstyle='round,pad=0.012,rounding_size=0.02',
                               facecolor=fc, edgecolor=ec, linewidth=lw)
        ax.add_patch(patch)
        ax.text(x + 0.018, y + h - 0.04, title, fontsize=7.8, fontweight='bold',
                va='top', ha='left', color=INK)
        ax.text(x + 0.018, y + h - 0.085, '\n'.join(lines), fontsize=6.63,
                va='top', ha='left', color=INK, linespacing=1.35)

    def arrow(a, b, color=INK, lw=1.4):
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle='-|>', mutation_scale=12,
                                     linewidth=lw, color=color))

    box(0.03, 0.64, 0.23, 0.24, 'Input cohorts',
        ['7 DKD transcriptome cohorts',
         'frozen 9,900-gene space',
         'patient-level audit',
         'ERCB other-CKD biopsies as matched controls'],
        fc='#f3efe5')
    box(0.03, 0.33, 0.23, 0.22, 'External resources',
        ['GSE175759, GSE162830, GSE20602',
         'GSE142153 and small-sample checks',
         'KPMP single-nucleus and proteomics',
         'ST003255 metabolomics, GWAS'],
        fc='#eef3f8')

    box(0.33, 0.73, 0.28, 0.14, 'Main analysis',
        ['harmonise cohorts -> LODO benchmark',
         'repeat over all 26 cohort subsets'],
        fc='#f8eadf', ec=ACCENT, lw=1.2)
    box(0.33, 0.53, 0.28, 0.14, 'Two benchmark axes',
        ['reproducibility: cross-fold Jaccard',
         'predictive performance: external AUROC',
         'reported separately, not interchangeably'],
        fc='#fff7f1')
    box(0.33, 0.30, 0.28, 0.17, 'Confounder diagnosis',
        ['immediate-early module across cohorts',
         'procurement-matched negative controls',
         'module-size mechanism and simulation'],
        fc='#fff1e8')
    box(0.33, 0.08, 0.28, 0.15, 'Correction benchmark',
        ['handling-score residualisation',
         'SVA, RUVg, ComBat',
         'remove confounder vs preserve AUROC'],
        fc='#fff7f1')

    box(0.69, 0.71, 0.27, 0.17, 'Method claim',
        ['cohort composition changes the',
         'reproducibility verdict',
         'absolute AUROC is weakly informative here'],
        fc='#e9f2fb', ec=BLUE, lw=1.2)
    box(0.69, 0.45, 0.27, 0.19, 'Disease-specific finding',
        ['the most reproducible signal tracks',
         'tissue procurement, not DKD biology',
         'the design problem persists across',
         'eight kidney diagnoses'],
        fc='#fbe9e2', ec=ACCENT, lw=1.2)
    box(0.69, 0.18, 0.27, 0.19, 'Candidate and extension layer',
        ['candidate list is stress-tested, not promoted',
         'cross-omics / GWAS / external datasets',
         'are corroboration or falsification layers'],
        fc='#edf6ee', ec=GREEN, lw=1.2)

    ax.add_patch(Rectangle((0.665, 0.08), 0.305, 0.33, facecolor='none',
                           edgecolor=MUTED, linewidth=1.0, linestyle='--'))
    ax.text(0.682, 0.392, 'Exploratory / stress-test branch', fontsize=6.63,
            color=MUTED, va='bottom')

    arrow((0.26, 0.76), (0.33, 0.80))
    arrow((0.26, 0.44), (0.33, 0.38))
    arrow((0.47, 0.73), (0.47, 0.67))
    arrow((0.47, 0.53), (0.47, 0.47))
    arrow((0.47, 0.30), (0.47, 0.23))
    arrow((0.61, 0.80), (0.69, 0.80), color=BLUE)
    arrow((0.61, 0.39), (0.69, 0.54), color=ACCENT)
    arrow((0.61, 0.16), (0.69, 0.27), color=GREEN)
    arrow((0.26, 0.44), (0.69, 0.27), color=MUTED, lw=1.2)

    ax.text(0.03, 0.95,
            'Study overview: the benchmark, the confounder controls, and the stress-test layers',
            fontsize=9.36, fontweight='bold', ha='left')
    ax.text(0.03, 0.91,
            'Solid boxes are manuscript-core analyses; the dashed band marks extensions that '
            'can corroborate or weaken candidate claims without changing the main design result.',
            fontsize=6.63, color=MUTED, ha='left')
    # 이 그림은 본문에서 빠졌다. overview_figure.py 의 P0_overview 가 같은 자리를 대신하고,
    # 그쪽은 숫자를 results/ 에서 읽으므로 본문과 어긋나지 않는다. 여기는 저장소 문서용으로
    # 남겨 둔다 — 자료 목록과 보조 분기를 한 장에 보여주는 용도로는 아직 쓸모가 있다.
    save(fig, 'P0_framework_overview', pdf=True)


# ---------------------------------------------------------------- P1
def f1_cohort_sensitivity():
    d = pd.read_csv('results/cohort_sensitivity/subsets.tsv', sep='\t')
    fig, axes = plt.subplots(1, 3, figsize=(W, 2.9),
                             gridspec_kw={'width_ratios': [1.15, 1.15, 1]})

    # (a) winner distribution by subset size, reproducibility
    #
    # 동률을 세면 안 된다. 두 코호트만 쓰면 LODO fold 가 둘뿐이라 cross-fold Jaccard 가
    # 네 방법 모두 0.000 으로 붙는 부분집합이 나오는데, idxmax() 는 그 경우 첫 방법을
    # 임의로 승자로 돌려준다. 그렇게 세면 size=2 막대 전체가 열 순서의 산물이 된다.
    # 본문 표는 이미 size=2 를 빼고 있었으므로 그림도 맞춘다.
    ax = axes[0]
    sizes = [s for s in sorted(d['size'].unique()) if s >= 3]
    cols = {'DEG_meta': GREY, 'WGCNA_hub': BLUE, 'RBS': ACCENT, 'RBS_orth': GREEN}
    n_dec = {}
    for s in sizes:
        piv = d[d['size'] == s].pivot(index='subset', columns='method',
                                      values='cross_fold_jaccard').dropna()
        top = piv.max(axis=1)
        n_dec[s] = piv[(piv.ge(top - 1e-9, axis=0)).sum(axis=1) == 1]
    bottom = np.zeros(len(sizes))
    for m in METHODS:
        frac = []
        for s in sizes:
            piv = n_dec[s]
            w = piv.idxmax(axis=1).value_counts()
            frac.append(100 * w.get(m, 0) / len(piv) if len(piv) else 0)
        ax.bar(range(len(sizes)), frac, bottom=bottom, color=cols[m], edgecolor=INK,
               linewidth=0.5, label=NICE[m])
        for i, (f, b) in enumerate(zip(frac, bottom)):
            if f >= 12:
                ax.text(i, b + f / 2, '%.0f' % f, ha='center', va='center', fontsize=6.24,
                        color='white' if m in ('WGCNA_hub', 'RBS') else INK)
        bottom += np.array(frac)
    ax.set_xticks(range(len(sizes)))
    ax.set_xticklabels(['%d\n(n=%d)' % (s, len(n_dec[s])) for s in sizes], fontsize=6.63)
    ax.set_xlabel('cohorts in the benchmark')
    ax.set_ylabel('% of subsets where the method wins')
    ax.set_ylim(0, 100)
    ax.set_title('(a) winner by reproducibility', loc='left', fontsize=7)
    # 범례를 (a) 축에 매달면 그 축의 x 제목과 겹치고, tight_layout 이 축을 위로 밀어
    # 아래가 비어 버린다. 범례는 네 패널 공통이므로 그림 전체의 하단에 한 줄로 둔다.
    legend_handles = ax.get_legend_handles_labels()

    # (b) the difference, every subset
    ax = axes[1]
    for i, s in enumerate(sizes):
        piv = d[d['size'] == s].pivot(index='subset', columns='method',
                                      values='cross_fold_jaccard').dropna()
        if piv.empty:
            continue
        v = (piv['RBS'] - piv['WGCNA_hub']).values
        x = np.full(len(v), i) + np.random.default_rng(0).normal(0, 0.06, len(v))
        ax.scatter(x, v, s=34, facecolor=ACCENT, edgecolor=INK, linewidth=0.4,
                   alpha=0.85, zorder=3)
        ax.plot([i - 0.24, i + 0.24], [v.mean()] * 2, color=INK, lw=2, zorder=4)
    ax.axhline(0, color=INK, lw=1, ls='--')
    ax.set_xticks(range(len(sizes)))
    ax.set_xticklabels(['%d' % s for s in sizes])
    ax.set_xlabel('cohorts in the benchmark')
    ax.set_ylabel('stability minus connectivity\n(cross-fold agreement)')
    ax.set_title('(b) every subset', loc='left', fontsize=7)
    # 크기 5 는 부분집합이 하나라 부호가 뒤집힐 수 없다. 캡션과 같은 말을 적는다.
    # 왼쪽 아래에는 -0.106 점이 있어 글과 겹치므로 오른쪽 아래에 둔다.
    ax.text(0.97, 0.04, 'sign flips at every\nsize except five',
            transform=ax.transAxes, fontsize=6.24, color=MUTED, va='bottom', ha='right')

    # (c) AUROC difference is stable and grows
    ax = axes[2]
    means, los, his = [], [], []
    for s in sizes:
        piv = d[d['size'] == s].pivot(index='subset', columns='method',
                                      values='external_auroc').dropna()
        v = (piv['RBS'] - piv['WGCNA_hub']).values
        means.append(v.mean()); los.append(v.min()); his.append(v.max())
    ax.errorbar(sizes, means, yerr=[np.array(means) - np.array(los),
                                    np.array(his) - np.array(means)],
                fmt='o-', color=BLUE, ecolor=MUTED, capsize=4, lw=1.6, ms=6)
    ax.axhline(0, color=INK, lw=1, ls='--')
    ax.set_xticks(sizes)
    ax.set_xlabel('cohorts in the benchmark')
    ax.set_ylabel('stability minus connectivity\n(external AUROC)')
    # 131mm 폭에서는 세 번째 패널 제목이 오른쪽으로 잘렸다. 줄인다.
    ax.set_title('(c) AUROC advantage', loc='left', fontsize=7)
    fig.legend(*legend_handles, fontsize=6.2, loc='lower center', ncol=4,
               columnspacing=1.4, handlelength=1.2, frameon=False,
               bbox_to_anchor=(0.5, 0.004))
    fig.tight_layout(rect=(0, 0.09, 1, 0.99))
    # 패널 글자를 기본 자리에 두면 아래쪽 범례와 같은 줄에 찍힌다. 범례 위로 올린다.
    panels_below(fig, axes, pad=0.21, letter_y=0.105)
    save(fig, 'P1_cohort_sensitivity' + SUF, pdf=True)


# ---------------------------------------------------------------- F3
def f3_procurement():
    d = pd.read_csv('results/ieg_heterogeneity.tsv', sep='\t')
    d['is_biopsy'] = (d['procurement'] == 'biopsy')
    per = []
    for ds, sub in d.groupby('dataset'):
        if sub['is_biopsy'].nunique() < 2:
            continue
        per.append((ds, sub[sub['is_biopsy']]['ieg_z'].mean(),
                    sub[~sub['is_biopsy']]['ieg_z'].mean()))
    per.sort(key=lambda r: r[1] - r[2])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W, 3.0),
                                   gridspec_kw={'width_ratios': [1.3, 1]})
    y = np.arange(len(per))
    for i, (ds, b, nb) in enumerate(per):
        ax1.plot([b, nb], [i, i], color=MUTED, lw=1.2, zorder=1)
    ax1.scatter([r[1] for r in per], y, s=52, facecolor=ACCENT, edgecolor=INK,
                linewidth=0.5, zorder=3, label='biopsy')
    ax1.scatter([r[2] for r in per], y, s=52, facecolor=BLUE, edgecolor=INK,
                linewidth=0.5, zorder=3, label='nephrectomy / donor')
    ax1.axvline(0, color=INK, lw=0.8, ls=':')
    ax1.set_yticks(y)
    ax1.set_yticklabels([r[0].replace('_', ' ') for r in per], fontsize=6.63)
    ax1.set_xlabel('immediate-early module, z within dataset')
    # 상관계수를 축 안에 놓을 자리가 없다. 아래쪽은 긴 막대가, 오른쪽은 범례가 쓴다.
    # 제목의 둘째 줄로 올리면 어디와도 겹치지 않는다.
    ax1.set_title('(a) nine kidney datasets, ordered by effect\n'
                  'Spearman(procurement, module) = -0.75, p < 0.001', loc='left')
    ax1.legend(fontsize=6.24, loc='upper right')

    # (b) kidney vs heart
    hs = pd.read_csv('results/second_disease/gse5406.tsv', sep='\t')
    hs = hs[hs['test'].isin(['disease vs donor', 'ischaemic vs idiopathic'])]
    labels = ['kidney (KPMP)\nDKD vs non-biopsy', 'kidney (KPMP)\nDKD vs other biopsy',
              'heart\ndisease vs donor', 'heart\nischaemic vs idiopathic']
    vals = [-0.843, 0.082,
            float(hs[hs['test'] == 'disease vs donor']['ieg_g'].iloc[0]),
            float(hs[hs['test'] == 'ischaemic vs idiopathic']['ieg_g'].iloc[0])]
    cols = [ACCENT, GREY, ACCENT, GREY]
    ax2.barh(range(4), vals, color=cols, edgecolor=INK, linewidth=0.6)
    ax2.set_yticks(range(4)); ax2.set_yticklabels(labels, fontsize=6.24)
    ax2.invert_yaxis()
    ax2.axvline(0, color=INK, lw=0.8)
    ax2.set_xlabel("immediate-early module, Hedges' g")
    # 제목이 패널보다 길면 오른쪽에서 잘린다. 짧게, 주장 대신 내용으로.
    ax2.set_title('(b) by contrast and organ', loc='left')
    for i, v in enumerate(vals):
        # 긴 막대의 값표를 막대 밖에 두면 축 왼쪽 끝의 눈금 이름과 겹친다. 길면 안쪽에,
        # 짧으면 바깥에 둔다.
        inside = abs(v) > 0.5
        pad = 0.03 if (v < 0) == inside else -0.03
        ax2.text(v + pad, i, '%+.2f' % v, va='center',
                 ha='left' if pad > 0 else 'right',
                 fontsize=6.24, color='white' if inside else INK)
    fig.tight_layout()
    panels_below(fig, (ax1, ax2))
    save(fig, 'P3_procurement' + SUF, pdf=True)


# ---------------------------------------------------------------- F5
def f5_block_size():
    d = pd.read_csv('results/module_size/real_rbs_sweep.tsv', sep='\t')
    cfg = d[['n_dis', 'n_conf']].drop_duplicates().values.tolist()
    cfg.sort(key=lambda r: r[0] / r[1])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(JOURNAL_W, 3.0))
    x = np.arange(len(cfg)); w = 0.27
    cols = {'effect_size': GREY, 'RBS_real': ACCENT, 'connectivity': BLUE}
    nice = {'effect_size': 'effect size', 'RBS_real': 'stability (RBS)',
            'connectivity': 'connectivity comparator'}
    for j, m in enumerate(['effect_size', 'RBS_real', 'connectivity']):
        v, e = [], []
        for nd, nc in cfg:
            s = d[(d['n_dis'] == nd) & (d['n_conf'] == nc) & (d['method'] == m)]['confounder_frac']
            v.append(100 * s.mean()); e.append(100 * 1.96 * s.std(ddof=1) / np.sqrt(len(s)))
        ax1.bar(x + (j - 1) * w, v, w, yerr=e, color=cols[m], edgecolor=INK, linewidth=0.5,
                ecolor=INK, capsize=2, label=nice[m])
    ax1.set_xticks(x)
    ax1.set_xticklabels(['%d / %d' % (a, b) for a, b in cfg], fontsize=6.63)
    ax1.set_xlabel('disease block / confounder block (genes)')
    ax1.set_ylabel('% of the top 50 that is confounder')
    ax1.set_title('(a) confounder pickup by block-size ratio', loc='left')
    ax1.legend(fontsize=6.24)

    diffs, los, his = [], [], []
    for nd, nc in cfg:
        sub = d[(d['n_dis'] == nd) & (d['n_conf'] == nc)]
        piv = sub.pivot(index='rep', columns='method', values='confounder_frac')
        v = (piv['connectivity'] - piv['RBS_real']).dropna().values
        m_, se = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
        diffs.append(m_); los.append(m_ - 1.96 * se); his.append(m_ + 1.96 * se)
    ax2.errorbar(x, diffs, yerr=[np.array(diffs) - np.array(los),
                                 np.array(his) - np.array(diffs)],
                 fmt='o', color=INK, ecolor=MUTED, capsize=4, ms=7)
    ax2.axhline(0, color=INK, lw=1, ls='--')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['%d / %d' % (a, b) for a, b in cfg], fontsize=6.63)
    ax2.set_xlabel('disease block / confounder block (genes)')
    ax2.set_ylabel('connectivity minus stability')
    ax2.set_title('(b) the sign flips with the ratio', loc='left')
    # 왼쪽 위에는 20/200 점과 그 오차막대가 있다. 왼쪽 아래가 비어 있으므로 그리로 옮긴다.
    ax2.text(0.03, 0.04, 'above 0: connectivity picks MORE confounder\n'
                         'below 0: connectivity picks LESS',
             transform=ax2.transAxes, va='bottom', fontsize=6, color=MUTED)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    save(fig, 'P5_block_size', pdf=True)


# ---------------------------------------------------------------- F6
def f6_corrections():
    d = pd.read_csv('results/ruv_benchmark/benchmark.tsv', sep='\t')
    bg = float(d['spec_background_rate'].iloc[0])
    agg = (d[d['K'] == 50].groupby('arm')
             .agg(proc=('procurement_frac', 'mean'), sp=('spec_pass_frac', 'mean'),
                  auc=('external_auroc', 'mean')).reset_index())
    order = ['none', 'ComBat', 'RUVg_empirical', 'RUVg_ieg', 'SVA_like', 'IEG_resid']
    lbl = {'none': 'no correction', 'ComBat': 'ComBat (batch = cohort)',
           'RUVg_empirical': 'RUVg, empirical controls', 'RUVg_ieg': 'RUVg, module controls',
           'SVA_like': 'surrogate variables', 'IEG_resid': 'handling-score residualisation'}
    agg = agg.set_index('arm').reindex(order).reset_index()

    fig, axes = plt.subplots(1, 3, figsize=(W, 2.9))
    y = np.arange(len(order))
    cols = [GREY if a in ('none', 'ComBat', 'RUVg_empirical') else ACCENT for a in order]

    axes[0].barh(y, 100 * agg['proc'], color=cols, edgecolor=INK, linewidth=0.6)
    # 372pt 를 셋으로 나누면 패널이 1.7in 이다. 제목이 그보다 길면 옆 패널을 덮는다.
    axes[0].set_title('(a) procurement-driven\ngenes in top 50', loc='left')
    axes[0].set_xlabel('% of top 50 flagged')
    for i, v in enumerate(100 * agg['proc']):
        axes[0].text(v + 0.8, i, '%.0f' % v, va='center', fontsize=6.24)

    axes[1].barh(y, agg['sp'] / bg, color=cols, edgecolor=INK, linewidth=0.6)
    axes[1].axvline(1.0, color=INK, lw=1, ls='--')
    axes[1].set_title('(b) specificity pass-rate\n(x %.1f%% background)'
                      % (100 * bg), loc='left')
    axes[1].set_xlabel('specificity, x background')
    for i, v in enumerate(agg['sp'] / bg):
        axes[1].text(v + 0.1, i, '%.1fx' % v, va='center', fontsize=6.24)

    axes[2].barh(y, agg['auc'], color=cols, edgecolor=INK, linewidth=0.6)
    axes[2].set_xlim(0.8, 0.96)
    axes[2].set_title('(c) external AUROC\n(cost of correcting)', loc='left', fontsize=7)
    axes[2].set_xlabel('external AUROC')
    for i, v in enumerate(agg['auc']):
        axes[2].text(v + 0.002, i, '%.3f' % v, va='center', fontsize=6.24)

    for ax in axes:
        ax.set_yticks(y)
        ax.set_yticklabels([lbl[a] for a in order], fontsize=6.63)
        ax.invert_yaxis()
    for ax in axes[1:]:
        ax.set_yticklabels([])
    # 데이터 좌표에 두면 축선 위로 내려앉아 눈금과 겹치고 옆 패널까지 넘어간다.
    # 짧은 막대 두 개 오른쪽이 비어 있으므로 그 자리에 축 좌표로 붙인다.
    axes[0].text(0.98, 0.04, 'ComBat matches\nno correction',
                 transform=axes[0].transAxes, fontsize=6.0, color=MUTED,
                 ha='right', va='bottom', linespacing=1.3)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    panels_below(fig, axes)
    save(fig, 'P6_corrections' + SUF, pdf=True)


# P2 and P4 are the exploratory figures unchanged, so they are copied rather than redrawn.
# They were previously copied by hand, which meant a clean checkout produced four of the six
# manuscript figures and silently omitted two.
ADOPTED = [
    ('F3_auroc_vs_null', 'P2_null_control'),
    ('F10_singlecell_control', 'P4_singlecell_control'),
]


def adopt():
    for src, dst in ADOPTED:
        for ext in ('.png', '.pdf'):
            a = os.path.join(OUT, src + ext)
            b = os.path.join(OUT, dst + ext)
            if not os.path.exists(a):
                log('  ! %s missing - run `python run_all.py --stage figures` first'
                    % (src + ext))
                continue
            shutil.copyfile(a, b)
            log('  wrote %s' % (dst + ext))


def main():
    global W, SUF, PANEL_BELOW
    os.makedirs(OUT, exist_ok=True)
    if '--bib' in sys.argv:
        # BiB 본문 폭으로 다시 그린다. 확대되지 않으니 글자가 본문과 같은 크기로 앉는다.
        W, SUF, PANEL_BELOW = BIB_W, '_bib', True
        for fn in (f1_cohort_sensitivity, f3_procurement, f6_corrections):
            fn()
        return
    for fn in (p0_overview, f1_cohort_sensitivity, f3_procurement, f5_block_size,
               f6_corrections):
        try:
            fn()
        except Exception as e:
            log('  ! %s failed: %s' % (fn.__name__, e))
    adopt()


if __name__ == '__main__':
    main()
