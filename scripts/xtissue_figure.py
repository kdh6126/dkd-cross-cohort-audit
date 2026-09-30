#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""세 조직 검정의 그림 한 장.

    a  코호트마다 즉시초기 모듈의 사례-대조 효과. 채취 비대칭 코호트와 대칭 코호트를 나누고,
       염증 점수를 뗀 값을 빈 표시로 겹친다. 비대칭 코호트는 조직을 가리지 않고 음수이고,
       대장의 양수는 절반쯤이 염증으로 설명된다는 것이 한 축에서 보인다. (P2, H5)
    b  보정 전후 RBS 상위 50 의 즉시초기 유전자 수. 교란이 공유된 신장에서만 크고, 보정이
       거기서만 걷어낸다. (P3)
    c  세 코호트 부분집합마다 RBS - WGCNA 재현성 차이. 신장과 간은 0 을 넘나들고 대장은
       넘지 않는다. 대장을 신장 크기로 줄인 반복이 있으면 옆에 둔다. (P1, H6)

값은 전부 results/xtissue 에서 읽는다.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import JOURNAL_W, BIB_W, save, INK, MUTED, ACCENT, BLUE, GREY  # noqa: E402

R = 'results/xtissue'
TNAME = {'kidney_dkd': 'Kidney, DKD', 'liver_masld': 'Liver, MASH', 'colon_uc': 'Colon, UC'}
ORDER = ['kidney_dkd', 'liver_masld', 'colon_uc']


def panel_a(ax):
    p2 = pd.read_csv(os.path.join(R, 'p2_ieg_by_cohort.tsv'), sep='\t')
    h5 = pd.read_csv(os.path.join(R, 'posthoc', 'h5_inflammation_adjustment.tsv'), sep='\t')
    d = p2.merge(h5[['cohort', 'case_coef_adjusted', 'case_coef_unadjusted']], on='cohort',
                 how='left')
    d = d[d['arms'] != 'unresolved']
    d['t'] = d['tissue'].map({t: i for i, t in enumerate(ORDER)})
    d = d.sort_values(['arms', 't', 'g_ieg'], ascending=[True, True, True])
    y = np.arange(len(d))[::-1]
    for yi, (_, r) in zip(y, d.iterrows()):
        col = ACCENT if r['arms'] == 'asymmetric' else BLUE
        ax.plot([r['ci_lo'], r['ci_hi']], [yi, yi], color=col, lw=0.8, alpha=0.6, zorder=1)
        ax.scatter(r['g_ieg'], yi, s=16, color=col, zorder=3)
        # 염증 점수를 뗀 값. 회귀 계수는 표준화 단위라 g 와 척도가 다르므로, 비율로 옮겨 적는다.
        if pd.notna(r['case_coef_adjusted']) and abs(r['case_coef_unadjusted']) > 1e-6:
            adj = r['g_ieg'] * r['case_coef_adjusted'] / r['case_coef_unadjusted']
            ax.scatter(adj, yi, s=16, facecolor='white', edgecolor=col, linewidth=0.8, zorder=4)
    ax.axvline(0, color=INK, lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(['%s  %s' % (TNAME[t].split(',')[0], c)
                        for t, c in zip(d['tissue'], d['cohort'])], fontsize=6.0)
    nasym = int((d['arms'] == 'asymmetric').sum())
    ax.axhline(y[nasym - 1] - 0.5, color=GREY, lw=0.7, ls='--')
    ax.text(0.02, 0.995, 'biopsy cases, surgical controls', transform=ax.transAxes,
            fontsize=6, color=ACCENT, va='top')
    # 구분선 바로 아래 왼쪽. 오른쪽에 두면 GSE48452 점과 겹쳤다.
    ax.set_xlim(-6.8, 6.8)
    ax.text(-6.6, y[nasym - 1] - 0.62, 'same procurement', fontsize=6, color=BLUE,
            va='top', ha='left')
    ax.set_xlabel("immediate-early module, Hedges' g (case minus control)"
                  + chr(10) + "hollow: rescaled by the inflammation-adjusted coefficient ratio")
    ax.set_title('a  immediate-early module by cohort', loc='left', fontsize=7,
                 fontweight='bold')
    # 좁은 폭에서는 범례가 자료 위로 올라온다. 채움은 측정값, 빈 표시는 환산값이라는
     # 설명은 축 라벨과 캡션에 있으므로 범례를 두지 않는다.


def panel_b(ax):
    vals = {}
    for t in ORDER:
        s = pd.read_csv(os.path.join(R, t, 'sweep', 'subsets.tsv'), sep='\t')
        f = s[s['size'] == s['size'].max()].set_index('method')['n_ieg_topK']
        vals[t] = (f['RBS'], f['RBS_orth'])
    x = np.arange(len(ORDER))
    ax.bar(x - 0.18, [vals[t][0] for t in ORDER], 0.34, color=ACCENT, label='uncorrected')
    ax.bar(x + 0.18, [vals[t][1] for t in ORDER], 0.34, color=GREY, label='handling score removed')
    for xi, t in zip(x, ORDER):
        ax.text(xi - 0.18, vals[t][0] + 0.15, '%.1f' % vals[t][0], ha='center', fontsize=6.0)
    ax.set_xticks(x)
    ax.set_xticklabels([TNAME[t] for t in ORDER], fontsize=6)
    ax.set_ylim(0, 10.5)
    ax.set_yticks(range(0, 9, 2))
    ax.set_ylabel('immediate-early genes in' + chr(10) + 'stability top 50')
    ax.set_title('b  in the stability top 50', loc='left', fontsize=7, fontweight='bold')
    ax.legend(fontsize=6.0, frameon=False)


def panel_c(ax):
    pts = []
    for t in ORDER:
        s = pd.read_csv(os.path.join(R, t, 'sweep', 'subsets.tsv'), sep='\t') \
            if t != 'kidney_dkd' else pd.read_csv('results/cohort_sensitivity/subsets.tsv', sep='\t')
        s = s[s['size'] == 3]
        pv = s.pivot(index='subset', columns='method', values='cross_fold_jaccard').dropna()
        pts.append((TNAME[t], (pv['RBS'] - pv['WGCNA_hub']).values))
    h6 = os.path.join(R, 'posthoc', 'h6', 'replicates.tsv')
    if os.path.exists(h6):
        h = pd.read_csv(h6, sep='\t')
        pts.append(('Colon' + chr(10) + 'subsampled', h['jaccard_diff'].values))
    conc = pd.read_csv(os.path.join(R, 'posthoc_effect_concordance.tsv'), sep='\t') \
        .set_index('tissue')['median_pairwise_rho_of_g']
    rng = np.random.default_rng(0)
    for i, (name, v) in enumerate(pts):
        col = [ACCENT if x < 0 else BLUE for x in v]
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=10, c=col, zorder=3,
                   edgecolor='white', linewidth=0.3)
    ax.axhline(0, color=INK, lw=0.7)
    ax.set_xticks(range(len(pts)))
    labels = [p[0] for p in pts]
    for i, t in enumerate(ORDER):
        labels[i] = TNAME[t].split(',')[0] + chr(10) + 'ρ %.2f' % conc[t]
    ax.set_xticklabels(labels, fontsize=6.0)
    ax.set_ylim(-0.14, 0.25)
    ax.set_ylabel('stability minus connectivity\n(cross-fold agreement)')
    ax.set_title('c  three-cohort verdicts', loc='left', fontsize=7, fontweight='bold')


def summary():
    """인쇄본용 요약판. 접근번호 24개 라벨 대신 조직 x 채취 설계 네 묶음으로 보이고 점을 키운다.

    코호트별 라벨이 붙은 원래 그림(P10_cross_tissue)은 보충 자료로 간다. 값은 같은 파일에서 읽는다.
    """
    p2 = pd.read_csv(os.path.join(R, 'p2_ieg_by_cohort.tsv'), sep='\t')
    h5 = pd.read_csv(os.path.join(R, 'posthoc', 'h5_inflammation_adjustment.tsv'), sep='\t')
    d = p2.merge(h5[['cohort', 'case_coef_adjusted', 'case_coef_unadjusted']], on='cohort', how='left')
    d = d[d['arms'] != 'unresolved']
    groups = [('kidney_dkd', 'asymmetric', 'Kidney\nbiopsy vs surgical'),
              ('liver_masld', 'asymmetric', 'Liver\nbiopsy vs surgical'),
              ('liver_masld', 'matched', 'Liver\nsame procurement'),
              ('colon_uc', 'matched', 'Colon\nsame procurement')]
    fig, axes = plt.subplots(1, 3, figsize=(BIB_W, 2.75),
                             gridspec_kw=dict(width_ratios=[1.3, 0.75, 1.25], wspace=0.6))
    ax = axes[0]
    rng = np.random.default_rng(1)
    labels = []
    for i, (t, arm, name) in enumerate(groups):
        g = d[(d['tissue'] == t) & (d['arms'] == arm)]
        y = len(groups) - 1 - i
        col = ACCENT if arm == 'asymmetric' else BLUE
        adj = g['g_ieg'] * g['case_coef_adjusted'] / g['case_coef_unadjusted']
        # 채운 점(측정값)과 빈 점(보정값)을 한 묶음 안에서 위아래 두 줄로 나눈다. 줄 간격이
        # 좁고 점이 크면 두 줄이 서로 닿아 겹쳐 보인다. 간격을 벌리고 흔들림과 점을 줄인다.
        jit = rng.uniform(-0.07, 0.07, len(g))
        ax.scatter(g['g_ieg'], y + 0.25 + jit, s=13, color=col, zorder=3, edgecolor='white',
                   linewidth=0.4, alpha=0.9)
        ax.scatter(adj, y - 0.25 + jit, s=13, facecolor='none', edgecolor=col, linewidth=0.8,
                   zorder=3, alpha=0.9)
        # 중앙값 막대는 점보다 아래(zorder 2)에 둔다. 위에 두면 점을 가린다.
        for yy, v in ((y + 0.25, g['g_ieg']), (y - 0.25, adj)):
            ax.plot([v.median()] * 2, [yy - 0.13, yy + 0.13], color=MUTED, lw=1.6, zorder=2)
        labels.append((y, '%s (%d)' % (name, len(g))))
    ax.plot([0, 0], [-0.6, len(groups) - 0.55], color=INK, lw=0.7)
    ax.set_yticks([y for y, _ in labels])
    ax.set_yticklabels([l for _, l in labels], fontsize=6.5)
    ax.set_ylim(-0.6, len(groups) + 0.05)
    ax.set_xlabel("immediate-early module, Hedges' g"
                  + chr(10) + "(hollow: rescaled by the inflammation-adjusted coefficient ratio)",
                  fontsize=6.5)
    # 패널 제목은 두지 않는다. 글자는 그림 아래, 설명은 캡션이 맡는다.
    ax.scatter([], [], s=22, color=MUTED, label="Hedges' g as measured")
    ax.scatter([], [], s=22, facecolor='white', edgecolor=MUTED,
               label='rescaled regression effect,' + chr(10) + 'inflammation-adjusted')
    # 범례를 축 안에 두면 맨 위 묶음(Kidney)과 붙는다. 축 위 바깥에 둔다.
    ax.legend(fontsize=6, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.0),
              ncol=2, handletextpad=0.2, columnspacing=0.8, borderaxespad=0.0)  # 6pt 이상 유지

    ax = axes[1]
    vals = {}
    for t in ORDER:
        s = pd.read_csv(os.path.join(R, t, 'sweep', 'subsets.tsv'), sep='\t')
        f = s[s['size'] == s['size'].max()].set_index('method')['n_ieg_topK']
        vals[t] = (f['RBS'], f['RBS_orth'])
    x = np.arange(len(ORDER))
    ax.bar(x - 0.19, [vals[t][0] for t in ORDER], 0.36, color=ACCENT, label='uncorrected')
    ax.bar(x + 0.19, [vals[t][1] for t in ORDER], 0.36, color=GREY, label='corrected')
    for xi, t in zip(x, ORDER):
        ax.text(xi - 0.19, vals[t][0] + 0.2, '%.1f' % vals[t][0], ha='center', fontsize=6.5)
    ax.set_xticks(x)
    ax.set_xticklabels([TNAME[t].split(',')[0] for t in ORDER], fontsize=6.5)
    ax.set_ylabel('immediate-early genes\nin stability top 50')

    ax.set_ylim(0, 11.5)
    ax.set_yticks(range(0, 9, 2))
    ax.legend(fontsize=6, frameon=False, loc='upper right', borderaxespad=0.1)

    ax = axes[2]
    pts = []
    for t in ORDER:
        s = pd.read_csv(os.path.join(R, t, 'sweep', 'subsets.tsv'), sep='\t') \
            if t != 'kidney_dkd' else pd.read_csv('results/cohort_sensitivity/subsets.tsv', sep='\t')
        s = s[s['size'] == 3]
        pv = s.pivot(index='subset', columns='method', values='cross_fold_jaccard').dropna()
        pts.append((t, (pv['RBS'] - pv['WGCNA_hub']).values))
    h = pd.read_csv(os.path.join(R, 'posthoc', 'h6', 'replicates.tsv'), sep='\t')
    pts.append(('colon_small', h['jaccard_diff'].values))
    conc = pd.read_csv(os.path.join(R, 'posthoc_effect_concordance.tsv'), sep='\t') \
        .set_index('tissue')['median_pairwise_rho_of_g']
    for i, (_, v) in enumerate(pts):
        w = 0.30 if len(v) > 20 else 0.18
        ax.scatter(i + rng.uniform(-w, w, len(v)), v, s=11 if len(v) > 20 else 15,
                   c=[ACCENT if z < 0 else BLUE for z in v], edgecolor='white', linewidth=0.25,
                   alpha=0.8, zorder=3)
    ax.axhline(0, color=INK, lw=0.7)
    ax.set_xticks(range(len(pts)))
    ax.set_xticklabels(['Kidney' + chr(10) + 'ρ %.2f' % conc['kidney_dkd'],
                        'Liver' + chr(10) + 'ρ %.2f' % conc['liver_masld'],
                        'Colon' + chr(10) + 'ρ %.2f' % conc['colon_uc'],
                        'Colon' + chr(10) + 'subsampled'], fontsize=6)
    ax.set_ylabel('stability minus connectivity\n(cross-fold agreement)')
    # 위쪽은 (a) 의 축 밖 범례 자리다.
    fig.subplots_adjust(left=0.19, right=0.975, top=0.86, bottom=0.30)
    for a, ch in zip(axes, 'abc'):
        bb = a.get_position()
        fig.text((bb.x0 + bb.x1) / 2, 0.025, ch, ha='center', va='bottom',
                 fontsize=8, fontweight='bold', color=INK)
    save(fig, 'P10_cross_tissue_summary', pdf=True)


def main():
    # 인쇄본 요약판(BiB 본문)과 코호트별 라벨판(BMC 본문, BiB 보충)을 한 번에 그린다
    summary()
    fig = plt.figure(figsize=(JOURNAL_W, 5.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1, 1],
                          hspace=0.55, wspace=0.55)
    panel_a(fig.add_subplot(gs[:, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[1, 1]))
    fig.subplots_adjust(left=0.17, right=0.985, top=0.95, bottom=0.10)
    save(fig, 'P10_cross_tissue', pdf=True)


if __name__ == '__main__':
    main()
