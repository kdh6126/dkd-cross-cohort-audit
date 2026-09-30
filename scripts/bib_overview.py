#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""BiB 판 Fig. 1 — 프로토콜이 먼저 보이는 개요.

BMC 판의 개요(P0_overview)는 DKD 분석 흐름을 그린다. BiB 판은 문제 해결 프로토콜로
읽혀야 하므로 네 단으로 나눈다.

    A  재사용 가능한 감사 절차      어느 질환에서나 같은 순서로 적용한다
    B  개발 사례: 당뇨병성 신장질환   각 검사가 무엇을 잡아냈는가
    C  적용 검정: 간 · 대장          데이터를 보기 전에 고정한 예측의 판정
    D  적용 조건과 한계              어디까지 옮겨 쓸 수 있는가

수치는 전부 results/ 에서 읽는다. 본문과 어긋나면 그림이 먼저 틀리기 때문이다.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import save, BIB_W, INK, MUTED, ACCENT, BLUE   # noqa: E402

HAIR = '#b8b3aa'
BAND = '#f4f2ee'
NL = chr(10)
MINUS = '−'
FS, FB, FH = 6.0, 6.4, 7.0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def numbers():
    n = {}
    coh = pd.read_csv('results/cohort_table.tsv', sep='\t')
    n['n_cohort'], n['n_disc'] = len(coh), int((coh['role'] == 'discovery').sum())
    sub = pd.read_csv('results/cohort_sensitivity/subsets.tsv', sep='\t')
    n['n_subsets'] = int(sub['subset'].nunique())
    p = sub[sub['size'] == 3].pivot(index='subset', columns='method', values='cross_fold_jaccard')
    d = p['RBS'] - p['WGCNA_hub']
    n['k_lo'], n['k_hi'] = float(d.min()), float(d.max())
    b = pd.read_csv('results/cohort_sensitivity/by_size.tsv', sep='\t')
    n['auroc5'] = float(b[b['size'] == 5]['auroc_diff_mean'].iloc[0])
    het = pd.read_csv('results/ieg_heterogeneity.tsv', sep='\t')
    n['n_datasets'] = int(het['dataset'].nunique())
    gen = pd.read_csv('results/ckd_generalization/per_diagnosis.tsv', sep='\t')
    n['n_dx'] = len(gen)
    n['shrink_lo'], n['shrink_hi'] = 100 * float(gen['shrink'].min()), 100 * float(gen['shrink'].max())
    ruv = pd.read_csv('results/ruv_benchmark/benchmark.tsv', sep='\t')
    r = ruv[ruv['K'] == 50].groupby('arm')[['procurement_frac']].mean()
    n['proc_pre'] = 100 * float(r.loc['none', 'procurement_frac'])
    n['proc_post'] = 100 * float(r.loc['IEG_resid', 'procurement_frac'])
    # 세 조직 검정
    ieg = {}
    for t in ('kidney_dkd', 'liver_masld', 'colon_uc'):
        s = pd.read_csv('results/xtissue/%s/sweep/subsets.tsv' % t, sep='\t')
        f = s[s['size'] == s['size'].max()].set_index('method')['n_ieg_topK']
        ieg[t] = float(f['RBS'])
    n['ieg'] = ieg
    lv = pd.read_csv('results/xtissue/liver_masld/sweep/subsets.tsv', sep='\t')
    lv = lv[lv['size'] == 3].pivot(index='subset', columns='method', values='cross_fold_jaccard')
    dl = lv['RBS'] - lv['WGCNA_hub']
    n['lv_lo'], n['lv_hi'] = float(dl.min()), float(dl.max())
    cl = pd.read_csv('results/xtissue/colon_uc/sweep/subsets.tsv', sep='\t')
    cl = cl[cl['size'] == 3].pivot(index='subset', columns='method', values='cross_fold_jaccard')
    dc = cl['RBS'] - cl['WGCNA_hub']
    n['cl_lo'], n['cl_hi'] = float(dc.min()), float(dc.max())
    conc = pd.read_csv('results/xtissue/posthoc_effect_concordance.tsv', sep='\t') \
        .set_index('tissue')['median_pairwise_rho_of_g']
    n['rho'] = {k: float(conc[k]) for k in ('kidney_dkd', 'liver_masld', 'colon_uc')}
    p2 = pd.read_csv('results/xtissue/p2_ieg_by_cohort.tsv', sep='\t')
    asym = p2[p2['arms'] == 'asymmetric']
    n['asym_neg'], n['asym_n'] = int((asym['g_ieg'] < 0).sum()), len(asym)
    col = p2[(p2['tissue'] == 'colon_uc')]
    n['colon_pos'], n['colon_n'] = int((col['g_ieg'] > 0).sum()), len(col)
    return n


def band(ax, y, h, letter, title):
    ax.add_patch(Rectangle((0.6, y), 98.8, h, facecolor=BAND, edgecolor=HAIR,
                           linewidth=0.6, zorder=1))
    ax.text(2.0, y + h - 1.6, letter, fontsize=FH, fontweight='bold', color=ACCENT,
            ha='left', va='top', zorder=4)
    ax.text(5.0, y + h - 1.6, title, fontsize=FH, fontweight='bold', color=INK,
            ha='left', va='top', zorder=4)


def box(ax, x, y, w, h, head, body, ec=INK, fc='#ffffff'):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, linewidth=0.7, zorder=3))
    ax.text(x + w / 2, y + h - 1.4, head, fontsize=FB, fontweight='bold', color=INK,
            ha='center', va='top', zorder=4, linespacing=1.35)
    ax.text(x + w / 2, y + h - 1.4 - 3.0 * (head.count(NL) + 1), body, fontsize=FS,
            color=MUTED, ha='center', va='top', zorder=4, linespacing=1.45)


def arrow(ax, x0, y0, x1, y1, colour=INK):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=7,
                                 linewidth=0.8, color=colour, shrinkA=0, shrinkB=0, zorder=2))


def main():
    n = numbers()
    fig = plt.figure(figsize=(BIB_W, 4.9))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')

    # ---------------------------------------------------------------- A 절차
    band(ax, 72.5, 26.0, 'A', 'Audit protocol, applied in this order to any cross-cohort benchmark')
    xs = [3.0, 27.0, 51.0, 75.0]
    w, h = 22.0, 15.0
    heads = ['1  Audit the input',
             '2  Vary the cohorts',
             '3  Vary the control',
             '4  Benchmark the corrections']
    bodies = ['patient overlap across' + NL + 'accessions; procurement' + NL + 'route of each arm',
              'repeat the benchmark on' + NL + 'every admissible subset;' + NL + 'report the interval',
              'replace the control with' + NL + 'one obtained by the same' + NL + 'procedure',
              'what each correction' + NL + 'removes and what it' + NL + 'costs in AUROC']
    for x, hd, bd in zip(xs, heads, bodies):
        box(ax, x, 74.0, w, h, hd, bd)
    for x in xs[:-1]:
        arrow(ax, x + w, 81.5, x + 24.0, 81.5)
    # 띠 A 와 B 사이의 빈 줄에 둔다. 띠 안에 두면 네 번째 상자의 밑변과 겹쳤다.
    ax.text(98.0, 71.2, 'output: reporting checklist (Table 4)', fontsize=FS, color=ACCENT,
            ha='right', va='center', zorder=4)

    # ---------------------------------------------------------------- B 개발 사례
    band(ax, 41.0, 29.0, 'B', 'Development case: diabetic kidney disease (%d cohorts, %d discovery)'
         % (n['n_cohort'], n['n_disc']))
    bb = [('%d subsets' % n['n_subsets'],
           'reproducibility verdict' + NL + '%s%.3f to +%.3f;' % (MINUS, abs(n['k_lo']), n['k_hi'])
           + NL + 'mean AUROC +%.3f at five' % n['auroc5']),
          ('procurement',
           'immediate-early module' + NL + 'separates biopsy from' + NL
           + 'surgical in 8 of %d sets' % n['n_datasets']),
          ('control swap',
           'a biopsy control retains' + NL + '%.0f–%.0f%% of the signal' % (n['shrink_lo'], n['shrink_hi'])
           + NL + 'in %d diagnoses' % n['n_dx']),
          ('corrections',
           'procurement-driven genes' + NL + '%.0f%% → %.0f%%; batch by' % (n['proc_pre'], n['proc_post'])
           + NL + 'cohort is inert')]
    for x, (hd, bd) in zip(xs, bb):
        box(ax, x, 42.5, w, 15.0, hd, bd, ec=BLUE)
    for x in xs[:-1]:
        arrow(ax, x + w, 50.0, x + 24.0, 50.0, colour=BLUE)

    # ---------------------------------------------------------------- C 적용 검정
    band(ax, 14.5, 24.5, 'C', 'Application test: liver and colon, predictions hash-fixed before the data were read')
    cc = [('P1  composition',
           'liver %s%.3f to +%.3f' % (MINUS, abs(n['lv_lo']), n['lv_hi']) + NL
           + 'colon +%.3f to +%.3f' % (n['cl_lo'], n['cl_hi']) + NL
           + 'held in liver, failed in colon'),
          ('P2  module by cohort',
           '%d of %d biopsy-vs-surgical' % (n['asym_neg'], n['asym_n']) + NL
           + 'cohorts negative; %d of %d colon' % (n['colon_pos'], n['colon_n']) + NL
           + 'cohorts positive: failed'),
          ('P3  reproducibility',
           'top 50 holds %.1f genes in' % n['ieg']['kidney_dkd'] + NL
           + 'kidney, %.1f in liver, %.1f in colon' % (n['ieg']['liver_masld'], n['ieg']['colon_uc'])
           + NL + 'held')]
    for x, (hd, bd) in zip([3.0, 36.0, 69.0], cc):
        box(ax, x, 16.0, 28.0, 15.0, hd, bd, ec=ACCENT)

    # ---------------------------------------------------------------- D 조건과 한계
    band(ax, 1.0, 12.0, 'D', 'Where the findings apply, and where they do not')
    dd = ['instability went with low cross-cohort effect' + NL
          + 'concordance (ρ %.2f colon, %.2f kidney, %.2f liver),' % (
              n['rho']['colon_uc'], n['rho']['kidney_dkd'], n['rho']['liver_masld']) + NL
          + 'not with cohort size',
          'the module also rises with inflammation, so' + NL
          + 'correcting on it is tissue-conditional; a control' + NL
          + 'by the same procedure is preferred where available',
          'three tissues, one disease per tissue; digests are' + NL
          + 'local, so they fix the files, not an external date;' + NL
          + 'no new biomarker is claimed']
    for x, t in zip([3.0, 36.0, 69.0], dd):
        ax.text(x, 9.2, t, fontsize=FS, color=INK, ha='left', va='top', zorder=4,
                linespacing=1.5)

    save(fig, 'P0_overview_bib', pdf=True)


if __name__ == '__main__':
    main()
