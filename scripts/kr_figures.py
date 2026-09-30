#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""설명자료용 그림 — 글로만 된 설명을 그림으로 바꾼다.

설명자료가 글이 많다는 지적을 받았습니다. 말로 세 문단이 필요한 것 중 그림 한 장이면
끝나는 것이 둘 있습니다.

    LODO 흐름   코호트 5개에서 하나를 빼고 나머지로 고르고 뺀 것으로 채점하는 절차.
                "발굴 코호트가 하나냐" 는 질문이 나온 것은 이 그림이 없어서입니다.
    AUROC       환자와 대조군 점수를 짝지어 세는 방식. 식만 적으면 안 읽힙니다.

한글 글꼴을 쓰므로 matplotlib 에 맞는 글꼴을 찾아 지정합니다. 없으면 네모로 깨집니다.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np

OUT = 'results/figures'
INK = '#1a1a1a'
MUTED = '#6f6f6f'
ACCENT = '#c1440e'
BLUE = '#2a5d8f'
GREY = '#c9c4bc'
PALE = '#f2efea'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def korean_font():
    """설치된 한글 글꼴 중 하나를 고른다. 없으면 글자가 네모로 깨진다."""
    names = {f.name for f in fm.fontManager.ttflist}
    for cand in ('Malgun Gothic', 'NanumGothic', 'AppleGothic', 'Noto Sans KR',
                 'Gulim', 'Batang'):
        if cand in names:
            return cand
    return None


def setup():
    f = korean_font()
    if f:
        plt.rcParams['font.family'] = f
        log('  한글 글꼴: %s' % f)
    else:
        log('  ! 한글 글꼴을 찾지 못했습니다. 글자가 깨질 수 있습니다.')
    plt.rcParams.update({'axes.unicode_minus': False, 'savefig.dpi': 300,
                         'savefig.bbox': None, 'pdf.fonttype': 42})
    return f


def box(ax, x, y, w, h, text, face=PALE, edge=GREY, fs=8.5, bold=False, tc=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle='round,pad=0.012,rounding_size=0.02',
                                facecolor=face, edgecolor=edge, linewidth=1.0, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs,
            color=tc, fontweight='bold' if bold else 'normal', zorder=3, linespacing=1.5)


def arrow(ax, p, q, colour=INK, lw=1.2):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=11,
                                 linewidth=lw, color=colour, shrinkA=2, shrinkB=2, zorder=1))


# ---------------------------------------------------------------- LODO
def fig_lodo():
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')

    ax.text(1, 97, 'LODO — 코호트 하나를 통째로 빼고 나머지로 고른다',
            fontsize=12, fontweight='bold', va='top', color=INK)
    ax.text(1, 92, '한 번 돌릴 때의 절차. 빼는 코호트를 바꿔가며 5번 반복한다.',
            fontsize=8.5, va='top', color=MUTED)

    names = ['GSE30528', 'GSE96804', 'GSE104948', 'GSE142025', 'GSE294519']
    for i, nm in enumerate(names):
        held = (i == 4)
        box(ax, 2, 74 - i * 7.6, 20, 6.4, nm,
            face='#fbe6dd' if held else PALE,
            edge=ACCENT if held else GREY,
            tc=ACCENT if held else INK, bold=held, fs=8)
    ax.text(2, 82, '발굴 코호트 5개', fontsize=9, fontweight='bold', color=INK)
    ax.text(12, 33, '↑ 이번 회차에 빼둔 것', fontsize=7.6, color=ACCENT, ha='center')

    box(ax, 30, 55, 22, 14, '학습용 4개를 합친다\n(환자·대조 라벨 사용)', fs=8.2)
    box(ax, 30, 36, 22, 14, '유전자 50개를 고른다\n네 방법이 각자 고름', fs=8.2,
        face='#fff', edge=BLUE, tc=BLUE)
    box(ax, 30, 17, 22, 14, 'L2 로지스틱 회귀\n학습', fs=8.2)
    # 학습 줄기는 위에서 아래로. 빼둔 코호트는 학습에 들어가지 않고 채점 상자로만 간다.
    for y0, y1 in ((55, 50), (36, 31)):
        arrow(ax, (41, y0), (41, y1))
    arrow(ax, (22, 60), (30, 62))

    box(ax, 60, 36, 24, 14, '빼둔 코호트에서\n점수를 낸다', fs=8.2,
        face='#fbe6dd', edge=ACCENT, tc=ACCENT)
    box(ax, 60, 17, 24, 13, 'AUROC\n= 라벨을 얼마나 맞혔나', fs=8.2,
        face='#e9f2fb', edge=BLUE, tc=BLUE)

    # 빼둔 코호트 -> 채점 상자 (주황). 학습 쪽으로 가면 안 된다.
    ax.add_patch(FancyArrowPatch((22, 39), (60, 43), connectionstyle='arc3,rad=-0.22',
                                 arrowstyle='-|>', mutation_scale=11, linewidth=1.2,
                                 color=ACCENT, shrinkA=2, shrinkB=2, zorder=1))
    # 학습된 분류기 -> 채점 상자
    arrow(ax, (52, 24), (60, 24))
    arrow(ax, (72, 36), (72, 30), BLUE)

    ax.text(88, 43, '5번\n반복', fontsize=8.6, color=MUTED, ha='center', va='center',
            linespacing=1.5)
    ax.add_patch(FancyArrowPatch((84, 50), (84, 26), connectionstyle='arc3,rad=-0.9',
                                 arrowstyle='-|>', mutation_scale=10, color=MUTED,
                                 linewidth=1.0))

    ax.add_patch(Rectangle((1, 2), 96, 10, facecolor='#f6f4f0', edgecolor=GREY, lw=0.8))
    ax.text(49, 9.5, '중요 — 빼둔 코호트의 라벨은 유전자 선택에도 학습에도 쓰지 않는다.',
            fontsize=8.6, ha='center', va='top', color=INK, fontweight='bold')
    ax.text(49, 5, '같은 코호트 안에서 나누는 교차검증과 다르다. 다른 병원 자료에서도 되는지를 '
            '보려면 코호트째 빼야 한다.',
            fontsize=8, ha='center', va='top', color=MUTED)

    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, 'K1_lodo_flow.' + ext))
    plt.close(fig)
    log('  K1_lodo_flow')


# ---------------------------------------------------------------- AUROC
def fig_auroc():
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.3),
                                  gridspec_kw={'width_ratios': [1.15, 1]})
    for a in (ax, ax2):
        a.set_xlim(0, 100); a.set_ylim(0, 100); a.axis('off')

    ax.text(0, 97, 'AUROC를 세는 법', fontsize=11.5, fontweight='bold', va='top')
    ax.text(0, 91, '환자 한 명과 대조군 한 명을 짝지어, 환자 점수가 더 높은 짝을 센다.',
            fontsize=8.4, va='top', color=MUTED)

    cases = [('환자 A', 0.9), ('환자 B', 0.4)]
    ctrls = [('대조 C', 0.6), ('대조 D', 0.3)]
    for i, (nm, sc) in enumerate(cases):
        box(ax, 2, 70 - i * 11, 30, 9, '%s   점수 %.1f' % (nm, sc),
            face='#fbe6dd', edge=ACCENT, tc=ACCENT, fs=8.4)
    for i, (nm, sc) in enumerate(ctrls):
        box(ax, 2, 42 - i * 11, 30, 9, '%s   점수 %.1f' % (nm, sc),
            face='#e9f2fb', edge=BLUE, tc=BLUE, fs=8.4)

    pairs = [('A vs C', 0.9, 0.6, True), ('A vs D', 0.9, 0.3, True),
             ('B vs C', 0.4, 0.6, False), ('B vs D', 0.4, 0.3, True)]
    ax.text(40, 80, '짝 4개', fontsize=9, fontweight='bold', color=INK)
    for i, (nm, a_, b_, win) in enumerate(pairs):
        y = 70 - i * 11
        box(ax, 40, y, 40, 9,
            '%s :  %.1f %s %.1f   →  %s' % (nm, a_, '>' if win else '<', b_,
                                            '1점' if win else '0점'),
            face='#fff', edge=ACCENT if win else GREY,
            tc=ACCENT if win else MUTED, fs=8.2)
    ax.text(40, 20, 'AUROC = 3 / 4 = 0.75', fontsize=11, fontweight='bold', color=INK)
    ax.text(40, 13, '이긴 짝 3개 ÷ 전체 짝 4개', fontsize=8.2, color=MUTED)
    ax.text(40, 5, '점수가 같으면 0.5점으로 센다.', fontsize=8, color=MUTED)

    # 오른쪽: ROC 곡선과 넓이
    ax2.text(0, 97, '왜 곡선 아래 넓이라고 부르나', fontsize=11.5, fontweight='bold', va='top')
    # 설명글과 겹치지 않도록 그래프를 위로 올린다. 설명은 그림 하단에 따로 놓는다.
    ins = fig.add_axes([0.62, 0.28, 0.28, 0.52])
    t = np.linspace(0, 1, 200)
    for p, lab, col, lw in ((1.0, '0.50  동전 던지기', GREY, 1.2),
                            (2.6, '0.75  이 예시', ACCENT, 1.8),
                            (7.0, '0.95  좋은 모형', BLUE, 1.4)):
        ins.plot(t, t ** (1 / p), color=col, lw=lw, label=lab)
    ins.fill_between(t, t ** (1 / 2.6), t, color=ACCENT, alpha=0.10)
    ins.set_xlabel('대조군을 환자로 잘못 본 비율', fontsize=7.6)
    ins.set_ylabel('환자를 환자로 맞춘 비율', fontsize=7.6)
    ins.tick_params(labelsize=7)
    ins.set_xlim(0, 1); ins.set_ylim(0, 1)
    ins.legend(fontsize=7, loc='lower right', frameon=False)
    for s in ('top', 'right'):
        ins.spines[s].set_visible(False)
    fig.text(0.615, 0.15, '점수 기준선을 옮기며 두 비율을 찍으면 곡선이 된다.\n'
             '그 아래 넓이가 왼쪽에서 센 짝의 비율과 정확히 같다.',
             fontsize=8.2, va='top', color=MUTED, linespacing=1.6)

    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, 'K2_auroc.' + ext))
    plt.close(fig)
    log('  K2_auroc')


def fig_proteome():
    """단백질 층에서 후보가 어떻게 되었는가 — 두 코호트를 한 장에.

    말로 하면 "코호트가 둘이고 기술이 다르고 각각은 유의하지 않은데 합치면 유의하다"
    가 됩니다. 네 문장입니다. 그림이면 왼쪽에 유전자별 결과, 오른쪽에 배경 대비를
    나란히 두는 것으로 끝납니다.
    """
    import pandas as pd
    so = pd.read_csv('results/proteome_validation/protein_dkd_vs_control.tsv', sep=chr(9))
    ms = pd.read_csv('results/proteome_ms/protein_ms_dkd_vs_control.tsv', sep=chr(9))
    meta = pd.read_csv('results/proteome_meta/summary.tsv', sep=chr(9)).iloc[0]
    # 분석과 같은 필터를 써야 배경 비율이 본문 숫자와 맞는다.
    so = so[so['on_panel_and_measured']].rename(columns={'is_candidate': 'is_cand'})

    fig = plt.figure(figsize=(11.6, 4.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.30,
                          left=0.085, right=0.975, top=0.86, bottom=0.215)

    # ---- 왼쪽: 유전자별 단백질 효과크기
    ax = fig.add_subplot(gs[0, 0])
    rows = []
    for lab, t in (('SOMAscan', so), ('질량분석', ms)):
        c = t[t['is_cand']]
        for _, r in c.iterrows():
            rows.append((r['gene'], lab, r['prot_g'], r['prot_q'], r['rna_g']))
    d = pd.DataFrame(rows, columns=['gene', 'src', 'g', 'q', 'rna'])
    order = (d.groupby('gene')['q'].min().sort_values().index.tolist())
    ypos = {g: len(order) - 1 - i for i, g in enumerate(order)}
    for _, r in d.iterrows():
        y = ypos[r['gene']] + (0.17 if r['src'] == 'SOMAscan' else -0.17)
        sig = r['q'] < 0.05
        ax.plot([0, r['g']], [y, y], color=GREY, lw=1.0, zorder=1)
        ax.scatter([r['g']], [y], s=64 if sig else 34,
                   marker='o' if r['src'] == 'SOMAscan' else 's',
                   color=(ACCENT if sig else 'white'),
                   edgecolor=(ACCENT if sig else MUTED), lw=1.3, zorder=3)
    ax.axvline(0, color=INK, lw=0.9)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order[::-1], fontsize=9)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel('단백질 효과크기 (Hedges g, DKD - 대조)', fontsize=9.5)
    ax.set_title('가  후보 12개의 단백질 수준 결과', fontsize=10.5, loc='left', color=INK)
    ax.tick_params(labelsize=8.5)
    for s_ in ('top', 'right'):
        ax.spines[s_].set_visible(False)
    h = [plt.Line2D([], [], marker='o', ls='', color=ACCENT, ms=7, label='SOMAscan · q<0.05'),
         plt.Line2D([], [], marker='o', ls='', mfc='white', mec=MUTED, ms=6,
                    label='SOMAscan · 유의하지 않음'),
         plt.Line2D([], [], marker='s', ls='', color=ACCENT, ms=7, label='질량분석 · q<0.05'),
         plt.Line2D([], [], marker='s', ls='', mfc='white', mec=MUTED, ms=6,
                    label='질량분석 · 유의하지 않음')]
    ax.legend(handles=h, fontsize=7.6, loc='lower right', frameon=False, ncol=1)

    # ---- 오른쪽: 배경 대비와 결합 검정
    ax2 = fig.add_subplot(gs[0, 1])
    labs, cand, bg = [], [], []
    for lab, t in (('SOMAscan', so), ('질량분석', ms)):
        c = t[t['is_cand']]
        b = t[~t['is_cand']]
        labs.append(lab)
        cand.append(100 * (c['prot_q'] < 0.05).mean())
        bg.append(100 * (b['prot_q'] < 0.05).mean())
    x = np.arange(len(labs))
    ax2.bar(x - 0.19, bg, 0.34, color=GREY, edgecolor=MUTED, lw=0.8, label='배경 단백질')
    ax2.bar(x + 0.19, cand, 0.34, color=ACCENT, edgecolor=ACCENT, lw=0.8, label='우리 후보')
    for xi, (b_, c_) in enumerate(zip(bg, cand)):
        ax2.text(xi - 0.19, b_ + 1.6, '%.0f%%' % b_, ha='center', fontsize=8.5, color=MUTED)
        ax2.text(xi + 0.19, c_ + 1.6, '%.0f%%' % c_, ha='center', fontsize=8.5, color=ACCENT)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labs, fontsize=9.5)
    ax2.set_ylabel('q<0.05 에 도달한 비율', fontsize=9.5)
    ax2.set_ylim(0, max(cand) * 1.42)
    ax2.set_title('나  배경보다 자주 유의한가', fontsize=10.5, loc='left', color=INK)
    ax2.tick_params(labelsize=8.5)
    for s_ in ('top', 'right'):
        ax2.spines[s_].set_visible(False)
    ax2.legend(fontsize=8, frameon=False, loc='upper left')
    ax2.text(0.5, -0.155,
             ('각각은 유의하지 않다 (p=%.2f, %.2f).' + chr(10)
              + '층화해 합치면 CMH p=%.3f, 순열 p=%.3f')
             % (pd.read_csv('results/proteome_validation/summary.tsv',
                            sep=chr(9)).loc[0, 'fisher_p'],
                pd.read_csv('results/proteome_ms/summary.tsv',
                            sep=chr(9)).loc[0, 'fisher_p'],
                meta['cmh_p'], meta['perm_p_sig']),
             transform=ax2.transAxes, ha='center', va='top', fontsize=8.4,
             color=INK, linespacing=1.5)

    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, 'K3_proteome.' + ext))
    plt.close(fig)
    log('  K3_proteome')


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)
    setup()
    fig_lodo()
    fig_auroc()
    fig_proteome()
    return 0


if __name__ == '__main__':
    sys.exit(main())
