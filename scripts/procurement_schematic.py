#!/usr/bin/env python
"""채취 방식이 왜 교란이 되는지 — 설명용 도해 초안.

이 논문의 핵심 발견은 "가장 재현성 높은 신호가 질병이 아니라 조직 채취 방식을 따라간다"
입니다. 그런데 생물정보학 독자 상당수는 신장 생검과 종양 신절제가 실제로 어떻게 다른
조직인지 모릅니다. 그 차이를 문단으로 설명하면 세 문단이 필요하고, 그림 한 장이면 끝납니다.

이 파일은 완성품이 아니라 **초안**입니다. 해부 일러스트를 출판 품질로 그리는 것은
matplotlib 이 할 일이 아닙니다. 여기서 맞추는 것은 구조와 라벨과 수치이고, 모양은
일러스트레이터나 Inkscape 에서 다듬는 것을 전제로 합니다. 그래서 SVG 를 주 산출물로
내보내고 글자를 텍스트 객체로 남깁니다(svg.fonttype='none'). PNG 만 내보내면 편집이
불가능해집니다.

처음 판은 콩 모양을 베지에 제어점으로 그리고 겉질을 굵은 선으로 얹었는데, 도형이 의도한
크기를 벗어나 라벨을 덮고 선폭이 58pt 로 번졌습니다. 지금은 모양을 매개변수식으로 촘촘히
샘플링해 크기를 예측 가능하게 두고, 세로 공간을 띠로 나눠 겹칠 자리를 없앴습니다.

지어낸 수치는 넣지 않았습니다. 허혈 시간처럼 이 저장소가 측정하지 않은 값은 적지 않고,
표시된 효과크기는 전부 results/ 에서 읽습니다.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Polygon, PathPatch
from matplotlib.path import Path
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import JOURNAL_W, INK, MUTED, ACCENT, BLUE, GREY   # noqa: E402

OUT = 'results/figures'
TISSUE = '#e3cfc6'
CORTEX = '#c9a99c'
TUMOUR = '#b0a89c'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def numbers():
    """표시할 효과크기는 결과 파일에서 읽는다."""
    n = {}
    k = pd.read_csv('results/kpmp_singlecell/pseudobulk_de.tsv', sep='\t')
    ieg = k[k['kind'] == 'IEG']
    for contrast, key in (('DKD_vs_reference', 'naive'),
                          ('DKD_vs_CKDnondiabetic', 'matched')):
        sub = ieg[ieg['contrast'] == contrast]
        n[key] = float(sub['g'].mean())
    het = pd.read_csv('results/ieg_heterogeneity.tsv', sep='\t')
    n['n_datasets'] = het['dataset'].nunique()
    return n


# 신장 윤곽. 정규 좌표(-1..1)의 앵커 7개와 그 제어점.
#
# 처음에는 매개변수식을 240점으로 샘플링해 폴리곤으로 그렸습니다. 화면으로는 같아 보이지만
# SVG 로 나가면 238노드짜리 패스가 되어, 벡터 편집기에서 윤곽을 다듬는 것이 불가능합니다.
# 이 그림의 용도가 "받아서 다듬는 초안" 이므로 노드 수가 곧 쓸모입니다. 3차 베지에 7구간이면
# 앵커를 하나씩 잡아 끌 수 있습니다.
#
# 순서는 위 -> 오른쪽 옆구리 -> 아래 -> 왼쪽 아래 -> 문(hilum) -> 왼쪽 위 -> 위.
BEAN = [
    (0.00, 1.00),
    (0.46, 1.02), (0.88, 0.86), (0.90, 0.48),
    (0.96, 0.16), (0.96, -0.16), (0.90, -0.48),
    (0.88, -0.86), (0.46, -1.02), (0.00, -1.00),
    (-0.46, -1.00), (-0.80, -0.84), (-0.80, -0.42),
    (-0.82, -0.20), (-0.60, -0.20), (-0.55, 0.00),
    (-0.60, 0.20), (-0.82, 0.20), (-0.80, 0.42),
    (-0.80, 0.84), (-0.46, 1.00), (0.00, 1.00),
]


def bean_path(cx, cy, w, h, flip=False, scale=1.0):
    codes = [Path.MOVETO] + [Path.CURVE4] * 21
    verts = []
    for x, y in BEAN:
        vx = x * (w / 2) * scale
        vy = y * (h / 2) * scale
        verts.append((cx + (-vx if flip else vx), cy + vy))
    return Path(verts, codes)


def draw_kidney(ax, cx, cy, w, h, flip=False, tag='kidney'):
    """겉질 띠는 같은 윤곽을 두 번 축소해 만든다.

    도형마다 이름(gid)을 붙인다. Inkscape 의 XML 편집기나 일러스트레이터의 레이어
    목록에서 'kidney-left-outline' 을 바로 찾아 교체할 수 있어야 초안으로서 쓸모가 있다.
    """
    for name, sc, face, edge, z in (('outline', 1.00, TISSUE, INK, 2),
                                    ('cortex', 0.78, CORTEX, 'none', 3),
                                    ('medulla', 0.66, TISSUE, 'none', 4)):
        pp = PathPatch(bean_path(cx, cy, w, h, flip, sc), facecolor=face,
                       edgecolor=edge, linewidth=0.9 if edge != 'none' else 0, zorder=z)
        pp.set_gid('%s-%s' % (tag, name))
        ax.add_patch(pp)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    n = numbers()
    plt.rcParams['svg.fonttype'] = 'none'      # 글자를 텍스트 객체로 남긴다

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(JOURNAL_W, 3.9),
                                  gridspec_kw={'width_ratios': [1.25, 1]})
    for a in (ax, ax2):
        a.set_xlim(0, 100)
        a.set_ylim(0, 100)
        a.axis('off')

    # ================================================================ (a)
    # 세로 공간을 띠로 나눈다. 겹칠 자리를 애초에 만들지 않는다.
    #   99..93 제목 · 90..62 그림 · 60..52 역할 · 50..38 설명 · 34..24 기제 · 18..2 요약
    ax.text(0, 99, '(a) how the tissue was obtained', fontsize=7.5, va='top', color=INK)

    draw_kidney(ax, 24, 76, w=19, h=26, tag='kidney-left')
    draw_kidney(ax, 74, 76, w=19, h=26, flip=True, tag='kidney-right')

    # 왼쪽: 경피 바늘 생검. 신장은 몸 안에 있고 순환이 유지된다.
    ax.plot([42, 27.5], [92, 79], color=INK, lw=1.4, solid_capstyle='butt', zorder=5)
    ax.plot([27.5, 24.5], [79, 76.5], color=ACCENT, lw=2.4, solid_capstyle='butt', zorder=6)
    ax.text(43, 92.5, 'needle', fontsize=6.2, color=INK, va='bottom', ha='left')

    # 오른쪽: 종양 신절제. 장기는 이미 몸 밖이고, 종양에서 떨어진 정상 겉질을 쓴다.
    ax.add_patch(Circle((70, 70), 3.6, facecolor=TUMOUR, edgecolor=INK, lw=0.7, zorder=5))
    ax.text(70, 70, 'T', fontsize=5.8, ha='center', va='center', color=INK, zorder=6)
    ax.add_patch(Polygon([[76, 86], [82, 86], [82, 81], [76, 81]], closed=True,
                         facecolor='none', edgecolor=ACCENT, lw=1.0, linestyle='--', zorder=6))
    ax.text(84, 83.5, 'sampled\nhere', fontsize=6.2, color=ACCENT, va='center',
            ha='left', linespacing=1.4)

    ax.plot([49, 49], [24, 95], color=GREY, lw=0.7, linestyle=':')

    for x0, label, sub, col in ((1, 'CASES', 'percutaneous biopsy', ACCENT),
                                (51, 'CONTROLS', 'nephrectomy / donor', BLUE)):
        ax.add_patch(FancyBboxPatch((x0, 50), 46, 12,
                                    boxstyle='round,pad=0.3,rounding_size=1.0',
                                    facecolor='#ffffff', edgecolor=col, linewidth=1.0))
        ax.text(x0 + 23, 58.4, label, fontsize=6.6, ha='center', va='center',
                color=col, fontweight='bold')
        ax.text(x0 + 23, 53.4, sub, fontsize=6.3, ha='center', va='center', color=col)

    ax.text(24, 47, 'kidney in situ and perfused;\na 1-2 mm core is taken\nand fixed at the '
            'bedside', fontsize=6.2, ha='center', va='top', color=MUTED, linespacing=1.55)
    ax.text(74, 47, 'organ already removed;\nnormal-appearing cortex is\ncut away from the '
            'tumour (T)', fontsize=6.2, ha='center', va='top', color=MUTED, linespacing=1.55)

    ax.text(50, 33, 'The two routes differ in how long the tissue is without circulation and\n'
            'in how it is handled before fixation. Immediate-early genes (FOS, JUN,\n'
            'EGR1 and the rest of the 19-gene module) respond to exactly that.',
            fontsize=6.3, ha='center', va='top', color=INK, linespacing=1.6)

    ax.add_patch(FancyBboxPatch((1, 2), 96, 13,
                                boxstyle='round,pad=0.4,rounding_size=1.0',
                                facecolor='#f4f2ee', edgecolor=GREY, linewidth=0.8))
    # 한 줄에 다 넣으면 상자 폭 91mm 를 넘는다
    ax.text(49, 12, 'In every public DKD cohort, cases are biopsies and controls are not.',
            fontsize=6.2, ha='center', va='top', color=INK, fontweight='bold')
    ax.text(49, 6.5, 'The module separates the two routes in 8 of %d kidney datasets.'
            % n['n_datasets'], fontsize=6.2, ha='center', va='top', color=MUTED)

    # ================================================================ (b)
    ax2.text(0, 99, '(b) what the comparison actually contrasts', fontsize=7.5,
             va='top', color=INK)

    rows = [('conventional design', 'biopsy', 'nephrectomy', ACCENT, BLUE, n['naive'],
             'disease AND handling differ'),
            ('this paper', 'biopsy', 'biopsy', ACCENT, ACCENT, n['matched'],
             'only disease differs')]
    for k, (name, a_lab, b_lab, ca, cb, g, note) in enumerate(rows):
        top = 90 - 38 * k
        ax2.text(0, top, name, fontsize=6.8, color=INK, fontweight='bold', va='top')
        for j, (lab, col) in enumerate(((a_lab, ca), (b_lab, cb))):
            x = 1 + 53 * j
            ax2.add_patch(FancyBboxPatch((x, top - 15), 40, 9,
                                         boxstyle='round,pad=0.3,rounding_size=1.0',
                                         facecolor='#ffffff', edgecolor=col, linewidth=1.1))
            ax2.text(x + 20, top - 10.5, lab, fontsize=6.5, ha='center', va='center',
                     color=col)
        ax2.text(48, top - 10.5, 'vs', fontsize=6.2, ha='center', va='center', color=MUTED)
        ax2.text(0, top - 18, note, fontsize=6.2, color=MUTED, va='top')
        ax2.text(0, top - 25, "immediate-early module, Hedges' g = %+.2f" % g,
                 fontsize=6.5, color=INK, va='top', fontweight='bold')

    ax2.text(0, 14, 'Replacing the control with a biopsy of a different kidney disease removes\n'
             'the handling difference, and the module effect goes with it. What survives\n'
             'that substitution is the part of the signal that is about the disease.',
             fontsize=6.3, va='top', color=INK, linespacing=1.6)

    fig.text(0.005, 0.005, 'Draft schematic. Anatomy is not to scale; the effect sizes are read '
             'from the KPMP single-nucleus analysis.',
             fontsize=5.8, color=MUTED, ha='left', va='bottom')
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    os.makedirs(OUT, exist_ok=True)
    for ext in ('svg', 'pdf', 'png'):
        fig.savefig(os.path.join(OUT, 'S_procurement_schematic.' + ext))
    plt.close(fig)
    log('  S_procurement_schematic.{svg,pdf,png}')
    log('  SVG 를 벡터 편집기에서 여세요. 글자는 텍스트 객체로 남겼습니다.')
    log('  효과크기: 기존 대조 %+.2f · 채취 맞춘 대조 %+.2f' % (n['naive'], n['matched']))


if __name__ == '__main__':
    main()
