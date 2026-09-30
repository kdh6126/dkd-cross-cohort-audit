#!/usr/bin/env python
"""후보의 증거를 한 장으로 — 본문에서 길게 설명하는 대신 보여준다.

원고는 후보를 문단으로 설명하고 있었습니다. 후보가 30개고 증거 축이 여섯이라 문단으로는
독자가 "어느 후보가 왜 그 등급인가"를 따라올 수 없습니다.

두 판을 만듭니다.

    본문용   2·3등급 12개만. 주장에 걸리는 것은 이 둘뿐입니다 — 데이터가 강하거나
             새롭거나, 둘 중 하나는 되는 후보들입니다. 4등급 18개는 둘 다 아닌 잔여이고,
             30줄을 본문 폭에 밀어 넣으면 유전자명이 읽히지 않습니다.
    부록용   30개 전부. 잘라낸 것을 감추지 않도록 같은 그림의 완전판을 부록으로 냅니다.

    (b) 2x2  새로움 x 데이터 강도. 여기는 두 판 모두 30개 전부를 셉니다. 본문 그림이
             12줄만 보인다고 해서 세는 대상까지 줄이면, 한 그림 안에서 두 패널이 서로
             다른 모집단을 말하게 됩니다.

등급 기준은 master_candidate_table.py 가 정합니다. 여기서 다시 계산하지 않고 그 결과를
읽기만 합니다. 두 곳에서 따로 판정하면 언젠가 어긋납니다.

(b) 의 제목은 한 번 낮췄습니다. 처음에 'the two axes never meet' 이라고 썼는데 이건 일반
법칙처럼 읽힙니다. 실제로 말할 수 있는 것은 '지금 기준과 지금 30개에서 그 칸이 비어 있다'
뿐입니다.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import save, JOURNAL_W, INK, MUTED, ACCENT, GREY   # noqa: E402

SRC = 'results/candidates_v2/master_candidate_table.tsv'
PALE = '#f0ede8'
MAIN_TIERS = (2, 3)      # 본문에 싣는 등급


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def evidence(d):
    """증거 축 여섯. 각각 참/거짓/판정불가 셋 중 하나."""
    e = pd.DataFrame(index=d.index)
    e['above null'] = d['above_perm_ceiling']
    e['cross-diagnosis consistent'] = d['pattern_class'].isin(['A', 'B'])
    e['not control-driven'] = ~d['donor_driven']
    # 단일세포와 단백체는 조회되지 않은 유전자가 있다. 없는 것을 '아니오'로 칠하면
    # 검사해서 떨어진 것과 애초에 못 잰 것이 구분되지 않는다.
    sn = d['sn_top_foldchange']
    e['single-nucleus support'] = np.where(sn.isna(), np.nan, (sn.abs() > 1.5))
    pr = d['protein_adjP']
    # This is localization in human kidney, not a DKD-versus-control protein result.
    e['kidney region (not disease)'] = np.where(pr.isna(), np.nan, (pr < 0.05))
    e['prior DKD literature'] = ~d['lit_novel']
    return e


def matrix_panel(ax, d, e):
    # 칸 사이를 띄우고 테두리를 그리면 히트맵이 아니라 가로 막대 그래프처럼 읽힌다.
    # 칸을 맞붙이고 흰 선으로만 나누는 편이 표로 읽힌다.
    n, m = len(d), e.shape[1]
    for j in range(m):
        for i in range(n):
            v = e.iloc[i, j]
            if v != v:                      # NaN — 조회되지 않음
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor=PALE,
                                       edgecolor='white', linewidth=1.1))
                ax.plot([j + 0.36, j + 0.64], [i + 0.5, i + 0.5], color=MUTED, lw=0.7)
            else:
                ax.add_patch(Rectangle((j, i), 1, 1,
                                       facecolor=ACCENT if v else '#efece7',
                                       edgecolor='white', linewidth=1.1))
    bounds, labels = [], []
    for _, sub in d.groupby('tier_n'):
        bounds.append(int(sub.index.max()) + 1)
        labels.append((float(np.mean(sub.index.values)) + 0.5, d.loc[sub.index[0], 'tier']))
    for b in bounds[:-1]:
        ax.plot([-0.02, m + 0.02], [b, b], color=INK, lw=0.8)

    ax.set_xlim(0, m)
    ax.set_ylim(n, 0)
    ax.set_xticks(np.arange(m) + 0.5)
    # 열 폭이 6mm 남짓이다. 이름을 그대로 쓰면 옆 칸을 덮고, 세로로 세우면 세로 공간을
    # 1.1in 먹어 행렬이 눌린다. 번호만 올리고 이름은 그림 아래 범례로 내린다.
    ax.set_xticklabels([str(i + 1) for i in range(m)], fontsize=6.5)
    ax.xaxis.set_ticks_position('top')
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_yticklabels(d['gene'], fontsize=6.4 if n <= 15 else 6.0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    groups = ((0, 3, 'selection and control'),
              (3, 5, 'kidney context'),
              (5, 6, 'literature'))
    # 세로 라벨이 차지하는 높이를 축 좌표로 환산해 그 위에 둔다
    for lo, hi, label in groups:
        ax.text((lo + hi) / 2, -0.55, label, ha='center', va='bottom', fontsize=6.0,
                fontweight='bold', color=MUTED, clip_on=False,
                transform=ax.get_xaxis_transform())
    for y, lab in labels:
        ax.text(m + 0.35, y, 'tier ' + lab[0], rotation=270, va='center', ha='center',
                fontsize=6, color=MUTED)


def grid_panel(ax2, d_all):
    """2x2 는 언제나 후보 30개 전부를 센다."""
    grid = np.zeros((2, 2), dtype=int)
    for _, r in d_all.iterrows():
        grid[int(bool(r['data_strong'])), int(bool(r['lit_novel']))] += 1
    names = [['reported,\nweaker', 'unreported,\nweaker'],
             ['reported,\nstronger', 'unreported,\nstronger']]
    for a in (0, 1):
        for b in (0, 1):
            empty = grid[a, b] == 0
            ax2.add_patch(Rectangle((b, a), 1, 1, facecolor='#ffffff' if empty else PALE,
                                    edgecolor=ACCENT if empty else GREY,
                                    linewidth=1.4 if empty else 0.7,
                                    linestyle='--' if empty else '-'))
            ax2.text(b + 0.5, a + 0.30, str(grid[a, b]), ha='center', va='center',
                     fontsize=14 if empty else 12,
                     color=ACCENT if empty else INK, fontweight='bold')
            ax2.text(b + 0.5, a + 0.72, names[a][b], ha='center', va='center', fontsize=6.0,
                     color=MUTED)
    ax2.set_xlim(-0.03, 2.03)
    ax2.set_ylim(-0.03, 2.03)
    ax2.set_xticks([0.5, 1.5])
    ax2.set_xticklabels(['prior DKD\nliterature', 'no prior DKD\nliterature'], fontsize=6.0,
                        linespacing=1.3)
    ax2.set_yticks([0.5, 1.5])
    ax2.set_yticklabels(['weaker\ncurrent evidence', 'stronger\ncurrent evidence'], fontsize=6.0,
                        linespacing=1.3)
    for s in ax2.spines.values():
        s.set_visible(False)
    ax2.tick_params(length=0)


def draw(d_all, tiers, name, height, headline):
    d = (d_all[d_all['tier_n'].isin(tiers)] if tiers else d_all).reset_index(drop=True)
    e = evidence(d)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(JOURNAL_W, height),
                                  gridspec_kw={'width_ratios': [2.3, 1]})
    matrix_panel(ax, d, e)
    grid_panel(ax2, d_all)

    cols = list(e.columns)
    key = ('   '.join('%d %s' % (i + 1, c) for i, c in enumerate(cols[:3]))
           + chr(10) + '   '.join('%d %s' % (i + 4, c) for i, c in enumerate(cols[3:])))
    fig.text(0.005, 0.005,
             key + chr(10)
             + 'orange, criterion met   |   pale, criterion not met   |   dash, not assessed'
             + chr(10) + 'kidney-region evidence is not DKD disease evidence',
             fontsize=6.0, color=MUTED, ha='left', va='bottom')
    fig.tight_layout(rect=(0, 0.13, 1, 0.95))
    # 제목을 축에 붙이면 (a) 는 기울인 열 이름 때문에 더 올라가고, pad 로 맞추면 축 높이가
    # 깎여 유전자 줄이 겹친다. 제목을 축에서 떼어 그림 좌표에 같은 높이로 놓는다.
    for a, t in ((ax, '(a) evidence lines, tiers 2 and 3'),
                 (ax2, '(b) tier summary, all 30')):
        fig.text(a.get_position().x0, 0.965, t, fontsize=7.5, ha='left', va='bottom',
                 color=INK)
    # 주장을 이 표 안으로 한정한다. 그 칸이 비어 있다는 것은 지금 기준과 지금 30개에
    # 대한 사실이지 일반 법칙이 아니다.
    # 축 아래에 붙이면 두 줄짜리 눈금 이름과 겹치고, 더 내리면 캔버스 밖으로 잘린다.
    # 왼쪽 범례와 같은 높이의 오른쪽 끝이 유일하게 비어 있는 자리다.
    save(fig, name, pdf=True)
    return len(d)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    d = pd.read_csv(SRC, sep='\t')
    d['tier_n'] = d['tier'].str.slice(0, 1).astype(int)
    d = d.sort_values(['tier_n', 'mean_rank']).reset_index(drop=True)

    n_main = draw(d, MAIN_TIERS, 'P7_candidates', 3.2,
                  'The candidates that could support a claim, and what stands behind each')
    n_full = draw(d, None, 'P7_candidates_full', 5.0,
                  'All %d specificity-gated candidates against the six evidence lines' % len(d))

    log('  본문용 %d개 (등급 %s) · 부록용 %d개 (전부)'
        % (n_main, '/'.join(str(t) for t in MAIN_TIERS), n_full))
    log('  등급 판정은 master_candidate_table.py 의 결과를 읽기만 했습니다.')
    log('  1등급 %d · 2등급 %d · 3등급 %d · 4등급 %d'
        % tuple(int((d['tier_n'] == t).sum()) for t in (1, 2, 3, 4)))


if __name__ == '__main__':
    main()
