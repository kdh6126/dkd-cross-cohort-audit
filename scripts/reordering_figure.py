#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""대조군을 바꾸면 순위의 꼭대기가 다시 짜인다 — 그림.

3.6절이 말로 하면 네 문단입니다. 효과크기가 줄어든다는 것, 그런데 고르게 줄지 않는다는
것, 상위 50개가 한 개만 겹친다는 것, 그리고 어느 쪽이 오르고 어느 쪽이 내리는지.
왼쪽에 한 대조의 순위, 오른쪽에 다른 대조의 순위를 두고 선으로 이으면 한 장이 됩니다.

    왼쪽   사례 대 비생검 대조(신절제·공여자).  조달 방식이 사례와 다르다.
    오른쪽 사례 대 다른 진단의 생검.            조달 방식이 사례와 같다.

같은 사례 집합을 두 대조군과 각각 비교한 것입니다. 같은 사람에게서 두 방식으로 채취해
짝지은 자료가 아닙니다. 제목에 그것을 적습니다.

순위는 로그로 놓습니다. 1위와 20위의 차이가 5000위와 5500위의 차이보다 훨씬 중요한데,
선형으로 그리면 꼭대기가 눌려 아무것도 안 보입니다.

그리는 유전자는 후보 목록(보정 전 30 + 보정 후 30)에 한 번이라도 오른 것들입니다.
9,900개를 다 그리면 잉크뿐입니다. 다만 상위 50 겹침은 전체 9,900개로 계산한 값을
오른쪽에 따로 적어, 고른 것과 계산한 것을 구분합니다.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                  # noqa: E402
import numpy as np                                               # noqa: E402
import pandas as pd                                              # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import save, JOURNAL_W, BIB_W, INK, MUTED, ACCENT, BLUE, GREY   # noqa: E402

SHIFT = 'results/comparator_reordering/candidate_rank_shift.tsv'
SUMM = 'results/comparator_reordering/summary.tsv'
TOPS = 'results/comparator_reordering/top_overlap.tsv'
NAME = 'P9_reordering'
LABEL_N = 6          # 위아래로 이름을 붙일 개수


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    d = pd.read_csv(SHIFT, sep='\t').dropna(subset=['gene'])
    m = pd.read_csv(SUMM, sep='\t').iloc[0]
    tops = pd.read_csv(TOPS, sep='\t')
    d = d.sort_values('moved', ascending=False).reset_index(drop=True)

    bib = '--bib' in sys.argv
    width = BIB_W if bib else JOURNAL_W
    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(width, 4.3), gridspec_kw=dict(width_ratios=[2.05, 1.0], wspace=0.42))

    # ---------------------------------------------------------------- 왼쪽: 순위 이동
    x0, x1 = 0.0, 1.0
    for _, r in d.iterrows():
        up = r['moved'] > 0
        ax.plot([x0, x1], [r['rank_naive'], r['rank_matched']],
                color=(BLUE if up else ACCENT), lw=0.9, alpha=0.55,
                solid_capstyle='round', zorder=2)
        ax.scatter([x0, x1], [r['rank_naive'], r['rank_matched']], s=7,
                   color=(BLUE if up else ACCENT), zorder=3, linewidths=0)

    def place(vals, gap):
        """로그 좌표에서 라벨이 겹치지 않게 최소 간격만큼 밀어낸다.

        값 위치에 그대로 찍으면 상위권에서 이름이 서로 덮어 읽히지 않는다. 순서는
        지키면서 간격만 벌린다.
        """
        out = list(vals)
        for i in range(1, len(out)):
            if out[i] - out[i - 1] < gap:
                out[i] = out[i - 1] + gap
        return out

    # 오르는 것은 오른쪽 끝에, 내리는 것은 왼쪽 끝에 이름을 붙인다. 각각 그쪽이
    # 눈에 띄는 자리이기 때문이다.
    for sub, xoff, ha, col, key in (
            (d.head(LABEL_N), x1 + 0.06, 'left', BLUE, 'rank_matched'),
            (d.tail(LABEL_N), x0 - 0.06, 'right', ACCENT, 'rank_naive')):
        sub = sub.sort_values(key)
        ys = place([np.log10(v) for v in sub[key]], 0.21)
        for (_, r), y in zip(sub.iterrows(), ys):
            ax.plot([xoff - (0.04 if ha == 'left' else -0.04), xoff],
                    [r[key], 10 ** y], color=col, lw=0.5, alpha=0.55, zorder=1)
            ax.text(xoff, 10 ** y, r['gene'], va='center', ha=ha,
                    fontsize=6.3, color=col)

    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_xlim(-0.50, 1.50)
    ax.set_ylim(9900, 1)
    ax.set_yticks([1, 10, 100, 1000, 9900])
    ax.set_yticklabels(['1', '10', '100', '1,000', '9,900'], fontsize=7)
    ax.set_ylabel('rank by $|g|$  (1 = strongest)', fontsize=7.6)
    ax.set_xticks([x0, x1])
    ax.set_xticklabels(['control obtained\ndifferently from cases\n(nephrectomy, donor)',
                        'control obtained\nthe same way as cases\n(other-diagnosis biopsy)'],
                       fontsize=7)
    ax.tick_params(axis='x', length=0, pad=6)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.set_title('a  the top of the ranking is rebuilt, not just shrunk',
                 fontsize=8.2, loc='left', color=INK, pad=8)
    ax.plot([], [], color=BLUE, lw=1.4, label='rises under the matched control')
    ax.plot([], [], color=ACCENT, lw=1.4, label='falls')
    ax.legend(fontsize=6.6, frameon=False, loc='lower center', ncol=2,
              bbox_to_anchor=(0.5, -0.30))

    # ---------------------------------------------------------------- 오른쪽: 겹침
    t = tops.sort_values('shared_top50')
    nice = {'dn': 'diabetic nephropathy', 'iga nephropathy': 'IgA nephropathy',
            'membranous glomerulonephropathy': 'membranous GN',
            'minimal change disease': 'minimal change',
            'hypertensive nephropathy': 'hypertensive',
            'anca associated vasculitis': 'ANCA vasculitis',
            'systemic lupus erythematosus': 'lupus nephritis'}
    lab = [nice.get(x, x) for x in t['diagnosis']]
    y = np.arange(len(t))
    col = [ACCENT if x == 'dn' else GREY for x in t['diagnosis']]
    ax2.barh(y, t['shared_top50'], color=col, edgecolor=MUTED, linewidth=0.5, height=0.62)
    for yi, v in zip(y, t['shared_top50']):
        ax2.text(v + 0.7, yi, '%d' % v, va='center', fontsize=6.8, color=INK)
    ax2.set_yticks(y)
    ax2.set_yticklabels(lab, fontsize=6.8)
    ax2.set_xlim(0, 50)
    ax2.set_xticks([0, 25, 50])
    ax2.set_xticklabels(['0', '25', '50'], fontsize=7)
    ax2.set_xlabel('genes shared', fontsize=7.4)
    ax2.invert_yaxis()
    for s in ('top', 'right'):
        ax2.spines[s].set_visible(False)
    ax2.set_title('b  shared out of 50',
                  fontsize=8.2, loc='left', color=INK, pad=8)
    # 막대 옆에 두면 lupus 막대와 그 숫자를 덮는다. 위쪽 빈 곳으로 옮긴다.
    ax2.text(0.98, 0.99,
             'rank correlation %.2f\n(all 9,900 genes,\ndiabetic nephropathy)'
             % m['dn_spearman'],
             transform=ax2.transAxes, ha='right', va='top', fontsize=6.3, color=MUTED)

    # 이 문장은 캡션이 그대로 담고 있다. 그림 안에 두 번 적지 않는다.
    fig.subplots_adjust(left=0.135, right=0.955, top=0.92, bottom=0.20)
    if bib:
        # 패널 글자를 아래로. 설명은 캡션이 맡는다.
        for a in (ax, ax2):
            for loc in ('left', 'center', 'right'):
                a.set_title('', loc=loc)
        fig.subplots_adjust(bottom=0.26)
        for a, ch in zip((ax, ax2), 'ab'):
            bb = a.get_position()
            fig.text((bb.x0 + bb.x1) / 2, 0.015, ch, ha='center', va='bottom',
                     fontsize=8, fontweight='bold', color=INK)
    save(fig, NAME + ('_bib' if bib else ''), pdf=True)
    log('  후보 %d개를 그렸습니다. 오르는 것 %d개, 내리는 것 %d개.'
        % (len(d), int((d['moved'] > 0).sum()), int((d['moved'] <= 0).sum())))
    return 0


if __name__ == '__main__':
    sys.exit(main())
