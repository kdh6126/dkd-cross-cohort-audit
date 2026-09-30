#!/usr/bin/env python
"""F12: what changed when the candidate list was re-derived with the corrected selector."""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import save, INK, MUTED, ACCENT, BLUE, GREY
from dkd_deconfound import IEG

OUT = 'results/figures'


def main():
    v1 = pd.read_csv('results/final_candidates.tsv', sep='\t')
    v1 = v1.rename(columns={v1.columns[0]: 'entrez_id'})
    v1['entrez_id'] = v1['entrez_id'].astype(str)
    v2 = pd.read_csv('results/candidates_v2/candidates_all.tsv', sep='\t')
    v2 = v2.rename(columns={v2.columns[0]: 'entrez_id'})
    v2['entrez_id'] = v2['entrez_id'].astype(str)
    m1 = v1.set_index('entrez_id'); m2 = v2.set_index('entrez_id')
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = {s: g for g, s in zip(gs['entrez_id'], gs['symbol']) if isinstance(s, str)}

    fig = plt.figure(figsize=(12.6, 4.2))
    gsp = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1])

    # ---- panel A: immediate-early rank before/after
    ax = fig.add_subplot(gsp[0])
    rows = []
    for s in IEG:
        g = sym2id.get(s)
        if g in m1.index and g in m2.index:
            rows.append((s, m1.loc[g, 'mean_rank'], m2.loc[g, 'mean_rank']))
    rows.sort(key=lambda r: r[1])
    rows = rows[:12]
    y = np.arange(len(rows))
    for i, (s, a, b) in enumerate(rows):
        ax.annotate('', xy=(b, i), xytext=(a, i),
                    arrowprops=dict(arrowstyle='->', color=MUTED, lw=1.1))
    ax.scatter([r[1] for r in rows], y, s=42, facecolor=ACCENT, edgecolor=INK,
               linewidth=0.5, zorder=3, label='uncorrected')
    ax.scatter([r[2] for r in rows], y, s=42, facecolor=BLUE, edgecolor=INK,
               linewidth=0.5, zorder=3, label='corrected')
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xscale('log')
    ax.set_xlabel('consensus rank among 9,900 genes  (log scale)')
    ax.set_title('Immediate-early module falls out of the ranking', loc='left')
    ax.legend(fontsize=8, loc='lower right')

    # ---- panel B: specificity pass-rate, by signature size
    # Deliberately by K and not over the top 200: the correction concentrates its benefit at
    # the very top of the ranking, which is where a signature is actually taken from. Averaged
    # over 200 genes the two selectors look identical (46% vs 45%) and the effect disappears.
    ax = fig.add_subplot(gsp[1])
    dec = pd.read_csv('results/deconfound.tsv', sep='\t')
    agg = dec.groupby(['arm', 'K'])['spec_pass_frac'].mean()
    Ks = sorted(dec['K'].unique())
    x = np.arange(len(Ks)); w = 0.35
    a = [100 * agg[('none', k)] for k in Ks]
    b = [100 * agg[('orth', k)] for k in Ks]
    ax.bar(x - w / 2, a, w, color=ACCENT, edgecolor=INK, linewidth=0.6, label='uncorrected')
    ax.bar(x + w / 2, b, w, color=BLUE, edgecolor=INK, linewidth=0.6, label='corrected')
    for xi, v in zip(x - w / 2, a):
        ax.text(xi, v + 1, '%.0f%%' % v, ha='center', fontsize=9)
    for xi, v in zip(x + w / 2, b):
        ax.text(xi, v + 1, '%.0f%%' % v, ha='center', fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(['K=%d' % k for k in Ks])
    ax.set_ylabel('% of the selected set passing\nthe DKD-specificity gate')
    ax.set_ylim(0, 55)
    ax.set_title('Specificity, by signature size', loc='left')
    p1 = 100 * m1.sort_values('mean_rank').head(200)['procurement_flag'].mean()
    p2 = 100 * m2.sort_values('mean_rank').head(200)['procurement_flag'].mean()
    ax.text(0.02, 0.97, 'procurement-flagged in the top 200:\n%.0f%% -> %.0f%%' % (p1, p2),
            transform=ax.transAxes, va='top', fontsize=7.5, color=MUTED)
    ax.legend(fontsize=8, loc='lower right')

    # ---- panel C: tier composition of the final 30
    ax = fig.add_subplot(gsp[2])
    t1 = pd.read_csv('results/master_candidate_table.tsv', sep='\t')['tier'].value_counts()
    t2 = pd.read_csv('results/candidates_v2/master_candidate_table.tsv',
                     sep='\t')['tier'].value_counts()
    tiers = ['1 - novel AND strong', '2 - strong, already reported',
             '3 - novel, but weak in data', '4 - reported and weak']
    short = ['1 novel+strong', '2 strong,\npublished', '3 novel,\nweak', '4 reported,\nweak']
    x = np.arange(len(tiers)); w = 0.35
    va = [int(t1.get(t, 0)) for t in tiers]
    vb = [int(t2.get(t, 0)) for t in tiers]
    ax.bar(x - w / 2, va, w, color=ACCENT, edgecolor=INK, linewidth=0.6, label='uncorrected')
    ax.bar(x + w / 2, vb, w, color=BLUE, edgecolor=INK, linewidth=0.6, label='corrected')
    for xi, v in zip(x - w / 2, va):
        ax.text(xi, v + 0.3, str(v), ha='center', fontsize=9)
    for xi, v in zip(x + w / 2, vb):
        ax.text(xi, v + 0.3, str(v), ha='center', fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(short, fontsize=8)
    ax.set_ylabel('genes in the final 30')
    ax.set_title('Tier composition', loc='left')
    ax.text(0.02, 0.93, 'Tier 1 goes to zero:\nthe two "novel" candidates\nwere selection artifacts',
            transform=ax.transAxes, va='top', fontsize=7.5, color=MUTED)
    ax.legend(fontsize=8, loc='upper right')

    fig.suptitle('Re-deriving the candidates with the confounder-aware selector',
                 x=0.02, ha='left', fontsize=11, y=1.04)
    fig.tight_layout()
    save(fig, 'F12_v1_vs_v2', pdf=True)


if __name__ == '__main__':
    main()
