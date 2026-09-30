#!/usr/bin/env python
"""Figures for the DKD cross-cohort study.

One file per figure under results/figures/, PNG at 200 dpi plus PDF for the ones likely to go
into a manuscript. Deliberately plain: no gridlines competing with data, no colour used to
carry information that a label could carry, and every axis states its units.
"""
import os, sys, glob, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = 'results/figures'
INK = '#1a1a1a'
MUTED = '#8a8a8a'
ACCENT = '#c1440e'
BLUE = '#2a5d8f'
GREY = '#d9d9d9'

# 저널 인쇄 폭은 170mm(6.69in)이다. 캔버스를 그 폭으로 잡고 글자를 절대 크기로 두면
# 저장된 PDF 가 인쇄될 크기 그대로여서, 편집부에서 축소하며 글자가 뭉개지는 일이 없다.
# BMC(sn-jnl) 본문 폭 372pt. 캔버스를 이 폭으로 그려야 삽입 배율이 1 이 되고,
# 코드에 적은 글자 크기가 지면 크기와 같아진다. 6.69in(170mm)로 그리던 동안에는
# 6pt 글자가 지면에서 4.6pt 였다.
JOURNAL_W = 372.0 / 72
# BiB(OUP 2단)에서 두 단을 가로지르는 그림의 폭 526pt. BiB 전용 그림만 이 폭으로 그린다.
BIB_W = 526.4 / 72

RULE = '#9a958c'      # 축선. 검정으로 그으면 데이터보다 축이 먼저 보인다.

plt.rcParams.update({
    'figure.dpi': 110, 'savefig.dpi': 400, 'savefig.bbox': None,
    # 저널은 그림 글꼴로 Arial 또는 Helvetica 를 지정한다. matplotlib 기본값인
    # DejaVu Sans 를 그대로 두면 어느 그림이든 '기본 설정으로 그렸다'는 티가 난다.
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    # 기본 Type 3 글꼴은 상당수 출판사가 거부한다. 42 는 TrueType 임베딩이다.
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5,
    'axes.edgecolor': RULE, 'axes.linewidth': 0.7, 'axes.spines.top': False,
    'axes.spines.right': False, 'text.color': INK, 'axes.labelcolor': INK,
    'xtick.color': RULE, 'ytick.color': RULE,
    'xtick.labelcolor': INK, 'ytick.labelcolor': INK,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'legend.frameon': False, 'axes.axisbelow': True,
    'lines.solid_capstyle': 'round', 'patch.linewidth': 0.6,
})


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def save(fig, name, pdf=False):
    os.makedirs(OUT, exist_ok=True)
    w = fig.get_size_inches()[0]
    if name.startswith('P') and not name.endswith('_bib') and w > JOURNAL_W + 0.01:
        log('  ! %s is %.2f in wide; the journal prints at %.2f in. It will be scaled '
            'down and the type will fall below 7 pt.' % (name, w, JOURNAL_W))
    fig.savefig(os.path.join(OUT, name + '.png'))
    if pdf:
        fig.savefig(os.path.join(OUT, name + '.pdf'))
    plt.close(fig)
    log('  wrote %s.png' % name)


# ------------------------------------------------------------------ F1 data inventory
def fig_inventory():
    d = pd.read_csv('results/data_inventory.tsv', sep='\t')
    order = ['used', 'excluded', 'unused', 'not_downloaded']
    labels = {'used': 'entered the analysis', 'excluded': 'excluded by the audit',
              'unused': 'downloaded, unused', 'not_downloaded': 'not downloaded'}
    cols = {'used': BLUE, 'excluded': ACCENT, 'unused': GREY, 'not_downloaded': '#f0f0f0'}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.6),
                                   gridspec_kw={'width_ratios': [1, 1.3]})
    n = [len(d[d['state'] == s]) for s in order]
    ax1.barh(range(len(order)), n, color=[cols[s] for s in order], edgecolor=INK, linewidth=0.6)
    ax1.set_yticks(range(len(order)))
    ax1.set_yticklabels([labels[s] for s in order])
    ax1.invert_yaxis()
    ax1.set_xlabel('number of datasets')
    ax1.set_title('Datasets by fate', loc='left')
    for i, v in enumerate(n):
        ax1.text(v + 0.15, i, str(v), va='center', fontsize=9)

    sub = d[d['state'].isin(['used', 'unused', 'excluded'])].sort_values('size_gb', ascending=True)
    sub = sub[sub['size_gb'] > 0.004]
    ax2.barh(range(len(sub)), sub['size_gb'],
             color=[cols[s] for s in sub['state']], edgecolor=INK, linewidth=0.6)
    ax2.set_yticks(range(len(sub)))
    ax2.set_yticklabels(sub['asset'], fontsize=7.5)
    ax2.set_xlabel('GB on disk')
    ax2.set_title('Volume is dominated by assets that were not used', loc='left')
    ax2.text(0.98, 0.04, '96% of bytes unused,\nbut only because single-cell\nfiles are large',
             transform=ax2.transAxes, ha='right', va='bottom', fontsize=7.5, color=MUTED)
    save(fig, 'F1_data_inventory', pdf=True)


# ------------------------------------------------------------------ F2 cross-fold stability
def fig_stability():
    g3 = pd.read_csv('results/comparison_glom3.tsv', sep='\t')
    a4 = pd.read_csv('results/comparison_all4.tsv', sep='\t')
    # NOT sharey: each panel is sorted by its own values, so each must carry its own tick
    # labels. Sharing the y axis silently relabels the right panel with the left panel's
    # method order, which is wrong whenever the two rankings differ - and they do.
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for ax, df, title in ((axes[0], g3, '3 glomerular cohorts (2 training per fold)'),
                          (axes[1], a4, '4 cohorts (3 training per fold)')):
        s = df[df['K'] == 50].sort_values('cross_fold_jaccard')
        cols = [ACCENT if m == 'rbs' else GREY for m in s['method']]
        ax.barh(range(len(s)), s['cross_fold_jaccard'], color=cols,
                edgecolor=INK, linewidth=0.6)
        ax.set_yticks(range(len(s)))
        ax.set_yticklabels(['RBS (proposed)' if m == 'rbs' else m for m in s['method']])
        ax.set_xlabel('cross-fold Jaccard at K = 50')
        ax.set_title(title, loc='left')
        for i, v in enumerate(s['cross_fold_jaccard']):
            ax.text(v + 0.006, i, '%.3f' % v, va='center', fontsize=8)
        ax.set_xlim(0, max(s['cross_fold_jaccard']) * 1.25)
    fig.suptitle('How much of the signature survives changing which cohort is held out',
                 x=0.02, ha='left', fontsize=11, y=0.995, va='top')
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, 'F2_cross_fold_stability', pdf=True)


# ------------------------------------------------------------------ F3 AUROC vs null
def fig_null():
    draws = pd.read_csv('results/null_control_draws.tsv', sep='\t')
    comp = pd.read_csv('results/comparison_all4.tsv', sep='\t')
    res = pd.concat([pd.read_csv('results/baselines_all4/results.tsv', sep='\t'),
                     pd.read_csv('results/proposed_all4/results.tsv', sep='\t')])
    folds = ['GSE30528', 'GSE96804', 'GSE104948', 'GSE142025']
    fig, axes = plt.subplots(1, 4, figsize=(JOURNAL_W, 2.5), sharey=True)
    for ax, f in zip(axes, folds):
        col = '%s|50' % f
        d = draws[col].dropna().values
        ax.hist(d, bins=26, color=GREY, edgecolor='white', linewidth=0.4)
        sub = res[(res['held_out'] == f) & (res['K'] == 50)]
        ymax = ax.get_ylim()[1]
        for _, r in sub.iterrows():
            is_rbs = r['method'] == 'rbs'
            ax.axvline(r['external_auroc'], color=ACCENT if is_rbs else BLUE,
                       lw=1.8 if is_rbs else 0.8, alpha=1.0 if is_rbs else 0.5)
        ax.axvline(np.median(d), color=INK, lw=0.9, ls=':')
        ax.set_title(f, loc='left', fontsize=7.5)
        ax.set_xlabel('external AUROC, K=50')
        ax.set_xlim(0.3, 1.02)
    axes[0].set_ylabel('random 50-gene signatures')
    # 첫 패널 안에 두면 막대와 겹친다. 네 패널에 공통으로 해당하는 설명이므로 그림
    # 하단에 한 줄로 둔다.
    fig.text(0.5, 0.005, 'grey, 300 random gene sets   |   red, stability (RBS)   |   '
             'blue, baselines   |   dotted, median of the null',
             ha='center', va='bottom', fontsize=6.2, color=MUTED)
    fig.tight_layout(rect=(0, 0.07, 1, 0.99))
    save(fig, 'F3_auroc_vs_null', pdf=True)


# ------------------------------------------------------------------ F4 procurement artifact
def fig_artifact():
    d = pd.read_csv('results/ieg_artifact.tsv', sep='\t')
    m = d.groupby('gene')[['dkd_vs_control', 'otherckd_vs_control']].mean()
    fig, ax = plt.subplots(figsize=(5.2, 5))
    lim = 3.1
    ax.plot([-lim, lim], [-lim, lim], color=MUTED, lw=0.8, ls='--')
    ax.axhline(0, color=INK, lw=0.6)
    ax.axvline(0, color=INK, lw=0.6)
    ax.scatter(m['dkd_vs_control'], m['otherckd_vs_control'], s=42,
               facecolor=ACCENT, edgecolor=INK, linewidth=0.5, zorder=3)
    # label only the genes far enough from the origin to matter, alternating the offset so
    # the dense cluster near (-0.7,-1.2) stays readable
    lab = m[(m.abs().max(axis=1) > 0.55)]
    for i, (g, r) in enumerate(lab.iterrows()):
        dx, dy = ((7, 4) if i % 2 == 0 else (-7, -9))
        ax.annotate(g, (r['dkd_vs_control'], r['otherckd_vs_control']),
                    textcoords='offset points', xytext=(dx, dy), fontsize=7,
                    ha='left' if dx > 0 else 'right')
    ax.set_xlim(-lim, 1.0); ax.set_ylim(-lim, 1.0)
    ax.set_xlabel("Hedges' g,  DKD vs control")
    ax.set_ylabel("Hedges' g,  non-diabetic CKD vs control")
    ax.set_title('Immediate-early genes: non-diabetic CKD moves as far or further\n'
                 'than DKD, so the module tracks procurement, not diabetes', loc='left')
    # Careful with the diagonal: every gene here is DOWN-regulated, so both axes are negative
    # and "below the diagonal" means the non-diabetic shift is the LARGER one. Describing it as
    # "DKD-specific" (the naive reading of below-the-line) inverts the conclusion.
    ax.text(0.03, 0.24,
            'a DKD-specific gene would sit near the x-axis\n'
            '(|y| much smaller than |x|)\n\n'
            'on the diagonal = non-diabetic CKD moves\n'
            'exactly as much as DKD\n'
            'beyond it = non-diabetic CKD moves MORE',
            transform=ax.transAxes, va='top', fontsize=7.5, color=MUTED)
    save(fig, 'F4_procurement_artifact', pdf=True)


# ------------------------------------------------------------------ F5 per-diagnosis heatmap
def fig_per_diagnosis():
    d = pd.read_csv('results/per_diagnosis_profile.tsv', sep='\t')
    d = d[d['cohort'] == 'GSE104948']
    keep = ['FMOD', 'MMP2', 'MOXD1', 'LUM', 'OLFML3', 'PRSS23', 'FOLR2', 'CFD', 'MS4A6A',
            'PCOLCE', 'COL1A2', 'S100A4', 'CALHM2', 'CCDC91', 'DPP6', 'C1orf21']
    cols = [c for c in d.columns if c.startswith('grp_')]
    order = (['grp_DN'] + [c for c in cols if c == 'grp_living_donor'] +
             sorted(c for c in cols if c not in ('grp_DN', 'grp_living_donor')))
    sub = d[d['gene'].isin(keep)].set_index('gene').reindex(keep).dropna(how='all')
    M = sub[order].values
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    v = np.nanmax(np.abs(M))
    im = ax.imshow(M, cmap='RdBu_r', vmin=-v, vmax=v, aspect='auto')
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([c[4:].replace('_', ' ')[:16] for c in order], rotation=40, ha='right', fontsize=8)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(sub.index, fontsize=8.5)
    for j, c in enumerate(order):
        if c == 'grp_DN':
            ax.add_patch(Rectangle((j - .5, -.5), 1, len(sub), fill=False,
                                   edgecolor=INK, lw=1.8, zorder=5))
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if abs(M[i, j]) > 0.9:
                ax.text(j, i, '%.1f' % M[i, j], ha='center', va='center', fontsize=6,
                        color='white' if abs(M[i, j]) > v * 0.6 else INK)
    fig.colorbar(im, ax=ax, shrink=0.7, label='mean z-score (glomeruli, ERCB)')
    ax.set_title('Candidate expression across nine diagnoses\n'
                 'a DKD-specific gene is high only in the boxed column', loc='left')
    save(fig, 'F5_per_diagnosis', pdf=True)


# ------------------------------------------------------------------ F6 stress test
def fig_stress():
    d = pd.read_csv('results/stress_test.tsv', sep='\t')
    kinds = ['subsample', 'noise', 'mask', 'batch']
    nice = {'subsample': 'subsampling', 'noise': 'gaussian noise',
            'mask': 'feature masking', 'batch': 'batch shift'}
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.1), sharey=True)
    for ax, k in zip(axes, kinds):
        g = d[d['perturbation'] == k].groupby('severity')
        s, r = g['signature_auroc'].mean(), g['random_auroc'].mean()
        ax.plot(s.index, s.values, '-o', color=ACCENT, ms=4, lw=1.6, label='signature')
        ax.plot(r.index, r.values, '-o', color=MUTED, ms=3, lw=1.1, label='random genes')
        ax.fill_between(s.index, r.values, s.values, color=ACCENT, alpha=0.12)
        ax.set_title(nice[k], loc='left')
        ax.set_xlabel('severity')
        ax.set_ylim(0.5, 1.0)
    axes[0].set_ylabel('external AUROC')
    axes[0].legend(fontsize=8, loc='lower left')
    fig.suptitle('The margin over random signatures is what stays constant, not the AUROC',
                 x=0.02, ha='left', fontsize=11, y=0.995, va='top')
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, 'F6_stress_test', pdf=True)


# ------------------------------------------------------------------ F7 candidate tiers
def fig_candidates():
    d = pd.read_csv('results/master_candidate_table.tsv', sep='\t')
    perm = pd.read_csv('results/specificity_null.tsv', sep='\t')
    ceil = float(np.percentile(perm['perm_max_abs_g'], 95))
    tiers = sorted(d['tier'].unique())
    marks = {t: m for t, m in zip(tiers, ['o', 's', '^', 'v'])}
    colr = {t: c for t, c in zip(tiers, [ACCENT, BLUE, '#7a9e3f', GREY])}
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    ax.axhline(ceil, color=INK, lw=0.9, ls='--')
    ax.axhline(-ceil, color=INK, lw=0.9, ls='--')
    ax.axhline(0, color=MUTED, lw=0.6)
    for t in tiers:
        s = d[d['tier'] == t]
        ax.scatter(s['mean_rank'], s['g_DKD_vs_otherCKD'], s=52, marker=marks[t],
                   facecolor=colr[t], edgecolor=INK, linewidth=0.5, label=t, zorder=3)
    for _, r in d.iterrows():
        if abs(r['g_DKD_vs_otherCKD']) > ceil or r['tier'].startswith('1'):
            ax.annotate(r['gene'], (r['mean_rank'], r['g_DKD_vs_otherCKD']),
                        textcoords='offset points', xytext=(6, 3), fontsize=8)
    ax.text(255, ceil + 0.04, '95th pct of label-permutation (|g| = %.2f)' % ceil,
            fontsize=7.5, color=MUTED, ha='right')
    ax.set_xlabel('RBS consensus rank (lower = selected more consistently)')
    ax.set_ylabel("Hedges' g,  DKD vs non-diabetic CKD")
    ax.set_title('Candidates: cross-cohort selection vs confounder-matched specificity', loc='left')
    ax.legend(fontsize=8, loc='lower right')
    save(fig, 'F7_candidates', pdf=True)


# ------------------------------------------------------------------ F8 gene space
def fig_genespace():
    import json
    gs = json.load(open('data/processed/genesets_per_cohort.json'))
    names = ['GSE142025', 'GSE96804', 'GSE30528', 'GSE104948']
    sizes = [len(gs[n]) for n in names]
    inter = len(set.intersection(*[set(gs[n]) for n in names]))
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    y = np.arange(len(names))
    ax.barh(y, sizes, color=GREY, edgecolor=INK, linewidth=0.6)
    ax.barh(y, [inter] * len(names), color=BLUE, edgecolor=INK, linewidth=0.6)
    ax.set_yticks(y); ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlabel('genes mapped to Entrez')
    for i, s in enumerate(sizes):
        ax.text(s + 200, i, '%d' % s, va='center', fontsize=8)
    ax.text(inter / 2, len(names) - 0.3, 'shared: %d' % inter, ha='center',
            fontsize=9, color='white', fontweight='bold')
    ax.set_title('Common gene space is bound by the ENTREZG custom-CDF cohorts', loc='left')
    save(fig, 'F8_gene_space', pdf=True)



# ------------------------------------------------------------------ F9 artifact replication
def fig_replication():
    """The procurement artifact in two independent cohorts, different platforms."""
    ercb = pd.read_csv('results/ieg_artifact.tsv', sep='	')
    e = ercb.groupby('gene')['otherckd_vs_control'].mean()
    rep = pd.read_csv('results/gse175759_validation_effects.tsv', sep='	')
    rep = rep.dropna(subset=['symbol']).set_index('symbol')['g_otherCKD_vs_control']
    genes = [g for g in e.index if g in rep.index]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2),
                                   gridspec_kw={'width_ratios': [1.15, 1]})

    y = np.arange(len(genes))
    ax1.barh(y - 0.2, e[genes].values, height=0.38, color=GREY, edgecolor=INK,
             linewidth=0.5, label='ERCB (array, glom+tub)')
    ax1.barh(y + 0.2, rep[genes].values, height=0.38, color=ACCENT, edgecolor=INK,
             linewidth=0.5, label='GSE175759 (RNA-seq, tubulo)')
    ax1.set_yticks(y); ax1.set_yticklabels(genes, fontsize=8)
    ax1.invert_yaxis()
    ax1.axvline(0, color=INK, lw=0.7)
    ax1.set_xlabel("Hedges' g,  non-diabetic CKD vs nephrectomy control")
    ax1.set_title('Immediate-early module, two independent cohorts', loc='left')
    ax1.legend(fontsize=8, loc='lower left')

    ieg = set(genes)
    allg = rep.dropna()
    other = allg[~allg.index.isin(ieg)]
    ax2.hist(other.values, bins=60, color=GREY, edgecolor='white', linewidth=0.3,
             label='all other genes')
    for v in rep[genes].values:
        ax2.axvline(v, color=ACCENT, lw=1.1, alpha=0.85)
    ax2.axvline(other.mean(), color=INK, lw=1.0, ls=':')
    ax2.set_xlabel("Hedges' g,  GSE175759")
    ax2.set_ylabel('genes')
    ax2.set_title('The module sits in the tail of the transcriptome\n'
                  'Mann-Whitney p = 5.4e-14', loc='left')
    ax2.legend(fontsize=8)
    fig.tight_layout()
    save(fig, 'F9_artifact_replication', pdf=True)


# ------------------------------------------------------------------ F10 single-cell control
def fig_singlecell():
    """KPMP single-nucleus: the artifact and its negative control side by side.

    The middle panel is the one that matters. DKD and non-diabetic CKD are BOTH biopsies in
    KPMP, so if the module tracked disease it would separate them; it does not (p = 0.21),
    while it separates DKD from non-biopsy controls at p = 6e-76.
    """
    d = pd.read_csv('results/kpmp_singlecell/pseudobulk_de.tsv', sep='\t')
    order = ['DKD_vs_reference', 'DKD_vs_CKDnondiabetic', 'reference_vs_nephrectomy']
    titles = ['DKD vs healthy reference\n(biopsy vs NON-biopsy)',
              'DKD vs non-diabetic CKD\n(biopsy vs biopsy)',
              'reference vs nephrectomy\n(non-biopsy vs non-biopsy)']
    kinds = ['random', 'candidate', 'IEG']
    cols = {'random': GREY, 'candidate': BLUE, 'IEG': ACCENT}
    fig, axes = plt.subplots(1, 3, figsize=(JOURNAL_W, 2.8), sharey=True)
    for ax, cname, ttl in zip(axes, order, titles):
        sub = d[d['contrast'] == cname]
        for i, k in enumerate(kinds):
            v = sub[sub['kind'] == k]['g'].dropna().values
            if not len(v):
                continue
            jitter = np.random.default_rng(0).normal(0, 0.07, len(v))
            ax.scatter(v, np.full(len(v), i) + jitter, s=7, alpha=0.45,
                       facecolor=cols[k], edgecolor='none')
            ax.plot([np.median(v)] * 2, [i - 0.28, i + 0.28], color=INK, lw=2, zorder=5)
            # 중앙값이 0 근처라 값표가 0 점선 위에 얹힌다. 흰 바탕을 깔아 선을 끊는다.
            ax.text(np.median(v), i + 0.36, '%+.2f' % np.mean(v), ha='center',
                    fontsize=6.4, color=INK, zorder=6,
                    bbox=dict(facecolor='white', edgecolor='none', pad=0.6))
        ax.axvline(0, color=INK, lw=0.8, ls=':')
        ax.set_yticks(range(len(kinds)))
        ax.set_yticklabels(['random genes', 'candidates', 'immediate-early'])
        # 세 패널이 같은 축이다. 각각에 긴 이름을 붙이면 인쇄 폭에서 서로 이어 붙는다.
        ax.set_title(ttl, loc='left', fontsize=7.2)
        ax.set_xlim(-3.2, 3.2)
    axes[0].text(0.02, 0.03, 'IEG vs random\np = 6e-76', transform=axes[0].transAxes,
                 fontsize=6.2, color=ACCENT, va='bottom')
    axes[1].text(0.02, 0.03, 'IEG vs random\np = 0.21  (n.s.)', transform=axes[1].transAxes,
                 fontsize=6.2, color=MUTED, va='bottom')
    axes[2].text(0.02, 0.03, 'IEG vs random\np = 0.89  (n.s.)', transform=axes[2].transAxes,
                 fontsize=6.2, color=MUTED, va='bottom')
    fig.supxlabel("Hedges' g per (gene, cell type), donor-level", fontsize=7, y=0.02)
    fig.suptitle('KPMP single-nucleus, 27 DKD-proxy / 10 non-diabetic CKD / 12 reference / 6 nephrectomy '
                 'donors', x=0.02, ha='left', fontsize=8.2, y=0.995, va='top')
    fig.tight_layout(rect=(0, 0.05, 1, 0.91))
    save(fig, 'F10_singlecell_control', pdf=True)


# ------------------------------------------------------------------ F11 deconfounded selection
def fig_deconfound():
    """Moving the procurement correction inside the selector: what it buys and what it costs."""
    d = pd.read_csv('results/deconfound.tsv', sep='\t')
    agg = (d.groupby(['arm', 'K'])
             .agg(auc=('external_auroc', 'mean'), sp=('spec_pass_frac', 'mean'),
                  pr=('procurement_frac', 'mean'), xf=('cross_fold_jaccard', 'mean'))
             .reset_index())
    arms = ['none', 'orth', 'penalty']
    labels = {'none': 'no correction', 'orth': 'residualise on\nhandling score',
              'penalty': 'penalise by\ncorrelation'}
    cols = {'none': GREY, 'orth': ACCENT, 'penalty': BLUE}
    Ks = sorted(agg['K'].unique())

    panels = [('pr', 'procurement-flagged genes\nin the selected set', True),
              ('sp', 'genes passing the\nDKD-specificity gate', False),
              ('auc', 'external AUROC', False),
              ('xf', 'cross-fold Jaccard', False)]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6))
    w = 0.26
    for ax, (key, ttl, lower_better) in zip(axes, panels):
        x = np.arange(len(Ks))
        for i, arm in enumerate(arms):
            v = [float(agg[(agg['arm'] == arm) & (agg['K'] == k)][key].iloc[0]) for k in Ks]
            if key in ('pr', 'sp'):
                v = [100 * z for z in v]
            ax.bar(x + (i - 1) * w, v, width=w, color=cols[arm], edgecolor=INK, linewidth=0.5,
                   label=labels[arm].replace('\n', ' ') if ax is axes[0] else None)
            for xi, vi in zip(x + (i - 1) * w, v):
                ax.text(xi, vi + (1.2 if key in ('pr', 'sp') else 0.012),
                        ('%.0f' % vi) if key in ('pr', 'sp') else ('%.2f' % vi),
                        ha='center', fontsize=7)
        ax.set_xticks(x); ax.set_xticklabels(['K=%d' % k for k in Ks])
        ax.set_title(ttl + ('   (lower is better)' if lower_better else ''),
                     loc='left', fontsize=9)
        if key in ('pr', 'sp'):
            ax.set_ylabel('%')
            ax.set_ylim(0, 68)
        elif key == 'auc':
            ax.set_ylim(0.6, 1.0)
        else:
            ax.set_ylim(0, 0.38)
    axes[0].legend(fontsize=8, loc='upper right')
    fig.suptitle('Residualising on the immediate-early handling score removes the procurement '
                 'axis from the selection, at a modest cost',
                 x=0.02, ha='left', fontsize=11, y=0.995, va='top')
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, 'F11_deconfound', pdf=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', default='')
    args = ap.parse_args()
    figs = {'inventory': fig_inventory, 'stability': fig_stability, 'null': fig_null,
            'artifact': fig_artifact, 'diagnosis': fig_per_diagnosis, 'stress': fig_stress,
            'candidates': fig_candidates, 'genespace': fig_genespace,
            'replication': fig_replication, 'singlecell': fig_singlecell,
            'deconfound': fig_deconfound}
    want = args.only.split(',') if args.only else list(figs)
    os.makedirs(OUT, exist_ok=True)
    for k in want:
        try:
            figs[k]()
        except Exception as e:
            log('  ! %s failed: %s' % (k, e))
    log('\nfigures in %s/' % OUT)


if __name__ == '__main__':
    main()
