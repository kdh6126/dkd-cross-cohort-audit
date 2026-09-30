#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""사후 분석 H4, H5, H6, S1. 계획은 docs/xtissue/POSTHOC_PLAN.md 에 실행 전에 고정했다.

    h4   간 GSE48452 안에서 대조군 정의 치환(종양 수술 대조 vs 비만수술 정상 간)
    h5   즉시초기 모듈의 사례 효과에서 염증 점수를 뗐을 때 무엇이 남는가
    h6   대장 코호트를 신장 크기로 줄였을 때 판정 불안정이 생기는가 (계산이 크다)
    s1   간 RNA-seq 발현 필터를 완화했을 때 P3 순서가 유지되는가
"""
import argparse
import gzip
import itertools
import os
import shutil
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_cohort, load_many                               # noqa: E402
from dkd_deconfound import IEG                                            # noqa: E402
from xtissue_config import COHORTS, arm_match                             # noqa: E402
from xtissue_harmonize import matrix_path, parse_matrix, probe_map        # noqa: E402

OUT = 'results/xtissue/posthoc'
INFL = ['PTPRC', 'CD14', 'CD68', 'FCGR3A', 'S100A8', 'S100A9', 'LCN2', 'CXCL8', 'IL1B', 'TNF',
        'SAA1', 'MMP3', 'CHI3L1']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


_SYM = None


def sym2id():
    global _SYM
    if _SYM is None:
        _SYM = {}
        with gzip.open('data/raw/annotation/Homo_sapiens.gene_info.gz', 'rt', encoding='utf-8',
                       errors='replace') as f:
            hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
            ix = {c: i for i, c in enumerate(hdr)}
            for line in f:
                p = line.rstrip('\n').split('\t')
                _SYM[p[ix['Symbol']]] = p[ix['GeneID']]
    return _SYM


def hedges_vec(A, B):
    n1, n0 = len(A), len(B)
    sp = np.sqrt(((n1 - 1) * A.var(0, ddof=1) + (n0 - 1) * B.var(0, ddof=1)) / (n1 + n0 - 2))
    return (A.mean(0) - B.mean(0)) / (sp + 1e-12) * (1 - 3 / (4 * (n1 + n0) - 9))


# ================================================================ H4
def h4():
    samples, df = parse_matrix(matrix_path('GSE48452', 'GPL11532'))
    pmap, _ = probe_map('GPL11532')
    grp = {}
    for s in samples:
        src = s['source'].lower()
        if 'after surgery' in src:
            continue
        g = s['chars'].get('group', '').lower()
        grp[s['gsm']] = {'nash': 'NASH', 'control': 'oncological_surgery_control',
                         'healthy obese': 'bariatric_normal_liver'}.get(g)
    keep = [k for k, v in grp.items() if v]
    df = df[keep]
    df = df.loc[[i for i in df.index if i in pmap]]
    df.insert(0, 'entrez', [pmap[i] for i in df.index])
    df['_m'] = df.drop(columns='entrez').mean(axis=1)
    df = df.sort_values('_m', ascending=False).drop_duplicates('entrez').drop(columns='_m')
    X = df.set_index('entrez')
    X = X.loc[X.notna().all(axis=1)]
    Z = ((X.T - X.T.mean()) / X.T.std(ddof=1)).fillna(0.0)          # samples x genes
    lab = pd.Series({k: grp[k] for k in Z.index})
    genes = [str(g) for g in Z.columns]
    ieg = [i for i, g in enumerate(genes) if g in {sym2id().get(s) for s in IEG}]
    V = Z.values
    case = V[(lab == 'NASH').values]
    rng = np.random.default_rng(20260914)
    rows, gvec = [], {}
    for ctrl in ('oncological_surgery_control', 'bariatric_normal_liver'):
        ctl = V[(lab == ctrl).values]
        g = hedges_vec(case, ctl)
        gvec[ctrl] = g
        h_case, h_ctl = case[:, ieg].mean(1), ctl[:, ieg].mean(1)
        gi = hedges_vec(h_case[:, None], h_ctl[:, None])[0]
        boot = []
        for _ in range(2000):
            a = rng.choice(len(h_case), len(h_case))
            b = rng.choice(len(h_ctl), len(h_ctl))
            boot.append(hedges_vec(h_case[a][:, None], h_ctl[b][:, None])[0])
        rows.append(dict(control=ctrl, n_case=len(case), n_control=len(ctl),
                         g_ieg=float(gi), ci_lo=float(np.percentile(boot, 2.5)),
                         ci_hi=float(np.percentile(boot, 97.5)),
                         genes_abs_g_ge_0_5=int((np.abs(g) >= 0.5).sum()),
                         median_abs_g=float(np.median(np.abs(g)))))
    a, b = gvec['oncological_surgery_control'], gvec['bariatric_normal_liver']
    top = lambda v: set(np.argsort(-np.abs(v))[:50])
    summary = dict(genes=len(genes), ieg_present=len(ieg),
                   top50_shared=len(top(a) & top(b)),
                   spearman_all_genes=float(spearmanr(a, b).statistic),
                   median_abs_g_ratio_matched_over_naive=float(np.median(np.abs(b)) /
                                                               np.median(np.abs(a))))
    os.makedirs(OUT, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'h4_gse48452_contrasts.tsv'), sep='\t',
                              index=False)
    pd.DataFrame([summary]).to_csv(os.path.join(OUT, 'h4_gse48452_summary.tsv'), sep='\t',
                                   index=False)
    log(pd.DataFrame(rows).round(3).to_string(index=False))
    log(summary)


# ================================================================ H5
def h5():
    ids = sym2id()
    ieg_ids = {ids[s] for s in IEG if s in ids}
    infl_ids = {ids[s] for s in INFL if s in ids} - ieg_ids
    p2 = pd.read_csv('results/xtissue/p2_ieg_by_cohort.tsv', sep='\t')
    rows = []
    for _, r in p2.iterrows():
        root = ('data/processed/harmonized' if r['tissue'] == 'kidney_dkd'
                else 'data/processed/xtissue/%s' % r['tissue'])
        Z, y, _ = load_cohort(r['cohort'], root=root)
        genes = [str(g) for g in Z.columns]
        ic = [i for i, g in enumerate(genes) if g in ieg_ids]
        fc = [i for i, g in enumerate(genes) if g in infl_ids]
        if len(fc) < 5:
            continue
        V = Z.values
        h = V[:, ic].mean(1)
        f = V[:, fc].mean(1)
        hz = (h - h.mean()) / h.std(ddof=1)
        fz = (f - f.mean()) / f.std(ddof=1)
        yc = y - y.mean()
        X0 = np.column_stack([np.ones(len(y)), yc])
        X1 = np.column_stack([np.ones(len(y)), yc, fz])
        b0 = np.linalg.lstsq(X0, hz, rcond=None)[0][1]
        b1 = np.linalg.lstsq(X1, hz, rcond=None)[0][1]
        g_inf = hedges_vec(fz[y == 1][:, None], fz[y == 0][:, None])[0]
        rows.append(dict(tissue=r['tissue'], cohort=r['cohort'], arms=r['arms'],
                         n_infl_genes=len(fc), g_inflammation=float(g_inf),
                         case_coef_unadjusted=float(b0), case_coef_adjusted=float(b1),
                         retained_fraction=float(b1 / b0) if abs(b0) > 1e-6 else np.nan,
                         corr_handling_inflammation=float(np.corrcoef(hz, fz)[0, 1])))
    df = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(os.path.join(OUT, 'h5_inflammation_adjustment.tsv'), sep='\t', index=False)
    log(df.round(2).to_string(index=False))
    for (t, a), g in df.groupby(['tissue', 'arms']):
        log('%-12s %-11s n=%d  median unadjusted %+.2f  adjusted %+.2f  retained %.2f'
            % (t, a, len(g), g['case_coef_unadjusted'].median(), g['case_coef_adjusted'].median(),
               g['retained_fraction'].median()))


# ================================================================ H6
KIDNEY_SIZES = {'GSE30528': (9, 13), 'GSE142025': (27, 9), 'GSE294519': (23, 13),
                'GSE104948': (12, 26), 'GSE96804': (41, 20)}


def h6(reps, b_inner):
    from bootstrap_ci import run_all_methods
    cohorts = [c for c, cfg in COHORTS['colon_uc'].items() if cfg['status'] == 'core']
    data, genes = load_many(cohorts, root='data/processed/xtissue/colon_uc')
    size_rank = sorted(cohorts, key=lambda c: len(data[c][1]))
    kid_rank = sorted(KIDNEY_SIZES, key=lambda c: sum(KIDNEY_SIZES[c]))
    target = {c: KIDNEY_SIZES[k] for c, k in zip(size_rank, kid_rank)}
    genes = [str(g) for g in genes]
    ids = sym2id()
    ieg_cols = [i for i, g in enumerate(genes) if g in {ids.get(s) for s in IEG}]
    rng = np.random.default_rng(20260915)
    out = os.path.join(OUT, 'h6')
    os.makedirs(out, exist_ok=True)
    rows = []
    for rep in range(reps):
        dat, sizes = {}, {}
        for c in cohorts:
            X, y = data[c][0].values, data[c][1]
            n1, n0 = target[c]
            i1 = np.flatnonzero(y == 1)
            i0 = np.flatnonzero(y == 0)
            k1, k0 = min(n1, len(i1)), min(n0, len(i0))
            idx = np.concatenate([rng.choice(i1, k1, replace=False),
                                  rng.choice(i0, k0, replace=False)])
            Xs, ys = X[idx], y[idx]
            Xs = (Xs - Xs.mean(0)) / (Xs.std(0, ddof=1) + 1e-12)   # 줄인 표본 안에서 다시 z
            dat[c] = (Xs, ys)
            sizes[c] = (k1, k0)
        G = np.vstack([hedges_vec(dat[c][0][dat[c][1] == 1], dat[c][0][dat[c][1] == 0])
                       for c in cohorts])
        conc = float(np.median([spearmanr(G[i], G[j]).statistic
                                for i, j in itertools.combinations(range(len(cohorts)), 2)]))
        for sub in itertools.combinations(cohorts, 3):
            res = run_all_methods({c: dat[c] for c in sub}, list(sub), ieg_cols, 50, b_inner,
                                  'relieff', 20260825 + rep, 6)
            rows.append(dict(rep=rep, subset='+'.join(sub), concordance=conc,
                             jaccard_diff=res['RBS']['cross_fold_jaccard'] -
                             res['WGCNA_hub']['cross_fold_jaccard'],
                             auroc_diff=res['RBS']['external_auroc'] -
                             res['WGCNA_hub']['external_auroc'],
                             sizes=';'.join('%s:%d/%d' % (c, *sizes[c]) for c in sub)))
            pd.DataFrame(rows).to_csv(os.path.join(out, 'replicates.tsv'), sep='\t', index=False)
        r = pd.DataFrame([x for x in rows if x['rep'] == rep])
        log('rep %d: concordance %.3f, Jaccard diff %+.3f..%+.3f (%d negative), AUROC diff '
            '%+.3f..%+.3f (%d negative)' % (rep, conc, r['jaccard_diff'].min(),
                                             r['jaccard_diff'].max(), (r['jaccard_diff'] < 0).sum(),
                                             r['auroc_diff'].min(), r['auroc_diff'].max(),
                                             (r['auroc_diff'] < 0).sum()))


# ================================================================ S1
def s1():
    """RNA-seq 필터를 원시 계수 1 이상, 표본 20% 이상으로 완화한 간 코호트 디렉터리를 만든다."""
    import xtissue_harmonize as xh
    src = 'data/processed/xtissue/liver_masld'
    dst = 'data/processed/xtissue/liver_masld_relaxed'
    os.makedirs(dst, exist_ok=True)
    cfgs = COHORTS['liver_masld']
    for c, cfg in cfgs.items():
        if cfg['status'] != 'core':
            continue
        if not cfg.get('counts'):
            for suf in ('_expr.tsv', '_pheno.tsv'):
                shutil.copy(os.path.join(src, c + suf), os.path.join(dst, c + suf))
            continue
        samples, _ = parse_matrix(matrix_path(c, cfg['platform']))
        cnt = xh.counts_matrix(c, cfg, samples)
        keep = (cnt >= 1).mean(axis=1) >= 0.2
        lib = cnt.sum(axis=0)
        logcpm = np.log2(cnt.loc[keep] / lib * 1e6 + 1)
        ids = [xh.ensembl2entrez(i) for i in logcpm.index]
        logcpm = logcpm.loc[[i is not None for i in ids]]
        logcpm.insert(0, 'entrez', [i for i in ids if i is not None])
        logcpm['_m'] = logcpm.drop(columns='entrez').mean(axis=1)
        logcpm = logcpm.sort_values('_m', ascending=False).drop_duplicates('entrez')
        expr = logcpm.drop(columns='_m').set_index('entrez')
        ph = pd.read_csv(os.path.join(src, c + '_pheno.tsv'), sep='\t')
        expr = expr[[s for s in ph['sample'] if s in expr.columns]]
        expr.to_csv(os.path.join(dst, c + '_expr.tsv'), sep='\t')
        shutil.copy(os.path.join(src, c + '_pheno.tsv'), os.path.join(dst, c + '_pheno.tsv'))
        log('%s relaxed filter: %d genes' % (c, expr.shape[0]))


# ================================================================ concordance
def conc():
    """조직마다 코어 코호트끼리 유전자별 g 의 순위 상관. P1 이 갈린 뒤 사후에 본 값이다."""
    from xtissue_sweep import core_cohorts, KIDNEY5
    rows = []
    for t in ('kidney_dkd', 'liver_masld', 'colon_uc'):
        cs, root = ((KIDNEY5, 'data/processed/harmonized') if t == 'kidney_dkd'
                    else (core_cohorts(t), 'data/processed/xtissue/%s' % t))
        d, genes = load_many(cs, root=root)
        G = np.vstack([hedges_vec(d[c][0].values[d[c][1] == 1], d[c][0].values[d[c][1] == 0])
                       for c in cs])
        rho = [spearmanr(G[i], G[j]).statistic
               for i, j in itertools.combinations(range(len(cs)), 2)]
        top = np.median([np.sort(np.abs(G[i]))[::-1][:100].mean() for i in range(len(cs))])
        rows.append(dict(tissue=t, cohorts=len(cs), genes=len(genes),
                         median_pairwise_rho_of_g=float(np.median(rho)), min_rho=float(np.min(rho)),
                         median_top100_abs_g=float(top)))
    df = pd.DataFrame(rows)
    df.to_csv('results/xtissue/posthoc_effect_concordance.tsv', sep='\t', index=False)
    log(df.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('which', choices=['h4', 'h5', 'h6', 's1', 'conc'])
    ap.add_argument('--reps', type=int, default=5)
    ap.add_argument('--B-inner', type=int, default=30)
    a = ap.parse_args()
    {'h4': h4, 'h5': h5, 's1': s1, 'conc': conc}.get(a.which, lambda: h6(a.reps, a.B_inner))()


if __name__ == '__main__':
    main()
