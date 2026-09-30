#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""두 번째 단백체 코호트 — 질량분석으로 잰 DKD 대 대조.

첫 번째 단백체(Mendeley 83k89shdx5, SOMAscan)는 표적 패널이라 우리 후보 30개 중
7개만 쟀습니다. 그것만으로는 집합 수준의 판단이 불가능했습니다(Fisher p=0.11).

PRIDE PXD041884 는 완전히 다른 자료입니다.

    기술    LC-MS/MS  (앱타머가 아니라 질량분석)
    시료    사람 신장 FFPE 절편, 피질 60% 이상
    대조    당뇨 결절성 경화증 5명  대  비당뇨 대조 7명
    출처    Johns Hopkins 병리 보관 케이스, 양쪽 모두 적출조직(explant)

설계에 생검 대 신절제 대비가 없다는 점이 중요합니다. 이 프로젝트가 반복해서 지적해 온
교란은 질병군을 생검으로, 대조군을 종양신절제나 공여자 신장으로 받는 데서 생깁니다.
여기서는 양쪽이 같은 기관의 같은 종류(보관 FFPE 적출조직)입니다.

다만 여기까지만 말할 수 있습니다. "조달이 맞춰졌다" 고 쓰면 허혈시간이나 고정 시간을
실제로 쟀다는 뜻이 되는데, 이 자료도 그것을 보고하지 않습니다. 전분석 노출이 측정되지
않은 것은 우리 전사체 코호트와 똑같습니다. 설계에서 그 대비가 없다는 것과 노출이
같다는 것은 다른 말입니다.

검정은 Welch t검정을 씁니다. Mann-Whitney 를 쓰면 안 됩니다. n=5 대 7 에서 양측
최소 p 가 2/C(12,5)=0.0025 이므로, 1,600여 개를 BH 보정하면 q<0.05 에 도달하는 것이
원리적으로 불가능합니다. 순위검정의 해상도가 표본 수에 막히는 경우입니다.

첫 코호트와 같은 분석 계획을 그대로 씁니다. 후보를 배경(패널에 있으나 후보가 아닌
단백질)과 비교하고, 유의 비율과 방향 일치 비율 두 가지를 봅니다.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

XLSX = ('data/raw/proteomics/pride_PXD041884/'
        'Secocnd_Extraction_GPQ_Norm_TS_All_Peptides_FFPE_Protein_Rollup_Workbook.xlsx')
MAP = 'data/raw/proteomics/pride_PXD041884/uniprot_to_gene.tsv'
CAND = 'results/candidates_v2/master_candidate_table.tsv'
GVAL = 'data/processed/effect_sizes_hedges_g.tsv'
GSPACE = 'data/processed/gene_space_all6.tsv'
OUT = 'results/proteome_ms'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)

    d = pd.read_excel(XLSX, sheet_name='All Protien Rolled Up')
    u2g = dict(pd.read_csv(MAP, sep='\t').values)
    ctl = [c for c in d.columns if c.endswith(': Log2') and '_C' in c]
    dkd = [c for c in d.columns if c.endswith(': Log2') and '_D' in c]
    log('질량분석 단백체: 단백질 %d개 · 대조 %d명 · DKD %d명' % (len(d), len(ctl), len(dkd)))

    A = d[dkd].apply(pd.to_numeric, errors='coerce').values
    B = d[ctl].apply(pd.to_numeric, errors='coerce').values
    acc = d['Accession #'].astype(str).str.strip()
    pep = pd.to_numeric(d['Peptide Count'], errors='coerce')

    rows = []
    for i in range(len(d)):
        a = A[i][np.isfinite(A[i])]
        b = B[i][np.isfinite(B[i])]
        if len(a) < 3 or len(b) < 3:
            continue
        n1, n0 = len(a), len(b)
        sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / (n1 + n0 - 2))
        if sp == 0 or not np.isfinite(sp):
            continue
        pv = stats.ttest_ind(a, b, equal_var=False).pvalue
        if not np.isfinite(pv):
            continue
        g = (a.mean() - b.mean()) / sp * (1 - 3 / (4 * (n1 + n0) - 9))
        rows.append((u2g.get(acc.iloc[i]), g, pv,
                     float(np.nanmean(np.concatenate([a, b]))),
                     float(pep.iloc[i]) if np.isfinite(pep.iloc[i]) else np.nan,
                     float(np.isnan(np.concatenate([A[i], B[i]])).mean())))

    r = pd.DataFrame(rows, columns=['gene', 'prot_g', 'prot_p', 'abundance',
                                    'peptides', 'miss_frac']).dropna(subset=['gene'])
    # 한 유전자에 여러 accession 이 붙는 경우가 있다. 효과는 평균, p 는 가장 작은 것.
    r = r.groupby('gene', as_index=False).agg(prot_g=('prot_g', 'mean'),
                                              prot_p=('prot_p', 'min'),
                                              abundance=('abundance', 'mean'),
                                              peptides=('peptides', 'max'),
                                              miss_frac=('miss_frac', 'mean'))
    r['prot_q'] = bh(r['prot_p'].values)

    gs = pd.read_csv(GSPACE, sep='\t', dtype=str)
    gt = pd.read_csv(GVAL, sep='\t', dtype={0: str}).set_index('entrez_id')
    rna = {s: float(gt.loc[e].mean())
           for s, e in zip(gs['symbol'], gs['entrez_id']) if e in gt.index}
    r['rna_g'] = r['gene'].map(rna)
    r = r[r['rna_g'].notna()].copy()          # 우리 유전자 공간 안에 있는 것만 비교한다

    cand = set(pd.read_csv(CAND, sep='\t').iloc[:, 0].astype(str))
    r['is_cand'] = r['gene'].isin(cand)
    r['conc'] = np.sign(r['prot_g']) == np.sign(r['rna_g'])
    r.to_csv(os.path.join(OUT, 'protein_ms_dkd_vs_control.tsv'), sep='\t', index=False)

    c = r[r['is_cand']].sort_values('prot_q')
    bg = r[~r['is_cand']]
    log('')
    log('우리 유전자 공간과 겹치는 단백질 %d개, 그중 후보 %d개' % (len(r), len(c)))
    log('  %-10s %9s %9s %8s %6s' % ('유전자', '전사체 g', '단백질 g', 'q', '방향'))
    for _, x in c.iterrows():
        log('  %-10s %+9.2f %+9.2f %8.3f %6s'
            % (x['gene'], x['rna_g'], x['prot_g'], x['prot_q'],
               '일치' if x['conc'] else '불일치'))

    n_sig = int((c['prot_q'] < 0.05).sum())
    n_conc = int(c['conc'].sum())
    log('')
    log('  후보  q<0.05 %d/%d (%.0f%%) · 방향일치 %d/%d (%.0f%%)'
        % (n_sig, len(c), 100 * n_sig / len(c), n_conc, len(c), 100 * n_conc / len(c)))
    log('  배경  q<0.05 %d/%d (%.0f%%) · 방향일치 %d/%d (%.0f%%)'
        % (int((bg['prot_q'] < 0.05).sum()), len(bg), 100 * (bg['prot_q'] < 0.05).mean(),
           int(bg['conc'].sum()), len(bg), 100 * bg['conc'].mean()))

    tab = [[n_sig, len(c) - n_sig],
           [int((bg['prot_q'] < 0.05).sum()), int((bg['prot_q'] >= 0.05).sum())]]
    pf = stats.fisher_exact(tab, alternative='greater')[1]
    pb = stats.binomtest(n_conc, len(c), float(bg['conc'].mean()),
                         alternative='greater').pvalue
    log('  Fisher p = %.4g · 방향 이항검정 p = %.4g' % (pf, pb))
    log('  이 코호트 하나로는 유의하지 않습니다. 합쳐야 합니다 — proteome_meta.py.')

    pd.DataFrame([dict(n_proteins=len(r), n_candidates_measured=len(c), n_sig=n_sig,
                       n_concordant=n_conc, bg_n=len(bg),
                       bg_sig_frac=float((bg['prot_q'] < 0.05).mean()),
                       bg_conc_frac=float(bg['conc'].mean()),
                       fisher_p=pf, binom_p=pb, n_dkd=len(dkd), n_control=len(ctl))]
                 ).to_csv(os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
