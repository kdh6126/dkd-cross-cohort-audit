#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""전사체 후보가 단백질 수준에서도 재현되는가 — 올바른 대조로.

이 프로젝트는 오랫동안 "DKD 에는 두 번째 오믹스 층이 없다" 고 적어 왔습니다. 근거는
KPMP 지역 단백체에 당뇨 환자가 0명이라는 것이었습니다. 그 사실 자체는 맞지만, 거기서
"공개 데이터 전체에 없다" 로 넘어간 것이 틀렸습니다.

Mendeley 83k89shdx5 (Diabetes 2024, AKR1A1 논문이 쓴 자료) 는 사람 신장 피질 생검
SOMAscan 단백체이고, DKD 23명과 정상 10명을 담고 있습니다. 즉 우리가 필요로 하던
'DKD 대 정상' 대조를 단백질 층에서 제공합니다.

여기서 하는 것은 교차 오믹스 검증입니다.

    전사체   우리 코호트 7개에서 고른 후보 30개
    단백체   전혀 다른 환자 33명에서 잰 단백질 1,305개
    질문     후보가 단백질에서도 같은 방향으로 움직이는가

같은 환자가 아니므로 멀티오믹스는 아닙니다. 그러나 대조가 같으므로, 지금까지 우리가
가진 어떤 축보다 강한 확인입니다. 단일세포는 층이 같았고, KPMP 단백체는 대조가 달랐습니다.

귀무 비교를 반드시 함께 합니다. 패널에 있는 1,305개 중 후보가 아닌 것들이 보이는
분포와 비교해야, 후보가 특별한지 아니면 아무 단백질이나 그런지 알 수 있습니다.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

XLSX = 'data/raw/proteomics/mendeley_83k89shdx5/SOMAscan.xlsx'
CAND = 'results/candidates_v2/master_candidate_table.tsv'
GVAL = 'data/processed/effect_sizes_hedges_g.tsv'
GSPACE = 'data/processed/gene_space_all6.tsv'
OUT = 'results/proteome_validation'
SEED = 0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    n1, n0 = len(a), len(b)
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / (n1 + n0 - 2))
    if sp == 0 or not np.isfinite(sp):
        return np.nan
    return (a.mean() - b.mean()) / sp * (1 - 3 / (4 * (n1 + n0) - 9))


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    v = p[ok]
    o = np.argsort(v)
    r = np.empty_like(o)
    r[o] = np.arange(1, len(v) + 1)
    qq = v * len(v) / r
    qq = np.minimum.accumulate(qq[o][::-1])[::-1]
    out = np.empty_like(v)
    out[o] = np.minimum(qq, 1)
    q[ok] = out
    return q


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)

    d = pd.read_excel(XLSX)
    con = [c for c in d.columns if str(c).startswith('Con_')]
    dkd = [c for c in d.columns if str(c).startswith('DKD_')]
    log('단백체: 단백질 %d개 · 대조 %d명 · DKD %d명' % (len(d), len(con), len(dkd)))

    # SOMAscan 값은 오른쪽으로 크게 치우친다. 로그를 취해야 t검정이 의미를 가진다.
    X = d[con + dkd].apply(pd.to_numeric, errors='coerce')
    X = np.log2(X.clip(lower=1))
    A = X[dkd].values
    B = X[con].values

    res = pd.DataFrame({
        'gene': d['Entrez Gene Symbol'].astype(str).str.strip(),
        'protein': d['Protein Name'],
        'uniprot': d['UniProt ID'],
    })
    res['prot_g'] = [hedges_g(A[i][np.isfinite(A[i])], B[i][np.isfinite(B[i])])
                     for i in range(len(d))]
    res['prot_p'] = [stats.mannwhitneyu(A[i][np.isfinite(A[i])], B[i][np.isfinite(B[i])],
                                        alternative='two-sided').pvalue
                     for i in range(len(d))]
    res['prot_q'] = bh(res['prot_p'].values)
    # 매칭 민감도 분석용 공변량. "후보가 원래 잘 잡히는 단백질이라 유의한 것 아닌가"
    # 라는 질문은 이 값들을 맞춘 배경과 비교해야 답할 수 있다.
    res['abundance'] = np.nanmean(np.concatenate([A, B], axis=1), axis=1)
    res['miss_frac'] = np.isnan(np.concatenate([A, B], axis=1)).mean(axis=1)
    res = res.replace([np.inf, -np.inf], np.nan)

    # ---------------------------------------------------------- 전사체 쪽
    cand = pd.read_csv(CAND, sep='\t')
    cand_genes = set(cand[cand.columns[0]].astype(str))
    gs = pd.read_csv(GSPACE, sep='\t', dtype=str)
    sym2id = dict(zip(gs['symbol'], gs['entrez_id']))
    gtab = pd.read_csv(GVAL, sep='\t', dtype={0: str}).set_index('entrez_id')
    rna_g = {}
    for sym, eid in sym2id.items():
        if eid in gtab.index:
            rna_g[sym] = float(gtab.loc[eid].mean())

    res['is_candidate'] = res['gene'].isin(cand_genes)
    res['rna_g'] = res['gene'].map(rna_g)
    res['on_panel_and_measured'] = res['rna_g'].notna()
    res.to_csv(os.path.join(OUT, 'protein_dkd_vs_control.tsv'), sep='\t', index=False)

    # ---------------------------------------------------------- 후보 성적
    c = res[res['is_candidate']].copy()
    log('')
    log('=' * 78)
    log('후보 30개 중 이 패널에 있는 것: %d개' % len(c))
    if not len(c):
        log('  패널에 후보가 없습니다.')
        return 1
    c['concordant'] = np.sign(c['prot_g']) == np.sign(c['rna_g'])
    c = c.sort_values('prot_q')
    log('  %-10s %8s %8s %10s %6s' % ('유전자', '전사체 g', '단백질 g', '단백질 q', '방향'))
    for _, r in c.iterrows():
        log('  %-10s %8.2f %8.2f %10.3g %6s'
            % (r['gene'], r['rna_g'] if pd.notna(r['rna_g']) else np.nan,
               r['prot_g'], r['prot_q'], '일치' if r['concordant'] else '불일치'))

    n_sig = int((c['prot_q'] < 0.05).sum())
    n_conc = int(c['concordant'].sum())
    log('')
    log('  단백질에서 q<0.05: %d / %d' % (n_sig, len(c)))
    log('  방향 일치        : %d / %d' % (n_conc, len(c)))

    # ---------------------------------------------------------- 귀무 비교
    bg = res[(~res['is_candidate']) & res['on_panel_and_measured']].copy()
    bg['concordant'] = np.sign(bg['prot_g']) == np.sign(bg['rna_g'])
    log('')
    log('  배경(후보가 아닌 패널 단백질 %d개)' % len(bg))
    log('    q<0.05  %d개 (%.0f%%)   vs 후보 %.0f%%'
        % (int((bg['prot_q'] < 0.05).sum()), 100 * (bg['prot_q'] < 0.05).mean(),
           100 * n_sig / len(c)))
    log('    방향 일치 %.0f%%   vs 후보 %.0f%%'
        % (100 * bg['concordant'].mean(), 100 * n_conc / len(c)))

    tab = [[int((c['prot_q'] < 0.05).sum()), int((c['prot_q'] >= 0.05).sum())],
           [int((bg['prot_q'] < 0.05).sum()), int((bg['prot_q'] >= 0.05).sum())]]
    odds, p = stats.fisher_exact(tab, alternative='greater')
    log('')
    log('  후보가 배경보다 단백질에서 더 자주 유의한가:  Fisher p = %.4g (odds %.2f)'
        % (p, odds))
    pc = stats.binomtest(n_conc, len(c), float(bg['concordant'].mean()),
                         alternative='greater').pvalue
    log('  방향 일치가 배경보다 높은가:                  이항검정 p = %.4g' % pc)

    pd.DataFrame([dict(n_candidates_on_panel=len(c), n_sig=n_sig, n_concordant=n_conc,
                       bg_n=len(bg), bg_sig_frac=float((bg['prot_q'] < 0.05).mean()),
                       bg_conc_frac=float(bg['concordant'].mean()),
                       fisher_p=p, binom_p=pc)]).to_csv(
        os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)
    log('')
    log('  같은 환자가 아닙니다. 교차 오믹스이지 멀티오믹스가 아닙니다.')
    log('  다만 대조가 DKD 대 정상으로 같으므로, 지금까지의 어떤 축보다 직접적입니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
