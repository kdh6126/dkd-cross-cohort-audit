# -*- coding: utf-8 -*-
"""과제 발표 장표가 쓰는 수치를 결과 파일에서 읽는다.

장표에 숫자를 손으로 적지 않습니다. 분석이 바뀌면 장표도 함께 바뀌어야 하고, 그러지
않으면 발표장에서 원고와 다른 숫자를 말하게 됩니다. 이 프로젝트에서 이미 한 번 그런
일이 있었습니다 — 현황 덱이 옛 FDR 13.7% 를 몇 주 동안 들고 있었습니다.
"""
import os
import sys

import numpy as np
import pandas as pd

R = 'results'


def _csv(p, **kw):
    return pd.read_csv(p, sep='\t', **kw)


def facts():
    f = {}

    # ---- 자료 확보 현황
    inv = _csv(os.path.join(R, 'data_inventory.tsv'))
    f['inv_total'] = len(inv)
    f['inv_used'] = int((inv['state'] == 'used').sum())
    f['inv_excluded'] = int((inv['state'] == 'excluded').sum())
    f['inv_gb'] = float(inv['size_gb'].sum())

    # ---- KPMP: 다층 설계가 가능한가
    dm = _csv(os.path.join(R, 'kpmp_overlap/donor_meta.tsv'))
    f['kpmp_donor'] = len(dm)
    f['kpmp_diab'] = int(dm['is_diabetic'].astype(bool).sum())
    f['kpmp_proxy'] = int(dm['dkd'].astype(bool).sum())

    # ---- 유전자 공간
    gs = _csv('data/processed/gene_space_all6.tsv', dtype=str)
    f['genes'] = len(gs)

    # ---- 코호트 구성
    rows = []
    import glob
    for p in sorted(glob.glob('data/processed/harmonized/*_pheno.tsv')):
        d = _csv(p, dtype=str)
        name = os.path.basename(p).replace('_pheno.tsv', '')
        lab = d['label'].astype(str)
        ct = d.loc[lab == '0', 'control_type'].dropna()
        rows.append(dict(cohort=name, case=int((lab == '1').sum()),
                         control=int((lab == '0').sum()),
                         compartment=str(d['compartment'].iloc[0]),
                         control_type=', '.join(sorted(set(ct))) if len(ct) else '-'))
    f['cohorts'] = pd.DataFrame(rows)
    f['n_case'] = int(f['cohorts']['case'].sum())
    f['n_ctrl'] = int(f['cohorts']['control'].sum())

    # ---- 무작위 서명 널
    nl = _csv(os.path.join(R, 'null_control.tsv'))
    n50 = nl[nl['K'] == 50]
    f['null_lo'] = float(n50['random_mean'].min())
    f['null_hi'] = float(n50['random_mean'].max())
    f['null_perm'] = float(nl['permuted_mean'].mean())
    f['null_max'] = float(n50['random_max'].max())

    # ---- LODO 벤치마크
    cm = _csv(os.path.join(R, 'comparison_all4.tsv'))
    c100 = cm[cm['K'] == 100].sort_values('cross_fold_jaccard', ascending=False)
    f['auc_lo'] = float(cm[cm['K'] == 50]['external'].min())
    f['auc_hi'] = float(cm[cm['K'] == 50]['external'].max())
    f['xf_best'] = float(c100['cross_fold_jaccard'].iloc[0])
    f['xf_next'] = float(c100['cross_fold_jaccard'].iloc[1])
    f['xf_method'] = str(c100['method'].iloc[0])

    # ---- 대조군 순위 재배열
    rs = _csv(os.path.join(R, 'comparator_reordering/candidate_rank_shift.tsv'))
    f['shift'] = rs.reindex(rs['moved'].abs().sort_values(ascending=False).index)
    ov = _csv(os.path.join(R, 'comparator_reordering/top_overlap.tsv'))
    f['overlap'] = ov

    # ---- 특이성 게이트
    gc = _csv(os.path.join(R, 'gate_characteristics.tsv'))
    f['bg'] = 100 * float(gc['background_rate'].iloc[0])
    f['fdr'] = 100 * float(gc['gate_fdr_estimate'].iloc[0])
    f['n_pass'] = int(gc['genes_passing'].iloc[0])

    # ---- 후보표
    mt = _csv(os.path.join(R, 'candidates_v2/master_candidate_table.tsv'), index_col=0)
    f['cand'] = mt
    f['tier'] = mt['tier'].value_counts().sort_index()
    f['tier2'] = mt[mt['tier'].astype(str).str.startswith('2')]

    # ---- 보정 전후 문헌 회복
    un = _csv(os.path.join(R, 'literature_review.tsv'))
    co = _csv(os.path.join(R, 'candidates_v2/literature_review.tsv'))
    mk = lambda d: int((d['verdict'] == 'ALREADY PROPOSED as DKD marker').sum())
    f['lit_un'], f['lit_co'] = mk(un), mk(co)
    f['lit_overlap'] = len(set(un['gene']) & set(co['gene']))

    # ---- 상위 대 하위 대 무작위 (문헌 계층)
    f['strata'] = {}
    for key, p in (('상위 30 (제안 후보)', 'candidates_v2/literature_review.tsv'),
                   ('하위 30', 'candidates_v2/literature_bottom30.tsv'),
                   ('무작위 30', 'candidates_v2/literature_random30.tsv')):
        fp = os.path.join(R, p)
        if not os.path.exists(fp):
            continue
        d = _csv(fp)
        f['strata'][key] = dict(
            n=len(d),
            dkd=int((d['n_dkd'] > 0).sum()),
            marker=int((d['n_dkd_biomarker'] > 0).sum()),
            med_dkd=float(d['n_dkd'].median()),
            none=int((d['n_kidney'] == 0).sum()))

    # 상위 대 하위·무작위의 Fisher 단측 검정. 이 층 비교가 결론 장표의 근거다.
    from scipy import stats
    f['strata_p'] = {}
    top = f['strata'].get('상위 30 (제안 후보)')
    if top:
        for k, v in f['strata'].items():
            if k == '상위 30 (제안 후보)':
                continue
            f['strata_p'][k] = float(stats.fisher_exact(
                [[top['marker'], top['n'] - top['marker']],
                 [v['marker'], v['n'] - v['marker']]], alternative='greater')[1])

    # ---- 보정 벤치마크
    f['deconf'] = _csv(os.path.join(R, 'deconfound.tsv'))

    # ---- 알고리즘별 비교표 (선택기 8종, K=50)
    f['sel_table'] = (cm[cm['K'] == 50]
                      .sort_values('cross_fold_jaccard', ascending=False)
                      .reset_index(drop=True))

    # ---- 교란 보정 비교표 (6종, K=50, LODO 평균)
    rb = _csv(os.path.join(R, 'ruv_benchmark/benchmark.tsv'))
    rb = rb[rb['K'] == 50]
    # spec_enrichment 는 propagate_gate_fdr 가 나중에 붙이는 열이다. 그 단계를 아직
    # 돌리지 않았거나 파일을 다시 만든 직후면 없을 수 있으므로, 없으면 직접 계산한다.
    # 열 하나 때문에 장표 전체가 안 만들어지는 일은 없어야 한다.
    if 'spec_enrichment' not in rb.columns:
        rb = rb.assign(spec_enrichment=rb['spec_pass_frac'] / (f['bg'] / 100.0))
    f['corr_table'] = (rb.groupby('arm')
                       .agg(auroc=('external_auroc', 'mean'),
                            ieg=('n_ieg', 'mean'),
                            proc=('procurement_frac', 'mean'),
                            spec=('spec_pass_frac', 'mean'),
                            enrich=('spec_enrichment', 'mean'),
                            xfold=('cross_fold_jaccard', 'mean'))
                       .sort_values('proc'))

    # ---- 기존 방식 비교 (DEG 메타분석 · WGCNA 허브)
    cp = _csv(os.path.join(R, 'conventional_pipeline/comparison.tsv'))
    f['conv'] = cp[cp['K'] == 50].set_index('pipeline')

    # ---- 생성모델 (VAE · diffusion)
    gr = _csv(os.path.join(R, 'generative_robustness.tsv'))
    f['gen'] = gr.groupby('generator').agg(
        sig=('signature_auroc', 'mean'), rnd=('random_auroc', 'mean'),
        margin=('margin', 'mean'))

    # ---- 섭동 스트레스 · 코호트 민감도
    st = _csv(os.path.join(R, 'stress_test.tsv'))
    f['stress'] = st.groupby('perturbation').agg(
        sig=('signature_auroc', 'mean'), rnd=('random_auroc', 'mean'),
        margin=('margin', 'mean'))
    f['sens_size'] = _csv(os.path.join(R, 'cohort_sensitivity/by_size.tsv')).set_index('size')

    # ---- WGCNA 32변형: IEG 가 상위 50에 든 변형 수
    wg = _csv(os.path.join(R, 'wgcna_validation/variants.tsv'))
    f['wgcna_n'] = len(wg)
    f['wgcna_ieg_top50'] = int((wg['ieg_in_top50'] > 0).sum())

    # ---- GWAS
    gw = _csv(os.path.join(R, 'gwas_layer.tsv'))
    f['gwas_n'] = len(gw)
    f['gwas_sig'] = int((gw['empirical_p'] < 0.05).sum())

    # ---- 민감도
    sp = os.path.join(R, 'sensitivity_no104948/summary.tsv')
    if os.path.exists(sp):
        s = _csv(sp)
        v = s[~s['arm'].str.startswith('main')].iloc[0]
        f['sens_keep'] = int(v['n_overlap'])
        f['sens_p'] = float(v['mwu_p'])

    # ---- 선택기별 보정 전후 sweep (있으면)
    sp2 = os.path.join(R, 'correction_sweep/sweep.tsv')
    f['sweep'] = _csv(sp2) if os.path.exists(sp2) else None

    # ---- 파이프라인 규모
    sys.path.insert(0, os.getcwd())
    import run_all
    f['n_stage'] = len(run_all.S)
    f['n_paper'] = len(run_all.PAPER)
    return f
