#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""공개 자료 목록을 하나의 조회 가능한 카탈로그로 만든다.

이 프로젝트는 "같은 환자에서 두 오믹스 층을 잴 수 있는가" 라는 질문에 세 번 답했고
두 번 틀렸습니다. 매번 원인이 같았습니다 — 층 하나만 세고 다른 층은 안 세었거나,
환자군만 세고 대조군을 안 세었습니다. 표를 만들어 두었으면 한눈에 보였을 것들입니다.

그래서 세 가지를 담습니다.

    dataset    무엇을 어디서 받았고 어떤 오믹스이며 몇 명인가
    subject    참여자 한 명이 어떤 층을 가지고 있는가, 질환 분류는 무엇인가
    view       질환 분류 x 층 쌍마다 사례와 대조가 각각 몇 명인가

세 번째가 이 카탈로그의 존재 이유입니다. 다중 오믹스가 가능하려면 **한 층 쌍에 대해
사례와 대조 양쪽에 사람이 있어야** 합니다. 어느 한쪽이 0명이면 짝은 있어도 비교가
성립하지 않습니다. 그 판정을 사람이 세지 않고 뷰가 하게 둡니다.

참여자 식별자는 KPMP 공개 아틀라스가 배포하는 비식별 번호입니다. 그래도 공개본에는
넣지 않습니다. 논문이 쓰는 것은 질환별 집계뿐이고, 개별 행을 배포할 이유가 없습니다.
"""
import os
import sqlite3
import sys

import numpy as np
import pandas as pd

OUT = 'db/omics_catalog.sqlite'
KEGG = 'results/enrichment/KEGG_2021_Human.tsv'
KPMP = 'results/kpmp_overlap/participant_layers.tsv'
DONOR = 'results/kpmp_overlap/donor_meta.tsv'
COHORTS = 'results/cohort_table.tsv'
PROV = 'results/dataset_provenance.tsv'
INV = 'results/data_inventory.tsv'
FEAS = 'results/multiomics_feasibility/layer_pair_contrast.tsv'
MAIN_DB = 'db/dkd.sqlite'
CAND = 'results/candidates_v2/master_candidate_table.tsv'
CAND_V1 = 'results/master_candidate_table.tsv'
COV = 'results/proteome_meta/candidate_protein_coverage.tsv'
SHIFT = 'results/comparator_reordering/candidate_rank_shift.tsv'
EFF = 'data/processed/effect_sizes_hedges_g.tsv'
GSPACE = 'data/processed/gene_space_all6.tsv'

# 군을 나누는 규칙. 카탈로그의 판정이 이 한 곳에서만 나온다.
#
# 처음에는 사례/대조/기타 셋뿐이었는데, 그러면 ERCB 의 다른 진단 생검 전부가 '기타' 로
# 묶였다. 그 생검들은 이 논문의 **조달 맞춘 비교군**이고 결론이 거기서 나온다. 이름이
# 없으면 카탈로그를 봐도 그 사실이 보이지 않으므로 comparator 를 따로 둔다.
CASE = {'diabetic nephropathy', 'ckd', 'aki', 'dkd', 'dmr',
        'diabetes mellitus with renal disease',
        # 등록 진단이 아니라 병력으로 고른 조작적 부분집합. 이름에 그 사실을 남긴다.
        'dkd proxy (ckd + diabetes history)'}
# 조직을 질병 때문이 아니라 다른 이유로 얻은 대조. 조달 경로가 사례와 다르다.
CONTROL = {'healthy reference tissue', 'healthy stone donor', 'tumor nephrectomy',
           'deceased donor', 'living donor', 'diabetes mellitus resilient', 'control'}

# 우리가 실제로 쓴 단백체 두 코호트. GEO 밖이라 dataset_provenance 에 없다.
EXTRA_DATASETS = [
    dict(accession='83k89shdx5', repository='Mendeley Data', omics='proteomics',
         assay='SOMAscan aptamer panel', tissue='kidney cortex',
         n_case=23, n_control=10, case_label='DKD', control_label='healthy participant',
         subject_ids_public=0,
         note='1,305 proteins. Source publication Hirohama 2023 (JASN).'),
    dict(accession='PXD041884', repository='PRIDE', omics='proteomics',
         assay='label-free LC-MS/MS', tissue='kidney cortex FFPE',
         n_case=5, n_control=7, case_label='DKD', control_label='non-diabetic case',
         subject_ids_public=1,
         note='Archival explant cases, one institution. Schwab 2024.'),
]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def omics_of(accession):
    """접근번호 접두사가 오믹스를 말해 준다.

    GSE 는 GEO 의 전사체, ST 는 Metabolomics Workbench 의 대사체, PXD 는 PRIDE 의
    단백체입니다. 한 값으로 몰아 넣으면 대사체 자료가 전사체로 둔갑합니다.
    """
    a = str(accession).upper()
    if a.startswith('GSE'):
        return 'transcriptomics'
    if a.startswith('ST'):
        return 'metabolomics'
    if a.startswith('PXD'):
        return 'proteomics'
    return 'other'


SYNONYM = {
    'anca-associated vasculitis': 'anca associated vasculitis',
    'focal segmental glomerular sclerosis': 'focal segmental glomerulosclerosis',
    'focal segmental glomerular sclerosis/minimal change disease':
        'focal segmental glomerulosclerosis and minimal change disease',
    'thin membrande disease': 'thin membrane disease',      # 원자료의 오타
}


def normalise_class(x):
    """질환 분류의 표기를 하나로 모은다.

    같은 병이 표기만 달라 여러 줄로 갈리면 집계가 나뉘고 판정이 틀어집니다. 실제로
    루푸스가 33 과 32 로, ANCA 혈관염이 22 와 21 로, 국소분절사구체경화증이 18 과 13 으로
    갈려 있었습니다. 대소문자와 하이픈만으로는 부족해 동의어 표까지 둡니다. 원자료의
    오타(thin membrande)도 여기서 고칩니다 — 고쳤다는 사실이 코드에 남습니다.
    """
    if x is None or (isinstance(x, float) and x != x) or str(x).strip() == '':
        return 'unstated'
    d = str(x).strip().lower().replace('_', ' ')
    d = d.replace('-', ' ')
    d = ' '.join(d.split())
    return SYNONYM.get(d, d)


def arm_of(disease):
    """질환 분류를 군으로 옮긴다. 입력은 normalise_class 를 거친 값이어야 한다."""
    d = normalise_class(disease)
    if d in CASE:
        return 'case'
    if d in CONTROL:
        return 'control'
    if d == 'unstated':
        return 'unstated'
    # 남는 것은 전부 다른 진단의 생검이다. 사례와 같은 방식으로 얻은 조직이므로
    # '기타' 가 아니라 비교군이다.
    return 'comparator'


def patient_deviation(c):
    """환자 한 명 한 명이 정상에서 얼마나 벗어났는지.

    z = (환자의 값 - 그 코호트 대조군 평균) / 대조군 표준편차. 코호트마다 따로
    계산합니다 — 플랫폼과 눈금이 달라 한꺼번에 재면 코호트 차이가 환자 차이처럼
    보입니다. 피처는 후보 30개로 한정합니다. 9,900개를 다 넣으면 250만 행이 되고
    화면에서 읽히지도 않습니다. 유전자를 묶는 '주소' 는 KEGG 경로이고, 어느
    경로에도 없는 것은 억지로 묶지 않고 미분류로 둡니다.
    """
    cand = list(pd.read_csv(CAND, sep='\t')['gene'])
    gs = pd.read_csv(GSPACE, sep='\t', dtype=str)
    sym2id = dict(zip(gs['symbol'], gs['entrez_id']))
    ids = {sym2id[g]: g for g in cand if g in sym2id}

    m = sqlite3.connect(MAIN_DB)
    smp = pd.read_sql(
        'select st.accession as cohort, s.accession as sample, s.label, s.compartment '
        'from sample s join study st on s.study_id = st.study_id '
        'where s.label is not null', m)
    m.close()

    rows = []
    for coh, g in smp.groupby('cohort'):
        path = 'db/columnar/%s_expr.parquet' % coh
        if not os.path.exists(path):
            continue
        X = pd.read_parquet(path)
        X['entrez_id'] = X['entrez_id'].astype(str)
        X = X[X['entrez_id'].isin(ids)].set_index('entrez_id')
        have = [s for s in g['sample'] if s in X.columns]
        if not have:
            continue
        sub = g[g['sample'].isin(have)]
        ctl = list(sub[sub['label'] == 0]['sample'])
        if len(ctl) < 3:
            continue
        mu = X[ctl].mean(axis=1)
        sd = X[ctl].std(axis=1).replace(0, np.nan)
        Z = X[have].sub(mu, axis=0).div(sd, axis=0)
        for gid, r in Z.iterrows():
            for s_, z in r.items():
                if z == z:
                    lab = sub.loc[sub['sample'] == s_, 'label'].iloc[0]
                    rows.append((coh, s_, {1.0: 'case', 0.0: 'control'}.get(lab, 'other'),
                                 ids[gid], float(z)))
    z = pd.DataFrame(rows, columns=['cohort', 'sample', 'arm', 'gene', 'z'])

    # 유전자의 '주소'. KEGG 경로 중 가장 유의한 것 하나를 붙인다.
    addr = {}
    if os.path.exists(KEGG):
        k = pd.read_csv(KEGG, sep='\t').sort_values('p_value')
        for _, r in k.iterrows():
            for gene in str(r['genes']).replace(';', ',').split(','):
                gene = gene.strip().upper()
                if gene and gene not in addr:
                    addr[gene] = r['term']
    z['pathway'] = [addr.get(g, '미분류') for g in z['gene']]

    z.to_sql('patient_deviation', c, index=False, if_exists='replace')
    c.execute('create index if not exists ix_dev_gene on patient_deviation(gene)')
    c.commit()
    print('환자별 편차 %d행 · 코호트 %d개 · 유전자 %d개 · 경로 %d개'
          % (len(z), z['cohort'].nunique(), z['gene'].nunique(), z['pathway'].nunique()))
    print(z.groupby('arm')['z'].agg(['count', 'mean']).round(2).to_string())


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs('db', exist_ok=True)
    if os.path.exists(OUT):
        os.remove(OUT)
    c = sqlite3.connect(OUT)

    # ------------------------------------------------------------------ dataset
    prov = pd.read_csv(PROV, sep='\t')
    coh = pd.read_csv(COHORTS, sep='\t')
    ds = prov.copy()
    ds['omics'] = [omics_of(a) for a in ds['accession']]
    ds = ds.merge(coh[['cohort', 'role', 'compartment', 'n_case', 'n_control',
                       'case_procurement', 'control_procurement']],
                  left_on='accession', right_on='cohort', how='left').drop(columns=['cohort'])
    for e in EXTRA_DATASETS:
        ds = pd.concat([ds, pd.DataFrame([dict(
            accession=e['accession'], repository=e['repository'], title=e['note'],
            platform=e['assay'], n_samples=e['n_case'] + e['n_control'],
            omics=e['omics'], role='proteomic corroboration',
            compartment=e['tissue'], n_case=e['n_case'], n_control=e['n_control'],
            case_procurement='biopsy or explant', control_procurement=e['control_label'])])],
            ignore_index=True)
    ds.to_sql('dataset', c, index=False)

    # ------------------------------------------------------------------ subject
    k = pd.read_csv(KPMP, sep='\t')
    # KPMP 의 등록 범주에는 DKD 항목이 없다. 범주만 쓰면 "KPMP 에 당뇨병성 신장질환
    # 참여자가 없다" 는 틀린 표가 만들어진다. h5ad 의 diabetes_history 를 보면 CKD 로
    # 등록된 공여자 중 27명이 당뇨 병력을 가진다. 등록 범주와 실제 병력은 다른 것이다.
    if os.path.exists(DONOR):
        dm = pd.read_csv(DONOR, sep='\t')
        dkd = set(dm.loc[dm['dkd'].astype(bool), 'donor_id'].astype(str))
        k['disease'] = [('DKD-proxy (CKD + diabetes history)' if str(p_) in dkd else d)
                        for p_, d in zip(k['participant'], k['disease'])]
    sub = pd.DataFrame(dict(
        participant=k['participant'].astype(str), resource='KPMP',
        disease_class=[normalise_class(x) for x in k['disease']],
        arm=[arm_of(x) for x in k['disease']],
        transcriptomics=k['has_transcriptomics'].astype(int),
        proteomics=k['has_proteomics'].astype(int),
        metabolomics=k['has_metabolomics'].astype(int)))
    sub['n_layers'] = sub[['transcriptomics', 'proteomics', 'metabolomics']].sum(axis=1)

    # 우리 전사체 코호트의 참여자는 GEO 시료 수준에서만 식별된다. 층은 전사체 하나뿐이
    # 므로 다중 오믹스 판정에는 기여하지 않지만, 카탈로그가 반쪽이 되지 않게 넣는다.
    if os.path.exists(MAIN_DB):
        m = sqlite3.connect(MAIN_DB)
        g = pd.read_sql(
            'select external_id as participant, source_cohort as resource, '
            'diagnosis as disease_class, is_diabetic from subject', m)
        m.close()
        g['disease_class'] = [normalise_class(x) for x in g['disease_class']]
        # 진단명이 있으면 그것으로 군을 정한다. is_diabetic 만 보면 다른 진단의 생검이
        # 전부 '기타' 가 되어, 이 논문의 비교군이 카탈로그에서 사라진다.
        g['arm'] = [arm_of(d) if d != 'unstated' else ('case' if x == 1 else 'unstated')
                    for d, x in zip(g['disease_class'], g['is_diabetic'].fillna(-1))]
        g = g.drop(columns=['is_diabetic'])
        g['transcriptomics'], g['proteomics'], g['metabolomics'] = 1, 0, 0
        g['n_layers'] = 1
        sub = pd.concat([sub, g], ignore_index=True)
    sub.to_sql('subject_layer', c, index=False)

    # ------------------------------------------------------------------ sample meta
    # 질병과 다른 변수가 **완전히 묶여 있는가**를 보려면 시료 수준 메타데이터가 필요하다.
    # 한 층(예: biopsy) 안에 사례와 대조가 함께 있어야 질병 효과를 그 변수와 분리해
    # 관찰할 수 있다. 어느 층에도 양쪽이 없으면 사후 보정으로도 나눌 수 없다 —
    # 나눌 정보 자체가 자료에 없기 때문이다.
    if os.path.exists(MAIN_DB):
        m = sqlite3.connect(MAIN_DB)
        sm = pd.read_sql(
            'select st.accession as cohort, s.label, s.procurement, s.compartment, '
            's.platform, s.disease_group from sample s '
            'join study st on s.study_id = st.study_id where s.label is not null', m)
        m.close()
        sm['arm'] = sm['label'].map({1.0: 'case', 0.0: 'control', -1.0: 'comparator'})
        sm.to_sql('sample_meta', c, index=False)

        c.executescript("""
        -- 코호트 x 축마다, 사례와 대조가 함께 있는 층이 몇 개인가.
        -- 0이면 완전 교란이다.
        create view v_confounding as
        with x(axis) as (values ('procurement'),('compartment'),('platform'))
        select cohort, 'procurement' as axis, procurement as stratum,
               sum(arm='case') as n_case, sum(arm='control') as n_control
        from sample_meta where arm in ('case','control') and procurement is not null
        group by cohort, procurement
        union all
        select cohort, 'compartment', compartment,
               sum(arm='case'), sum(arm='control')
        from sample_meta where arm in ('case','control') and compartment is not null
        group by cohort, compartment
        union all
        select cohort, 'platform', platform,
               sum(arm='case'), sum(arm='control')
        from sample_meta where arm in ('case','control') and platform is not null
        group by cohort, platform;

        -- 축마다 한 줄로 요약. separable_strata 가 0 이면 그 축에서 질병 효과를
        -- 떼어낼 수 없다.
        create view v_confounding_summary as
        select cohort, axis,
               count(*) as n_strata,
               sum(case when n_case>0 and n_control>0 then 1 else 0 end)
                   as separable_strata,
               sum(n_case) as n_case, sum(n_control) as n_control
        from v_confounding group by cohort, axis;
        """)
        c.commit()

    # ------------------------------------------------------------------ biomarker
    # "무엇이 나왔는가" 가 카탈로그에 없으면 반쪽이다. 후보와 그 근거를 함께 담는다.
    cand = pd.read_csv(CAND, sep='\t')
    bm = cand[['gene', 'tier', 'verdict', 'g_DKD_vs_control', 'g_DKD_vs_otherCKD',
               'above_perm_ceiling', 'pattern_class', 'donor_driven',
               'n_dkd', 'n_dkd_biomarker', 'lit_novel', 'data_strong',
               'sn_top_celltype', 'sn_top_foldchange',
               'protein_glom_vs_TI_FC', 'protein_adjP']].copy()
    bm['selector'] = 'corrected (handling-score residualised)'
    if os.path.exists(CAND_V1):
        v1 = set(pd.read_csv(CAND_V1, sep='\t')['gene'])
        bm['in_uncorrected_list'] = bm['gene'].isin(v1).astype(int)
    if os.path.exists(COV):
        cov = pd.read_csv(COV, sep='\t')
        bm = bm.merge(cov, on='gene', how='left')
    if os.path.exists(SHIFT):
        sh = pd.read_csv(SHIFT, sep='\t')[['gene', 'rank_naive', 'rank_matched']]
        bm = bm.merge(sh, on='gene', how='left')
    bm.to_sql('biomarker', c, index=False)

    # ------------------------------------------------------------------ model result
    # "그래서 모델을 돌리면 어떻게 되는가" 가 없으면 카탈로그가 자료 목록에서 끝난다.
    # 교란을 보정했을 때 무엇이 달라지는지가 이 프로젝트의 결론이므로 함께 담는다.
    for name, path in (('model_method', 'results/comparison_all4.tsv'),
                       ('model_correction', 'results/ruv_benchmark/benchmark.tsv'),
                       ('model_residualisation', 'results/deconfound_v2.tsv'),
                       ('model_null', 'results/null_control.tsv')):
        if os.path.exists(path):
            pd.read_csv(path, sep='	').to_sql(name, c, index=False)

    c.executescript("""
    -- 보정 방법마다 예측력과 교란 유래 유전자 비율이 어떻게 되는가.
    -- 좋은 보정은 procurement_frac 을 낮추면서 external_auroc 를 지킨다.
    create view v_correction_effect as
    select arm as method, K,
           round(avg(external_auroc), 3) as auroc,
           round(avg(procurement_frac), 3) as procurement_frac,
           round(avg(cross_fold_jaccard), 3) as reproducibility,
           count(*) as n_folds
    from model_correction group by arm, K;

    -- 무작위 유전자로도 나오는 AUROC. 절대값을 믿을 수 없다는 근거다.
    create view v_null_auroc as
    select K,
           round(avg(random_mean), 3) as random_mean,
           round(min(random_p05), 3) as random_p05,
           round(max(random_p95), 3) as random_p95,
           round(avg(permuted_mean), 3) as permuted_mean
    from model_null group by K;
    """)
    c.commit()

    # ------------------------------------------------------------------ gene effect
    # 측정값 자체(발현 행렬)는 db/columnar 의 Parquet 에 있고 46 MB 라 여기 넣지 않는다.
    # 조회에 쓸모 있는 것은 코호트별 효과크기이므로 그것을 담는다.
    eff = pd.read_csv(EFF, sep='\t', dtype={'entrez_id': str})
    gs = pd.read_csv(GSPACE, sep='\t', dtype=str)
    eff = eff.merge(gs[['entrez_id', 'symbol']], on='entrez_id', how='left')
    eff = eff.melt(id_vars=['entrez_id', 'symbol'], var_name='cohort',
                   value_name='hedges_g').dropna(subset=['hedges_g'])
    eff.to_sql('gene_effect', c, index=False)
    c.execute('create index ix_gene_effect_symbol on gene_effect(symbol)')

    # 성별과 연령대는 KPMP 에만 있고 GEO 코호트에는 없다. 측정 정의가 맞지 않으므로
    # 핵심 교란 보정 변수로 올리지 않고, KPMP 안의 탐색적 공변량으로 따로 둔다.
    if os.path.exists(DONOR):
        dm = pd.read_csv(DONOR, sep='	')
        CTL = {'Healthy_reference_tissue', 'Healthy_stone_donor', 'Tumor_nephrectomy',
               'Deceased_donor', 'Diabetes_mellitus_resilient'}
        dm['arm'] = ['case' if d else ('control' if c in CTL else 'other case')
                     for d, c in zip(dm['dkd'], dm['disease_category'])]
        dm.to_sql('kpmp_donor_covariate', c, index=False, if_exists='replace')
        c.executescript("""
        -- KPMP 안에서만 보는 탐색적 공변량. 사례와 대조가 함께 있는 층이 0이면
        -- 그 변수와 질병을 나눌 수 없다. GEO 코호트에는 이 값이 없으므로 합치지 않는다.
        create view v_kpmp_covariate as
        select 'sex' as axis, sex as stratum,
               sum(arm='case') as n_case, sum(arm='control') as n_control
        from kpmp_donor_covariate where arm in ('case','control') group by sex
        union all
        select 'age_stage', age_stage,
               sum(arm='case'), sum(arm='control')
        from kpmp_donor_covariate where arm in ('case','control') group by age_stage;
        """)
        c.commit()

    patient_deviation(c)

    # ------------------------------------------------------------------ views
    c.executescript("""
    -- 층 쌍마다 사례와 대조가 각각 몇 명인가. 어느 한쪽이 0이면 비교가 성립하지 않는다.
    create view v_multiomics_feasibility as
    with pair(layer_a, layer_b) as (
        values ('transcriptomics','proteomics'),
               ('transcriptomics','metabolomics'),
               ('proteomics','metabolomics'))
    select p.layer_a, p.layer_b, s.resource,
           sum(case when s.arm='case' and
                    ((p.layer_a='transcriptomics' and s.transcriptomics=1) or
                     (p.layer_a='proteomics'      and s.proteomics=1)) and
                    ((p.layer_b='proteomics'      and s.proteomics=1) or
                     (p.layer_b='metabolomics'    and s.metabolomics=1))
               then 1 else 0 end) as n_case,
           sum(case when s.arm='control' and
                    ((p.layer_a='transcriptomics' and s.transcriptomics=1) or
                     (p.layer_a='proteomics'      and s.proteomics=1)) and
                    ((p.layer_b='proteomics'      and s.proteomics=1) or
                     (p.layer_b='metabolomics'    and s.metabolomics=1))
               then 1 else 0 end) as n_control
    from pair p, subject_layer s
    group by p.layer_a, p.layer_b, s.resource;

    -- 질병 하나를 딱 집어 묻는 판정. 사례를 합쳐 세면 "DKD 에서 되는가" 를 알 수 없다.
    create view v_multiomics_by_disease as
    with pair(layer_a, layer_b) as (
        values ('transcriptomics','proteomics'),
               ('transcriptomics','metabolomics'))
    select p.layer_a, p.layer_b, s.resource, s.disease_class,
           sum(case when ((p.layer_a='transcriptomics' and s.transcriptomics=1)) and
                    ((p.layer_b='proteomics'   and s.proteomics=1) or
                     (p.layer_b='metabolomics' and s.metabolomics=1))
               then 1 else 0 end) as n_with_both,
           count(*) as n_in_class
    from pair p, subject_layer s
    where s.arm in ('case','control')
    group by p.layer_a, p.layer_b, s.resource, s.disease_class;

    -- 질환 분류별로 어떤 층을 몇 명이 가지고 있는가
    create view v_layer_by_class as
    select resource, disease_class, arm, count(*) as n_subjects,
           sum(transcriptomics) as transcriptomics,
           sum(proteomics) as proteomics,
           sum(metabolomics) as metabolomics,
           sum(case when n_layers > 1 then 1 else 0 end) as n_multilayer
    from subject_layer group by resource, disease_class, arm;

    -- 오믹스별 자료 목록
    create view v_dataset_by_omics as
    select omics, repository, count(*) as n_datasets,
           sum(coalesce(n_case,0)) as n_case, sum(coalesce(n_control,0)) as n_control
    from dataset group by omics, repository;
    """)
    c.commit()

    log('=' * 78)
    log('공개 자료 카탈로그  %s' % OUT)
    log('=' * 78)
    for t in ('dataset', 'subject_layer'):
        n = c.execute('select count(*) from "%s"' % t).fetchone()[0]
        log('  %-22s %5d행' % (t, n))
    log('')
    log('  다중 오믹스가 성립하는가 (사례와 대조 양쪽에 사람이 있어야 한다)')
    f = pd.read_sql('select * from v_multiomics_feasibility '
                    "where n_case>0 or n_control>0", c)
    for _, r in f.iterrows():
        ok = r['n_case'] > 0 and r['n_control'] > 0
        log('    %-14s + %-14s %-6s 사례 %2d · 대조 %2d   %s'
            % (r['layer_a'], r['layer_b'], r['resource'], r['n_case'], r['n_control'],
               '가능' if ok else '불가'))
    c.close()
    log('')
    log('  db/omics_catalog.sqlite 를 열거나, omics_catalog_html.py 로 브라우저에서 봅니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
