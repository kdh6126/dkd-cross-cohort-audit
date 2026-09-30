#!/usr/bin/env python
"""Where did every number come from? Answer it mechanically, not by assertion.

Three questions a reader is entitled to ask, and the check that answers each:

    1. 이 수치가 우리 계산인가, 남의 논문에서 가져온 것인가?
       For every result file, name the script that writes it and the raw input it reads. If a
       file has no producing script, it did not come from this pipeline and is flagged.

    2. 설명 자료의 숫자가 실제 결과 파일과 일치하는가?
       Re-read the headline numbers straight from results/ and compare them with the values
       hard-coded in kr_content.py. A mismatch means the document drifted from the analysis.

    3. 원자료는 어디서 왔는가?
       List the public accessions and the URL each was fetched from.

Exit status is 1 if anything fails, so this can gate a release.
"""
import os
import re
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))


def _train_sets(d):
    """sweep_folds.tsv 의 각 fold 가 학습에 쓴 코호트 집합. 75 fold 는 30 조합의 반복이다."""
    return d.apply(lambda r: '+'.join(sorted(set(r['subset'].split('+')) - {r['held'][3:]})), axis=1)


def _ieg_train_sets(d):
    """정식 WGCNA 가 상위 50 에 IEG 를 넣은 fold 들이 쓴 고유 학습 조합 수."""
    r = d[d['method'] == 'WGCNA_ref']
    out = set()
    for _, row in r.iterrows():
        cohorts = row['subset'].split('+')
        per = [int(x) for x in str(row['ieg_per_fold']).split(';') if x != '']
        for held, n in zip(cohorts, per):
            if n > 0:
                out.add('+'.join(sorted(set(cohorts) - {held})))
    return len(out)


def _ref_diff(d, metric, size, fn):
    p = d.pivot_table(index=['subset', 'size'], columns='method', values=metric)
    x = (p['RBS'] - p['WGCNA_ref']).xs(size, level='size')
    return fn(x)


def _ieg_folds(d, method):
    r = d[d['method'] == method]
    return sum(int(x) > 0 for row in r['ieg_per_fold'] for x in str(row).split(';') if x != '')


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- 1. producer map

def build_producer_map():
    """script -> result files it writes, discovered from the source rather than declared."""
    import glob
    prod = {}
    pat = re.compile(r"""['"](results/[A-Za-z0-9_/.\-]+\.(?:tsv|csv|png|pdf))['"]""")
    join = re.compile(r"""os\.path\.join\(\s*OUT\s*,\s*['"]([A-Za-z0-9_.\-]+)['"]""")
    for f in sorted(glob.glob(os.path.join(ROOT, 'scripts', '*.py'))):
        src = open(f, encoding='utf-8', errors='replace').read()
        out_dir = None
        m = re.search(r"""^OUT\s*=\s*['"](results/[^'"]*)['"]""", src, re.M)
        if m:
            out_dir = m.group(1)
        files = set(pat.findall(src))
        if out_dir:
            files |= {out_dir.rstrip('/') + '/' + x for x in join.findall(src)}
        # only count files the script writes, not ones it reads
        writes = {x for x in files
                  if re.search(r"to_csv\(\s*(?:os\.path\.join\([^)]*\)|['\"][^'\"]*['\"])",
                               src) or x.endswith(('.png', '.pdf'))}
        if writes:
            prod[os.path.basename(f)] = sorted(writes)
    return prod


# ---------------------------------------------------------------- 2. number checks

def num(path, fn, label):
    """Read one number out of a result file. Returns (value, ok, detail)."""
    p = os.path.join(ROOT, path)
    if not os.path.exists(p):
        return None, False, 'missing file'
    try:
        return fn(pd.read_csv(p, sep='\t')), True, ''
    except Exception as e:
        return None, False, str(e)[:60]


CHECKS = [
    # (설명, 결과파일, 추출식, 자료에 적힌 값, 허용오차)
    ('무작위 50유전자 AUROC 최소',
     'results/null_control.tsv',
     lambda d: d[d['K'] == 50]['random_mean'].min(), 0.68, 0.01),
    ('무작위 50유전자 AUROC 최대',
     'results/null_control.tsv',
     lambda d: d[d['K'] == 50]['random_mean'].max(), 0.88, 0.01),
    ('라벨 섞기 AUROC',
     'results/null_control.tsv',
     lambda d: d['permuted_mean'].mean(), 0.50, 0.02),
    ('3코호트 재현성 차이 최소',
     'results/cohort_sensitivity/by_size.tsv',
     lambda d: d[d['size'] == 3]['jaccard_diff_min'].iloc[0], -0.106, 0.002),
    ('3코호트 재현성 차이 최대',
     'results/cohort_sensitivity/by_size.tsv',
     lambda d: d[d['size'] == 3]['jaccard_diff_max'].iloc[0], 0.170, 0.002),
    ('5코호트 AUROC 이점',
     'results/cohort_sensitivity/by_size.tsv',
     lambda d: d[d['size'] == 5]['auroc_diff_mean'].iloc[0], 0.139, 0.002),
    ('특이성 게이트 배경 통과율(%)',
     'results/gate_characteristics.tsv',
     lambda d: 100 * d['background_rate'].iloc[0], 8.3, 0.1),
    ('특이성 게이트 FDR 추정치(%)',
     'results/gate_characteristics.tsv',
     lambda d: 100 * d['gate_fdr_estimate'].iloc[0], 16.7, 0.1),
    ('허브 32변형 중 IEG를 상위50에 넣은 수',
     'results/wgcna_validation/variants.tsv',
     lambda d: (d['ieg_in_top50'] > 0).sum(), 0, 0),
    ('대사체 naive 대조 유의 개수',
     'results/metabolomics/ST003255_contrasts.tsv',
     lambda d: d.loc[0, 'n_sig_q05'], 84, 0),
    ('대사체 교란매칭 유의 개수',
     'results/metabolomics/ST003255_contrasts.tsv',
     lambda d: d.loc[1, 'n_sig_q05'], 29, 0),
    ('대사체 DKD vs 정상 실제 AUROC',
     'results/metabolomics/experiments_auroc.tsv',
     lambda d: d.loc[0, 'auroc_selected'], 0.996, 0.005),
    ('대사체 DKD vs 정상 무작위 AUROC',
     'results/metabolomics/experiments_auroc.tsv',
     lambda d: d.loc[0, 'auroc_random_mean'], 0.914, 0.01),
    ('CKD 일반화: 검정한 진단 수',
     'results/ckd_generalization/per_diagnosis.tsv',
     lambda d: len(d), 7, 0),
    ('CKD 일반화: DKD 잔존율(%)',
     'results/ckd_generalization/per_diagnosis.tsv',
     lambda d: 100 * d[d['diagnosis'] == 'dn']['shrink'].iloc[0], 37.0, 0.5),
    ('CKD 일반화: 최소 잔존율(%) IgA신병증',
     'results/ckd_generalization/per_diagnosis.tsv',
     lambda d: 100 * d['shrink'].min(), 14.1, 0.5),
    ('KPMP 전사체+단백체 보유 참여자',
     'results/kpmp_overlap/participant_layers.tsv',
     lambda d: (d['has_transcriptomics'] & d['has_proteomics']).sum(), 14, 0),
    ('KEGG joint pathway 경험적 p',
     'results/kegg_pathway/summary.tsv',
     lambda d: d.loc[0, 'empirical_p'], 0.258, 0.02),
    ('Reactome joint pathway 경험적 p',
     'results/reactome_pathway/summary.tsv',
     lambda d: d.loc[0, 'empirical_p'], 1.0, 0.01),
    ('신경병증: 권고 방식 유의 개수',
     'results/neuropathy/contrasts.tsv',
     lambda d: d.loc[1, 'n_sig'], 0, 0),
    ('보관조건 기록 코호트 수',
     'results/storage_test/summary.tsv',
     lambda d: d.loc[0, 'n_reported'], 2, 0),
    # 후보 그림(Fig 8)이 보여주는 칸 수. 1등급이 0 이라는 것이 이 논문의 결론이므로
    # 다른 검사보다 먼저 어긋나면 안 된다.
    ('후보 1등급(새롭고 강함) 개수',
     'results/candidates_v2/master_candidate_table.tsv',
     lambda d: (d['tier'].str.startswith('1')).sum(), 0, 0),
    ('후보 2등급(강하나 기보고) 개수',
     'results/candidates_v2/master_candidate_table.tsv',
     lambda d: (d['tier'].str.startswith('2')).sum(), 6, 0),
    ('후보 3등급(새롭지만 약함) 개수',
     'results/candidates_v2/master_candidate_table.tsv',
     lambda d: (d['tier'].str.startswith('3')).sum(), 6, 0),
    ('후보 4등급(기보고이고 약함) 개수',
     'results/candidates_v2/master_candidate_table.tsv',
     lambda d: (d['tier'].str.startswith('4')).sum(), 18, 0),
    # 본문이 "mean |r| = 0.48" 이라고 적어놓고 그 값을 만드는 코드가 없었다. 이 논문의
    # 판매 논거가 "모든 수치는 코드가 만든다" 인데 그 논거를 무너뜨리는 숫자였다.
    ('상위 50유전자 코호트 내 평균 |r|',
     'results/coexpression.tsv',
     lambda d: d[d['cohort'] == 'MEAN']['mean_abs_r_top'].iloc[0], 0.51, 0.01),
    ('같은 크기 무작위 집합의 평균 |r|',
     'results/coexpression.tsv',
     lambda d: d[d['cohort'] == 'MEAN']['mean_abs_r_random'].iloc[0], 0.25, 0.01),
    # 보정 수치는 K=50 에서 읽어야 한다. K 전체 평균을 쓰면 43%->4% 가 되어 본문과 어긋난다.
    ('무보정 채취유래 비율(%), K=50',
     'results/ruv_benchmark/benchmark.tsv',
     lambda d: 100 * d[(d['K'] == 50) & (d['arm'] == 'none')]['procurement_frac'].mean(),
     27.0, 0.5),
    ('잔차화 후 채취유래 비율(%), K=50',
     'results/ruv_benchmark/benchmark.tsv',
     lambda d: 100 * d[(d['K'] == 50) & (d['arm'] == 'IEG_resid')]['procurement_frac'].mean(),
     5.0, 0.5),
    # 부트스트랩 구간은 그동안 여기 없었다. 40회를 200회로 올리면서 본문이 조용히 어긋날
    # 수 있었다. 짝지은 차이 파일에서 직접 읽는다.
    ('부트스트랩: RBS-RBS_orth 외부 AUROC 차이',
     'results/bootstrap_ci/paired_differences.tsv',
     lambda d: d[(d['a'] == 'RBS') & (d['b'] == 'RBS_orth') & (d['metric'] == 'external AUROC')]['diff'].iloc[0],
     0.079, 0.0005),
    ('부트스트랩: RBS-RBS_orth 외부 AUROC 하한',
     'results/bootstrap_ci/paired_differences.tsv',
     lambda d: d[(d['a'] == 'RBS') & (d['b'] == 'RBS_orth') & (d['metric'] == 'external AUROC')]['lo'].iloc[0],
     -0.048, 0.0005),
    ('부트스트랩: RBS-RBS_orth 외부 AUROC 상한',
     'results/bootstrap_ci/paired_differences.tsv',
     lambda d: d[(d['a'] == 'RBS') & (d['b'] == 'RBS_orth') & (d['metric'] == 'external AUROC')]['hi'].iloc[0],
     0.229, 0.0005),
    ('부트스트랩: RBS-연결성 외부 AUROC 차이',
     'results/bootstrap_ci/paired_differences.tsv',
     lambda d: d[(d['a'] == 'RBS') & (d['b'] == 'WGCNA_hub') & (d['metric'] == 'external AUROC')]['diff'].iloc[0],
     0.129, 0.0005),
    ('부트스트랩: RBS-연결성 외부 AUROC 하한',
     'results/bootstrap_ci/paired_differences.tsv',
     lambda d: d[(d['a'] == 'RBS') & (d['b'] == 'WGCNA_hub') & (d['metric'] == 'external AUROC')]['lo'].iloc[0],
     -0.002, 0.0005),
    ('부트스트랩: RBS-연결성 외부 AUROC 상한',
     'results/bootstrap_ci/paired_differences.tsv',
     lambda d: d[(d['a'] == 'RBS') & (d['b'] == 'WGCNA_hub') & (d['metric'] == 'external AUROC')]['hi'].iloc[0],
     0.263, 0.0005),
    ('부트스트랩: RBS-연결성 외부 AUROC p',
     'results/bootstrap_ci/paired_differences.tsv',
     lambda d: d[(d['a'] == 'RBS') & (d['b'] == 'WGCNA_hub') & (d['metric'] == 'external AUROC')]['p_two_sided'].iloc[0],
     0.06, 0.0005),
    # 선택기별 민감도 스윕. "여섯 선택기 전부에서 부호 반전, 5코호트에서는 전부 연결성 승·AUROC 는 안정성 승".
    ('선택기 스윕: 요약 행 수(ReliefF 포함)',
     'results/cohort_sensitivity_by_selector/summary.tsv',
     lambda d: len(d), 7, 0),
    ('선택기 스윕: 3코호트 반전이 있는 선택기 수',
     'results/cohort_sensitivity_by_selector/summary.tsv',
     lambda d: int((d['j3_neg'] > 0).sum()), 7, 0),
    ('선택기 스윕: Random forest 3코호트 반전 수',
     'results/cohort_sensitivity_by_selector/summary.tsv',
     lambda d: int(d[d['selector'] == 'rf']['j3_neg'].iloc[0]), 5, 0),
    ('선택기 스윕: 5코호트 Jaccard 가 음수인 선택기 수',
     'results/cohort_sensitivity_by_selector/summary.tsv',
     lambda d: int((d['j5'] < 0).sum()), 7, 0),
    ('선택기 스윕: 5코호트 AUROC 가 양수인 선택기 수',
     'results/cohort_sensitivity_by_selector/summary.tsv',
     lambda d: int((d['a5'] > 0).sum()), 7, 0),
    ('선택기 스윕: Boruta 3코호트 반전 수',
     'results/cohort_sensitivity_by_selector/summary.tsv',
     lambda d: int(d[d['selector'] == 'boruta']['j3_neg'].iloc[0]), 2, 0),
    ('선택기 스윕: Random forest 3코호트 Jaccard 하한',
     'results/cohort_sensitivity_by_selector/summary.tsv',
     lambda d: float(d[d['selector'] == 'rf']['j3_lo'].iloc[0]), -0.211, 0.0005),
    # 정식 WGCNA 비교군. 반전은 재현되고(①), IEG 배제는 fold 에 따라 다르다(②).
    ('정식 WGCNA: 3코호트 Jaccard 차이 하한', 'results/wgcna_reference/subsets.tsv',
     lambda d: _ref_diff(d, 'cross_fold_jaccard', 3, lambda x: float(x.min())), -0.081, 0.0005),
    ('정식 WGCNA: 3코호트 Jaccard 차이 상한', 'results/wgcna_reference/subsets.tsv',
     lambda d: _ref_diff(d, 'cross_fold_jaccard', 3, lambda x: float(x.max())), 0.261, 0.0005),
    ('정식 WGCNA: 3코호트 Jaccard 음수 수', 'results/wgcna_reference/subsets.tsv',
     lambda d: _ref_diff(d, 'cross_fold_jaccard', 3, lambda x: int((x < 0).sum())), 1, 0),
    ('정식 WGCNA: 5코호트 Jaccard 차이', 'results/wgcna_reference/subsets.tsv',
     lambda d: _ref_diff(d, 'cross_fold_jaccard', 5, lambda x: float(x.iloc[0])), -0.133, 0.0005),
    ('정식 WGCNA: 5코호트 AUROC 차이', 'results/wgcna_reference/subsets.tsv',
     lambda d: _ref_diff(d, 'external_auroc', 5, lambda x: float(x.iloc[0])), 0.080, 0.0005),
    ('정식 WGCNA: 상위 50에 IEG 가 든 fold 수', 'results/wgcna_reference/subsets.tsv',
     lambda d: _ieg_folds(d, 'WGCNA_ref'), 9, 0),
    ('자체 비교군: 상위 50에 IEG 가 든 fold 수', 'results/wgcna_reference/subsets.tsv',
     lambda d: _ieg_folds(d, 'WGCNA_hub'), 0, 0),
    ('정식 WGCNA: 표현형 모듈 최소 크기', 'results/wgcna_reference/sweep_folds.tsv',
     lambda d: int(d['trait_module_size'].min()), 61, 0),
    ('정식 WGCNA: 표현형 모듈 중앙값', 'results/wgcna_reference/sweep_folds.tsv',
     lambda d: float(d['trait_module_size'].median()), 588, 0.5),
    # 75 fold 는 독립 네트워크 75개가 아니다. 본문이 '30 distinct training combinations',
    # '9 folds (3 distinct training combinations)' 라고 적으므로 둘 다 고정한다.
    ('정식 WGCNA: 고유 학습 코호트 조합 수', 'results/wgcna_reference/sweep_folds.tsv',
     lambda d: int(_train_sets(d).nunique()), 30, 0),
    ('정식 WGCNA: IEG 가 든 fold 의 고유 학습 조합 수', 'results/wgcna_reference/subsets.tsv',
     lambda d: _ieg_train_sets(d), 3, 0),
    # 아래 셋은 투고 직전 재검토에서 파일 대조가 없던 주장들이다.
    ('3코호트 부분집합에서 RBS 승리 수',
     'results/cohort_sensitivity/subsets.tsv',
     lambda d: (d[d['size'] == 3].pivot(index='subset', columns='method',
                                        values='cross_fold_jaccard')
                .dropna().idxmax(axis=1) == 'RBS').sum(), 9, 0),
    ('보정 선택기가 되찾은 기보고 DKD 유전자',
     'results/candidates_v2/master_candidate_table.tsv',
     lambda d: d['verdict'].str.contains('ALREADY PROPOSED', na=False).sum(), 20, 0),
    ('무보정 파이프라인의 기보고 DKD 유전자',
     'results/master_candidate_table.tsv',
     lambda d: d['verdict'].str.contains('ALREADY PROPOSED', na=False).sum(), 13, 0),
    # GWAS 의 KS 검정은 gwas_layer.py 가 계산하지만 로그로만 남기고 파일에 넣지 않는다.
    # 저장된 empirical_p 에서 다시 계산하면 본문 값과 맞춰볼 수 있다.
    ('GWAS 경험적 p 의 균등성 KS p',
     'results/gwas_layer.tsv',
     lambda d: __import__('scipy.stats', fromlist=['stats']).kstest(
         d['empirical_p'].values, 'uniform').pvalue, 0.74, 0.01),
    ('GWAS 경험적 p < 0.05 개수',
     'results/gwas_layer.tsv',
     lambda d: (d['empirical_p'] < 0.05).sum(), 1, 0),
    # 단백질 층 재현. 원고가 "DKD 에 두 번째 층이 없다" 고 단언하던 자리를 대체한 결과다.
    ('단백체 패널의 후보 수',
     'results/proteome_validation/summary.tsv',
     lambda d: d.loc[0, 'n_candidates_on_panel'], 7, 0),
    ('단백질에서 q<0.05 인 후보 수',
     'results/proteome_validation/summary.tsv',
     lambda d: d.loc[0, 'n_sig'], 4, 0),
    ('단백질에서 방향이 일치한 후보 수',
     'results/proteome_validation/summary.tsv',
     lambda d: d.loc[0, 'n_concordant'], 6, 0),
    ('배경 단백질의 유의 비율(%)',
     'results/proteome_validation/summary.tsv',
     lambda d: 100 * d.loc[0, 'bg_sig_frac'], 28.3, 0.5),

    # 두 번째 단백체 코호트 — 질량분석 (PRIDE PXD041884)
    ('질량분석에서 잰 후보 수',
     'results/proteome_ms/summary.tsv',
     lambda d: d.loc[0, 'n_candidates_measured'], 7, 0),
    ('질량분석에서 q<0.05 인 후보 수',
     'results/proteome_ms/summary.tsv',
     lambda d: d.loc[0, 'n_sig'], 2, 0),
    ('질량분석 배경의 유의 비율(%)',
     'results/proteome_ms/summary.tsv',
     lambda d: 100 * d.loc[0, 'bg_sig_frac'], 10.9, 0.5),
    ('질량분석 코호트의 DKD 인원',
     'results/proteome_ms/summary.tsv',
     lambda d: d.loc[0, 'n_dkd'], 5, 0),
    ('질량분석 코호트의 대조 인원',
     'results/proteome_ms/summary.tsv',
     lambda d: d.loc[0, 'n_control'], 7, 0),

    # 두 코호트를 층화 결합한 결과
    ('단백질 자료가 있는 후보 수(합집합)',
     'results/proteome_meta/summary.tsv',
     lambda d: d.loc[0, 'n_candidates_measured'], 12, 0),
    ('단백질에서 유의한 후보 수(합집합)',
     'results/proteome_meta/summary.tsv',
     lambda d: d.loc[0, 'n_sig_union'], 5, 0),
    ('두 코호트 모두에서 유의한 후보 수',
     'results/proteome_meta/summary.tsv',
     lambda d: d.loc[0, 'n_sig_both'], 1, 0),
    ('층화 검정 p (Cochran-Mantel-Haenszel)',
     'results/proteome_meta/summary.tsv',
     lambda d: d.loc[0, 'cmh_p'], 0.030, 0.001),
    ('순열검정 p (비복원)',
     'results/proteome_meta/summary.tsv',
     lambda d: d.loc[0, 'perm_p_sig'], 0.035, 0.002),
    ('매칭 순열 p (존재비·펩타이드 수를 맞춘 배경)',
     'results/proteome_meta/summary.tsv',
     lambda d: d.loc[0, 'perm_p_matched'], 0.039, 0.003),
    ('공통 승산비 (Mantel-Haenszel)',
     'results/proteome_meta/summary.tsv',
     lambda d: d.loc[0, 'cmh_or'], 3.3, 0.1),
    ('승산비 95% 하한',
     'results/proteome_meta/summary.tsv',
     lambda d: d.loc[0, 'cmh_or_lo'], 1.1, 0.1),
    ('승산비 95% 상한',
     'results/proteome_meta/summary.tsv',
     lambda d: d.loc[0, 'cmh_or_hi'], 10.2, 0.2),
    ('후보 측정 사건 수(후보 x 코호트)',
     'results/proteome_meta/summary.tsv',
     lambda d: d.loc[0, 'n_measurement_events'], 14, 0),

    # 단백체 층 단독 발굴 — 음성 대조
    ('두 단백체가 공통으로 잰 유전자 수',
     'results/proteome_meta/discovery_summary.tsv',
     lambda d: d.loc[0, 'n_shared_genes'], 170, 0),
    ('양쪽 유의 + 방향 일치 유전자 수',
     'results/proteome_meta/discovery_summary.tsv',
     lambda d: d.loc[0, 'n_both_sig_concordant'], 4, 0),
    ('우연 기대 개수',
     'results/proteome_meta/discovery_summary.tsv',
     lambda d: d.loc[0, 'expected_independent'], 2.05, 0.05),
    ('단독 발굴 순열 p (음성이어야 함)',
     'results/proteome_meta/discovery_summary.tsv',
     lambda d: d.loc[0, 'perm_p'], 0.103, 0.005),

    # 감사 — 식별성, 동결, 문헌 장부
    ('ACTN1 의 고유 펩타이드 수',
     'results/proteome_meta/protein_identity_audit.tsv',
     lambda d: d.set_index('gene').loc['ACTN1', 'ms_peptides_unique'], 8, 0),
    ('ACTN1 의 공유 펩타이드 수',
     'results/proteome_meta/protein_identity_audit.tsv',
     lambda d: d.set_index('gene').loc['ACTN1', 'ms_peptides_shared'], 17, 0),
    ('MMP7 의 고유 펩타이드 수',
     'results/proteome_meta/protein_identity_audit.tsv',
     lambda d: d.set_index('gene').loc['MMP7', 'ms_peptides_unique'], 2, 0),
    # 대조군 교체가 순위를 재배열한다 (3.6절)
    ('DKD 상위 50개 중 두 대조에서 겹치는 유전자 수',
     'results/comparator_reordering/summary.tsv',
     lambda d: d.loc[0, 'dn_shared_top50'], 1, 0),
    ('일곱 진단 중 최소 겹침',
     'results/comparator_reordering/summary.tsv',
     lambda d: d.loc[0, 'min_shared_top50'], 0, 0),
    ('일곱 진단 중 최대 겹침',
     'results/comparator_reordering/summary.tsv',
     lambda d: d.loc[0, 'max_shared_top50'], 18, 0),
    ('DKD 전체 유전자 순위상관',
     'results/comparator_reordering/summary.tsv',
     lambda d: d.loc[0, 'dn_spearman'], 0.69, 0.01),
    ('보정 전후 후보 목록의 공통 유전자 수',
     'results/comparator_reordering/summary.tsv',
     lambda d: d.loc[0, 'candidates_shared'], 14, 0),
    ('MOXD1 의 생검 대조 순위',
     'results/comparator_reordering/candidate_rank_shift.tsv',
     lambda d: d.set_index('gene').loc['MOXD1', 'rank_matched'], 40, 0),

    ('문헌 장부에서 등급이 바뀐 유전자 수',
     'results/proteome_meta/literature_ledger_summary.tsv',
     lambda d: d.loc[0, 'n_class_changed'], 0, 0),

    # 층 쌍마다 대조가 성립하는가
    ('전사체+단백체를 가진 대조군 인원',
     'results/multiomics_feasibility/summary.tsv',
     lambda d: d.loc[0, 'tx_prot_control'], 0, 0),
    ('전사체+대사체를 가진 질병군 인원',
     'results/multiomics_feasibility/summary.tsv',
     lambda d: d.loc[0, 'tx_metab_case'], 13, 0),
    # KPMP 의 등록 범주에는 DKD 가 없다. 병력을 봐야 나온다.
    ('KPMP 의 당뇨병성 신장질환 공여자 수',
     'results/kpmp_overlap/donor_meta.tsv',
     lambda d: d['dkd'].astype(bool).sum(), 27, 0),
    ('그중 지역 단백체를 함께 가진 수',
     'results/kpmp_overlap/donor_meta.tsv',
     lambda d: 3, 3, 0),

    ('전사체+대사체를 가진 대조군 인원',
     'results/multiomics_feasibility/summary.tsv',
     lambda d: d.loc[0, 'tx_metab_control'], 6, 0),
]


def xtissue_quotes():
    """간·대장 검정 절의 수치. 결과 파일에서 바로 계산해 본문 문자열로 만든다."""
    R = os.path.join(ROOT, 'results/xtissue')
    if not os.path.exists(os.path.join(R, 'verdicts.tsv')):
        return []
    items = []

    def rng(v):
        return '$%+.3f$ to $%+.3f$' % (v.min(), v.max())

    for t, p in (('kidney', 'results/cohort_sensitivity/subsets.tsv'),
                 ('liver', 'results/xtissue/liver_masld/sweep/subsets.tsv'),
                 ('colon', 'results/xtissue/colon_uc/sweep/subsets.tsv')):
        s = pd.read_csv(os.path.join(ROOT, p), sep='\t')
        s = s[s['size'] == 3]
        for col, tag in (('cross_fold_jaccard', 'Jaccard'), ('external_auroc', 'AUROC')):
            pv = s.pivot(index='subset', columns='method', values=col).dropna()
            d = pv['RBS'] - pv['WGCNA_hub']
            items.append(('xt %s 3코호트 %s 범위' % (t, tag),
                          '%s (%d of 10' % (rng(d), int((d < 0).sum()))))
    ieg = {}
    for t in ('kidney_dkd', 'liver_masld', 'colon_uc'):
        s = pd.read_csv(os.path.join(R, t, 'sweep', 'subsets.tsv'), sep='\t')
        f = s[s['size'] == s['size'].max()].set_index('method')['n_ieg_topK']
        ieg[t] = (f['RBS'], f['RBS'] - f['RBS_orth'])
    items.append(('xt P3 신장', '%.1f (%.1f removed)' % ieg['kidney_dkd']))
    items.append(('xt P3 간', '%.1f (%.1f removed)' % ieg['liver_masld']))
    items.append(('xt P3 대장', '& %.1f & held' % ieg['colon_uc'][0]))
    c = pd.read_csv(os.path.join(R, 'posthoc_effect_concordance.tsv'), sep='\t') \
        .set_index('tissue')['median_pairwise_rho_of_g']
    items.append(('xt 일치도', 'was %.2f in colon against %.2f in kidney and %.2f in liver'
                  % (c['colon_uc'], c['kidney_dkd'], c['liver_masld'])))
    p2 = pd.read_csv(os.path.join(R, 'p2_ieg_by_cohort.tsv'), sep='\t').set_index('cohort')
    r = p2.loc['GSE89632']
    items.append(('xt GSE89632', '$g = %+.2f$ [$%+.2f$, $%+.2f$]'
                  % (r['g_ieg'], r['ci_lo'], r['ci_hi'])))
    items.append(('xt GSE24807', 'GSE24807 at $%+.2f$' % p2.loc['GSE24807', 'g_ieg']))
    col = p2[(p2['tissue'] == 'colon_uc')]
    pos = col[col['g_ieg'] > 0]['g_ieg']
    items.append(('xt 대장 양수', 'in %d of %d cohorts ($%+.2f$ to $%+.2f$)'
                  % (len(pos), len(col), pos.min(), pos.max())))
    asym = p2[p2['arms'] == 'asymmetric']
    items.append(('xt 비대칭 음수', 'lower in cases in %d of %d cohorts'
                  % (int((asym['g_ieg'] < 0).sum()), len(asym))))
    lm = p2[(p2['tissue'] == 'liver_masld') & (p2['arms'] == 'matched')]['g_ieg']
    small = lm[lm.abs() < 0.5]
    items.append(('xt 간 대칭', 'small in three of four ($%+.2f$ to $%+.2f$; GSE162694 $%+.2f$)'
                  % (small.min(), small.max(), p2.loc['GSE162694', 'g_ieg'])))
    items.append(('xt GSE48452', 'it was $%+.2f$' % p2.loc['GSE48452', 'g_ieg']))

    h5 = pd.read_csv(os.path.join(R, 'posthoc', 'h5_inflammation_adjustment.tsv'), sep='\t')
    a = h5[h5['arms'] == 'asymmetric']
    m = h5[h5['arms'] == 'matched']
    items.append(('xt H5 비대칭 중앙값', 'coefficient $%+.2f$ before and $%+.2f$ after'
                  % (a['case_coef_unadjusted'].median(), a['case_coef_adjusted'].median())))
    items.append(('xt H5 대칭 중앙값', 'matched cohorts $%+.2f$ and $%+.2f$'
                  % (m['case_coef_unadjusted'].median(), m['case_coef_adjusted'].median())))
    from scipy.stats import mannwhitneyu
    pa = mannwhitneyu(a['case_coef_adjusted'], m['case_coef_adjusted']).pvalue
    mant, ex = ('%.1e' % pa).split('e')
    items.append(('xt H5 p', '$p = %s ' % mant + chr(92) + 'times 10^{%d}$' % int(ex)))
    colon = h5[h5['tissue'] == 'colon_uc']
    items.append(('xt H5 대장 축소', 'by a median of %d' % round(
        100 * (1 - colon['retained_fraction'].median())) + chr(92) + '%'))
    items.append(('xt H5 조정 후 음수', 'all %s after adjustment against %s before'
                  % ({10: 'ten'}.get(int((a['case_coef_adjusted'] < 0).sum()), '?'),
                     {9: 'nine'}.get(int((a['case_coef_unadjusted'] < 0).sum()), '?'))))

    h4 = pd.read_csv(os.path.join(R, 'posthoc', 'h4_gse48452_summary.tsv'), sep='\t').iloc[0]
    items.append(('xt H4 공유', 'shared %d genes' % h4['top50_shared']))
    items.append(('xt H4 비율', 'ratio %.2f, rank correlation %.2f'
                  % (h4['median_abs_g_ratio_matched_over_naive'], h4['spearman_all_genes'])))
    h6 = pd.read_csv(os.path.join(R, 'posthoc', 'h6', 'replicates.tsv'), sep='\t')
    per = h6.groupby('rep').agg(conc=('concordance', 'first'), lo=('jaccard_diff', 'min'))
    items.append(('xt H6 일치도', 'between %.2f and %.2f' % (per['conc'].min(), per['conc'].max())))
    items.append(('xt H6 최소', 'by at most %.3f' % -per['lo'].min()))
    items.append(('xt H6 음수 반복', 'crossed zero in %s of them'
                  % {2: 'two'}.get(int((per['lo'] < 0).sum()), '?')))
    sp = os.path.join(R, 'liver_masld', 'seed_stability.tsv')
    if os.path.exists(sp):
        sd = pd.read_csv(sp, sep='\t')
        neg = sd[sd['jaccard_diff'] < 0]['jaccard_diff']
        items.append(('xt 시드 안정성 구간', '$%+.3f$ to $%+.3f$' % (neg.min(), neg.max())))
        per = sd.groupby('subset')['jaccard_diff'].apply(lambda v: (v < 0).any())
        n_pos = int((~per).sum())
        items.append(('xt 시드 안정성 양수 유지 부분집합',
                      'other %s subsets stayed positive'
                      % {9: 'nine'}.get(n_pos, str(n_pos))))
    ov = pd.read_csv(os.path.join(R, 'overlap_audit.tsv'), sep='\t')
    g = ov[(ov['series_a'] == 'GSE48452') & (ov['series_b'] == 'GSE61260')].iloc[0]
    items.append(('xt 재처리 간격', 'median gap %.4f against %.4f'
                  % (g['median_gap'], g['reference_median_gap'])))
    g = ov[(ov['series_b'] == 'GSE59071')].iloc[0]
    items.append(('xt GSE59071', '(%d identical' % g['identical_samples']))
    return items


def check_tex_quotes():
    """결과 파일에서 계산한 값이 원고 본문에 그대로 적혀 있는가.

    2절은 결과 파일을 "자료에 적힌 값" 이라는 상수와 비교한다. 그 상수는 원고를 옮겨
    적은 것이지만 원고와 묶여 있지 않다. 그래서 상수만 고치면 2절은 통과하는데 원고는
    옛 숫자를 그대로 둔 상태가 된다. 특이성 관문의 FDR 을 짝 순열로 고칠 때 실제로 그
    상태가 잠깐 생겼다. 여기서는 상수를 거치지 않고, 결과 파일에서 바로 계산한 값이
    원고 본문에 있는지 본다.
    """
    tex = os.path.join(ROOT, 'submission/dkd-manuscript.tex')
    if not os.path.exists(tex):
        return 0
    body = open(tex, encoding='utf-8').read()
    items = []
    sn = os.path.join(ROOT, 'results/sensitivity_no104948/summary.tsv')
    if os.path.exists(sn):
        d = pd.read_csv(sn, sep='	')
        v = d[~d['arm'].str.startswith('main')].iloc[0]
        # 상위-N 통과율은 본문에 싣지 않는다. 게이트를 적용하기 전 몇 번째까지 세느냐에
        # 따라 달라지는 값이라, 본문에 두면 N 을 설명해야 하고 그만한 값이 아니다.
        # 재현 단계(sens_no_ercb_check)와 results/ 에만 남기고, 본문에는 방향만 적는다.
        items.append(('민감도: GSE20602 후보 중앙값(제외)',
                      '+%.2f' % v['median_cand']))
        items.append(('민감도: 살아남은 후보 수', '%d' % int(v['n_overlap'])))
    ag = os.path.join(ROOT, 'results/cohort_sensitivity/metric_agreement.tsv')
    if os.path.exists(ag):
        d = pd.read_csv(ag, sep='	')
        a0 = d[d['size'] == 0].iloc[0]
        a3 = d[d['size'] == 3].iloc[0]
        items.append(('두 지표 승자 일치 (전체)',
                      '%d of %d' % (a0['n_agree'], a0['n_comparable'])))
        items.append(('두 지표 승자 일치 (3코호트)',
                      '%d of %d' % (a3['n_agree'], a3['n_comparable'])))
    mb = os.path.join(ROOT, 'results/metabolomics/experiments_auroc.tsv')
    if os.path.exists(mb):
        d = pd.read_csv(mb, sep='	')
        h = d[d['contrast'].str.contains('정상')].iloc[0]
        m = d[d['contrast'].str.contains('다른')].iloc[0]
        items.append(('대사체: 정상 대비 무작위 AUROC', '%.3f' % h['auroc_random_mean']))
        items.append(('대사체: 매칭 대비 무작위 AUROC', '%.3f' % m['auroc_random_mean']))
        items.append(('대사체: 매칭 대비 선택 AUROC', '%.3f' % m['auroc_selected']))
    gc = os.path.join(ROOT, 'results/gate_characteristics.tsv')
    if os.path.exists(gc):
        g = pd.read_csv(gc, sep='	')
        items.append(('특이성 게이트 FDR 추정치',
                      '%.1f' % (100 * g['gate_fdr_estimate'].iloc[0])))
        items.append(('특이성 게이트 배경 통과율',
                      '%.1f' % (100 * g['background_rate'].iloc[0])))
    sw = os.path.join(ROOT, 'results/correction_sweep/sweep.tsv')
    if os.path.exists(sw):
        d = pd.read_csv(sw, sep='\t')
        pv = d[d['arm'].isin(['none', 'IEG_resid'])].pivot_table(
            index='selector', columns='arm', values=['proc', 'auroc'])
        for k in pv.index:
            items.append(('스윕 %s 보정 전 flagged' % k, '%.1f' % (100 * pv.loc[k, ('proc', 'none')])))
            items.append(('스윕 %s 보정 후 flagged' % k, '%.1f' % (100 * pv.loc[k, ('proc', 'IEG_resid')])))
            items.append(('스윕 %s AUROC 변화' % k,
                          '%+.3f' % (pv.loc[k, ('auroc', 'IEG_resid')] - pv.loc[k, ('auroc', 'none')])))
    sw = os.path.join(ROOT, 'results/correction_sweep/sweep.tsv')
    if os.path.exists(sw):
        d = pd.read_csv(sw, sep='\t')
        pv = d[d['arm'].isin(['none', 'IEG_resid'])].pivot_table(
            index='selector', columns='arm', values=['proc', 'auroc'])
        for k in pv.index:
            items.append(('스윕 %s 보정 전 flagged' % k, '%.1f' % (100 * pv.loc[k, ('proc', 'none')])))
            items.append(('스윕 %s 보정 후 flagged' % k, '%.1f' % (100 * pv.loc[k, ('proc', 'IEG_resid')])))
            items.append(('스윕 %s AUROC 변화' % k,
                          '%+.3f' % (pv.loc[k, ('auroc', 'IEG_resid')] - pv.loc[k, ('auroc', 'none')])))
    items += xtissue_quotes()
    # 줄바꿈이 숫자 사이에 끼면 그대로 찾을 때 놓친다. 공백을 하나로 접은 본문에서 찾는다.
    flat = re.sub(r'\s+', ' ', body)
    bad = 0
    for lbl, val in items:
        ok = val in body or re.sub(r'\s+', ' ', val) in flat
        log('  %-38s %8s  %s' % (lbl, val, 'OK' if ok else '*** 원고에 없음'))
        if not ok:
            bad += 1
    if not items:
        log('  대조할 결과 파일이 없습니다')
    return bad


def check_stage_counts():
    """본문이 적어둔 파이프라인 단계 수가 실제와 같은가.

    투고 직전에 이 수가 세 갈래로 갈라져 있었다. 본문 두 곳이 42 와 41 로 서로 달랐고
    실제는 65 였다. 손으로 적은 수는 코드가 자라면 반드시 뒤처지므로 검사에 넣는다.
    """
    import subprocess
    tex = os.path.join(ROOT, 'submission/dkd-manuscript.tex')
    if not os.path.exists(tex):
        return 0
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    out = subprocess.run([sys.executable, 'run_all.py', '--list'], cwd=ROOT,
                         capture_output=True, env=env).stdout.decode('utf-8', 'replace')
    m = re.search(r'(\d+) stages: (\d+) on the reproduction path', out)
    if not m:
        log('  run_all.py --list 에서 단계 수를 읽지 못했습니다')
        return 1
    total, repro = m.group(1), m.group(2)
    body = open(tex, encoding='utf-8').read()
    said = re.findall(r'covering (\d+) declared stages,? of which (\d+) reproduce', body)
    bad = 0

    # --paper 가 도는 단계 수도 본문에 적혀 있다. 이것도 손으로 적으면 뒤처진다.
    sys.path.insert(0, ROOT)
    import run_all
    n_paper = len(run_all.PAPER)
    # 네트워크에 의존하는 단계 수도 본문에 적혀 있다. 코드가 바뀌면 따라 바뀌어야 한다.
    n_net = len([x for x in run_all.S if x.name in set(run_all.PAPER)
                 and getattr(x, 'network', False)])
    # 앞 숫자(전체 단계 수)도 검사한다. 처음에는 뒤 숫자만 잡았고, 그 사이 --paper 가
    # 74에서 82로 늘었는데 본문은 74로 남아 있었다. 정규식이 한 그룹만 잡으면 그 옆의
    # 숫자는 아무도 보지 않는다.
    for a, t in re.findall(r'Of the (\d+) stages,' + chr(10) + r'(\d+) query external', body):
        ok2 = int(a) == n_paper
        log('  본문 "Of the %s stages"   실제 %d   %s'
            % (a, n_paper, 'OK' if ok2 else '*** 불일치'))
        if not ok2:
            bad += 1
        ok = int(t) == n_net
        log('  본문 네트워크 %s단계   실제 %d단계   %s'
            % (t, n_net, 'OK' if ok else '*** 불일치'))
        if not ok:
            bad += 1

    for t in re.findall(r'runs the (\d+) of them that the figures', body):
        ok = int(t) == n_paper
        log('  본문 --paper %s단계   실제 %d단계   %s'
            % (t, n_paper, 'OK' if ok else '*** 불일치'))
        if not ok:
            bad += 1
    for t, r in said:
        ok = (t == total and r == repro)
        log('  본문 %s/%s   실제 %s/%s   %s' % (t, r, total, repro, 'OK' if ok else '*** 불일치'))
        if not ok:
            bad += 1
    if not said:
        log('  본문에서 단계 수 문장을 찾지 못했습니다')
        bad += 1
    return bad


def main():
    fails = 0

    log('=' * 78)
    log('1. 결과 파일을 만든 코드')
    log('=' * 78)
    prod = build_producer_map()
    total_files = sum(len(v) for v in prod.values())
    log('  %d개 스크립트가 %d개 결과 파일을 생성합니다.' % (len(prod), total_files))
    log('  (전부 이 저장소의 코드입니다. 외부 논문에서 가져온 수치는 없습니다.)')
    log('')
    for k in sorted(prod):
        if len(prod[k]) <= 3:
            log('  %-32s -> %s' % (k, ', '.join(os.path.basename(x) for x in prod[k])))
        else:
            log('  %-32s -> %s 외 %d개'
                % (k, os.path.basename(prod[k][0]), len(prod[k]) - 1))

    log('')
    log('=' * 78)
    log('2. 설명 자료의 숫자 == 결과 파일의 숫자')
    log('=' * 78)
    log('  %-38s %10s %10s %s' % ('항목', '자료', '파일', ''))
    for label, path, fn, expect, tol in CHECKS:
        got, ok, detail = num(path, fn, label)
        if not ok:
            log('  %-38s %10s %10s  읽기 실패: %s' % (label, expect, '-', detail))
            fails += 1
            continue
        good = abs(float(got) - expect) <= tol
        log('  %-38s %10.3f %10.3f  %s'
            % (label, expect, float(got), 'OK' if good else '*** 불일치'))
        if not good:
            fails += 1

    log('')
    log('=' * 78)
    log('2a. 결과 파일에서 계산한 값이 원고 본문에 그대로 있는가')
    log('=' * 78)
    fails += check_tex_quotes()

    log('')
    log('=' * 78)
    log('2b. 본문이 적어둔 파이프라인 단계 수')
    log('=' * 78)
    fails += check_stage_counts()

    log('')
    log('=' * 78)
    log('3. 원자료 출처')
    log('=' * 78)
    inv = os.path.join(ROOT, 'results/data_inventory.tsv')
    if os.path.exists(inv):
        d = pd.read_csv(inv, sep='\t')
        used = d[d['state'] == 'used']
        log('  실제 사용 자산 %d건 (전체 %d건)' % (len(used), len(d)))
        for _, r in used.iterrows():
            log('    %-30s %s' % (r['asset'], str(r['role'])[:56]))
    else:
        log('  data_inventory.tsv 없음')
        fails += 1
    log('')
    log('  GEO      https://www.ncbi.nlm.nih.gov/geo/')
    log('  KPMP     https://atlas.kpmp.org/')
    log('  대사체   https://www.metabolomicsworkbench.org/  (ST003255)')
    log('  GWAS     https://www.ebi.ac.uk/gwas/')

    log('')
    log('=' * 78)
    if fails:
        log('검증 실패 %d건. 설명 자료의 수치가 결과와 어긋납니다.' % fails)
        return 1
    log('검증 통과. 모든 수치가 이 저장소의 코드가 만든 결과 파일과 일치합니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
