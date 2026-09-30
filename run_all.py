#!/usr/bin/env python
"""Single entry point for the whole analysis.

The project accumulated 51 scripts in dependency order that only existed in the author's head.
This file makes that order explicit and executable, so a reader can reproduce any stage or the
whole thing without knowing which script feeds which.

    python run_all.py --list                 show every stage and whether its outputs exist
    python run_all.py --stage harmonise      run one stage
    python run_all.py --from benchmark       run from a stage onward
    python run_all.py --all                  everything, in order
    python run_all.py --paper                only what the paper's figures and tables need
    python run_all.py --dry-run --all        print the commands without running them

Stages declare their outputs, so a stage whose outputs already exist is skipped unless --force
is given. Download stages are separated from analysis stages because they need network access
and take far longer than everything else combined.
"""
import argparse
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable


class Stage:
    def __init__(self, name, desc, cmd, outputs, group='analysis', minutes=1, network=False):
        self.name = name
        self.desc = desc
        self.cmd = cmd
        self.outputs = outputs
        self.group = group
        self.minutes = minutes
        self.network = network

    def done(self):
        # A stage that declares no outputs cannot be verified as complete, so it is never
        # skipped. all([]) is True, which would silently mark such a stage done forever.
        if not self.outputs:
            return False
        return all(os.path.exists(os.path.join(ROOT, o)) for o in self.outputs)


S = [
    # ------------------------------------------------------------------ acquisition
    Stage('download', 'fetch GEO series, KPMP atlas, annotations, GWAS',
          [PY, 'scripts/download_geo.py', 'GSE142025', 'GSE30528', 'GSE96804', 'GSE104948',
           'GSE30529', 'GSE104954', 'GSE294519', 'GSE162830', 'GSE175759', 'GSE5406'],
          ['data/raw/geo/GSE30528/GSE30528_series_matrix.txt.gz'],
          group='acquire', minutes=25, network=True),
    Stage('inventory', 'catalogue what was acquired and what was used',
          [PY, 'scripts/data_inventory.py'], ['results/data_inventory.tsv'],
          group='acquire', minutes=1),

    # ------------------------------------------------------------------ audit and harmonisation
    Stage('audit', 'parse sample metadata; detect patient reuse across accessions',
          [PY, 'scripts/parse_series_matrix.py'], ['data/metadata/sample_inventory.tsv'],
          minutes=1),
    Stage('genespace', 'probe to Entrez maps and the frozen 9,900-gene intersection',
          [PY, 'scripts/build_gene_space.py'],
          ['data/processed/gene_space_all6.tsv', 'data/processed/probe2gene.tsv'], minutes=3),
    Stage('harmonise', 'collapse probes to genes; write expression and phenotype tables',
          [PY, 'scripts/harmonize.py'],
          ['data/processed/harmonized/GSE30528_expr.tsv'], minutes=4),
    Stage('add_cohort5', 'map GSE294519 into the frozen gene space',
          [PY, 'scripts/add_gse294519.py'],
          ['data/processed/harmonized/GSE294519_expr.tsv'], minutes=1),
    Stage('feasibility', 'per-cohort effect sizes and cross-cohort concordance',
          [PY, 'scripts/sanity_effect_sizes.py'],
          ['data/processed/effect_sizes_hedges_g.tsv'], minutes=2),

    # ------------------------------------------------------------------ controls
    Stage('null', 'random-signature and permuted-label nulls per fold and K',
          [PY, 'scripts/null_control.py', '--n-rep', '300'],
          ['results/null_control_draws.tsv', 'results/null_control.tsv'], minutes=12),

    # ------------------------------------------------------------------ benchmark
    Stage('baselines', 'seven selectors under LODO, bootstrap stability selection',
          [PY, 'scripts/run_baselines.py', '--B', '200', '--k', '50',
           '--cohorts', 'GSE30528,GSE96804,GSE104948,GSE142025',
           '--out', 'results/baselines_all4'],
          ['results/baselines_all4/results.tsv'], minutes=120),
    Stage('proposed', 'the cross-cohort stability selector (RBS)',
          [PY, 'scripts/run_proposed.py', '--B', '200', '--base', 'relieff', '--agg', 'geomean',
           '--cohorts', 'GSE30528,GSE96804,GSE104948,GSE142025',
           '--out', 'results/proposed_all4'],
          ['results/proposed_all4/results.tsv'], minutes=4),
    Stage('compare', 'merge benchmark tables and normalise AUROC against the null',
          [PY, 'scripts/compare_methods.py', '--baselines', 'results/baselines_all4',
           '--proposed', 'results/proposed_all4',
           '--cohorts', 'GSE30528,GSE96804,GSE104948,GSE142025',
           '--out', 'results/comparison_all4.tsv'],
          ['results/comparison_all4.tsv'], minutes=1),
    Stage('baselines_glom3', 'the same seven selectors on the three glomerular cohorts',
          [PY, 'scripts/run_baselines.py', '--B', '200', '--k', '50',
           '--cohorts', 'GSE30528,GSE96804,GSE104948',
           '--out', 'results/baselines'],
          ['results/baselines/results.tsv'], minutes=90),
    Stage('proposed_glom3', 'RBS on the three glomerular cohorts',
          [PY, 'scripts/run_proposed.py', '--B', '200', '--base', 'relieff', '--agg', 'geomean',
           '--cohorts', 'GSE30528,GSE96804,GSE104948',
           '--out', 'results/proposed'],
          ['results/proposed/results.tsv'], minutes=4),
    Stage('null_glom3', 'the random-signature null for the three-cohort folds',
          [PY, 'scripts/null_control.py', '--n-rep', '300',
           '--cohorts', 'GSE30528,GSE96804,GSE104948',
           '--out', 'results/null_control_glom3.tsv'],
          ['results/null_control_glom3.tsv',
           'results/null_control_glom3_draws.tsv'], minutes=9),
    Stage('compare_glom3', 'merge the three-cohort tables against their own null',
          [PY, 'scripts/compare_methods.py', '--baselines', 'results/baselines',
           '--proposed', 'results/proposed',
           '--cohorts', 'GSE30528,GSE96804,GSE104948',
           '--null', 'results/null_control_glom3_draws.tsv',
           '--out', 'results/comparison_glom3.tsv'],
          ['results/comparison_glom3.tsv'], minutes=1),
    Stage('conventional', 'DEG meta-analysis and WGCNA-style hub ranking',
          [PY, 'scripts/conventional_pipeline.py'],
          ['results/conventional_pipeline/comparison.tsv'], minutes=6),

    # ------------------------------------------------------------------ uncertainty
    Stage('bootstrap_ci', 'paired patient-level bootstrap intervals on every comparison',
          [PY, 'scripts/bootstrap_ci.py', '--n-boot', '200', '--B-inner', '30',
           '--cohorts', 'GSE30528,GSE96804,GSE104948,GSE142025,GSE294519'],
          ['results/bootstrap_ci/paired_differences.tsv'], minutes=60),
    Stage('cohort_sensitivity', 'repeat the benchmark over all 26 cohort subsets',
          [PY, 'scripts/cohort_sensitivity.py', '--min-size', '2'],
          ['results/cohort_sensitivity/subsets.tsv',
           'results/cohort_sensitivity/by_size.tsv'], minutes=25),
    # 같은 스윕을 다른 내부 선택기로. 불안정성이 선택기 한 쌍의 성질인지 묻는다.
    # 선택기마다 편차가 크다(univariate 1분, boruta 12시간 이상).
    Stage('sens_selectors', 'the 26-subset sweep under each inner selector of stability selection',
          [PY, 'scripts/selector_sensitivity.py'],
          ['results/cohort_sensitivity_by_selector/summary.tsv'], minutes=900),
    # 비교군을 정식 WGCNA 패키지로 바꿔 같은 두 질문을 묻는다. R 4.6 + WGCNA 1.74 가 필요하다.
    Stage('wgcna_ref', 'the reference WGCNA package as comparator: reversals and IEG exclusion',
          [PY, 'scripts/wgcna_reference.py', '--headline', '--sweep'],
          ['results/wgcna_reference/subsets.tsv', 'results/wgcna_reference/headline.tsv',
           'results/wgcna_reference/sweep_folds.tsv'], minutes=240),

    # ------------------------------------------------------------------ the confounder
    Stage('artifact', 'immediate-early module: DKD vs non-diabetic CKD in ERCB',
          [PY, 'scripts/ieg_artifact_test.py'], ['results/ieg_artifact.tsv'], minutes=2),
    Stage('artifact_across', 'module level by procurement across all datasets',
          [PY, 'scripts/ieg_heterogeneity.py'], ['results/ieg_heterogeneity.tsv'], minutes=8),
    Stage('artifact_rnaseq', 'independent replication in GSE175759',
          [PY, 'scripts/validate_gse175759.py'],
          ['results/gse175759_validation.tsv',
           'results/gse175759_validation_effects.tsv'], minutes=6),
    Stage('artifact_matched', 'histology-matched contrast in GSE162830',
          [PY, 'scripts/validate_gse162830.py'], ['results/gse162830_validation.tsv'], minutes=2),
    Stage('artifact_singlecell', 'KPMP single-nucleus with a biopsy-vs-biopsy control',
          [PY, 'scripts/kpmp_singlecell.py'],
          ['results/kpmp_singlecell/pseudobulk_de.tsv',
           'results/kpmp_singlecell/expr_cache.npz'], minutes=20),
    Stage('artifact_heart', 'external organ check: heart failure GSE5406',
          [PY, 'scripts/second_disease.py'], ['results/second_disease/gse5406.tsv'], minutes=4),

    # ------------------------------------------------------------------ mechanism
    Stage('wgcna_validation', 'does S3.4 survive a canonical WGCNA comparator?',
          [PY, 'scripts/wgcna_validation.py'],
          ['results/wgcna_validation/variants.tsv',
           'results/wgcna_validation/soft_threshold.tsv'], minutes=8),
    Stage('gse163603', 'procurement effect where storage is documented and uniform',
          [PY, 'scripts/gse163603_procurement.py'],
          ['results/gse163603/summary.tsv'], minutes=2),
    Stage('storage_test', 'do cohorts that report storage show a smaller procurement effect?',
          [PY, 'scripts/storage_metadata_test.py'],
          ['results/storage_test/per_dataset.tsv',
           'results/storage_test/summary.tsv'], minutes=1),
    Stage('module_size', 'is it a boundary accident or module-size bias?',
          [PY, 'scripts/module_size_test.py', '--reps', '8'],
          ['results/module_size/real_module_membership.tsv'], minutes=6),
    Stage('sim_blocks', 'block-size sweep with intervals',
          [PY, 'scripts/sim_block_size.py', '--reps', '40'],
          ['results/module_size/block_size_ci.tsv',
           'results/module_size/block_size_paired.tsv'], minutes=12),
    Stage('sim_real', 'the same sweep driven by the actual selector',
          [PY, 'scripts/sim_real_rbs.py', '--reps', '25', '--B', '40'],
          ['results/module_size/real_rbs_sweep.tsv'], minutes=25),

    # ------------------------------------------------------------------ corrections
    Stage('deconfound', 'residualisation arms and their cost',
          [PY, 'scripts/dkd_deconfound.py', '--B', '200'],
          ['results/deconfound.tsv', 'results/deconfound_v2.tsv'], minutes=45),
    Stage('ruv', 'RUVg, SVA and ComBat against the handling-score correction',
          [PY, 'scripts/ruv_benchmark.py', '--B', '150'],
          ['results/ruv_benchmark/summary.tsv',
           'results/ruv_benchmark/benchmark.tsv'], minutes=30),
    # 보정이 relieff 에만 드는지 묻는 스윕. 선택기 일곱 개 x 보정 여섯 개를 B=40 으로 돈다.
    # Boruta 하나가 몇 시간을 쓰므로 --paper 에서 가장 긴 단계다. 값의 정밀도가 아니라
    # 보정 전후의 방향을 묻는 것이라 B 를 낮춘 것이고, 본문도 그렇게 적는다.
    Stage('corr_sweep', 'the handling-score correction under each of the seven selectors',
          [PY, 'scripts/correction_sweep.py'],
          ['results/correction_sweep/sweep.tsv'], minutes=360),

    # ------------------------------------------------------------------ candidates
    Stage('candidates', 'specificity-gated candidate list, uncorrected selector',
          [PY, 'scripts/final_candidates.py'], ['results/final_candidates.tsv', 'results/final_candidates_gated.tsv'], minutes=2),
    Stage('candidates_v2', 're-derive with the corrected selector',
          [PY, 'scripts/rederive_candidates.py', '--arm', 'orth', '--B', '200'],
          ['results/candidates_v2/candidates_gated.tsv',
           'results/candidates_v2/candidates_all.tsv'], minutes=5),
    Stage('spec_null', 'permutation false-discovery rate of the specificity gate',
          [PY, 'scripts/specificity_null.py', '--n-perm', '500'],
          ['results/specificity_null.tsv'], minutes=4),
    Stage('per_diagnosis', 'candidate levels across nine diagnoses',
          [PY, 'scripts/per_diagnosis_profile.py'],
          ['results/per_diagnosis_profile.tsv'], minutes=2),
    Stage('literature', 'PubMed novelty counts with ambiguity flags',
          [PY, 'scripts/literature_review.py'], ['results/literature_review.tsv'],
          minutes=3, network=True),
    Stage('kpmp_evidence', 'KPMP cell type and regional proteomics per candidate',
          [PY, 'scripts/kpmp_corroborate.py'], ['results/kpmp_corroboration.tsv'],
          minutes=4, network=True),
    Stage('master_table', 'join every line of evidence and tier the candidates',
          [PY, 'scripts/master_candidate_table.py'],
          ['results/master_candidate_table.tsv'], minutes=1),
    # --- 보정 선택기(v2) 쪽 증거 결합.
    # 원고의 후보 30개는 v2 목록에서 나온다. 그런데 오랫동안 이 네 파일을 손으로 한 번
    # 만들어 두고 그대로 썼고, 파이프라인에는 선언이 없었다. 새로 복제하면 하류가
    # 전부 실패하는 상태였다. 명령은 기존 파일을 그대로 재현하는 것으로 확인했다.
    # PubMed 집계는 시간이 지나면 늘어난다. 원고가 보고하는 것은 동결된 파일이고,
    # 검색일은 literature_ledger.py 가 남긴다. run_all 은 산출물이 있으면 건너뛰므로
    # 이 저장소에서는 동결본이 덮이지 않는다.
    Stage('literature_v2', 'PubMed novelty counts for the corrected selector list',
          [PY, 'scripts/literature_review.py',
           '--candidates', 'results/candidates_v2/candidates_gated.tsv',
           '--out', 'results/candidates_v2/literature_review.tsv'],
          ['results/candidates_v2/literature_review.tsv'], minutes=3, network=True),
    Stage('kpmp_evidence_v2', 'KPMP evidence for the corrected selector list',
          [PY, 'scripts/kpmp_corroborate.py',
           '--candidates', 'results/candidates_v2/candidates_gated.tsv',
           '--out', 'results/candidates_v2/kpmp_corroboration.tsv'],
          ['results/candidates_v2/kpmp_corroboration.tsv'], minutes=4, network=True),
    Stage('per_diagnosis_v2', 'per-diagnosis pattern for the corrected selector list',
          [PY, 'scripts/per_diagnosis_profile.py',
           '--candidates', 'results/candidates_v2/candidates_gated.tsv',
           '--out', 'results/candidates_v2/per_diagnosis_profile.tsv'],
          ['results/candidates_v2/per_diagnosis_profile.tsv'], minutes=2),
    Stage('master_table_v2', 'the candidate table the manuscript reports',
          [PY, 'scripts/master_candidate_table.py',
           '--gated', 'results/candidates_v2/candidates_gated.tsv',
           '--lit', 'results/candidates_v2/literature_review.tsv',
           '--kpmp', 'results/candidates_v2/kpmp_corroboration.tsv',
           '--pattern', 'results/candidates_v2/per_diagnosis_profile.tsv',
           '--out', 'results/candidates_v2/master_candidate_table.tsv'],
          ['results/candidates_v2/master_candidate_table.tsv'], minutes=1),
    Stage('metabolomics', 'second-modality feasibility: DKD vs other CKD in plasma metabolites',
          [PY, 'scripts/metabolomics_feasibility.py'],
          ['results/metabolomics/ST003255_contrasts.tsv',
           'results/metabolomics/ST003255_DKD_vs_other_kidney_disease.tsv'], minutes=2),
    Stage('metabolomics_exp', 'do the transcriptomic findings hold in the second modality?',
          [PY, 'scripts/metabolomics_experiments.py'],
          ['results/metabolomics/experiments_auroc.tsv',
           'results/metabolomics/experiments_k_sensitivity.tsv'], minutes=3),
    Stage('metabolomics_rob', 'two challenges to the second-modality result',
          [PY, 'scripts/metabolomics_robustness.py'],
          ['results/metabolomics/robust_per_disease.tsv'], minutes=2),
    Stage('method_vs_design', 'how much of the outcome is method choice vs cohort choice?',
          [PY, 'scripts/method_vs_design.py'],
          ['results/method_vs_design/variance.tsv'], minutes=1),
    Stage('spec_sweep', 'does our own specificity gate depend on which comparators we pooled?',
          [PY, 'scripts/specificity_comparator_sweep.py'],
          ['results/specificity_sweep/candidate_consistency.tsv'], minutes=4),
    Stage('sens_no_ercb', 'candidates re-derived with GSE104948 dropped from discovery',
          [PY, 'scripts/rederive_candidates.py', '--arm', 'orth', '--B', '200',
           '--cohorts', 'GSE30528,GSE96804,GSE142025',
           '--out', 'results/sensitivity_no104948'],
          ['results/sensitivity_no104948/candidates_gated.tsv',
           'results/sensitivity_no104948/candidates_all.tsv'], minutes=6),
    Stage('sens_no_ercb_check', 'does the specificity conclusion survive that removal',
          [PY, 'scripts/sensitivity_drop_ercb.py'],
          ['results/sensitivity_no104948/summary.tsv'], minutes=2),
    Stage('external_val', 'stress-test the candidates in three never-used datasets',
          [PY, 'scripts/candidate_external_validation.py'],
          ['results/external_validation/summary.tsv'], minutes=3),
    Stage('design_survey', 'how often do kidney metabolomics studies hold disease constant?',
          [PY, 'scripts/comparator_design_survey.py'],
          ['results/design_survey/kidney_metabolomics_designs.tsv'], minutes=4, network=True),
    Stage('kpmp_overlap', 'do the same KPMP participants carry more than one omics layer?',
          [PY, 'scripts/kpmp_modality_overlap.py'],
          ['results/kpmp_overlap/participant_layers.tsv',
           'results/kpmp_overlap/layer_summary.tsv'], minutes=3, network=True),
    Stage('ckd_general', 'is the control-group problem specific to DKD or general to kidney biopsy?',
          [PY, 'scripts/ckd_generalization.py'],
          ['results/ckd_generalization/per_diagnosis.tsv'], minutes=4),
    Stage('neuropathy', 'does the confounder-matched design generalise to another complication?',
          [PY, 'scripts/neuropathy_generalization.py'],
          ['results/neuropathy/contrasts.tsv'], minutes=3),
    Stage('supp_stress', 'small-sample stress tests, kept as supplementary by design',
          [PY, 'scripts/supplementary_stress.py'],
          ['results/external_validation/supplementary.tsv'], minutes=2),
    Stage('dossier', 'evidence-tiered candidate record for the R&D log',
          [PY, 'scripts/candidate_dossier.py'],
          ['results/candidates_v2/candidate_dossier.tsv'], minutes=1),
    Stage('kegg_pathway', 'rule-based joint pathway test, replacing the ad hoc collagen link',
          [PY, 'scripts/kegg_joint_pathway.py'],
          ['results/kegg_pathway/summary.tsv'], minutes=6, network=True),
    Stage('reactome_pathway', 'the same joint-pathway test under a second ontology',
          [PY, 'scripts/reactome_joint_pathway.py'],
          ['results/reactome_pathway/summary.tsv'], minutes=8, network=True),
    Stage('crossomics', 'pathway-level link between transcriptomic candidates and metabolites',
          [PY, 'scripts/crossomics_link.py'],
          ['results/metabolomics/crossomics_collagen.tsv'], minutes=1),
    Stage('reordering', 'does swapping the control comparator rebuild the top of the ranking',
          [PY, 'scripts/comparator_reordering.py'],
          ['results/comparator_reordering/top_overlap.tsv',
           'results/comparator_reordering/summary.tsv',
           'results/comparator_reordering/candidate_rank_shift.tsv'], minutes=2),
    Stage('gwas', 'germline test with SNP-count-matched nulls',
          [PY, 'scripts/gwas_layer.py'], ['results/gwas_layer.tsv'], minutes=6),

    # ------------------------------------------------------------------ robustness
    Stage('stress', 'perturbation stress test of the finished signature',
          [PY, 'scripts/stress_test.py'], ['results/stress_test.tsv'], minutes=8),
    Stage('generative', 'VAE and diffusion synthetic populations',
          [PY, 'scripts/generative_robustness.py'],
          ['results/generative_robustness.tsv'], minutes=20),
    Stage('per_dataset', 'single-dataset feature catalogue over eleven datasets',
          [PY, 'scripts/per_dataset_features.py', '--B', '60'],
          ['results/per_dataset/consensus_across_datasets.tsv'], minutes=15),

    # ------------------------------------------------------------------ outputs
    Stage('database', 'build the hybrid relational + columnar database',
          [PY, 'scripts/build_db.py', '--reset'], ['db/dkd.sqlite'], minutes=3),
    Stage('gate_fdr', 'attach background rate and gate FDR to every pass-rate table',
          [PY, 'scripts/propagate_gate_fdr.py'],
          ['results/gate_characteristics.tsv'], minutes=1),
    Stage('overview_fig', 'the framework overview figure (Figure 1)',
          [PY, 'scripts/overview_figure.py'],
          ['results/figures/P0_overview.pdf'], minutes=1),
    Stage('cohort_table', 'the main-text table of cohorts used (Table 1)',
          [PY, 'scripts/cohort_table.py'],
          ['results/cohort_table.tsv', 'results/subject_overlap.tsv',
           'submission/table_cohorts.tex',
           'submission/table_overlap.tex'], minutes=1),
    Stage('proteome_val', 'do the candidates replicate at protein level (Mendeley 83k89shdx5)',
          [PY, 'scripts/proteome_validation.py'],
          ['results/proteome_validation/summary.tsv',
           'results/proteome_validation/protein_dkd_vs_control.tsv'], minutes=2),
    Stage('pride_fetch', 'fetch the second proteomic cohort (PRIDE PXD041884)',
          [PY, 'scripts/fetch_pride_proteome.py'],
          ['data/raw/proteomics/pride_PXD041884/uniprot_to_gene.tsv'],
          group='acquire', minutes=4, network=True),
    Stage('proteome_ms', 'do the candidates replicate by mass spectrometry too',
          [PY, 'scripts/proteome_ms_validation.py'],
          ['results/proteome_ms/summary.tsv',
           'results/proteome_ms/protein_ms_dkd_vs_control.tsv'], minutes=2),
    Stage('proteome_meta', 'pool the two proteomic cohorts within strata',
          [PY, 'scripts/proteome_meta.py'],
          ['results/proteome_meta/summary.tsv',
           'results/proteome_meta/matching_balance.tsv',
           'results/proteome_meta/candidate_protein_coverage.tsv'], minutes=2),
    Stage('proteome_disc', 'can the protein layer discover on its own (negative control)',
          [PY, 'scripts/proteome_discovery_test.py'],
          ['results/proteome_meta/discovery_summary.tsv',
           'results/proteome_meta/discovery_attempt.tsv'], minutes=2),
    Stage('proteome_audit', 'what each proteomic source states about its cohorts',
          [PY, 'scripts/proteome_cohort_audit.py'],
          ['results/proteome_meta/cohort_audit.tsv'], minutes=1),
    Stage('multiomics_feas', 'which layer pairs can form a case/control contrast at all',
          [PY, 'scripts/multiomics_feasibility.py'],
          ['results/multiomics_feasibility/summary.tsv'], minutes=1),
    Stage('proteome_identity', 'did we measure the protein of that gene (aptamer, peptides)',
          [PY, 'scripts/proteome_identity_audit.py'],
          ['results/proteome_meta/protein_identity_audit.tsv'], minutes=2),
    Stage('freeze_audit', 'was the candidate list fixed before the proteomic results',
          [PY, 'scripts/candidate_freeze_audit.py'],
          ['results/proteome_meta/freeze_summary.tsv'], minutes=1),
    Stage('provenance_gap', 'is every result file something a stage actually produces',
          [PY, 'scripts/provenance_gap_audit.py'],
          ['results/proteome_meta/provenance_gap.tsv'], minutes=1),
    Stage('paper_closure', 'does --paper cover everything the manuscript depends on',
          [PY, 'scripts/paper_closure_audit.py'], [], minutes=1),
    Stage('prose_audit', 'sentence-level check of the typeset manuscript',
          [PY, 'scripts/prose_audit.py'], [], minutes=1),
    Stage('gene_claims', 'do named gene claims in the text match the candidate table',
          [PY, 'scripts/gene_claim_audit.py'], [], minutes=1),
    Stage('lit_ledger', 'per-candidate PubMed query, date and PMIDs behind each tier',
          [PY, 'scripts/literature_ledger.py'],
          ['results/proteome_meta/literature_ledger.tsv',
           'results/proteome_meta/literature_ledger_summary.tsv'], minutes=3, network=True),
    Stage('proteome_table', 'the main-text table of protein-level replication',
          [PY, 'scripts/proteome_table.py'], ['submission/table_proteome.tex'], minutes=1),
    Stage('coexpression', 'how correlated the top-relevance genes are (backs the null result)',
          [PY, 'scripts/coexpression_check.py'], ['results/coexpression.tsv'], minutes=2),
    Stage('comparator_table', 'the WGCNA correspondence table (main text)',
          [PY, 'scripts/comparator_table.py'], ['submission/table_comparator.tex'], minutes=1),
    Stage('bench_protocol', 'winner frequency and metric agreement across cohort subsets',
          [PY, 'scripts/benchmark_protocol.py'],
          ['results/cohort_sensitivity/winner_frequency.tsv',
           'results/cohort_sensitivity/metric_agreement.tsv'], minutes=1),
    Stage('reordering_fig', 'the comparator-reordering figure',
          [PY, 'scripts/reordering_figure.py'],
          ['results/figures/P9_reordering.pdf'], minutes=1),
    Stage('candidate_fig', 'the candidate evidence figure',
          [PY, 'scripts/candidate_figure.py'], ['results/figures/P7_candidates.pdf',
           'results/figures/P7_candidates_full.pdf'], minutes=1),
    # ------------------------------------------------------------------ cross-tissue extension
    # 간 MASH 와 대장 UC 로 논문의 두 주장을 가르는 사전 등록 검정. 예측과 사후 계획은
    # docs/xtissue 에 실행 전에 고정했고, xt_registration 이 그 해시를 매번 확인한다.
    Stage('xt_registration', 'pre-registered predictions and post hoc plan are unchanged',
          [PY, 'scripts/xtissue_registration_check.py'], [], minutes=1),
    Stage('xt_candidates', 'series-level metadata for candidate liver and colon cohorts',
          [PY, 'scripts/xtissue_candidates.py'],
          ['results/xtissue/candidates/liver_masld.json',
           'results/xtissue/candidates/colon_uc.json'], minutes=4, network=True),
    Stage('xt_harmonise_colon', 'download and harmonise the colon UC cohorts',
          [PY, 'scripts/xtissue_harmonize.py', '--domain', 'colon_uc'],
          ['results/xtissue/colon_uc/cohorts.tsv'], minutes=12, network=True),
    Stage('xt_harmonise_liver', 'download and harmonise the liver MASH cohorts',
          [PY, 'scripts/xtissue_harmonize.py', '--domain', 'liver_masld'],
          ['results/xtissue/liver_masld/cohorts.tsv'], minutes=10, network=True),
    Stage('xt_overlap', 'value-vector audit for shared samples across same-platform series',
          [PY, 'scripts/xtissue_overlap_audit.py'],
          ['results/xtissue/overlap_audit.tsv'], minutes=6, network=True),
    Stage('xt_p2', 'prediction P2: immediate-early effect by cohort and procurement',
          [PY, 'scripts/xtissue_p2.py'],
          ['results/xtissue/p2_ieg_by_cohort.tsv', 'results/xtissue/p2_verdict.tsv'], minutes=4),
    Stage('xt_sweep_colon', 'prediction P1/P3: cohort-composition sweep, colon',
          [PY, 'scripts/xtissue_sweep.py', '--tissue', 'colon_uc'],
          ['results/xtissue/colon_uc/sweep/by_size.tsv',
           'results/xtissue/colon_uc/sweep/subsets.tsv'], minutes=30),
    Stage('xt_sweep_liver', 'prediction P1/P3: cohort-composition sweep, liver',
          [PY, 'scripts/xtissue_sweep.py', '--tissue', 'liver_masld'],
          ['results/xtissue/liver_masld/sweep/by_size.tsv',
           'results/xtissue/liver_masld/sweep/subsets.tsv'], minutes=25),
    Stage('xt_sweep_kidney', 'prediction P3: kidney reference with the same code',
          [PY, 'scripts/xtissue_sweep.py', '--tissue', 'kidney_dkd', '--full-only'],
          ['results/xtissue/kidney_dkd/sweep/subsets.tsv'], minutes=5),
    Stage('xt_verdicts', 'verdicts on P1-P3 against the pre-registered wording',
          [PY, 'scripts/xtissue_verdicts.py'], ['results/xtissue/verdicts.tsv'], minutes=1),
    Stage('xt_h4', 'post hoc H4: control-definition substitution in liver GSE48452',
          [PY, 'scripts/xtissue_posthoc.py', 'h4'],
          ['results/xtissue/posthoc/h4_gse48452_summary.tsv'], minutes=2),
    Stage('xt_h5', 'post hoc H5: inflammation adjustment of the immediate-early effect',
          [PY, 'scripts/xtissue_posthoc.py', 'h5'],
          ['results/xtissue/posthoc/h5_inflammation_adjustment.tsv'], minutes=2),
    Stage('xt_h6', 'post hoc H6: colon downsampled to kidney cohort sizes',
          [PY, 'scripts/xtissue_posthoc.py', 'h6', '--reps', '5'],
          ['results/xtissue/posthoc/h6/replicates.tsv'], minutes=50),
    Stage('xt_s1', 'post hoc S1: relaxed RNA-seq filter for liver',
          [PY, 'scripts/xtissue_posthoc.py', 's1'],
          ['data/processed/xtissue/liver_masld_relaxed/GSE126848_expr.tsv'], minutes=2),
    Stage('xt_s1_sweep', 'post hoc S1: liver P3 count under the relaxed filter',
          [PY, 'scripts/xtissue_sweep.py', '--tissue', 'liver_masld', '--full-only', '--root',
           'data/processed/xtissue/liver_masld_relaxed', '--tag', 'sweep_relaxed'],
          ['results/xtissue/liver_masld/sweep_relaxed/subsets.tsv'], minutes=5),
    Stage('xt_concordance', 'post hoc: cross-cohort effect concordance per tissue',
          [PY, 'scripts/xtissue_posthoc.py', 'conc'],
          ['results/xtissue/posthoc_effect_concordance.tsv'], minutes=3),
    Stage('xt_seed_stability', 'post hoc: is the small liver sign reversal seed-dependent',
          [PY, 'scripts/xtissue_seed_stability.py', '--tissue', 'liver_masld', '--seeds', '5'],
          ['results/xtissue/liver_masld/seed_stability.tsv',
           'results/xtissue/liver_masld/seed_stability_summary.tsv'], minutes=55),
    Stage('xt_figure', 'the cross-tissue figure',
          [PY, 'scripts/xtissue_figure.py'], ['results/figures/P10_cross_tissue.pdf',
           'results/figures/P10_cross_tissue_summary.pdf'], minutes=1),

    Stage('stage_figures', 'renumber figures to citation order and stage them as Fig1..N',
          [PY, 'scripts/stage_figures.py'],
          ['submission/Fig%d.pdf' % k for k in range(1, 11)], minutes=1),
    Stage('additional_files', 'the Additional files the manuscript promises',
          [PY, 'scripts/additional_files.py'],
          ['submission/Additional_file_1.pdf', 'submission/Additional_file_2.pdf',
           'submission/Additional_file_3.xlsx', 'submission/Additional_file_4.xlsx',
           'submission/Additional_file_5.pdf', 'submission/Additional_file_6.pdf',
           'submission/Additional_file_7.xlsx'],
          minutes=2),
    Stage('verify_refs', 'every reference checked against CrossRef (no fabricated citations)',
          [PY, 'scripts/verify_references.py'], [], minutes=2),
    Stage('submission_audit', 'BMC format, numbering, references and supplementary linkage',
          [PY, 'scripts/submission_audit.py'], [], minutes=1),
    Stage('paper_figures', 'the manuscript figures P0-P6',
          [PY, 'scripts/paper_figures.py'],
          ['results/figures/P0_framework_overview.png'], minutes=3),
    Stage('figures', 'exploratory figure set F1-F12 (not in the manuscript)',
          [PY, 'scripts/make_figures.py'],
          ['results/figures/F2_cross_fold_stability.png'], minutes=3),
    Stage('kpmp_donor_meta', 'donor-level diabetes history, sex and age from the atlas',
          [PY, 'scripts/kpmp_donor_meta.py'],
          ['results/kpmp_overlap/donor_meta.tsv'], minutes=3),
    Stage('omics_catalog', 'browsable catalogue of what each resource measures',
          [PY, 'scripts/build_omics_catalog.py'], ['db/omics_catalog.sqlite'], minutes=1),
    Stage('multiomics_db', 'multi-omics schema extension plus the Excel template',
          [PY, 'scripts/build_multiomics_db.py'], ['db/multiomics_schema.xlsx'], minutes=1),
    Stage('kr_docs', 'Korean explainer deck and Word document',
          [PY, 'scripts/make_kr_docs.py'], ['DKD_설명자료.pptx'], minutes=1),
    Stage('project_deck', '과제 발표 장표 (DKD발표.pptx 템플릿 서식)',
          [PY, 'scripts/make_project_deck.py'], ['DKD발표.pptx'],
          group='analysis', minutes=1),
    Stage('deck', 'status deck',
          [PY, 'scripts/make_ppt.py'], ['DKD_status.pptx'], minutes=2),

    # ------------------------------------------------------------------ exploratory
    # Run during the work and kept reproducible, but no manuscript claim rests on them.
    # Excluded from --all and --paper; reachable with --stage.
    Stage('x_ablation', 'RBS component ablation',
          [PY, 'scripts/ablation_rbs.py'], ['results/ablation/ablation.tsv'],
          group='extra', minutes=10),
    Stage('x_enrich', 'GO and KEGG over-representation of the candidate list',
          [PY, 'scripts/enrich.py'], ['results/enrichment/KEGG_2021_Human.tsv'],
          group='extra', minutes=3, network=True),
    Stage('x_sim_first', 'the first simulation, superseded by sim_blocks',
          [PY, 'scripts/simulation.py'], ['results/simulation/simulation.tsv'],
          group='extra', minutes=10),
    Stage('x_validate_cand', 'early candidate validation, superseded by master_table',
          [PY, 'scripts/validate_candidates.py'], [], group='extra', minutes=4),
    Stage('x_fig_v1_v2', 'uncorrected vs corrected selector comparison figure',
          [PY, 'scripts/fig_v1_v2.py'], [], group='extra', minutes=2),
    Stage('x_fill', 'write submission_config.py into the manuscript and cover letter',
          [PY, 'scripts/fill_submission.py', '--check'], [], group='extra', minutes=1),
    Stage('x_verify', 'trace every quoted number back to the code that produced it',
          [PY, 'scripts/verify_provenance.py'], [], group='extra', minutes=1),
    Stage('x_figtype', 'figure text sizes in print, from the insertion scale in each document',
          [PY, 'scripts/figure_typography.py', '--min', '6'], [], group='extra', minutes=1),
    Stage('x_preflight', 'scan the submission package for placeholders and internal notes',
          [PY, 'scripts/preflight.py'], [], group='extra', minutes=1),
    Stage('x_bib', 'Briefings in Bioinformatics version: condensed text, supplement, number trace',
          [PY, 'scripts/bib_build.py'], ['submission_bib/bib-manuscript.pdf',
           'submission_bib/bib-supplement.pdf'], group='extra', minutes=2),
    Stage('x_release', 'stage the public release tree under release/ (publishes nothing)',
          [PY, 'scripts/make_release.py'], ['release/CITATION.cff'], group='extra', minutes=2),
]

# 원고 산출물과 원고 수치의 의존성 폐쇄집합. 손으로 고치지 말고
# scripts/paper_closure_audit.py 가 시키는 대로 맞춘다.
PAPER = ['inventory', 'audit', 'genespace', 'harmonise', 'add_cohort5',
         'null', 'baselines', 'proposed', 'compare', 'conventional',
         'baselines_glom3', 'proposed_glom3', 'null_glom3', 'compare_glom3', 'stress',
         'bootstrap_ci', 'cohort_sensitivity', 'sens_selectors', 'wgcna_ref', 'artifact', 'artifact_across', 'artifact_rnaseq',
         'artifact_matched', 'artifact_singlecell', 'artifact_heart', 'wgcna_validation', 'gse163603',
         'storage_test', 'module_size', 'sim_blocks', 'sim_real', 'deconfound',
         'ruv', 'corr_sweep', 'candidates', 'candidates_v2', 'spec_null', 'per_diagnosis',
         'literature', 'kpmp_evidence', 'master_table', 'literature_v2', 'kpmp_evidence_v2',
         'per_diagnosis_v2', 'master_table_v2', 'metabolomics', 'metabolomics_exp', 'external_val',
         'sens_no_ercb', 'sens_no_ercb_check',
         'kpmp_overlap', 'kpmp_donor_meta', 'ckd_general', 'reordering', 'reordering_fig', 'neuropathy', 'kegg_pathway', 'reactome_pathway',
         'gwas', 'gate_fdr', 'bench_protocol', 'overview_fig', 'cohort_table', 'proteome_val',
         'proteome_ms', 'proteome_meta', 'proteome_disc', 'proteome_audit', 'multiomics_feas',
         'proteome_identity', 'freeze_audit', 'gene_claims',
         'prose_audit', 'paper_closure', 'lit_ledger', 'proteome_table',
         'coexpression', 'comparator_table', 'candidate_fig',
         'xt_registration', 'xt_candidates', 'xt_harmonise_colon', 'xt_harmonise_liver',
         'xt_overlap', 'xt_p2', 'xt_sweep_colon', 'xt_sweep_liver', 'xt_sweep_kidney',
         'xt_verdicts', 'xt_h4', 'xt_h5', 'xt_h6', 'xt_s1', 'xt_s1_sweep', 'xt_concordance',
         'xt_figure', 'xt_seed_stability', 'stage_figures', 'additional_files',
         'verify_refs', 'submission_audit', 'paper_figures', 'figures', 'x_verify']


# BMC 형식 원고(submission/dkd-manuscript.tex)를 읽는 단계. 공개 릴리스에는 투고본(BiB)만
# 들어가므로 그 파일이 없고, 이 단계들은 실패가 아니라 건너뛴다. 결과·그림 생성과는 무관하다.
NEEDS_MANUSCRIPT = {'prose_audit', 'gene_claims', 'stage_figures', 'additional_files',
                    'submission_audit', 'x_verify', 'x_fill', 'x_figtype', 'x_bib', 'x_release'}


def run(stage, dry, force):
    if stage.done() and not force:
        print('  SKIP  %-20s outputs present' % stage.name)
        return True
    if stage.name in NEEDS_MANUSCRIPT and not os.path.exists(
            os.path.join(ROOT, 'submission', 'dkd-manuscript.tex')):
        print('  SKIP  %-20s manuscript audit; the public release carries the submitted '
              'version only' % stage.name)
        return True
    print('  RUN   %-20s %s  (~%d min)' % (stage.name, stage.desc, stage.minutes))
    if dry:
        print('        %s' % ' '.join(stage.cmd))
        return True
    t0 = time.time()
    # Stages print their diagnostics to stderr and are long-running, so the console stays the
    # primary view; the log is a copy kept for the stages that scroll past.
    logdir = os.path.join(ROOT, 'logs')
    os.makedirs(logdir, exist_ok=True)
    logpath = os.path.join(logdir, '%s.log' % stage.name)
    with open(logpath, 'w', encoding='utf-8') as fh:
        fh.write('$ %s\n\n' % ' '.join(stage.cmd))
        fh.flush()
        r = subprocess.run(stage.cmd, cwd=ROOT, stdout=fh,
                           stderr=subprocess.STDOUT)
    dt = (time.time() - t0) / 60
    if r.returncode != 0:
        print('  FAIL  %-20s exit %d after %.1f min -- see logs/%s.log'
              % (stage.name, r.returncode, dt, stage.name))
        return False
    print('  OK    %-20s %.1f min' % (stage.name, dt))
    return True




# ---------------------------------------------------------------------------- 제출 고정
# 새 검사를 만들지 않는다. 있는 것을 정해진 순서로 돌리고 해시를 남긴다.
SUBMIT_STEPS = [
    ('BiB 판 빌드', [PY, 'scripts/bib_build.py']),
    ('BMC 제출 점검', [PY, 'scripts/preflight.py', '--target', 'bmc']),
    ('BiB 제출 점검', [PY, 'scripts/preflight.py', '--target', 'bib']),
    ('그림 글자 크기', [PY, 'scripts/figure_typography.py', '--min', '6']),
    ('쪽 배치', [PY, 'scripts/page_layout_audit.py']),
    ('BMC 형식', [PY, 'scripts/submission_audit.py']),
    ('수치 출처', [PY, 'scripts/verify_provenance.py']),
    ('문장 감사', [PY, 'scripts/prose_audit.py']),
    ('재생성 폐쇄집합', [PY, 'scripts/paper_closure_audit.py']),
]
SUBMIT_FILES = ['submission/dkd-manuscript.pdf', 'submission/dkd-manuscript.tex',
                'submission/dkd-references.bib', 'submission/COVER_LETTER.md',
                'submission_bib/bib-manuscript.pdf', 'submission_bib/bib-manuscript.tex',
                'submission_bib/bib-supplement.pdf', 'submission_bib/bib-supplement.tex',
                'submission_bib/COVER_LETTER_BiB.md']


def submit_check():
    """전체 절차를 돌리고, 전부 통과하면 제출 파일의 해시를 기록한다."""
    import glob
    import hashlib
    import time
    bad, blockers = [], []
    for label, cmd in SUBMIT_STEPS:
        r = subprocess.run(cmd, env=dict(os.environ, PYTHONIOENCODING='utf-8'),
                           capture_output=True)
        out = (r.stdout + r.stderr).decode('utf-8', 'replace')
        tail = [l for l in out.splitlines() if l.strip()][-1:] or ['(출력 없음)']
        print('%-16s %-4s %s' % (label, 'OK' if r.returncode == 0 else '***', tail[0][:96]))
        if r.returncode != 0:
            bad.append(label)
            # 무엇이 막는지 한 줄씩 보여 준다. 마지막 줄만 찍으면 저자 기입 항목만
            # 막는 것처럼 읽힌다. 폴더에 남은 낯선 파일도 막는다.
            # 파일별 요약 줄("... 3 internal marker(s) ...")은 항목이 아니다. 그것까지 세면
            # 출력 줄 수와 실제 항목 수가 어긋난다(16 대 12).
            keys = ('placeholder  ', 'INTERNAL', 'COMMENT ', 'unexpected')
            for line in out.splitlines():
                if any(k in line for k in keys) and 'NOT READY' not in line                         and 'placeholder(s)' not in line:
                    blockers.append('%s: %s' % (label, line.strip()[:96]))
    files = [f for f in SUBMIT_FILES if os.path.exists(f)]
    files += sorted(glob.glob('submission/Fig*.pdf') + glob.glob('submission/Additional_file_*')
                    + glob.glob('submission_bib/figures/*.pdf')
                    + glob.glob('submission_bib/Supplementary_File_*'))
    state = '최종 (모든 검사 통과)' if not bad else '잠정 (검사 미통과 - 제출본이 아니다)'
    lines = ['# 제출본 해시 %s' % time.strftime('%Y-%m-%dT%H:%M:%S'),
             '# 상태: ' + state,
             '# 이 해시는 파일이 같은지만 확인한다. 저자 기여·승인·투고 적합성을 증명하지',
             '# 않는다. 파일을 한 자라도 고치면 다시 기록해야 한다.', '']
    for f in files:
        h = hashlib.sha256(open(f, 'rb').read()).hexdigest()
        lines.append('%s  %s' % (h, f.replace(chr(92), '/')))
    os.makedirs('reports', exist_ok=True)
    final = 'reports/submission_manifest.txt'
    path = final if not bad else 'reports/submission_manifest_provisional.txt'
    open(path, 'w', encoding='utf-8').write(chr(10).join(lines) + chr(10))
    print('')
    if bad:
        print('제출 고정 못 함: %s' % ', '.join(bad))
        print('막는 항목 %d건:' % len(blockers))
        for b in blockers[:24]:
            print('  - %s' % b)
        print('해시는 %s 에 잠정으로만 적었습니다. 제출본이 아닙니다.' % path)
        if os.path.exists(final):
            print('이전 %s 는 지금 파일과 다를 수 있으므로 그대로 쓰지 마십시오.' % final)
        return 1
    prov = 'reports/submission_manifest_provisional.txt'
    if os.path.exists(prov):
        os.remove(prov)
    print('%d개 파일의 해시를 %s 에 적었습니다. 이 버전을 제출본으로 고정합니다.'
          % (len(files), path))
    print('해시는 파일 동일성만 보증합니다. 저자 기여와 승인은 저자가 확인할 몫입니다.')
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--stage', action='append', default=[])
    ap.add_argument('--from', dest='from_', default=None)
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--paper', action='store_true')
    ap.add_argument('--skip-network', action='store_true')
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--submit', action='store_true',
                    help='제출 직전 점검을 순서대로 돌리고 제출 파일 해시를 적는다')
    args = ap.parse_args()

    if args.submit:
        return submit_check()

    by_name = {s.name: s for s in S}

    if args.list or not (args.all or args.paper or args.stage or args.from_ or args.submit):
        print('%-20s %-9s %-6s %-5s %s' % ('stage', 'group', 'status', 'min', 'description'))
        total = 0
        for s in S:
            if s.group != 'extra':
                total += 0 if s.done() else s.minutes
            print('%-20s %-9s %-6s %-5d %s'
                  % (s.name, s.group, 'done' if s.done() else '-', s.minutes, s.desc))
        n_spine = sum(1 for x in S if x.group != 'extra')
        print('\n%d stages: %d on the reproduction path, %d exploratory. Roughly %d minutes outstanding.'
              % (len(S), n_spine, len(S) - n_spine, total))
        print('Network needed for: %s' % ', '.join(s.name for s in S if s.network))
        return

    if args.all:
        # 'extra' stages are superseded or side investigations; --all is the reproduction
        # path for the reported results, so they are opt-in via --stage.
        todo = [s for s in S if s.group != 'extra']
    elif args.paper:
        # PAPER 목록의 순서대로 돌리면 목록에 단계를 넣을 때마다 사람이 자리를
        # 맞춰야 하고, 틀리면 의존이 어긋난 채로 돈다. 선언 순서(S)가 이미 의존
        # 순서이므로 그것을 따르고, PAPER 은 포함 여부만 정한다.
        want = set(PAPER)
        todo = [st for st in S if st.name in want]
        unknown = want - {st.name for st in S}
        if unknown:
            print('S 에 없는 PAPER 단계 이름: %s' % ', '.join(sorted(unknown)))
            return 1
        # PAPER 은 코드에 적힌 정적 목록이다. 그것이 원고 산출물의 의존성 폐쇄집합과
        # 같은지 여기서, 아무것도 돌리기 전에 확인한다. 감사를 목록 안에 넣어만 두면
        # 58번째에서야 걸리고, 그 전에 57단계가 이미 돈다.
        chk = subprocess.run([PY, 'scripts/paper_closure_audit.py'],
                             cwd=ROOT, capture_output=True)
        if chk.returncode != 0:
            sys.stderr.write(chk.stderr.decode('utf-8', 'replace'))
            print('PAPER 이 폐쇄집합과 어긋납니다. 위 감사 결과대로 맞추고 다시 도세요.')
            return 1
    elif args.from_:
        spine = [s for s in S if s.group != 'extra']
        names = [s.name for s in spine]
        if args.from_ not in names:
            sys.exit('unknown stage: %s' % args.from_)
        todo = spine[names.index(args.from_):]
    else:
        for n in args.stage:
            if n not in by_name:
                sys.exit('unknown stage: %s' % n)
        todo = [by_name[n] for n in args.stage]

    if args.skip_network:
        todo = [s for s in todo if not s.network]

    print('running %d stages\n' % len(todo))
    failed = []
    for s in todo:
        if not run(s, args.dry_run, args.force):
            failed.append(s.name)
    print('\n%d/%d succeeded' % (len(todo) - len(failed), len(todo)))
    if failed:
        print('failed: %s' % ', '.join(failed))
        sys.exit(1)


if __name__ == '__main__':
    # main() 안의 return 1 이 종료 코드로 이어지지 않아, 폐쇄집합 검증이 막아도
    # 셸에서는 성공으로 보였다. 자동화가 그것을 믿으면 어긋난 채로 계속 간다.
    sys.exit(main() or 0)
