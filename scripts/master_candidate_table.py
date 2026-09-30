#!/usr/bin/env python
"""Join every line of evidence into one candidate table, and answer the deliverable question.

Columns come from six independent analyses:
  RBS rank / components         results/proposed_all4/  (cross-cohort selection)
  specificity                   results/final_candidates.tsv (DKD vs 310 non-diabetic CKD)
  above permutation ceiling     results/specificity_null.tsv (|g| > 0.92)
  per-diagnosis pattern         results/per_diagnosis_profile.tsv (A unique / B dominant / C donor-driven)
  literature verdict            results/literature_review.tsv (PubMed counts, verified titles)
  KPMP cross-omics              results/kpmp_corroboration.tsv (cell type + regional proteomics)

The deliverable asks for a NOVEL candidate. That requires a gene to be strong in the data AND
absent from the DKD literature, so the table is sorted to make that intersection visible rather
than asserting it.
"""
import os, sys, argparse
import numpy as np
import pandas as pd


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='results/master_candidate_table.tsv')
    ap.add_argument('--gated', default='results/final_candidates_gated.tsv')
    ap.add_argument('--lit', default='results/literature_review.tsv')
    ap.add_argument('--kpmp', default='results/kpmp_corroboration.tsv')
    ap.add_argument('--pattern', default='results/per_diagnosis_profile.tsv')
    args = ap.parse_args()

    gated = pd.read_csv(args.gated, sep='\t')
    gated = gated.rename(columns={gated.columns[0]: 'entrez_id'})
    lit = pd.read_csv(args.lit, sep='\t')
    kpmp = pd.read_csv(args.kpmp, sep='\t')
    perm = pd.read_csv('results/specificity_null.tsv', sep='\t')
    ceiling = float(np.percentile(perm['perm_max_abs_g'], 95))

    pat = pd.read_csv(args.pattern, sep='\t')
    pat_glom = pat[pat['cohort'] == 'GSE104948'].set_index('gene')['pattern']
    pat_tub = pat[pat['cohort'] == 'GSE104954'].set_index('gene')['pattern']

    df = gated[['symbol', 'mean_rank', 'mean_S', 'mean_R', 'mean_P',
                'g_DKD_vs_control', 'g_DKD_vs_otherCKD']].copy()
    df = df.rename(columns={'symbol': 'gene'})
    df['above_perm_ceiling'] = df['g_DKD_vs_otherCKD'].abs() > ceiling
    df['pattern_glom'] = df['gene'].map(pat_glom)
    df['pattern_tubulo'] = df['gene'].map(pat_tub)
    df['donor_driven'] = (df['pattern_glom'].fillna('').str.contains('donor') |
                          df['pattern_tubulo'].fillna('').str.contains('donor'))
    df['pattern_class'] = df['pattern_glom'].fillna('').str.slice(0, 1)

    df = df.merge(lit[['gene', 'n_dkd', 'n_dkd_biomarker', 'verdict', 'ambiguity_note']],
                  on='gene', how='left')
    df = df.merge(kpmp[['symbol', 'sn_top_celltype', 'sn_top_foldchange',
                        'protein_glom_vs_TI_FC', 'protein_adjP']].rename(columns={'symbol': 'gene'}),
                  on='gene', how='left')

    df['lit_novel'] = df['verdict'].isin(['NOVEL - no kidney literature',
                                          'kidney-studied, NOT in DKD'])
    df['data_strong'] = (df['pattern_class'].isin(['A', 'B'])) & (~df['donor_driven'])

    # tier
    def tier(r):
        if r['data_strong'] and r['lit_novel']:
            return '1 - novel AND strong'
        if r['data_strong']:
            return '2 - strong, already reported'
        if r['lit_novel']:
            return '3 - novel, but weak in data'
        return '4 - reported and weak'
    df['tier'] = df.apply(tier, axis=1)

    # 대조군 정의를 바꿨을 때 이 후보가 어디로 움직였는가. 지금까지 이 값은
    # comparator_reordering 쪽에만 있었고 후보표에는 없었다. 대장이 후보의 발견·선택·
    # 보강·반증 경로를 한 줄로 보여 주려면 이 열이 있어야 한다 — 순위가 대조군 정의에
    # 얼마나 의존하는지가 그 후보를 얼마나 믿을지에 직접 들어가기 때문이다.
    rs_path = 'results/comparator_reordering/candidate_rank_shift.tsv'
    if os.path.exists(rs_path):
        rs = pd.read_csv(rs_path, sep='\t').set_index('gene')
        df['rank_naive_control'] = df['gene'].map(rs['rank_naive'])
        df['rank_matched_control'] = df['gene'].map(rs['rank_matched'])
        df['rank_shift_on_substitution'] = -df['gene'].map(rs['moved'])
    else:
        for c in ('rank_naive_control', 'rank_matched_control',
                  'rank_shift_on_substitution'):
            df[c] = np.nan
    df = df.sort_values(['tier', 'mean_rank'])
    df.to_csv(args.out, sep='\t', index=False)

    log('permutation ceiling for |g| vs other-CKD: %.2f\n' % ceiling)
    for t, d in df.groupby('tier'):
        log('=== TIER %s  (n=%d) ===' % (t, len(d)))
        log('  %-9s %-7s %-9s %-9s %-6s %-5s %s'
            % ('gene', 'RBSrnk', 'g vs oth', '>ceiling', 'patt', 'DKDp', 'literature verdict'))
        for _, r in d.iterrows():
            log('  %-9s %-7.0f %+9.2f %-9s %-6s %-5s %s'
                % (r['gene'], r['mean_rank'], r['g_DKD_vs_otherCKD'],
                   'yes' if r['above_perm_ceiling'] else '',
                   (r['pattern_glom'] or '')[:5], int(r['n_dkd']) if pd.notna(r['n_dkd']) else '?',
                   r['verdict']))
        log('')

    log('=' * 90)
    t1 = df[df['tier'].str.startswith('1')]
    if len(t1):
        log('DELIVERABLE: %d candidate(s) are both novel in the literature and strong in the data:'
            % len(t1))
        for _, r in t1.iterrows():
            log('  %s' % r['gene'])
    else:
        log('DELIVERABLE: NO candidate is simultaneously (a) absent from the DKD literature and')
        log('(b) pattern A/B and not donor-driven. The genes that are strongest in the data are')
        log('already published DKD genes; the genes with no DKD literature are the ones whose')
        log('signal is largely living-donor-versus-biopsy. This must be reported as it stands.')
        log('')
        log('Closest to the requirement, and what to pursue:')
        best = df[df['lit_novel']].sort_values('mean_rank').head(5)
        for _, r in best.iterrows():
            log('  %-9s RBSrank=%-5.0f g vs otherCKD=%+.2f  pattern=%s  %s'
                % (r['gene'], r['mean_rank'], r['g_DKD_vs_otherCKD'],
                   (r['pattern_glom'] or '?'), r['verdict']))
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
