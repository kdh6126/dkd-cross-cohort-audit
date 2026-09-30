#!/usr/bin/env python
"""Attach the specificity gate's uncertainty to every table that quotes a pass-rate.

The gate keeps genes with |Hedges g| >= 0.5 between DKD and non-diabetic CKD in ERCB, from 12
and 17 DKD patients. Two numbers characterise it. They are both percentages and they are NOT the
same kind of ratio - the denominators differ - so they must never be set side by side as if one
could be compared with the other:

    BACKGROUND PASS RATE          denominator = all 9,900 genes in the frozen space.
                                  820 pass, so a random gene set passes at 8.3%. A selected set
                                  passing at 42% is 5.0x that, which is the number a reader
                                  needs; 42% on its own sounds like a coin flip.
    PERMUTATION-BASED FDR         denominator = the 820 observed passers.
      ESTIMATE                    under 500 patient-level label permutations 137 pass on
                                  average, so 16.7% of the observed passers are expected by
                                  chance. Patient-level because the two ERCB blocks profile the
                                  same patients in two compartments; permuting each block
                                  independently gives 13.7% and understates the gate.

Both are added as columns and as an enrichment ratio, so no table quotes a bare pass-rate.
"""
import os
import sys
import numpy as np
import pandas as pd

TABLES = [
    ('results/deconfound.tsv', 'spec_pass_frac', ['arm', 'K']),
    ('results/deconfound_v2.tsv', 'spec_pass_frac', ['arm', 'K']),
    ('results/ruv_benchmark/benchmark.tsv', 'spec_pass_frac', ['arm', 'K']),
    ('results/conventional_pipeline/comparison.tsv', 'spec_pass_frac', ['pipeline', 'K']),
]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    spec = pd.read_csv('results/final_candidates.tsv', sep='\t')
    spec = spec.rename(columns={spec.columns[0]: 'entrez_id'})
    n_pass = int(spec['passes_specificity'].sum())
    n_tot = len(spec)
    background = n_pass / n_tot

    perm = pd.read_csv('results/specificity_null.tsv', sep='\t')
    perm_mean = float(perm['perm_count'].mean())
    fdr = perm_mean / max(n_pass, 1)
    ceiling = float(np.percentile(perm['perm_max_abs_g'], 95))

    log('=== specificity gate characteristics ===')
    log('  genes passing            : %d of %d  -> background rate %.1f%%'
        % (n_pass, n_tot, 100 * background))
    log('  passing under permutation: %.0f  -> FDR estimate %.1f%%' % (perm_mean, 100 * fdr))
    log('  |g| reachable by chance  : 95th percentile %.2f' % ceiling)
    log('')

    for path, col, keys in TABLES:
        if not os.path.exists(path):
            log('  (missing) %s' % path)
            continue
        d = pd.read_csv(path, sep='\t')
        if col not in d.columns:
            log('  (no %s column) %s' % (col, path))
            continue
        # 이 단계는 표에 열을 덧붙인다. 예전 이름으로 붙여 둔 열이 남아 있으면 같은
        # 값이 두 이름으로 실려 나가므로, 이 단계가 만드는 열은 먼저 전부 지우고
        # 다시 붙인다. 그래야 다시 돌려도 결과가 같다.
        d = d.drop(columns=[c for c in ('spec_background_rate', 'spec_enrichment',
                                        'gate_fdr_estimate', 'gate_empirical_fdr')
                            if c in d.columns])
        d['spec_background_rate'] = background
        d['spec_enrichment'] = d[col] / background
        d['gate_fdr_estimate'] = fdr
        d.to_csv(path, sep='\t', index=False)
        agg = d.groupby(keys)[col].mean().reset_index()
        agg['enrichment_over_random'] = agg[col] / background
        log('=== %s ===' % path)
        log('  %-18s %-6s %14s %14s' % (keys[0], keys[1], 'pass rate', 'x background'))
        for _, r in agg.iterrows():
            log('  %-18s %-6s %13.0f%% %13.1fx'
                % (r[keys[0]], r[keys[1]], 100 * r[col], r['enrichment_over_random']))
        log('')

    pd.DataFrame([{
        'genes_total': n_tot,
        'genes_passing': n_pass,
        'background_rate': background,
        'mean_passing_under_permutation': perm_mean,
        'gate_fdr_estimate': fdr,
        'perm_abs_g_95th': ceiling,
    }]).to_csv('results/gate_characteristics.tsv', sep='	', index=False)
    log('wrote results/gate_characteristics.tsv')
    log('')
    log('reporting rule: wherever a pass-rate appears, quote the enrichment over the %.1f%%'
        % (100 * background))
    log('BACKGROUND PASS RATE (out of all 9,900 genes), and state the %.1f%% PERMUTATION-BASED'
        % (100 * fdr))
    log('FDR ESTIMATE (out of the observed passers) separately. They are not comparable rates.')


if __name__ == '__main__':
    main()
