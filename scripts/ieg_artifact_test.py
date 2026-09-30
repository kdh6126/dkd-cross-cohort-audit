#!/usr/bin/env python
"""Is the immediate-early gene module DKD biology, or a tissue-procurement artifact?

The proposed selector ranks FOS, JUN, JUNB, FOSB, EGR1, ATF3, DUSP1, ZFP36, NR4A2, BTG2 near
the top: reproducible across every cohort, stable under perturbation, and enriched for
"sequence-specific DNA binding". They are also the genes the DKD-specificity test rejects.

These are canonical immediate-early genes - they respond within minutes to ischaemia and
tissue handling. In every cohort here, cases are needle BIOPSIES while controls are tumour
NEPHRECTOMY or living-donor tissue, so warm-ischaemia time is confounded with case status
in exactly the same way in every cohort. A procurement effect would therefore look perfectly
"cross-cohort reproducible" - which is precisely what a robustness score rewards.

The decisive test: ERCB also contains NON-diabetic CKD biopsies, procured like the DKD
biopsies. If the module is DKD biology, other-CKD should sit with the controls. If it is a
procurement effect, other-CKD should sit with the DKD cases.
"""
import os, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_cohort

IEG = ['FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'ATF3', 'DUSP1', 'ZFP36', 'NR4A2', 'BTG2',
       'EGR2', 'EGR3', 'IER2', 'KLF2', 'KLF4', 'SOCS3', 'CYR61', 'DUSP2', 'JUND', 'NR4A1']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    n1, n0 = len(a), len(b)
    if n1 < 2 or n0 < 2:
        return np.full(a.shape[1], np.nan)
    sp = np.sqrt(((n1 - 1) * a.var(0, ddof=1) + (n0 - 1) * b.var(0, ddof=1)) / (n1 + n0 - 2))
    sp[sp == 0] = np.nan
    return (a.mean(0) - b.mean(0)) / sp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='results/ieg_artifact.tsv')
    args = ap.parse_args()

    sym2id = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = sym2id.dropna(subset=['symbol']).set_index('symbol')['entrez_id'].to_dict()
    ids = {s: sym2id[s] for s in IEG if s in sym2id}
    log('immediate-early module: %d of %d genes present in the 9,900-gene space\n' % (len(ids), len(IEG)))

    rows = []
    for c in ['GSE104948', 'GSE104954']:
        X, y, p = load_cohort(c, labelled_only=False)
        X.columns = X.columns.astype(str)
        p = p.reset_index(drop=True)
        cols = [v for v in ids.values() if v in X.columns]
        M = X[cols].values
        dkd = M[(p['label'] == 1).values]
        ctl = M[(p['label'] == 0).values]
        oth = M[(p['label'] == -1).values]

        g_dkd = hedges_g(dkd, ctl)
        g_oth = hedges_g(oth, ctl)
        log('=== %s  (DKD=%d, control=%d, other-CKD=%d) ===' % (c, len(dkd), len(ctl), len(oth)))
        log('  %-8s %-16s %-18s' % ('gene', 'DKD vs control', 'other-CKD vs control'))
        inv = {v: k for k, v in ids.items()}
        for j, col in enumerate(cols):
            log('  %-8s %+8.2f         %+8.2f' % (inv[col], g_dkd[j], g_oth[j]))
            rows.append(dict(cohort=c, gene=inv[col], dkd_vs_control=g_dkd[j],
                             otherckd_vs_control=g_oth[j]))
        log('  %-8s %+8.2f         %+8.2f   <- module mean'
            % ('MEAN', np.nanmean(g_dkd), np.nanmean(g_oth)))
        r = np.corrcoef(g_dkd[~np.isnan(g_dkd) & ~np.isnan(g_oth)],
                        g_oth[~np.isnan(g_dkd) & ~np.isnan(g_oth)])[0, 1]
        log('  correlation between the two effect vectors: r = %.3f\n' % r)

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)

    m_dkd = df['dkd_vs_control'].mean()
    m_oth = df['otherckd_vs_control'].mean()
    ratio = m_oth / m_dkd if m_dkd else np.nan
    log('=' * 78)
    log('module mean, DKD vs control      : %+.2f' % m_dkd)
    log('module mean, other-CKD vs control: %+.2f  (%.0f%% of the DKD shift)' % (m_oth, 100 * ratio))
    log('')
    if ratio > 0.6:
        log('VERDICT: non-diabetic CKD biopsies show the SAME shift as DKD biopsies. The module')
        log('tracks whether tissue arrived as a biopsy or as nephrectomy/donor tissue, not')
        log('whether the patient had diabetes. It is a PROCUREMENT ARTIFACT and must not be')
        log('reported as a DKD biomarker module, however reproducible it looks.')
    elif ratio < 0.3:
        log('VERDICT: the shift is specific to DKD - other CKD sits with the controls. The')
        log('module is disease biology after all.')
    else:
        log('VERDICT: partial overlap - the module is shared across kidney disease but not')
        log('purely procedural. Report as non-specific rather than artifactual.')
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
