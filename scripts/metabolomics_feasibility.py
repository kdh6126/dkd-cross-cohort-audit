#!/usr/bin/env python
"""Can a second modality carry the confounder-matched contrast that the transcriptomics could not?

The manuscript's Limitations say no genuine multi-omics contrast was constructible: KPMP
proteomics has no diabetic participants, and the DKD metabolomics depositions we found were
either raw spectra or had no comparison group. That conclusion came from searching the
Metabolomics Workbench for "diabetic nephropathy" and "diabetic kidney", which misses studies
that carry DKD as one arm of a broader design and are titled after the broader disease.

ST003255 (Seoul National University) is one of those: 286 human plasma samples, 187 named and
quantified metabolites, with DKD, three other kidney diseases and healthy controls in the same
deposition. It supplies exactly the design Section 3.3 argues for, in a modality where the
procurement confounder cannot arise, because a blood draw is a blood draw whether the patient
is a case or a control.

This script tests whether it is usable:

    naive contrast        DKD vs healthy. The design every public DKD study uses.
    matched contrast      DKD vs other kidney disease. What the paper argues is necessary.
    specificity           how much of the naive signal survives the matched contrast.

If the second number is a small fraction of the first, the metabolomics reproduces the
transcriptomic finding in a modality where handling cannot explain it.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

SRC = 'data/raw/metabolomics/ST003255/AN005337_datatable.tsv'
OUT = 'results/metabolomics'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    d = (a.mean() - b.mean()) / (sp + 1e-12)
    return d * (1 - 3 / (4 * (na + nb) - 9))


def main():
    os.makedirs(OUT, exist_ok=True)
    d = pd.read_csv(SRC, sep='\t')
    d = d.rename(columns={d.columns[0]: 'sample', d.columns[1]: 'group'})
    d['group'] = d['group'].str.replace('Disease:', '', regex=False)
    mets = [c for c in d.columns if c not in ('sample', 'group')]

    X = d[mets].apply(pd.to_numeric, errors='coerce')
    keep = X.notna().mean() >= 0.8
    X = X.loc[:, keep]
    mets = list(X.columns)
    X = X.fillna(X.median())
    # metabolomics intensities are strongly right-skewed; log then z-score per metabolite
    X = np.log2(X.clip(lower=X[X > 0].min().min() / 2))
    Z = (X - X.mean()) / (X.std(ddof=1) + 1e-12)

    g = d['group'].values
    log('=== ST003255 ===')
    for k, v in pd.Series(g).value_counts().items():
        log('  %-10s %3d' % (k, v))
    log('  metabolites quantified in >=80%% of samples: %d of %d' % (len(mets), len(keep)))
    log('')

    dkd = g == 'DKD'
    healthy = g == 'Normal'
    other = np.isin(g, ['IgAN', 'MN', 'HN'])

    rows = []
    for name, mask_a, mask_b in (
            ('DKD vs healthy', dkd, healthy),
            ('DKD vs other kidney disease', dkd, other),
            ('other kidney disease vs healthy', other, healthy)):
        gs, ps = [], []
        for m in mets:
            a, b = Z.loc[mask_a, m].values, Z.loc[mask_b, m].values
            gs.append(hedges_g(a, b))
            ps.append(stats.ttest_ind(a, b, equal_var=False).pvalue)
        gs, ps = np.array(gs), np.array(ps)
        # A constant metabolite yields a NaN p, and np.minimum.accumulate propagates that NaN
        # across the whole array, silently zeroing the significant count. Treat it as p = 1.
        ps = np.where(np.isfinite(ps), ps, 1.0)
        order = np.argsort(ps)
        m_ = len(ps)
        # Benjamini-Hochberg
        q = ps[order] * m_ / (np.arange(m_) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        qq = np.empty(m_)
        qq[order] = np.clip(q, 0, 1)
        rows.append(dict(contrast=name, n_a=int(mask_a.sum()), n_b=int(mask_b.sum()),
                         n_sig_q05=int((qq < 0.05).sum()),
                         n_large_g=int((np.abs(gs) >= 0.5).sum()),
                         max_abs_g=float(np.abs(gs).max())))
        pd.DataFrame({'metabolite': mets, 'hedges_g': gs, 'p': ps, 'q': qq}).to_csv(
            os.path.join(OUT, 'ST003255_%s.tsv' % name.replace(' ', '_')), sep='\t', index=False)

    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'ST003255_contrasts.tsv'), sep='\t', index=False)
    log('  %-32s %5s %5s %10s %10s %8s' % ('contrast', 'n_a', 'n_b', 'q<0.05', '|g|>=0.5', 'max|g|'))
    for _, r in t.iterrows():
        log('  %-32s %5d %5d %10d %10d %8.2f'
            % (r['contrast'], r['n_a'], r['n_b'], r['n_sig_q05'], r['n_large_g'], r['max_abs_g']))

    naive = t.loc[0, 'n_sig_q05']
    matched = t.loc[1, 'n_sig_q05']
    log('')
    log('=== what this means ===')
    log('  naive contrast   : %d of %d metabolites at q<0.05' % (naive, len(mets)))
    log('  matched contrast : %d of %d' % (matched, len(mets)))
    if naive:
        log('  survival rate    : %.0f%% of the naive signal is DKD-specific'
            % (100 * matched / naive))
    log('')
    log('  A blood draw is a blood draw, so the biopsy-versus-nephrectomy confounder of')
    log('  Section 3.3 cannot arise here. Whatever the matched contrast loses is disease')
    log('  non-specificity, not tissue handling.')


if __name__ == '__main__':
    main()
