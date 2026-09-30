#!/usr/bin/env python
"""Per-diagnosis expression profile of the final candidates.

The DKD-vs-otherCKD gate collapses nine diagnoses into one "other" group, which hides whether
a gene is uniquely DN or merely highest in DN. It also hides whether the control contrast is
being driven by living-donor tissue, which is procured differently from every biopsy group.

Printing the mean z-score per diagnosis separates three patterns that the gate cannot:
  A  uniquely elevated in DN, every other disease flat        -> disease-specific
  B  elevated across kidney disease, but far more in DN       -> quantitative, not qualitative
  C  the contrast is mostly living-donor vs everything else   -> procurement-adjacent, treat with care
"""
import os, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_cohort

ERCB = ['GSE104948', 'GSE104954']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', default='results/final_candidates_gated.tsv')
    ap.add_argument('--topn', type=int, default=30)
    ap.add_argument('--min-n', type=int, default=10)
    ap.add_argument('--out', default='results/per_diagnosis_profile.tsv')
    args = ap.parse_args()

    sym2id = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = sym2id.dropna(subset=['symbol']).set_index('symbol')['entrez_id'].to_dict()
    cand = pd.read_csv(args.candidates, sep='\t').head(args.topn)
    genes = [g for g in cand['symbol'].dropna().astype(str) if g in sym2id]

    rows = []
    for c in ERCB:
        X, y, p = load_cohort(c, labelled_only=False)
        X.columns = X.columns.astype(str)
        p = p.reset_index(drop=True)
        grp = p['group'].values
        keep = [g for g in pd.unique(grp) if (grp == g).sum() >= args.min_n]
        keep = ([g for g in keep if g == 'DN'] + [g for g in keep if g == 'living_donor'] +
                sorted(g for g in keep if g not in ('DN', 'living_donor')))
        log('\n=== %s : mean z-score by diagnosis ===' % c)
        log('%-9s' % 'gene' + ''.join('%-10s' % g[:9] for g in keep) + '  pattern')
        for s in genes:
            gid = sym2id[s]
            if gid not in X.columns:
                continue
            v = X[gid].values
            means = {g: float(v[grp == g].mean()) for g in keep}
            dn = means.get('DN', np.nan)
            ld = means.get('living_donor', np.nan)
            others = [means[g] for g in keep if g not in ('DN', 'living_donor')]
            oth_max = max(np.abs(others)) if others else np.nan
            # classify
            if abs(dn) >= 3 * oth_max:
                pat = 'A unique-to-DN'
            elif abs(dn) >= 1.8 * oth_max:
                pat = 'B DN-dominant'
            else:
                pat = 'C shared'
            if not np.isnan(ld) and abs(ld) > abs(dn) and abs(ld) > 1.5 * oth_max:
                pat += ' /donor-driven'
            log('%-9s' % s + ''.join('%-10.2f' % means[g] for g in keep) + '  ' + pat)
            rows.append(dict(cohort=c, gene=s, pattern=pat, dn=dn, living_donor=ld,
                             max_other_disease=oth_max, **{('grp_' + g): means[g] for g in keep}))
        log('%-9s' % 'n:' + ''.join('%-10d' % (grp == g).sum() for g in keep))

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)

    log('\n=== consensus pattern across both compartments ===')
    for s in genes:
        pats = df[df['gene'] == s]['pattern'].tolist()
        if len(pats) == 2 and pats[0].split()[0] == pats[1].split()[0]:
            log('  %-9s %s' % (s, pats[0]))
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
