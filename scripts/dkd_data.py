#!/usr/bin/env python
"""Data access for the DKD cross-cohort experiments.

Design notes
------------
* Cohorts live in data/processed/harmonized/ as gene x sample matrices on a single
  9,900-gene Entrez space (see docs/worklog/GENE_SPACE_AND_FEASIBILITY.md).
* Value scales differ fundamentally between cohorts (GSE30528/GSE30529 are mean-centred and
  go negative; the rest are positive log-scale), so every cross-cohort operation runs on
  **per-cohort gene-wise z-scores**.
* Z-scoring a cohort with its own statistics uses NO label information, so applying it to a
  held-out cohort is not leakage - it is the unsupervised domain-alignment step that makes
  LODO possible at all. Any *supervised* step (feature selection, model fitting) sees only
  training cohorts.
"""
import os
import numpy as np
import pandas as pd

H = 'data/processed/harmonized'

GLOM = ['GSE30528', 'GSE96804', 'GSE104948']          # primary LODO members, one compartment
CORTEX = ['GSE142025']                                 # cross-compartment generalisation test
TUBULO = ['GSE30529', 'GSE104954']                     # compartment reproducibility
ALL = GLOM + CORTEX + TUBULO

# subjects shared between a glomerular cohort and its tubulointerstitial counterpart, so the
# compartment analysis can be restricted to genuinely independent patients (see the audit).
SHARED_SUBJECTS = {
    'GSE30529': {'62', '67', '164', '168', '178', '76', '77', '81', '82'},
    'GSE104954': set(),   # matched on ERCB ID below instead
}


def load_cohort(name, labelled_only=True, root=H):
    """Return (X, y, pheno). X is samples x genes, gene-wise z-scored within this cohort.

    root lets the cross-tissue extension read its own harmonised directory with the same code.
    """
    e = pd.read_csv(os.path.join(root, '%s_expr.tsv' % name), sep='\t', index_col=0)
    p = pd.read_csv(os.path.join(root, '%s_pheno.tsv' % name), sep='\t', dtype={'sample': str})
    e.columns = e.columns.astype(str)
    if labelled_only:
        p = p[p['label'].isin([0, 1])].copy()
    e = e[p['sample'].values]

    X = e.T.astype(float)                              # samples x genes
    X = X.loc[:, ~X.isna().all(axis=0)]                # drop all-NaN genes
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1).replace(0, np.nan)
    Z = (X - mu) / sd
    Z = Z.fillna(0.0)                                  # constant genes carry no information
    return Z, p['label'].values.astype(int), p.reset_index(drop=True)


def load_many(names, labelled_only=True, root=H):
    """Load several cohorts onto the shared gene space. Returns dict name -> (X, y, pheno)."""
    out = {}
    genes = None
    for n in names:
        X, y, p = load_cohort(n, labelled_only, root)
        out[n] = (X, y, p)
        genes = X.columns if genes is None else genes.intersection(X.columns)
    for n in out:
        X, y, p = out[n]
        out[n] = (X[genes], y, p)
    return out, list(genes)


def lodo_folds(data, members=None):
    """Yield (test_name, X_train, y_train, cohort_train, X_test, y_test).

    Training data is the row-concatenation of every other cohort, each already z-scored
    within itself. cohort_train labels each training row with its source cohort, so a
    selector can weight or stratify by cohort if it wants to.
    """
    members = members or list(data)
    for held in members:
        tr = [m for m in members if m != held]
        Xtr = pd.concat([data[m][0] for m in tr], axis=0)
        ytr = np.concatenate([data[m][1] for m in tr])
        ctr = np.concatenate([[m] * len(data[m][1]) for m in tr])
        Xte, yte, _ = data[held]
        yield held, Xtr, ytr, ctr, Xte, yte


def summarise(data):
    rows = []
    for n, (X, y, p) in data.items():
        rows.append(dict(cohort=n, samples=X.shape[0], genes=X.shape[1],
                         case=int((y == 1).sum()), control=int((y == 0).sum()),
                         compartment=p['compartment'].iloc[0]))
    return pd.DataFrame(rows)


if __name__ == '__main__':
    data, genes = load_many(ALL)
    print('shared gene space: %d' % len(genes))
    print(summarise(data).to_string(index=False))
    print()
    for held, Xtr, ytr, ctr, Xte, yte in lodo_folds(data, GLOM):
        print('hold out %-11s train=%3d (%2d case) from %-28s test=%3d (%2d case)'
              % (held, len(ytr), ytr.sum(), '+'.join(sorted(set(ctr))), len(yte), yte.sum()))
