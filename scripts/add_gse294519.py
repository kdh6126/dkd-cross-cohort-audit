#!/usr/bin/env python
"""Add GSE294519 as a fifth cohort, without disturbing the frozen gene space.

23 DKD vs 13 normal (paraneoplastic tissue), renal tubules, NovaSeq. The natural temptation is
to re-derive the common gene space now that a fifth cohort exists. That must not happen: every
bootstrap selection frequency already computed is defined on the 9,900-gene space, and changing
the space silently invalidates all of them. So the new cohort is mapped INTO the frozen space
and its coverage is reported instead.

Same conventions as the other cohorts: HGNC symbol -> Entrez via NCBI gene_info, max-mean
collapse, log2 of the supplied normalised values, per-cohort z-scoring applied downstream.
"""
import os, sys, gzip, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

EXPR = 'data/raw/geo/GSE294519/GSE294519_Normalization_data_allsamples.txt.gz'
H = 'data/processed/harmonized'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def sym_maps():
    sym, syn = {}, {}
    with gzip.open('data/raw/annotation/Homo_sapiens.gene_info.gz', 'rt',
                   encoding='utf-8', errors='replace') as f:
        hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip('\n').split('\t')
            sym[p[ix['Symbol']]] = p[ix['GeneID']]
            for a in p[ix['Synonyms']].split('|'):
                if a and a != '-' and a not in syn:
                    syn[a] = p[ix['GeneID']]
    return sym, syn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--space', default='data/processed/gene_space_all6.tsv')
    args = ap.parse_args()

    gs = pd.read_csv(args.space, sep='\t', dtype=str)
    space = list(gs['entrez_id'])
    space_set = set(space)

    d = pd.read_csv(EXPR, sep='\t', index_col=0)
    log('raw: %d genes x %d samples' % d.shape)
    lg = np.log2(d.clip(lower=0) + 1)

    sym, syn = sym_maps()
    ent = pd.Series([sym.get(i) or syn.get(i) for i in lg.index], index=lg.index)
    lg = lg[ent.notna().values]
    ent = ent[ent.notna()]
    lg = lg.assign(__e=ent.values)
    means = lg.drop(columns='__e').mean(axis=1)
    order = np.lexsort((-means.values, lg['__e'].values))
    lg = lg.iloc[order]
    lg = lg[~lg['__e'].duplicated()]                       # max-mean collapse
    X = lg.drop(columns='__e')
    X.index = lg['__e'].values
    log('  mapped and collapsed: %d Entrez genes' % len(X))

    present = [g for g in space if g in set(X.index)]
    cov = len(present) / len(space)
    log('  coverage of the frozen %d-gene space: %d (%.1f%%)' % (len(space), len(present), 100 * cov))
    if cov < 0.85:
        log('  ! coverage below 85%% - this cohort would distort a joint analysis')

    # reindex onto the frozen space; genes absent here are left as NaN and dropped by callers
    E = X.reindex(space)
    E.index.name = 'entrez_id'

    labels = [1 if c.upper().startswith('DKD') else 0 for c in E.columns]
    pheno = pd.DataFrame({
        'cohort': 'GSE294519',
        'sample': list(E.columns),
        'title': list(E.columns),
        'subject_id': list(E.columns),
        'group': ['DKD' if l else 'control' for l in labels],
        'label': labels,
        'stage': '',
        'control_type': ['' if l else 'tumor_nephrectomy' for l in labels],
        'platform': 'GPL24676',
        'compartment': 'tubulointerstitium',
    })
    E.to_csv(os.path.join(H, 'GSE294519_expr.tsv'), sep='\t')
    pheno.to_csv(os.path.join(H, 'GSE294519_pheno.tsv'), sep='\t', index=False)
    log('  wrote %s/GSE294519_{expr,pheno}.tsv  (%d DKD / %d control)'
        % (H, sum(labels), len(labels) - sum(labels)))

    # sanity: does it recover the expected biology?
    Z = E.loc[present].T
    Z = ((Z - Z.mean()) / Z.std(ddof=1).replace(0, np.nan)).fillna(0.0)
    y = np.array(labels)
    a, b = Z.values[y == 1], Z.values[y == 0]
    sp = np.sqrt(((len(a) - 1) * a.var(0, ddof=1) + (len(b) - 1) * b.var(0, ddof=1))
                 / (len(a) + len(b) - 2))
    sp[sp == 0] = np.nan
    g = pd.Series((a.mean(0) - b.mean(0)) / sp, index=present)
    sym_by_id = dict(zip(gs['entrez_id'], gs['symbol']))
    log('\n  sanity check - strongest genes in this cohort:')
    for gid in g.abs().sort_values(ascending=False).head(12).index:
        log('    %-10s g = %+.2f' % (sym_by_id.get(gid, gid), g[gid]))
    known = ['FMOD', 'LUM', 'MMP2', 'COL1A2', 'THBS2', 'NPHS1', 'NPHS2', 'DUSP1', 'ZFP36']
    s2i = {v: k for k, v in sym_by_id.items() if isinstance(v, str)}
    log('\n  reference genes:')
    for s in known:
        gid = s2i.get(s)
        if gid in g.index:
            log('    %-10s g = %+.2f' % (s, g[gid]))


if __name__ == '__main__':
    main()
