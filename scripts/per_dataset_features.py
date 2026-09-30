#!/usr/bin/env python
"""Per-dataset feature importance catalogue - the first rung of the ladder.

The cross-cohort work needs several cohorts and therefore silently discards every dataset too
small or too odd to join a LODO rotation. That throws away information: what a single dataset
considers important is worth recording even when it cannot be pooled, because
  (a) it is a reference other people can use,
  (b) it makes dataset-specific features explicit rather than leaving them as an abstraction
      that the cross-cohort method is claimed to remove, and
  (c) the single-dataset -> cross-cohort -> cross-modality progression is only measurable if
      the single-dataset rung actually exists.

Every dataset with a usable two-group contrast is run on its own here: the four LODO cohorts,
the two compartment-validation cohorts, and five that the cross-cohort analysis could not use
(GSE1009, GSE111154, GSE20602, GSE142153 blood, GSE175759 non-diabetic-CKD-vs-control).

Each gets bootstrap stability selection with four selectors, and the output is a per-gene
selection frequency per dataset plus a cross-dataset agreement matrix.
"""
import os, re, sys, gzip, glob, json, argparse
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_selectors import SELECTORS
from dkd_stability import mean_pairwise_jaccard

ANN = 'data/raw/annotation'
H = 'data/processed/harmonized'
OUT = 'results/per_dataset'
SELECTORS_USED = ['univariate', 'lasso', 'relieff', 'rf']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- annotation
def annot_map(gpl):
    """probe -> entrez, dropping multi-gene probes.

    GPL17586 (HTA-2.0) has no .annot file on the GEO FTP - only the full platform table, whose
    gene assignment is a '//'-delimited blob. Handled separately.
    """
    out = {}
    if gpl == 'GPL17586':
        path = os.path.join(ANN, 'GPL17586.platform.txt')
        with open(path, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('ID\t'):
                    ix = {c: i for i, c in enumerate(line.rstrip('\n').split('\t'))}
                    break
            else:
                return out
            col = ix['gene_assignment']
            for line in f:
                if line.startswith('!'):
                    break
                p = line.rstrip('\n').split('\t')
                if len(p) <= col or p[col] in ('', '---'):
                    continue
                ids = set()
                for rec in p[col].split('///'):
                    fld = [x.strip() for x in rec.split('//')]
                    if len(fld) >= 5 and fld[4].isdigit():
                        ids.add(fld[4])
                if len(ids) == 1:
                    out[p[0]] = ids.pop()
        return out
    with gzip.open(os.path.join(ANN, '%s.annot.gz' % gpl), 'rt',
                   encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('ID\t'):
                ix = {c: i for i, c in enumerate(line.rstrip('\n').split('\t'))}
                break
        else:
            return out
        col = ix['Gene ID']
        for line in f:
            if line.startswith('!'):
                break
            p = line.rstrip('\n').split('\t')
            if len(p) <= col:
                continue
            ids = {g for g in p[col].split('///') if g.strip().isdigit()}
            if len(ids) == 1:
                out[p[0]] = ids.pop()
    return out


def sym_maps():
    sym, syn = {}, {}
    with gzip.open(os.path.join(ANN, 'Homo_sapiens.gene_info.gz'), 'rt',
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


def read_series_matrix(path):
    hdr, rows, cols = {}, [], None
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!'):
                p = line.rstrip('\n').split('\t')
                hdr.setdefault(p[0][1:], []).append([x.strip('"') for x in p[1:]])
        cols = [x.strip('"') for x in f.readline().rstrip('\n').split('\t')][1:]
        for line in f:
            if line.startswith('!series_matrix_table_end'):
                break
            p = line.rstrip('\n').split('\t')
            rows.append((p[0].strip('"'), p[1:]))
    data = np.array([[np.nan if v in ('', 'null', 'NA') else float(v) for v in r[1]]
                     for r in rows], dtype=float)
    return pd.DataFrame(data, index=[r[0] for r in rows], columns=cols), hdr


def collapse(df, pmap, space):
    keep = [p for p in df.index if pmap.get(p) in space]
    sub = df.loc[keep]
    genes = np.array([pmap[p] for p in keep])
    means = sub.mean(axis=1, skipna=True).values
    order = np.lexsort((-means, genes))
    seen, sel = set(), []
    for i in order:
        if genes[i] not in seen:
            seen.add(genes[i]); sel.append(i)
    out = sub.iloc[sel]
    out.index = genes[sel]
    return out


# ---------------------------------------------------------------- dataset loaders
def from_harmonised(name, space):
    e = pd.read_csv(os.path.join(H, '%s_expr.tsv' % name), sep='\t', index_col=0)
    p = pd.read_csv(os.path.join(H, '%s_pheno.tsv' % name), sep='\t', dtype={'sample': str})
    p = p[p['label'].isin([0, 1])]
    e = e[p['sample'].values]
    e.index = e.index.astype(str)
    return e, p['label'].values.astype(int), 'DKD vs control'


def from_series(acc, gpl, space, label_fn, contrast):
    df, hdr = read_series_matrix('data/raw/geo/%s/%s_series_matrix.txt.gz' % (acc, acc))
    pmap = annot_map(gpl)
    e = collapse(df, pmap, space)
    src = hdr.get('Sample_source_name_ch1', [[]])[0]
    ttl = hdr.get('Sample_title', [[]])[0]
    chs = hdr.get('Sample_characteristics_ch1', [])
    y = []
    for i in range(df.shape[1]):
        ch = ' | '.join(c[i] for c in chs if i < len(c))
        y.append(label_fn(ttl[i] if i < len(ttl) else '', src[i] if i < len(src) else '', ch))
    y = np.array(y)
    ok = y >= 0
    return e.loc[:, e.columns[ok]], y[ok].astype(int), contrast


def load_gse142153(space, sym, syn):
    """Blood PBMC. Use the submitter's gene-symbol matrix; columns encode the group."""
    d = pd.read_csv('data/raw/geo/GSE142153/GSE142153_normalized_data_with_genename.txt.gz',
                    sep='\t', index_col=0)
    d = d.drop(columns=[c for c in d.columns if c.lower() == 'probename'])
    ent = pd.Series([sym.get(i) or syn.get(i) for i in d.index], index=d.index)
    d = d[ent.notna().values]; ent = ent[ent.notna()]
    d = d[[c for c in d.columns]]
    d['__e'] = ent.values
    means = d.drop(columns='__e').mean(axis=1)
    order = np.lexsort((-means.values, d['__e'].values))
    d = d.iloc[order]
    d = d[~d['__e'].duplicated()]
    e = d.drop(columns='__e'); e.index = d['__e'].values
    e = e.loc[[g for g in e.index if g in space]]
    y = []
    for c in e.columns:
        cu = c.upper()
        if cu.startswith('HC'):
            y.append(0)
        elif cu.startswith('DM_NO') or cu.startswith('DM_MI') or cu.startswith('DM_MA'):
            y.append(1)          # diabetic with nephropathy (normo/micro/macroalbuminuria)
        else:
            y.append(-1)         # DM_ESRD - end stage, excluded from the two-group contrast
    y = np.array(y)
    ok = y >= 0
    return e.loc[:, e.columns[ok]], y[ok].astype(int), 'DN vs healthy control (PBMC)'


def load_gse175759(space):
    """Non-diabetic CKD vs nephrectomy control - the contrast this cohort can actually power."""
    from validate_gse175759 import ensembl_to_entrez, load_pheno, build_matrix
    X = build_matrix(ensembl_to_entrez())
    p = load_pheno().set_index('gsm')
    X = X[[c for c in X.columns if c in p.index]]
    p = p.loc[X.columns]
    keep = ~p['outlier'].fillna(False)
    X, p = X.loc[:, keep.values], p[keep]
    X = X.loc[[g for g in X.index if g in space]]
    y = np.where(p['diagnosis'].str.lower().str.contains('control'), 0, 1)
    return X, y.astype(int), 'any kidney disease vs nephrectomy control'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--B', type=int, default=200)
    ap.add_argument('--k', type=int, default=50)
    ap.add_argument('--seed', type=int, default=20260824)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    space = set(gs['entrez_id'])
    symname = dict(zip(gs['entrez_id'], gs['symbol']))
    sym, syn = sym_maps()

    datasets = {}
    for nm in ['GSE142025', 'GSE30528', 'GSE96804', 'GSE104948', 'GSE30529', 'GSE104954']:
        datasets[nm] = from_harmonised(nm, space)
    # datasets the cross-cohort analysis could not use
    datasets['GSE1009'] = from_series(
        'GSE1009', 'GPL8300', space,
        lambda t, s, c: 1 if 'diabet' in (t + s + c).lower() else 0, 'DN vs control')
    datasets['GSE111154'] = from_series(
        'GSE111154', 'GPL17586', space,
        lambda t, s, c: 1 if 'diabet' in (t + s + c).lower() and 'non-diabetic' not in
        (t + s + c).lower() else 0, 'early DN vs control')
    datasets['GSE20602'] = from_series(
        'GSE20602', 'GPL96', space,
        lambda t, s, c: 0 if 'nephrectomy' in (t + s + c).lower() else 1,
        'nephrosclerosis vs control')
    datasets['GSE142153'] = load_gse142153(space, sym, syn)
    datasets['GSE175759'] = load_gse175759(space)

    log('\n%-11s %6s %6s %5s %5s  %s' % ('dataset', 'genes', 'samp', 'case', 'ctrl', 'contrast'))
    for nm, (e, y, contrast) in datasets.items():
        log('%-11s %6d %6d %5d %5d  %s'
            % (nm, e.shape[0], e.shape[1], int((y == 1).sum()), int((y == 0).sum()), contrast))

    # ---------------------------------------------------------------- per-dataset selection
    freqs, tops = {}, {}
    for nm, (e, y, contrast) in datasets.items():
        if min((y == 1).sum(), (y == 0).sum()) < 3:
            log('\n  skip %s - fewer than 3 in a group' % nm)
            continue
        X = e.T.astype(float)
        Z = ((X - X.mean()) / X.std(ddof=1).replace(0, np.nan)).fillna(0.0).values
        genes = list(e.index)
        rng = np.random.default_rng(args.seed)
        acc = np.zeros((len(SELECTORS_USED), len(genes)))
        sets_all = []
        for si, sname in enumerate(SELECTORS_USED):
            fn = SELECTORS[sname]
            cnt = np.zeros(len(genes)); nb = 0; sets = []
            for b in range(args.B):
                idx = np.concatenate([rng.choice(np.flatnonzero(y == lab),
                                                 (y == lab).sum(), replace=True)
                                      for lab in (0, 1)])
                if len(np.unique(y[idx])) < 2:
                    continue
                sel, _ = fn(Z[idx], y[idx], k=args.k, rng=rng)
                sel = np.asarray(sel, int)
                cnt[sel] += 1; nb += 1; sets.append(sel)
            acc[si] = cnt / max(nb, 1)
            sets_all.append(sets)
        cons = acc.mean(0)                       # consensus across the four selectors
        freqs[nm] = pd.Series(cons, index=genes)
        order = np.argsort(-cons)
        tops[nm] = [genes[i] for i in order[:100]]
        wj = np.mean([mean_pairwise_jaccard(s) for s in sets_all if len(s) > 1])
        df = pd.DataFrame({'entrez_id': genes,
                           'symbol': [symname.get(g, '') for g in genes],
                           'consensus_freq': cons})
        for si, sname in enumerate(SELECTORS_USED):
            df[sname] = acc[si]
        df.sort_values('consensus_freq', ascending=False).to_csv(
            os.path.join(OUT, 'features_%s.tsv' % nm), sep='\t', index=False)
        log('\n  %-11s within-dataset bootstrap Jaccard = %.3f' % (nm, wj))
        log('     top 12: %s' % ', '.join(symname.get(g, g) for g in tops[nm][:12]))

    # ---------------------------------------------------------------- cross-dataset agreement
    names = list(freqs)
    log('\n=== cross-dataset agreement of the top-50 feature sets (Jaccard) ===')
    log('%-11s %s' % ('', ' '.join('%9s' % n[3:] for n in names)))
    M = pd.DataFrame(index=names, columns=names, dtype=float)
    for a in names:
        cells = []
        for b in names:
            if a == b:
                M.loc[a, b] = np.nan; cells.append('%9s' % '-'); continue
            sa, sb = set(tops[a][:50]), set(tops[b][:50])
            j = len(sa & sb) / len(sa | sb)
            M.loc[a, b] = j
            cells.append('%9.3f' % j)
        log('%-11s %s' % (a, ' '.join(cells)))
    M.to_csv(os.path.join(OUT, 'cross_dataset_jaccard.tsv'), sep='\t')

    F = pd.DataFrame(freqs)
    F.insert(0, 'symbol', [symname.get(g, '') for g in F.index])
    F['n_datasets_top50'] = [sum(1 for n in names if g in set(tops[n][:50])) for g in F.index]
    F.sort_values('n_datasets_top50', ascending=False).to_csv(
        os.path.join(OUT, 'consensus_across_datasets.tsv'), sep='\t')

    log('\n=== genes in the top-50 of the most datasets ===')
    top = F.sort_values(['n_datasets_top50'], ascending=False).head(15)
    for g, r in top.iterrows():
        log('  %-9s in %d/%d datasets' % (r['symbol'] or g, r['n_datasets_top50'], len(names)))
    log('\n  mean pairwise agreement across all %d datasets = %.3f'
        % (len(names), np.nanmean(M.values.astype(float))))
    log('\nwrote %s/' % OUT)


if __name__ == '__main__':
    main()
