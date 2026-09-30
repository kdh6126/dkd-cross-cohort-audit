#!/usr/bin/env python
"""Independent validation in GSE175759 - and it validates a different thing than expected.

GSE175759 is 90 manually microdissected tubulointerstitial RNA-seq samples: 46 IgAN, 22
nephrectomy controls, 9 minimal change, 4 membranous, 3 FSGS, 3 lupus, and only **3 diabetic
nephropathy**. Three samples are flagged technical outliers and are dropped here.

Three DN cases cannot test a DKD-vs-control claim. What this cohort *can* do is better:

  1. REPLICATE THE PROCUREMENT ARTIFACT independently. The stated weakness of this project was
     that the confounder-matched contrast existed only in ERCB. GSE175759 has non-diabetic CKD
     biopsies and nephrectomy controls from a different centre, a different platform (RNA-seq,
     not array) and a different compartment handling. If the immediate-early module shifts in
     non-diabetic CKD versus controls here too, the artifact conclusion replicates.

  2. TEST THE CANDIDATES AGAINST eGFR. This cohort carries CKD-EPI eGFR per sample. Correlating
     a candidate with eGFR needs no DN label at all, and within the biopsy group it is
     procurement-independent - the single most useful thing available here. FMOD/LUM/MMP2 serve
     as positive controls: published DKD ECM genes should track function loss.

  3. DN vs other CKD for the candidates, reported as underpowered (n=3) rather than omitted.
"""
import os, re, sys, gzip, glob, argparse
import numpy as np
import pandas as pd
from scipy import stats

SAMPLES = 'data/raw/geo/GSE175759/samples'
MATRIX = 'data/raw/geo/GSE175759/GSE175759_series_matrix.txt.gz'
IEG = ['FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'ATF3', 'DUSP1', 'ZFP36', 'NR4A2', 'BTG2',
       'EGR2', 'EGR3', 'IER2', 'KLF2', 'KLF4', 'SOCS3', 'DUSP2', 'JUND', 'NR4A1']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    n1, n0 = len(a), len(b)
    if n1 < 2 or n0 < 2:
        return np.full(a.shape[1], np.nan)
    sp = np.sqrt(((n1 - 1) * a.var(0, ddof=1) + (n0 - 1) * b.var(0, ddof=1)) / (n1 + n0 - 2))
    sp[sp == 0] = np.nan
    J = 1 - 3 / (4 * (n1 + n0) - 9)
    return (a.mean(0) - b.mean(0)) / sp * J


def ensembl_to_entrez():
    """From NCBI gene_info dbXrefs, which already lists Ensembl gene IDs - no extra download."""
    m = {}
    with gzip.open('data/raw/annotation/Homo_sapiens.gene_info.gz', 'rt',
                   encoding='utf-8', errors='replace') as f:
        hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip('\n').split('\t')
            for x in p[ix['dbXrefs']].split('|'):
                if x.startswith('Ensembl:'):
                    m[x.split(':', 1)[1]] = p[ix['GeneID']]
    return m


def load_pheno():
    hdr = {}
    with gzip.open(MATRIX, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!'):
                p = line.rstrip('\n').split('\t')
                hdr.setdefault(p[0][1:], []).append([x.strip('"') for x in p[1:]])
    gsm = hdr['Sample_geo_accession'][0]
    title = hdr['Sample_title'][0]
    chars = hdr.get('Sample_characteristics_ch1', [])
    rows = []
    for i, g in enumerate(gsm):
        rec = {'gsm': g, 'title': title[i]}
        for c in chars:
            v = c[i] if i < len(c) else ''
            if v.startswith('diagnosis:'):
                rec['diagnosis'] = v.split(':', 1)[1].strip()
            elif 'gfr' in v.lower():
                try:
                    rec['egfr'] = float(v.split(':', 1)[1])
                except ValueError:
                    rec['egfr'] = np.nan
            elif v.startswith('technical outlier:'):
                rec['outlier'] = 'Outlier' in v
        rows.append(rec)
    p = pd.DataFrame(rows)
    p['sample_col'] = p['title'].str.replace(r'^[A-Za-z]+_', '', regex=True)
    return p


def build_matrix(e2g):
    """Per-sample count files -> CPM, log2, collapsed to Entrez by max-mean."""
    files = sorted(glob.glob(os.path.join(SAMPLES, '*.txt.gz')))
    cols = {}
    for f in files:
        gsm = os.path.basename(f).split('_')[0]
        d = pd.read_csv(f, sep='\t', index_col=0)
        cols[gsm] = d.iloc[:, 0]
    cnt = pd.DataFrame(cols)
    cnt.index = [i.split('.')[0] for i in cnt.index]           # strip Ensembl version
    cnt = cnt.groupby(level=0).sum()
    log('  counts: %d ensembl genes x %d samples' % cnt.shape)

    cpm = cnt / cnt.sum(axis=0) * 1e6
    lg = np.log2(cpm + 1)

    ent = pd.Series([e2g.get(i) for i in lg.index], index=lg.index)
    lg = lg[ent.notna().values]
    ent = ent[ent.notna()]
    lg['__entrez'] = ent.values
    means = lg.drop(columns='__entrez').mean(axis=1)
    order = np.lexsort((-means.values, lg['__entrez'].values))
    lg = lg.iloc[order]
    lg = lg[~lg['__entrez'].duplicated(keep='first')]           # max-mean collapse
    out = lg.drop(columns='__entrez')
    out.index = lg['__entrez'].values
    log('  mapped + collapsed: %d entrez genes' % len(out))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='results/gse175759_validation.tsv')
    ap.add_argument('--master', default='results/master_candidate_table.tsv')
    args = ap.parse_args()

    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym = dict(zip(gs['entrez_id'], gs['symbol']))
    sym2id = {s: g for g, s in sym.items() if isinstance(s, str)}
    space = set(gs['entrez_id'])

    log('loading GSE175759 ...')
    X = build_matrix(ensembl_to_entrez())
    p = load_pheno().set_index('gsm')
    X = X[[c for c in X.columns if c in p.index]]
    p = p.loc[X.columns]

    keep = ~p['outlier'].fillna(False)
    X, p = X.loc[:, keep.values], p[keep]
    inter = [g for g in X.index if g in space]
    X = X.loc[inter]
    log('  after outlier removal and gene-space intersection: %d genes x %d samples'
        % (X.shape[0], X.shape[1]))

    # z-score within this cohort, as everywhere else in the project
    Z = ((X.T - X.T.mean()) / X.T.std(ddof=1).replace(0, np.nan)).fillna(0.0)   # samples x genes

    p['group'] = np.where(p['diagnosis'].str.lower().str.contains('control'), 'control',
                          np.where(p['diagnosis'].str.lower().str.contains('diabetic'), 'DN', 'otherCKD'))
    log('\n  groups: %s' % dict(p['group'].value_counts()))
    log('  eGFR available for %d/%d samples' % (p['egfr'].notna().sum(), len(p)))

    ctrl = Z[(p['group'] == 'control').values].values
    oth = Z[(p['group'] == 'otherCKD').values].values
    dn = Z[(p['group'] == 'DN').values].values

    # ---------------------------------------------------------------- 1. artifact replication
    log('\n=== 1. procurement artifact - independent replication ===')
    g_oth = pd.Series(hedges_g(oth, ctrl), index=Z.columns)
    ieg_ids = [sym2id[s] for s in IEG if s in sym2id and sym2id[s] in Z.columns]
    ieg_g = g_oth[ieg_ids]
    rest = g_oth.drop(index=ieg_ids)
    log('  non-diabetic CKD (n=%d) vs nephrectomy control (n=%d), tubulointerstitium, RNA-seq'
        % (len(oth), len(ctrl)))
    log('  immediate-early module mean g = %+.2f   (all other genes %+.2f)'
        % (ieg_g.mean(), rest.mean()))
    t = stats.mannwhitneyu(ieg_g.dropna(), rest.dropna(), alternative='two-sided')
    log('  module vs rest of transcriptome: Mann-Whitney p = %.2e' % t.pvalue)
    for s in ['DUSP1', 'ZFP36', 'FOS', 'EGR1', 'JUN', 'FOSB', 'ATF3', 'NR4A2', 'BTG2']:
        gid = sym2id.get(s)
        if gid in g_oth.index:
            log('     %-7s g = %+.2f' % (s, g_oth[gid]))
    verdict = ('REPLICATED - the module moves in non-diabetic CKD in an independent cohort, '
               'platform and centre' if abs(ieg_g.mean()) > 0.4 and
               abs(ieg_g.mean()) > abs(rest.mean()) + 0.2 else 'NOT replicated')
    log('  -> %s' % verdict)

    # ---------------------------------------------------------------- 2. eGFR association
    log('\n=== 2. candidates vs eGFR (no DN label needed) ===')
    mt = pd.read_csv(args.master, sep='\t')
    tier1 = mt[mt['tier'].str.startswith('1')]['gene'].tolist()
    tier2 = mt[mt['tier'].str.startswith('2')]['gene'].tolist()

    biopsy = (p['group'] != 'control').values          # procurement-matched subset
    rows = []
    for scope, mask in (('biopsies only', biopsy), ('all samples', np.ones(len(p), bool))):
        eg = p['egfr'].values[mask]
        ok = ~np.isnan(eg)
        sub = Z.values[mask][ok]
        eg = eg[ok]
        log('\n  -- %s (n=%d with eGFR) --' % (scope, ok.sum()))
        # null: distribution of |rho| over all genes, for context
        allrho = np.array([stats.spearmanr(sub[:, j], eg).statistic for j in range(sub.shape[1])])
        thr = np.nanpercentile(np.abs(allrho), 95)
        log('     genome-wide |rho| 95th percentile = %.3f  (context for the numbers below)' % thr)
        for tier, genes in (('T1', tier1), ('T2', tier2)):
            for s in genes:
                gid = sym2id.get(s)
                if gid not in Z.columns:
                    continue
                j = list(Z.columns).index(gid)
                r = stats.spearmanr(sub[:, j], eg)
                flag = '*' if abs(r.statistic) > thr else ' '
                log('     %s %-8s rho = %+.3f  p = %.3g %s' % (tier, s, r.statistic, r.pvalue, flag))
                rows.append(dict(scope=scope, tier=tier, gene=s, rho=r.statistic,
                                 p=r.pvalue, above_95pct=abs(r.statistic) > thr,
                                 n=int(ok.sum())))

    # ---------------------------------------------------------------- 3. DN vs other CKD
    log('\n=== 3. DN vs other CKD in this cohort (n=%d vs %d - underpowered) ===' % (len(dn), len(oth)))
    g_dn = pd.Series(hedges_g(dn, oth), index=Z.columns)
    for tier, genes in (('T1', tier1), ('T2', tier2)):
        for s in genes:
            gid = sym2id.get(s)
            if gid in g_dn.index:
                log('  %s %-8s g = %+.2f' % (tier, s, g_dn[gid]))
    log('  With 3 DN samples these are point estimates with no useful precision; reported for')
    log('  direction only, and they should not be used to promote or retire a candidate.')

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep='\t', index=False)
    pd.DataFrame({'entrez_id': g_oth.index, 'symbol': [sym.get(i, '') for i in g_oth.index],
                  'g_otherCKD_vs_control': g_oth.values,
                  'g_DN_vs_otherCKD': g_dn.reindex(g_oth.index).values}).to_csv(
        args.out.replace('.tsv', '_effects.tsv'), sep='\t', index=False)
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
