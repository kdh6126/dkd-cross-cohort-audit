#!/usr/bin/env python
"""Germline layer: is there GWAS signal near the candidate genes?

Data: GCST90179152, diabetic kidney disease in a Korean cohort (2,532 cases / 31,347 controls,
PMID 36627639). The file carries p_value, chromosome, position and rsID only - no effect sizes -
so this can be a positional gene-level test and nothing more.

Two things this test gets wrong if done naively, both handled:

  1. **min-p is biased by window size.** A gene with 400 SNPs in its window has a smaller
     minimum p than a gene with 40, with no biology involved. The null here is built from
     random genes **matched on SNP count** (nearest-neighbour matching in log SNP count), not
     from random genes generally.
  2. **LD makes SNPs non-independent**, so a window's SNPs are not 400 independent tests.
     min-p with a matched empirical null absorbs this to first order, because the null windows
     have the same LD structure on average. It is not MAGMA and is not presented as such.

Expected outcome, stated before looking: 2,532 cases is small for a GWAS. Nothing should reach
genome-wide significance, and the honest question is whether candidates rank above matched
random genes at all.
"""
import os, sys, gzip, argparse
import numpy as np
import pandas as pd
from scipy import stats

GWAS = 'data/raw/gwas/GCST90179152_KoreanDKD_GRCh37.tsv'
REFGENE = 'data/raw/annotation/refGene_hg19.txt.gz'
OUT = 'results/gwas_layer.tsv'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load_gene_coords():
    """refGene hg19 -> one canonical span per symbol (widest transcript on a primary chrom)."""
    cols = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
            'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat',
            'cdsEndStat', 'exonFrames']
    d = pd.read_csv(REFGENE, sep='\t', names=cols, compression='gzip', low_memory=False)
    d = d[d['chrom'].str.match(r'^chr(\d+|X|Y)$', na=False)]
    d['chr'] = d['chrom'].str.replace('chr', '', regex=False)
    g = (d.groupby(['name2', 'chr'])
           .agg(start=('txStart', 'min'), end=('txEnd', 'max')).reset_index())
    g['span'] = g['end'] - g['start']
    g = g.sort_values('span', ascending=False).drop_duplicates('name2')
    return g.set_index('name2')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--window', type=int, default=100000, help='bp added each side of the gene')
    ap.add_argument('--n-null', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=20260824)
    ap.add_argument('--gwas', default=GWAS)
    ap.add_argument('--master', default='results/master_candidate_table.tsv')
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args()

    log('loading gene coordinates ...')
    coords = load_gene_coords()
    log('  %d gene symbols with hg19 spans' % len(coords))

    log('loading GWAS summary statistics ...')
    op = gzip.open if args.gwas.endswith('.gz') else open
    gw = pd.read_csv(args.gwas, sep='\t', usecols=['p_value', 'chromosome', 'base_pair_location'],
                     dtype={'chromosome': str}, low_memory=False)
    gw = gw.dropna()
    gw['chromosome'] = gw['chromosome'].astype(str)
    log('  %d variants, min p = %.3g' % (len(gw), gw['p_value'].min()))
    log('  genome-wide significant (p<5e-8): %d' % int((gw['p_value'] < 5e-8).sum()))

    bychr = {c: d.sort_values('base_pair_location') for c, d in gw.groupby('chromosome')}
    pos = {c: d['base_pair_location'].values for c, d in bychr.items()}
    pval = {c: d['p_value'].values for c, d in bychr.items()}

    def window_stat(sym):
        if sym not in coords.index:
            return None
        r = coords.loc[sym]
        c = str(r['chr'])
        if c not in pos:
            return None
        lo, hi = r['start'] - args.window, r['end'] + args.window
        i0, i1 = np.searchsorted(pos[c], [lo, hi])
        if i1 <= i0:
            return None
        p = pval[c][i0:i1]
        return dict(gene=sym, chr=c, n_snps=int(i1 - i0), min_p=float(p.min()),
                    n_p05=int((p < 0.05).sum()), frac_p05=float((p < 0.05).mean()))

    # ---------------------------------------------------------------- candidates
    mt = pd.read_csv(args.master, sep='\t')
    cands = [g for g in mt['gene'].dropna().astype(str)]
    tier = dict(zip(mt['gene'], mt['tier']))

    log('\ncomputing window statistics ...')
    cand_stats = [s for s in (window_stat(g) for g in cands) if s]
    log('  %d of %d candidates have coordinates and variants in window' % (len(cand_stats), len(cands)))

    # ---------------------------------------------------------------- SNP-count-matched null
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    universe = [s for s in gs['symbol'].dropna().astype(str)
                if s in coords.index and s not in set(cands)]
    rng = np.random.default_rng(args.seed)
    sample = list(rng.choice(universe, size=min(args.n_null, len(universe)), replace=False))
    log('  building null from %d background genes ...' % len(sample))
    null = [s for s in (window_stat(g) for g in sample) if s]
    N = pd.DataFrame(null)
    log('  %d usable background genes' % len(N))

    logn = np.log10(N['n_snps'].values)
    rows = []
    for s in cand_stats:
        # match on SNP count: background genes within +-15% log10 SNP count
        target = np.log10(s['n_snps'])
        sel = N[np.abs(logn - target) < 0.06]
        if len(sel) < 30:
            k = np.argsort(np.abs(logn - target))[:120]
            sel = N.iloc[k]
        pct = float((sel['min_p'].values < s['min_p']).mean())     # lower = stronger signal
        rows.append(dict(gene=s['gene'], tier=tier.get(s['gene'], ''), chr=s['chr'],
                         n_snps=s['n_snps'], min_p=s['min_p'], frac_p05=s['frac_p05'],
                         null_n=len(sel), null_median_min_p=float(np.median(sel['min_p'])),
                         empirical_p=pct))
    res = pd.DataFrame(rows).sort_values('empirical_p')
    res.to_csv(args.out, sep='\t', index=False)

    log('\n=== candidate windows vs SNP-count-matched background ===')
    log('  %-9s %-22s %6s %10s %12s %s'
        % ('gene', 'tier', 'nSNP', 'min p', 'null median', 'empirical p'))
    for _, r in res.iterrows():
        flag = ' *' if r['empirical_p'] < 0.05 else ''
        log('  %-9s %-22s %6d %10.3g %12.3g %.3f%s'
            % (r['gene'], str(r['tier'])[:22], r['n_snps'], r['min_p'],
               r['null_median_min_p'], r['empirical_p'], flag))

    log('\n=== summary ===')
    log('  candidates with empirical p < 0.05 : %d / %d (expected by chance: %.1f)'
        % (int((res['empirical_p'] < 0.05).sum()), len(res), 0.05 * len(res)))
    log('  candidates reaching genome-wide significance (min p < 5e-8): %d'
        % int((res['min_p'] < 5e-8).sum()))
    log('  binomial test on the 0.05 tail: p = %.3f'
        % stats.binomtest(int((res['empirical_p'] < 0.05).sum()), len(res), 0.05,
                          alternative='greater').pvalue)
    ks = stats.kstest(res['empirical_p'].values, 'uniform')
    log('  KS test of empirical p against uniform: D = %.3f, p = %.3f' % (ks.statistic, ks.pvalue))
    log('\nwrote %s' % OUT)


if __name__ == '__main__':
    main()
