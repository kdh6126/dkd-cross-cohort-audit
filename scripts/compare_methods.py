#!/usr/bin/env python
"""Merge baseline and proposed results into the paper's comparison table.

Adds the metric the baselines' own runner cannot produce on its own: CROSS-FOLD stability -
how much of the top-K signature survives changing which cohort is held out. That is the
direct operationalisation of "is this gene dataset-specific?", and it is computed identically
for every method from the per-fold selection-frequency tables.
"""
import os, sys, glob, argparse
import numpy as np
import pandas as pd

TOPK = [5, 10, 20, 50, 100]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def cross_fold_from_freq(path, cohorts):
    """freq_<method>.tsv has one selection-frequency column per held-out fold."""
    df = pd.read_csv(path, sep='\t', index_col=0)
    cols = [c for c in cohorts if c in df.columns]
    if len(cols) < 2:
        return {}
    out = {}
    for K in TOPK:
        sets = [set(df[c].sort_values(ascending=False).head(K).index) for c in cols]
        js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
        out[K] = float(np.mean(js))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--baselines', default='results/baselines')
    ap.add_argument('--proposed', nargs='*', default=['results/proposed'])
    ap.add_argument('--cohorts', default='GSE30528,GSE96804,GSE104948')
    ap.add_argument('--out', default='results/comparison.tsv')
    ap.add_argument('--null', default='results/null_control_draws.tsv')
    args = ap.parse_args()
    cohorts = args.cohorts.split(',')

    frames = []
    b = os.path.join(args.baselines, 'results.tsv')
    if os.path.exists(b):
        frames.append(pd.read_csv(b, sep='\t'))
    for d in args.proposed:
        p = os.path.join(d, 'results.tsv')
        if os.path.exists(p):
            frames.append(pd.read_csv(p, sep='\t'))
    if not frames:
        sys.exit('no results found')
    res = pd.concat(frames, ignore_index=True)

    # ---- normalise every AUROC against the fold's random-signature null.
    # This is not optional here: random 50-gene sets already reach AUROC 0.71-0.87 on these
    # folds, because the DKD signature is so broadly co-expressed. A raw AUROC of 0.94 can sit
    # at the 95th percentile of chance. Percentile-vs-null is the only honest way to read it.
    if os.path.exists(args.null):
        draws = pd.read_csv(args.null, sep='	')
        def pct(r):
            col = '%s|%d' % (r['held_out'], r['K'])
            if col not in draws.columns or np.isnan(r['external_auroc']):
                return np.nan
            d = draws[col].dropna().values
            return float((d < r['external_auroc']).mean())
        def zed(r):
            col = '%s|%d' % (r['held_out'], r['K'])
            if col not in draws.columns or np.isnan(r['external_auroc']):
                return np.nan
            d = draws[col].dropna().values
            sd = d.std(ddof=1)
            return float((r['external_auroc'] - d.mean()) / sd) if sd > 0 else np.nan
        res['null_percentile'] = res.apply(pct, axis=1)
        res['null_z'] = res.apply(zed, axis=1)
    else:
        res['null_percentile'] = np.nan
        res['null_z'] = np.nan

    agg = (res.groupby(['method', 'K'])
              .agg(external=('external_auroc', 'mean'),
                   external_sd=('external_auroc', 'std'),
                   internal=('nested_cv_auroc', 'mean'),
                   boot_jaccard=('jaccard', 'mean'),
                   kuncheva=('kuncheva', 'mean'),
                   null_pct=('null_percentile', 'mean'),
                   null_z=('null_z', 'mean'))
              .reset_index())

    # ---- cross-fold stability
    xf = {}
    for f in glob.glob(os.path.join(args.baselines, 'freq_*.tsv')):
        m = os.path.basename(f)[5:-4]
        d = cross_fold_from_freq(f, cohorts)
        for K, v in d.items():
            xf[(m, K)] = v
    for d in args.proposed:
        p = os.path.join(d, 'cross_fold_stability.tsv')
        if os.path.exists(p):
            for _, r in pd.read_csv(p, sep='\t').iterrows():
                xf[(r['method'], int(r['K']))] = r['cross_fold_jaccard']
    agg['cross_fold_jaccard'] = [xf.get((m, K), np.nan) for m, K in zip(agg['method'], agg['K'])]
    agg['gap'] = agg['internal'] - agg['external']

    agg.to_csv(args.out, sep='\t', index=False)

    for K in TOPK:
        sub = agg[agg['K'] == K].sort_values('cross_fold_jaccard', ascending=False)
        if sub.empty:
            continue
        log('\n=== K = %d ===' % K)
        log('  %-12s %-9s %-11s %-9s %-11s %-11s %s'
            % ('method', 'external', 'vs null', 'internal', 'boot Jacc', 'cross-fold', 'gap'))
        for _, r in sub.iterrows():
            nz = '-' if np.isnan(r['null_z']) else '%+.2f sd (p%02d)' % (r['null_z'], round(100 * r['null_pct']))
            log('  %-12s %.3f     %-11s %.3f     %.3f       %-11s %+.3f'
                % (r['method'], r['external'], nz, r['internal'], r['boot_jaccard'],
                   '-' if np.isnan(r['cross_fold_jaccard']) else '%.3f' % r['cross_fold_jaccard'],
                   -r['gap']))
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
