#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Does the connectivity-comparator finding survive the reference WGCNA implementation?

The manuscript's comparator (conventional_pipeline.wgcna_hub) is a reimplementation and says
so. Two claims rest on it and they are separate:

  (2) IEG exclusion   the immediate-early module does not reach the comparator's top 50,
                      which the text explains by block size, not robustness.
  (1) verdict reversal the RBS-minus-comparator reproducibility difference changes sign
                      across cohort subsets.

This script re-runs both with the WGCNA R package itself (scripts/wgcna_reference.R: tutorial
defaults, soft power by scale-free fit, TOM + dynamic tree cut, module most correlated with
the trait, hubs by |kME|). Nothing about the reimplementation is assumed.

    python scripts/wgcna_reference.py --headline   # (2) on the four headline cohorts, LODO
    python scripts/wgcna_reference.py --sweep      # (1) on all 26 subsets, joined with the
                                                   #     stored RBS values of the same folds
    python scripts/wgcna_reference.py --smoke      # one fold, prints the R log

RBS is not recomputed: the stored sweep (results/cohort_sensitivity/subsets.tsv) is
deterministic under its seed, so its per-subset RBS rows are the same folds. What this script
adds is the reference comparator's top 50 on each fold: its cross-fold Jaccard, its external
AUROC with the same L2 logistic model, and how many immediate-early genes it holds.
"""
import argparse
import hashlib
import itertools
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM                      # noqa: E402
from cohort_sensitivity import ALL5                       # noqa: E402
from bootstrap_ci import auc_of                           # noqa: E402
from conventional_pipeline import wgcna_hub               # noqa: E402
from dkd_deconfound import IEG                            # noqa: E402

RSCRIPT = os.environ.get('RSCRIPT', r'C:\Program Files\R\R-4.6.1\bin\Rscript.exe')
OUT = 'results/wgcna_reference'
CACHE = os.path.join(OUT, 'cache')
K = 50
HEADLINE = GLOM + ['GSE142025']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def reference_scores(Xtr, ytr, genes, tag):
    """Per-gene |kME| inside the trait module, 0 elsewhere; cached per training set."""
    os.makedirs(CACHE, exist_ok=True)
    # 캐시 열쇠는 학습 코호트 집합만으로. headline 과 sweep 이 같은 fold 를 두 번 돌리지 않게.
    trainset = tag.split(':')[-1]
    key = hashlib.sha1((trainset + ':' + ','.join(map(str, genes))).encode()).hexdigest()[:16]
    path = os.path.join(CACHE, key + '.tsv')
    if not os.path.exists(path):
        with tempfile.TemporaryDirectory() as td:
            e, yp = os.path.join(td, 'expr.tsv'), os.path.join(td, 'y.tsv')
            pd.DataFrame(Xtr, columns=[str(g) for g in genes]).to_csv(e, sep='\t', index=False)
            np.savetxt(yp, ytr, fmt='%d')
            r = subprocess.run([RSCRIPT, 'scripts/wgcna_reference.R', e, yp, path],
                               capture_output=True, text=True, encoding='utf-8', errors='replace')
            if r.returncode != 0 or not os.path.exists(path):
                raise RuntimeError('Rscript failed for %s:\n%s' % (tag, r.stderr[-1500:]))
            log('    R: ' + r.stdout.strip().splitlines()[-1])
    t = pd.read_csv(path, sep='\t', dtype={'gene': str})
    score = pd.Series(0.0, index=[str(g) for g in genes])
    inmod = t[t['in_trait_module'] == 1]
    score.loc[inmod['gene']] = inmod['kME'].abs().values
    meta = t.iloc[0][['power', 'sft_r2', 'trait_module_cor', 'trait_module_size', 'n_modules']].to_dict()
    return score.values, meta


def fold_eval(data, genes, cohorts, ieg_idx, tag):
    """One LODO pass on `cohorts` for the reference and the custom comparator."""
    orders = {'WGCNA_ref': {}, 'WGCNA_hub': {}}
    aucs = {m: [] for m in orders}
    iegs = {m: [] for m in orders}
    metas = []
    for held in cohorts:
        tr = [c for c in cohorts if c != held]
        Xtr = np.vstack([data[c][0].values for c in tr])
        ytr = np.concatenate([data[c][1] for c in tr]).astype(int)
        Xte, yte = data[held][0].values, data[held][1]
        s_ref, meta = reference_scores(Xtr, ytr, genes, tag + ':' + '+'.join(tr))
        meta.update(held=held); metas.append(meta)
        s_hub = wgcna_hub({c: (data[c][0], data[c][1]) for c in tr}, tr)[0]
        for m, s in (('WGCNA_ref', s_ref), ('WGCNA_hub', s_hub)):
            o = np.argsort(-np.nan_to_num(s, nan=-np.inf))
            orders[m][held] = o
            aucs[m].append(auc_of(Xtr, ytr, Xte, yte, o[:K]))
            iegs[m].append(int(len(set(o[:K]) & set(ieg_idx))))
    res = {}
    for m in orders:
        sets = [set(orders[m][h][:K]) for h in cohorts]
        js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
        res[m] = dict(cross_fold_jaccard=float(np.mean(js)), external_auroc=float(np.nanmean(aucs[m])),
                      ieg_top50=float(np.mean(iegs[m])), ieg_per_fold=';'.join(map(str, iegs[m])))
    return res, metas


def ieg_index(genes):
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = {s: g for g, s in zip(gs['entrez_id'], gs['symbol']) if isinstance(s, str)}
    gpos = {str(g): i for i, g in enumerate(genes)}
    return [gpos[sym2id[s]] for s in IEG if s in sym2id and sym2id[s] in gpos]


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    ap = argparse.ArgumentParser()
    ap.add_argument('--headline', action='store_true')
    ap.add_argument('--sweep', action='store_true')
    ap.add_argument('--smoke', action='store_true')
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    data, genes = load_many(ALL5)
    ieg_idx = ieg_index(genes)
    log('%d genes, %d immediate-early genes on the space' % (len(genes), len(ieg_idx)))

    if args.smoke:
        tr = HEADLINE[:3]
        Xtr = np.vstack([data[c][0].values for c in tr])
        ytr = np.concatenate([data[c][1] for c in tr]).astype(int)
        s, meta = reference_scores(Xtr, ytr, genes, 'smoke:' + '+'.join(tr))
        top = np.argsort(-s)[:K]
        log('  trait module %d genes, top-50 IEG count %d, meta %s'
            % (int(meta['trait_module_size']), len(set(top) & set(ieg_idx)), meta))
        return

    if args.headline:
        res, metas = fold_eval(data, genes, HEADLINE, ieg_idx, 'headline')
        pd.DataFrame([dict(method=m, **v) for m, v in res.items()]).to_csv(
            os.path.join(OUT, 'headline.tsv'), sep='\t', index=False)
        pd.DataFrame(metas).to_csv(os.path.join(OUT, 'headline_folds.tsv'), sep='\t', index=False)
        for m, v in res.items():
            log('  %-10s Jaccard %.3f  AUROC %.3f  IEG in top 50: %.1f (%s)'
                % (m, v['cross_fold_jaccard'], v['external_auroc'], v['ieg_top50'], v['ieg_per_fold']))

    if args.sweep:
        stored = pd.read_csv('results/cohort_sensitivity/subsets.tsv', sep='\t')
        rows, metas = [], []
        subsets = [list(c) for r in range(2, 6) for c in itertools.combinations(ALL5, r)]
        for i, sub in enumerate(subsets, 1):
            res, mt = fold_eval(data, genes, sub, ieg_idx, 'sweep')
            name = '+'.join(s[3:] for s in sub)
            for m, v in res.items():
                rows.append(dict(subset=name, size=len(sub), method=m, **v))
            rbs = stored[(stored['subset'] == name) & (stored['method'] == 'RBS')]
            if len(rbs):
                rows.append(dict(subset=name, size=len(sub), method='RBS',
                                 cross_fold_jaccard=float(rbs['cross_fold_jaccard'].iloc[0]),
                                 external_auroc=float(rbs['external_auroc'].iloc[0]),
                                 ieg_top50=np.nan, ieg_per_fold=''))
            for x in mt:
                x.update(subset=name)
            metas += mt
            log('  %d / %d subsets' % (i, len(subsets)))
        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(OUT, 'subsets.tsv'), sep='\t', index=False)
        pd.DataFrame(metas).to_csv(os.path.join(OUT, 'sweep_folds.tsv'), sep='\t', index=False)
        p = df.pivot_table(index=['subset', 'size'], columns='method', values='cross_fold_jaccard')
        a = df.pivot_table(index=['subset', 'size'], columns='method', values='external_auroc')
        for comp in ('WGCNA_ref', 'WGCNA_hub'):
            dj = (p['RBS'] - p[comp]).xs(3, level='size')
            da = (a['RBS'] - a[comp]).xs(3, level='size')
            j5 = (p['RBS'] - p[comp]).xs(5, level='size').iloc[0]
            a5 = (a['RBS'] - a[comp]).xs(5, level='size').iloc[0]
            log('  RBS - %-9s three-cohort Jaccard %+.3f to %+.3f (%d/%d negative), AUROC %+.3f to %+.3f (%d/%d); five: %+.3f / %+.3f'
                % (comp, dj.min(), dj.max(), (dj < 0).sum(), len(dj), da.min(), da.max(), (da < 0).sum(), len(da), j5, a5))
    log('\nwrote %s/' % OUT)


if __name__ == '__main__':
    main()
