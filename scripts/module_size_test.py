#!/usr/bin/env python
"""Why did WGCNA avoid the immediate-early module when stability selection did not?

Two candidate explanations, both testable:

  (A) MODULE BOUNDARY ACCIDENT - the IEG genes never entered the disease-correlated module, so
      connectivity was never computed for them. Nothing methodological, just where the cut fell.

  (B) SIZE BIAS - the IEG genes did enter the module, but intramodular connectivity is a SUM of
      adjacencies over the whole module, so it rewards membership of the LARGEST coherent block.
      The ECM/immune programme is hundreds of genes; the IEG module is ~19. Stability selection
      has no such size term: a small block that is perfectly reproducible scores as high as a
      large one.

(B), if true, is a general principle and explains the earlier simulation failure: that
simulation gave the confounder and disease blocks the SAME size (30 vs 30), which is exactly
the case where a size-biased criterion has no advantage.

Part 1 measures which explanation holds in the real data.
Part 2 re-runs the simulation with unequal block sizes to see whether (B) reproduces it.
"""
import os, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM
from dkd_deconfound import IEG
from simulation import one_run

OUT = 'results/module_size'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def real_data_part(power=6, module_frac=0.10, K=50):
    cohorts = GLOM + ['GSE142025']
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym = dict(zip(gs['entrez_id'], gs['symbol']))
    sym2id = {s: g for g, s in sym.items() if isinstance(s, str)}
    ieg_ids = [sym2id[s] for s in IEG if s in sym2id]
    ieg_pos = np.array([int(np.where(genes == g)[0][0]) for g in ieg_ids if (genes == g).any()])

    X = np.vstack([data[c][0].values for c in cohorts])
    y = np.concatenate([data[c][1] for c in cohorts]).astype(float)
    Xz = (X - X.mean(0)) / (X.std(0, ddof=1) + 1e-12)
    n = len(y)
    yz = (y - y.mean()) / (y.std(ddof=1) + 1e-12)
    gsig = (Xz.T @ yz) / (n - 1)

    m = max(50, int(module_frac * X.shape[1]))
    mod = np.argsort(-np.abs(gsig))[:m]
    in_mod = np.isin(ieg_pos, mod)
    log('=== part 1: the real data ===')
    log('  disease-correlated module: %d genes of %d' % (len(mod), X.shape[1]))
    log('  immediate-early genes inside the module: %d / %d' % (in_mod.sum(), len(ieg_pos)))
    log('  |gene-trait correlation|: IEG median %.3f, module threshold %.3f, all-gene median %.3f'
        % (np.median(np.abs(gsig[ieg_pos])), np.abs(gsig[mod]).min(),
           np.median(np.abs(gsig))))

    S = Xz[:, mod]
    R = (S.T @ S) / (n - 1)
    A = np.abs(R) ** power
    np.fill_diagonal(A, 0.0)
    kIN = A.sum(1)
    rank_in_mod = np.argsort(np.argsort(-kIN)) + 1
    pos_in_mod = {g: i for i, g in enumerate(mod)}

    log('\n  connectivity rank inside the module (1 = most connected of %d):' % len(mod))
    rows = []
    for g, p_ in zip([sym.get(genes[i], genes[i]) for i in ieg_pos], ieg_pos):
        if p_ in pos_in_mod:
            i = pos_in_mod[p_]
            log('    %-8s |gs|=%.3f  kIN=%.1f  rank %d / %d'
                % (g, abs(gsig[p_]), kIN[i], rank_in_mod[i], len(mod)))
            rows.append(dict(gene=g, in_module=True, gs=abs(gsig[p_]), kIN=kIN[i],
                             rank=int(rank_in_mod[i])))
        else:
            log('    %-8s |gs|=%.3f  NOT in module' % (g, abs(gsig[p_])))
            rows.append(dict(gene=g, in_module=False, gs=abs(gsig[p_])))

    # how big is the block each gene belongs to?
    log('\n  block size around each gene (genes in the module with |r| > 0.5 to it):')
    Rabs = np.abs(R)
    blk = (Rabs > 0.5).sum(1)
    ieg_in = [pos_in_mod[p_] for p_ in ieg_pos if p_ in pos_in_mod]
    top_conn = np.argsort(-kIN)[:50]
    log('    immediate-early genes : median %d' % int(np.median(blk[ieg_in])) if ieg_in else '    none')
    log('    top-50 by connectivity: median %d' % int(np.median(blk[top_conn])))
    log('    whole module          : median %d' % int(np.median(blk)))

    verdict = ('B - size bias: the IEG genes ARE in the module but sit in a small block'
               if in_mod.sum() >= len(ieg_pos) * 0.6
               else 'A - boundary accident: most IEG genes never entered the module')
    log('\n  verdict: %s' % verdict)
    return pd.DataFrame(rows), verdict


def sim_part(reps=8, K=50):
    log('\n=== part 2: simulation with UNEQUAL block sizes ===')
    log('  the earlier sweep used 30 disease genes and 30 confounder genes, which is exactly')
    log('  the case where a size-biased criterion cannot separate them.')
    rng = np.random.default_rng(20260825)
    rows = []
    log('\n  %-8s %-8s | %s' % ('n_dis', 'n_conf',
                                ' '.join('%-26s' % m for m in
                                         ('effect_size', 'stability', 'connectivity'))))
    for n_dis, n_conf in ((30, 30), (100, 30), (200, 20), (300, 20), (20, 200)):
        acc = {m: [0.0, 0.0] for m in ('effect_size', 'stability', 'connectivity')}
        for r in range(reps):
            o = one_run(rng, 4, 0.95, 0.4, 0.05, 2000, 40, n_dis, n_conf,
                        1.0, 1.0, 0.35, K)
            for m, v in o.items():
                acc[m][0] += v['confounder_frac']
                acc[m][1] += v['disease_precision']
        log('  %-8d %-8d | ' % (n_dis, n_conf) + ' '.join(
            'conf %.2f / dis %.2f       ' % (acc[m][0] / reps, acc[m][1] / reps)
            for m in ('effect_size', 'stability', 'connectivity')))
        for m in acc:
            rows.append(dict(n_dis=n_dis, n_conf=n_conf, method=m,
                             confounder_frac=acc[m][0] / reps,
                             disease_precision=acc[m][1] / reps))
    log('\n  if size bias is the mechanism, connectivity should lose interest in the confounder')
    log('  as the disease block grows, while stability should not.')
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reps', type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    real, verdict = real_data_part()
    real.to_csv(os.path.join(OUT, 'real_module_membership.tsv'), sep='\t', index=False)
    sim = sim_part(reps=args.reps)
    sim.to_csv(os.path.join(OUT, 'block_size_sweep.tsv'), sep='\t', index=False)
    log('\nwrote %s/' % OUT)


if __name__ == '__main__':
    main()
