#!/usr/bin/env python
"""Downstream validation of a candidate signature.

Three questions, all answerable from data already on disk:

  1. COMPARTMENT   does a signature found in glomeruli reproduce in tubulointerstitium?
                   Reported twice - once per cohort pair - because the two pairs disagree
                   sharply (rho 0.014 in Woroniecka vs 0.511 in ERCB), and once restricted
                   to patients NOT shared between the paired cohorts, since GSE30528/GSE30529
                   and GSE104948/GSE104954 overlap by 5 DKD subjects each.

  2. SPECIFICITY   is the gene DKD-specific, or a generic CKD/fibrosis marker? ERCB carries
                   310 non-diabetic CKD biopsies (SLE, IgA, FSGS, MCD, MGN, RPGN, hypertensive)
                   that make this testable at no extra data cost. A gene that separates DKD
                   from healthy but NOT from other CKD is a kidney-damage marker, not a DKD
                   biomarker, and should be labelled as such rather than quietly promoted.

  3. CONSENSUS     how much does the signature change depending on which cohort was held out?
"""
import os, sys, glob, argparse
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_cohort

PAIRS = [('GSE30528', 'GSE30529', 'Woroniecka'), ('GSE104948', 'GSE104954', 'ERCB')]
OTHER_CKD_COHORTS = ['GSE104948', 'GSE104954']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(X, y):
    x1, x0 = X[y == 1], X[y == 0]
    n1, n0 = len(x1), len(x0)
    if n1 < 2 or n0 < 2:
        return np.full(X.shape[1], np.nan)
    v1, v0 = x1.var(0, ddof=1), x0.var(0, ddof=1)
    sp = np.sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2))
    sp[sp == 0] = np.nan
    d = (x1.mean(0) - x0.mean(0)) / sp
    J = 1 - 3 / (4 * (n1 + n0) - 9)
    return d * J


def consensus_signature(components_dir, topk):
    """Mean RBS rank across LODO folds -> one signature, plus per-fold agreement."""
    files = sorted(glob.glob(os.path.join(components_dir, 'components_holdout_*.tsv')))
    if not files:
        sys.exit('no components_holdout_*.tsv in %s' % components_dir)
    ranks, scores = {}, {}
    for f in files:
        held = os.path.basename(f)[len('components_holdout_'):-len('.tsv')]
        d = pd.read_csv(f, sep='\t', dtype={'entrez_id': str}).set_index('entrez_id')
        ranks[held] = d['RBS'].rank(ascending=False)
        scores[held] = d['RBS']
    R = pd.DataFrame(ranks)
    S = pd.DataFrame(scores)
    out = pd.DataFrame({'mean_rank': R.mean(axis=1), 'worst_rank': R.max(axis=1),
                        'mean_RBS': S.mean(axis=1), 'min_RBS': S.min(axis=1)})
    sym = pd.read_csv(files[0], sep='\t', dtype={'entrez_id': str}).set_index('entrez_id')['symbol']
    out['symbol'] = sym
    out = out.sort_values('mean_rank')
    return out, out.head(topk).index.tolist(), R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--components', default='results/proposed_glom')
    ap.add_argument('--topk', type=int, default=50)
    ap.add_argument('--out', default='results/validation')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    cons, cand, R = consensus_signature(args.components, args.topk)
    cons.to_csv(os.path.join(args.out, 'consensus_signature.tsv'), sep='\t')
    log('consensus signature: top %d of %d genes' % (len(cand), len(cons)))
    log('  top 20: %s' % ', '.join(cons.head(20)['symbol'].fillna('?').tolist()))

    # ---------------------------------------------------------------- 1. compartment
    log('\n=== 1. compartment reproducibility (glomerulus -> tubulointerstitium) ===')
    rows = []
    for glom, tub, lab in PAIRS:
        Xg, yg, pg = load_cohort(glom)
        Xt, yt, pt = load_cohort(tub)
        genes = Xg.columns.intersection(Xt.columns)
        gg = pd.Series(hedges_g(Xg[genes].values, yg), index=genes.astype(str))
        gt = pd.Series(hedges_g(Xt[genes].values, yt), index=genes.astype(str))
        sub = [g for g in cand if g in gg.index]
        rho_all = stats.spearmanr(gg, gt, nan_policy='omit').statistic
        rho_sig = stats.spearmanr(gg[sub], gt[sub], nan_policy='omit').statistic
        agree = float(np.mean(np.sign(gg[sub]) == np.sign(gt[sub])))
        log('  %-11s %s -> %s | signature: direction agrees %.0f%%, rho=%.3f  (all genes rho=%.3f)'
            % (lab, glom, tub, 100 * agree, rho_sig, rho_all))
        rows.append(dict(pair=lab, glom=glom, tub=tub, n_signature=len(sub),
                         direction_agreement=agree, rho_signature=rho_sig, rho_all_genes=rho_all))
    pd.DataFrame(rows).to_csv(os.path.join(args.out, 'compartment.tsv'), sep='\t', index=False)

    # ---------------------------------------------------------------- 2. specificity
    log('\n=== 2. DKD specificity: DKD vs OTHER (non-diabetic) CKD, ERCB only ===')
    spec = {}
    for c in OTHER_CKD_COHORTS:
        X, y, p = load_cohort(c, labelled_only=False)
        p = p.reset_index(drop=True)
        is_dkd = (p['label'] == 1).values
        is_ctrl = (p['label'] == 0).values
        is_other = (p['label'] == -1).values
        if is_other.sum() < 10:
            continue
        Xv = X.values
        g_vs_ctrl = hedges_g(Xv[is_dkd | is_ctrl], is_dkd[is_dkd | is_ctrl].astype(int))
        g_vs_other = hedges_g(Xv[is_dkd | is_other], is_dkd[is_dkd | is_other].astype(int))
        spec[c] = pd.DataFrame({'vs_control': g_vs_ctrl, 'vs_other_ckd': g_vs_other},
                               index=X.columns.astype(str))
        log('  %-11s DKD=%d control=%d other-CKD=%d (%s)'
            % (c, is_dkd.sum(), is_ctrl.sum(), is_other.sum(),
               ', '.join(sorted(p.loc[is_other, 'group'].unique())[:5])))

    if spec:
        tab = pd.DataFrame(index=[g for g in cand if g in list(spec.values())[0].index])
        tab['symbol'] = cons.loc[tab.index, 'symbol']
        for c, d in spec.items():
            tab['%s_vs_control' % c] = d['vs_control'].reindex(tab.index)
            tab['%s_vs_otherCKD' % c] = d['vs_other_ckd'].reindex(tab.index)
        ctrl_cols = [c for c in tab.columns if c.endswith('_vs_control')]
        other_cols = [c for c in tab.columns if c.endswith('_vs_otherCKD')]
        tab['mean_vs_control'] = tab[ctrl_cols].mean(axis=1)
        tab['mean_vs_otherCKD'] = tab[other_cols].mean(axis=1)
        # DKD-specific = still separates from other CKD, in the SAME direction
        tab['dkd_specific'] = ((tab['mean_vs_otherCKD'].abs() > 0.5) &
                               (np.sign(tab['mean_vs_otherCKD']) == np.sign(tab['mean_vs_control'])))
        tab = tab.sort_values('mean_vs_otherCKD', key=lambda s: -s.abs())
        tab.to_csv(os.path.join(args.out, 'specificity.tsv'), sep='\t')
        n_spec = int(tab['dkd_specific'].sum())
        log('\n  of the top %d candidates, %d (%.0f%%) still separate DKD from OTHER CKD '
            '(|g|>0.5, same direction) -> DKD-specific rather than generic kidney damage'
            % (len(tab), n_spec, 100 * n_spec / max(len(tab), 1)))
        log('\n  most DKD-specific:')
        for gid, r in tab.head(12).iterrows():
            log('    %-10s vs control g=%+.2f   vs other CKD g=%+.2f   %s'
                % (r['symbol'], r['mean_vs_control'], r['mean_vs_otherCKD'],
                   'DKD-SPECIFIC' if r['dkd_specific'] else 'generic CKD'))
        log('\n  least specific (strong vs control, weak vs other CKD - generic damage markers):')
        for gid, r in tab.tail(8).iterrows():
            log('    %-10s vs control g=%+.2f   vs other CKD g=%+.2f'
                % (r['symbol'], r['mean_vs_control'], r['mean_vs_otherCKD']))

    # ---------------------------------------------------------------- 3. consensus
    log('\n=== 3. how much does the signature depend on which cohort was held out? ===')
    for K in (10, 20, 50, 100):
        sets = [set(R[c].sort_values().head(K).index) for c in R.columns]
        js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
        log('  top-%-4d mean pairwise Jaccard between folds = %.3f' % (K, np.mean(js)))
    log('\nwrote %s/' % args.out)


if __name__ == '__main__':
    main()
