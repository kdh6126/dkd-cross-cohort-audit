#!/usr/bin/env python
"""How many genes pass the DKD-specificity gate by chance?

The gate keeps genes with |Hedges g| >= 0.5 between DKD and non-diabetic CKD in ERCB. DKD
there is only 12 (glomeruli) and 17 (tubulointerstitium) patients, so that statistic is noisy
and a gate applied to 9,900 genes will let some through on noise alone. Reporting the final
candidate list without knowing that rate would be unwarranted.

Two nulls:
  PERMUTED   shuffle the DKD / other-CKD labels, re-run the gate -> false-positive rate
  OBSERVED   the real gate

The comparison gives an empirical FDR for the gate, and the permutation distribution also
tells us whether the *top* candidates (|g| ~ 1.0+) are reachable by chance at all.

The permutation has to respect that the two ERCB blocks are the same patients
--------------------------------------------------------------------------------
GSE104948 (glomeruli) and GSE104954 (tubulointerstitium) profile the same ERCB patients in two
compartments. Of the samples entering this gate, 54 patients appear in both blocks and all 54
carry the same diagnosis in both. Shuffling the labels independently inside each block destroys
that dependence: the two per-gene g vectors become independent noise, averaging them cancels
more of it than it does in the real data, fewer genes clear the threshold, and the empirical
FDR comes out too low. Measured, it was 13.7% independent against 16.9% paired.

So the labels are permuted at the PATIENT level. The shared patients get one label each, used
in both blocks; the block-specific patients are permuted within their block. The case counts of
each block are preserved exactly (5 of the 54 shared patients are cases, plus 7 and 12 among the
patients unique to each block, giving the observed 12 and 17). The independent permutation is
still computed and reported alongside, so the difference stays visible rather than asserted.
"""
import os, re, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_cohort

ERCB = ['GSE104948', 'GSE104954']


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--thresh', type=float, default=0.5)
    ap.add_argument('--n-perm', type=int, default=500)
    ap.add_argument('--seed', type=int, default=20260824)
    ap.add_argument('--candidates', default='results/final_candidates_gated.tsv')
    ap.add_argument('--out', default='results/specificity_null.tsv')
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    blocks, obs_list, perm_counts = {}, [], None

    for c in ERCB:
        X, y, p = load_cohort(c, labelled_only=False)
        X.columns = X.columns.astype(str)
        p = p.reset_index(drop=True)
        M = X.values
        is_dkd = (p['label'] == 1).values
        is_oth = (p['label'] == -1).values
        sub = M[is_dkd | is_oth]
        lab = is_dkd[is_dkd | is_oth]
        sid = np.array([re.sub('-(Glom|Tub)-', '-', str(v))
                        for v in p['subject_id'].values[is_dkd | is_oth]])
        blocks[c] = (sub, lab, X.columns, sid)
        log('  %-11s DKD=%d vs other-CKD=%d' % (c, lab.sum(), (~lab).sum()))

    # 환자 단위 짝 순열의 뼈대. 두 블록에 모두 있는 환자는 라벨을 한 번만 받는다.
    i0, i1 = blocks[ERCB[0]][3], blocks[ERCB[1]][3]
    lab0, lab1 = dict(zip(i0, blocks[ERCB[0]][1])), dict(zip(i1, blocks[ERCB[1]][1]))
    shared = sorted(set(i0) & set(i1))
    uniq = {ERCB[0]: sorted(set(i0) - set(shared)), ERCB[1]: sorted(set(i1) - set(shared))}
    k_shared = sum(1 for v in shared if lab0[v])
    agree = sum(1 for v in shared if lab0[v] == lab1[v])
    n_case = {c: int(blocks[c][1].sum()) for c in ERCB}
    log('  두 블록에 모두 있는 환자 %d명 (라벨 일치 %d명), 그중 사례 %d명'
        % (len(shared), agree, k_shared))

    def paired_labels(rg):
        """환자 단위로 섞되 각 블록의 사례 수는 관측값 그대로 둔다."""
        pick = set(rg.choice(len(shared), k_shared, replace=False).tolist())
        asg = {v: (j in pick) for j, v in enumerate(shared)}
        for c in ERCB:
            u = uniq[c]
            need = n_case[c] - k_shared
            sel = set(rg.choice(len(u), need, replace=False).tolist()) if need > 0 else set()
            for j, v in enumerate(u):
                asg[v] = (j in sel)
        return {c: np.array([asg[v] for v in blocks[c][3]]) for c in ERCB}

    genes = blocks[ERCB[0]][2]
    obs = np.nanmean(np.vstack([hedges_g(blocks[c][0][blocks[c][1]],
                                         blocks[c][0][~blocks[c][1]]) for c in ERCB]), axis=0)
    n_obs = int(np.nansum(np.abs(obs) >= args.thresh))

    def permute(paired):
        rg = np.random.default_rng(args.seed)
        cnt, mx = [], []
        for _ in range(args.n_perm):
            pls = (paired_labels(rg) if paired
                   else {c: rg.permutation(blocks[c][1]) for c in ERCB})
            gs = [hedges_g(blocks[c][0][pls[c]], blocks[c][0][~pls[c]]) for c in ERCB]
            gp = np.nanmean(np.vstack(gs), axis=0)
            cnt.append(int(np.nansum(np.abs(gp) >= args.thresh)))
            mx.append(float(np.nanmax(np.abs(gp))))
        return np.array(cnt), np.array(mx)

    ind_counts, ind_maxes = permute(False)
    counts, maxes = permute(True)

    fdr = counts.mean() / max(n_obs, 1)
    log('\n=== specificity gate at |g| >= %.2f, %d genes tested ===' % (args.thresh, len(genes)))
    log('  observed passing        : %d (%.1f%%)' % (n_obs, 100 * n_obs / len(genes)))
    log('  permuted passing (mean) : %.0f  [5-95%%: %.0f-%.0f]'
        % (counts.mean(), np.percentile(counts, 5), np.percentile(counts, 95)))
    log('  empirical FDR of the gate: %.1f%%   <- patient-paired permutation' % (100 * fdr))
    log('  the same, shuffling each block independently: %.1f%%'
        % (100 * ind_counts.mean() / max(n_obs, 1)))
    log('  the blocks are the same patients, so independent shuffling cancels noise that')
    log('  the real data keeps, and reports a gate cleaner than it is')
    log('  largest |g| ever reached under permutation: %.2f  [95th pct %.2f]'
        % (maxes.max(), np.percentile(maxes, 95)))

    cand = pd.read_csv(args.candidates, sep='\t')
    if 'g_DKD_vs_otherCKD' in cand.columns:
        v = cand['g_DKD_vs_otherCKD'].abs()
        above = int((v > np.percentile(maxes, 95)).sum())
        log('\n  of the %d final candidates, %d exceed the 95th percentile of the permutation'
            % (len(cand), above))
        log('  maximum (%.2f) -> those cannot be explained by label noise in this contrast'
            % np.percentile(maxes, 95))
        log('  candidate |g| range: %.2f - %.2f' % (v.min(), v.max()))

    pd.DataFrame(dict(perm_count=counts, perm_max_abs_g=maxes,
                      indep_perm_count=ind_counts,
                      indep_perm_max_abs_g=ind_maxes)).to_csv(
        args.out, sep='\t', index=False)
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
