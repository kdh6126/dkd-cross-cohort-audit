#!/usr/bin/env python
"""Is the IEG-score correction anything more than RUVg with a hand-picked control set?

This has to be asked directly. RUVg (Risso et al. 2014) estimates factors of unwanted variation
from a set of NEGATIVE CONTROL GENES - genes assumed unaffected by the biology - and removes
them. Our correction takes the immediate-early module, forms a per-sample score, and
residualises on it. Structurally that is RUVg with k=1 and a hand-chosen control set, so the
novelty claim rests entirely on whether choosing THAT control set does something the standard
recipes do not.

Arms:
  none            no correction
  IEG_resid       this project's correction: residualise on the mean z of the IEG module
  RUVg_ieg        RUVg with the IEG genes as negative controls, k factors from SVD
  RUVg_empirical  RUVg with EMPIRICAL controls - the genes least associated with case status,
                  which is what one uses when no control set is known. The fair comparison.
  SVA_like        surrogate variables: SVD of residuals after projecting out the design matrix,
                  which is the protect-the-variable-of-interest idea behind SVA
  ComBat          empirical-Bayes location/scale adjustment by cohort, the standard batch fix

ComBat is included knowing it should fail: the confounder is collinear with case status INSIDE
each cohort, and adjusting by cohort cannot touch a within-cohort confounder. Showing that is
part of the point - it is the correction most people would reach for first.
"""
import os, sys, argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM
from dkd_rbs import rbs_fit
from dkd_deconfound import handling_score, IEG

OUT = 'results/ruv_benchmark'
TOPK = [10, 20, 50]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- corrections
def ruvg(Z, control_cols, k=1):
    """RUVg: SVD of the control-gene submatrix gives W; regress every gene on W and remove.

    Risso et al. 2014, in its simplest form. Controls are centred first so W captures variation
    shared by genes that should carry none of the signal of interest.
    """
    C = Z[:, control_cols]
    C = C - C.mean(0)
    U, S, Vt = np.linalg.svd(C, full_matrices=False)
    W = U[:, :k] * S[:k]                       # sample-level factors of unwanted variation
    W = W - W.mean(0)
    G = W.T @ W
    if np.linalg.matrix_rank(G) < W.shape[1]:
        return Z.copy()
    alpha = np.linalg.solve(G, W.T @ Z)
    return Z - W @ alpha


def sva_like(Z, y, k=2):
    """Surrogate variables: remove the design, take the leading PCs of what is left, then
    project those out of the data while leaving the design untouched."""
    X = np.column_stack([np.ones(len(y)), y.astype(float)])
    beta = np.linalg.lstsq(X, Z, rcond=None)[0]
    R = Z - X @ beta
    U, S, Vt = np.linalg.svd(R - R.mean(0), full_matrices=False)
    W = U[:, :k] * S[:k]
    # orthogonalise the surrogate variables against the design so the biology is protected
    W = W - X @ np.linalg.lstsq(X, W, rcond=None)[0]
    G = W.T @ W
    if np.linalg.matrix_rank(G) < W.shape[1]:
        return Z.copy()
    return Z - W @ np.linalg.solve(G, W.T @ Z)


def combat(Zs, ys):
    """Empirical-Bayes location/scale adjustment with cohort as the batch."""
    allZ = np.vstack(Zs)
    gmean, gsd = allZ.mean(0), allZ.std(0, ddof=1) + 1e-12
    out = []
    for Z in Zs:
        m, s = Z.mean(0), Z.std(0, ddof=1) + 1e-12
        # shrink each batch's parameters toward the grand values
        lam = len(Z) / (len(Z) + 10.0)
        m_ = lam * m + (1 - lam) * gmean
        s_ = lam * s + (1 - lam) * gsd
        out.append((Z - m_) / s_ * gsd + gmean)
    return out


def empirical_controls(Zs, ys, n=200):
    """Genes least associated with case status pooled over cohorts - the standard way to get
    negative controls when no curated set exists."""
    X = np.vstack(Zs); y = np.concatenate(ys)
    x1, x0 = X[y == 1], X[y == 0]
    sp = np.sqrt(x1.var(0, ddof=1) / len(x1) + x0.var(0, ddof=1) / len(x0)) + 1e-12
    t = np.abs(x1.mean(0) - x0.mean(0)) / sp
    return np.argsort(t)[:n]


def ext_auc(Xtr, ytr, Xte, yte, cols):
    if len(np.unique(yte)) < 2:
        return np.nan
    m = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000,
                           class_weight='balanced').fit(Xtr[:, cols], ytr)
    return roc_auc_score(yte, m.predict_proba(Xte[:, cols])[:, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--B', type=int, default=150)
    ap.add_argument('--k', type=int, default=50)
    ap.add_argument('--base', default='relieff')
    ap.add_argument('--ruv-k', type=int, default=2)
    ap.add_argument('--cohorts', default=','.join(GLOM) + ',GSE142025')
    ap.add_argument('--seed', type=int, default=20260825)
    # 선택기를 바꿔 가며 같은 보정 비교를 돌릴 수 있게 출력 경로를 인자로 받는다.
    # 기본값은 기존 경로라 이미 등록된 단계의 동작은 그대로다.
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args()
    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])
    gpos = {g: i for i, g in enumerate(genes)}
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = {s: g for g, s in zip(gs['entrez_id'], gs['symbol']) if isinstance(s, str)}
    ieg_cols = [gpos[sym2id[s]] for s in IEG if s in sym2id and sym2id[s] in gpos]

    spec = pd.read_csv('results/final_candidates.tsv', sep='\t')
    spec = spec.rename(columns={spec.columns[0]: 'entrez_id'})
    spec['entrez_id'] = spec['entrez_id'].astype(str)
    spec = spec.set_index('entrez_id')
    ieg_ids = set(genes[i] for i in ieg_cols)

    arms = ['none', 'IEG_resid', 'RUVg_ieg', 'RUVg_empirical', 'SVA_like', 'ComBat']
    rows, orders = [], {a: {} for a in arms}
    for arm in arms:
        log('arm %s ...' % arm)
        for held in cohorts:
            tr = [c for c in cohorts if c != held]
            Zs = [data[c][0].values.copy() for c in tr]
            ys = [data[c][1] for c in tr]

            if arm == 'ComBat':
                Zs = combat(Zs, ys)
            else:
                emp = empirical_controls(Zs, ys) if arm == 'RUVg_empirical' else None
                new = []
                for Z, yc in zip(Zs, ys):
                    if arm == 'IEG_resid':
                        s = handling_score(Z, ieg_cols)
                        d = float(s @ s)
                        Z = Z - np.outer(s, (s @ Z) / d) if d else Z
                    elif arm == 'RUVg_ieg':
                        Z = ruvg(Z, ieg_cols, k=min(args.ruv_k, len(ieg_cols) - 1))
                    elif arm == 'RUVg_empirical':
                        Z = ruvg(Z, emp, k=args.ruv_k)
                    elif arm == 'SVA_like':
                        Z = sva_like(Z, yc, k=args.ruv_k)
                    new.append(Z)
                Zs = new

            sub = {c: (Z, y) for c, Z, y in zip(tr, Zs, ys)}
            r = rbs_fit(sub, k=args.k, B=args.B, base=args.base, seed=args.seed, agg='geomean')
            order = np.argsort(-r['score'])
            orders[arm][held] = order

            Xtr = np.vstack([data[c][0].values for c in tr])
            ytr = np.concatenate([data[c][1] for c in tr])
            Xte, yte = data[held][0].values, data[held][1]
            for K in TOPK:
                top = genes[order[:K]]
                rows.append(dict(arm=arm, held_out=held, K=K,
                                 external_auroc=ext_auc(Xtr, ytr, Xte, yte, order[:K]),
                                 n_ieg=sum(1 for g in top if g in ieg_ids),
                                 procurement_frac=float(spec['procurement_flag']
                                                        .reindex(top).fillna(False).mean()),
                                 spec_pass_frac=float(spec['passes_specificity']
                                                      .reindex(top).fillna(False).mean())))

    df = pd.DataFrame(rows)
    for arm in arms:
        for K in TOPK:
            sets = [set(orders[arm][h][:K]) for h in cohorts]
            js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
            df.loc[(df['arm'] == arm) & (df['K'] == K), 'cross_fold_jaccard'] = float(np.mean(js))
    df.to_csv(os.path.join(out_dir, 'benchmark.tsv'), sep='\t', index=False)

    agg = (df.groupby(['arm', 'K'])
             .agg(auc=('external_auroc', 'mean'), ieg=('n_ieg', 'mean'),
                  proc=('procurement_frac', 'mean'), sp=('spec_pass_frac', 'mean'),
                  xf=('cross_fold_jaccard', 'mean')).reset_index())
    log('\n=== correction benchmark ===')
    for K in TOPK:
        log('\n  K = %d' % K)
        log('  %-16s %8s %10s %12s %10s %s'
            % ('arm', 'IEG', 'procure %', 'specificity %', 'AUROC', 'cross-fold J'))
        for _, r in agg[agg['K'] == K].sort_values('proc').iterrows():
            log('  %-16s %8.1f %9.0f%% %11.0f%% %10.3f %12.3f'
                % (r['arm'], r['ieg'], 100 * r['proc'], 100 * r['sp'], r['auc'], r['xf']))
    agg.to_csv(os.path.join(out_dir, 'summary.tsv'), sep='\t', index=False)
    log('\nwrote %s/' % out_dir)


if __name__ == '__main__':
    main()
