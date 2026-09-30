#!/usr/bin/env python
"""Do the transcriptomic findings hold in a second modality?

metabolomics_feasibility.py established that ST003255 supplies the confounder-matched contrast
the transcriptomics could not, and that only 35% of the naive signal survives it. This script
runs the rest of the paper's controls against the same data, so the comparison is like for like
rather than a single headline number:

    experiment 1   random-signature null. Section 3.1 showed random 50-gene sets reach AUROC
                   0.68-0.88 in transcriptomics, making the raw number uninformative. Does the
                   same hold for random metabolite sets?
    experiment 2   classification. How well does each contrast separate, measured honestly with
                   cross-validation and read against the null from experiment 1?
    experiment 3   shared damage. How much of the DKD signal is simply "kidney is failing"?
                   Measured as the correlation between the DKD-vs-healthy and
                   otherCKD-vs-healthy effect vectors.
    experiment 4   what survives. The metabolites that clear the matched contrast, named, so a
                   clinician can judge whether they are plausible.

One cohort, so leave-one-dataset-out is not available; stratified cross-validation is used and
labelled as such. Feature selection runs inside each fold, never on the full data.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

SRC = 'data/raw/metabolomics/ST003255/AN005337_datatable.tsv'
OUT = 'results/metabolomics'
SEED = 0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load():
    """원시 강도만 돌려준다. 대치·로그·정규화는 하지 않는다.

    처음에는 여기서 전체 표본에 대해 중앙값 대치와 z-정규화를 하고 그 뒤에 fold 를
    나눴다. 선택은 fold 안에서 했지만 전처리가 test fold 의 중앙값·평균·표준편차를
    이미 본 상태였다. 표본이 작아 영향은 크지 않을 수 있으나, 그 상태로 "cross-validated
    AUROC" 라고 쓸 수는 없다. 전처리도 학습 fold 에서 적합해 test fold 에 적용한다.
    """
    d = pd.read_csv(SRC, sep='\t')
    d = d.rename(columns={d.columns[0]: 'sample', d.columns[1]: 'group'})
    d['group'] = d['group'].str.replace('Disease:', '', regex=False)
    mets = [c for c in d.columns if c not in ('sample', 'group')]
    X = d[mets].apply(pd.to_numeric, errors='coerce')
    X = X.loc[:, X.notna().mean() >= 0.8]
    return X.values, d['group'].values, list(X.columns)


def prep_fit(Xtr):
    """학습 fold 에서만 전처리 모수를 뽑는다."""
    med = np.nanmedian(Xtr, axis=0)
    med = np.where(np.isnan(med), 0.0, med)
    f = Xtr[np.isfinite(Xtr) & (Xtr > 0)]
    floor = (f.min() / 2.0) if f.size else 1e-6
    A = np.log2(np.clip(np.where(np.isnan(Xtr), med, Xtr), floor, None))
    mu = A.mean(0)
    sd = A.std(0, ddof=1)
    return med, floor, mu, sd


def prep_apply(X, par):
    med, floor, mu, sd = par
    A = np.log2(np.clip(np.where(np.isnan(X), med, X), floor, None))
    return (A - mu) / (sd + 1e-12)


def prep_all(X):
    """CV 가 아닌 곳(효과크기)에서 쓰는 전역 전처리.

    Hedges g 는 특징별 선형변환에 불변이므로 여기서는 전역으로 해도 누수가 아니다.
    학습/검증 분리가 없는 계산이기 때문이다.
    """
    return prep_apply(X, prep_fit(X))


def cv_auc(X, y, feats=None, folds=5, seed=SEED, k=None):
    """AUROC by stratified CV. 전처리와 특징 선택 모두 학습 fold 안에서만 적합한다."""
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    aucs = []
    for tr, te in skf.split(X, y):
        par = prep_fit(X[tr])
        Ztr, Zte = prep_apply(X[tr], par), prep_apply(X[te], par)
        cols = feats
        if k is not None:
            t, _ = stats.ttest_ind(Ztr[y[tr] == 1], Ztr[y[tr] == 0], equal_var=False)
            cols = np.argsort(-np.abs(np.nan_to_num(t)))[:k]
        A, B = (Ztr[:, cols], Zte[:, cols]) if cols is not None else (Ztr, Zte)
        m = LogisticRegression(max_iter=2000, C=1.0)
        m.fit(A, y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(B)[:, 1]))
    return float(np.mean(aucs))


def main():
    os.makedirs(OUT, exist_ok=True)
    Xraw, g, mets = load()
    Z = prep_all(Xraw)          # 효과크기 전용. CV 는 fold 안에서 따로 전처리한다.
    rng = np.random.default_rng(SEED)
    p_ = Xraw.shape[1]
    K = 20                      # 미리 고정한 서명 크기. 아래에서 민감도를 함께 낸다.

    contrasts = {
        'DKD vs 정상': (g == 'DKD', g == 'Normal'),
        'DKD vs 다른 신장병': (g == 'DKD', np.isin(g, ['IgAN', 'MN', 'HN'])),
        '다른 신장병 vs 정상': (np.isin(g, ['IgAN', 'MN', 'HN']), g == 'Normal'),
    }

    # ---------------------------------------------------------- 1 + 2
    log('=== 실험 1·2: 무작위 대조와 분류 성능 ===')
    log('  %-22s %6s %8s %8s %10s %10s'
        % ('대조', 'n', '실제', '무작위', '무작위95%', '라벨섞기'))
    rows = []
    for name, (ma, mb) in contrasts.items():
        keep = ma | mb
        Zc, yc = Xraw[keep], ma[keep].astype(int)
        real = cv_auc(Zc, yc, k=K)
        rand = [cv_auc(Zc, yc, feats=rng.choice(p_, K, replace=False)) for _ in range(60)]
        perm = []
        for _ in range(30):
            yp = rng.permutation(yc)
            perm.append(cv_auc(Zc, yp, k=K))
        rows.append(dict(contrast=name, n=int(keep.sum()), auroc_selected=real,
                         auroc_random_mean=float(np.mean(rand)),
                         auroc_random_p95=float(np.percentile(rand, 95)),
                         auroc_permuted=float(np.mean(perm)),
                         beats_random=bool(real > np.percentile(rand, 95))))
        log('  %-22s %6d %8.3f %8.3f %10.3f %10.3f'
            % (name, keep.sum(), real, np.mean(rand), np.percentile(rand, 95), np.mean(perm)))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'experiments_auroc.tsv'), sep='\t', index=False)

    # ---------------------------------------------------------- 2b
    # K=20 은 미리 고정한 값이다. 그 선택이 결론을 만드는 것이 아님을 보이려고
    # 다른 크기에서도 같은 대조를 돌린다. 묻는 것은 절대 AUROC 가 아니라
    # "선택이 무작위를 넘는가" 이므로, 그 판정이 K 에 따라 뒤집히는지만 본다.
    log('')
    log('=== 실험 2b: 서명 크기 K 민감도 ===')
    log('  %-22s %5s %8s %8s %10s %s' % ('대조', 'K', '실제', '무작위', '무작위95%', '판정'))
    krows = []
    for name, (ma, mb) in contrasts.items():
        keep = ma | mb
        Zc, yc = Xraw[keep], ma[keep].astype(int)
        for k_ in (10, 20, 50):
            real = cv_auc(Zc, yc, k=k_)
            rand = [cv_auc(Zc, yc, feats=rng.choice(p_, k_, replace=False))
                    for _ in range(40)]
            p95 = float(np.percentile(rand, 95))
            krows.append(dict(contrast=name, K=k_, auroc_selected=real,
                              auroc_random_mean=float(np.mean(rand)),
                              auroc_random_p95=p95, beats_random=bool(real > p95)))
            log('  %-22s %5d %8.3f %8.3f %10.3f %s'
                % (name, k_, real, np.mean(rand), p95,
                   '무작위 초과' if real > p95 else '구분 안 됨'))
    pd.DataFrame(krows).to_csv(os.path.join(OUT, 'experiments_k_sensitivity.tsv'),
                               sep='\t', index=False)

    # ---------------------------------------------------------- 3
    log('')
    log('=== 실험 3: DKD 신호 중 "신장이 나빠서" 생긴 부분 ===')

    def effect(ma, mb):
        return np.array([(Z[ma, j].mean() - Z[mb, j].mean()) for j in range(p_)])

    e_dkd = effect(g == 'DKD', g == 'Normal')
    e_oth = effect(np.isin(g, ['IgAN', 'MN', 'HN']), g == 'Normal')
    r, pval = stats.pearsonr(e_dkd, e_oth)
    shared = r ** 2
    log('  DKD-vs-정상 효과와 다른신장병-vs-정상 효과의 상관: r = %.3f (p = %.1e)' % (r, pval))
    log('  즉 DKD 신호 분산의 %.0f%%는 신장병 일반과 공유됩니다.' % (100 * shared))
    pd.DataFrame({'metabolite': mets, 'g_DKD_vs_normal': e_dkd,
                  'g_otherCKD_vs_normal': e_oth}).to_csv(
        os.path.join(OUT, 'experiments_shared_effect.tsv'), sep='\t', index=False)

    # ---------------------------------------------------------- 4
    log('')
    log('=== 실험 4: 교란 매칭을 통과한 대사물질 ===')
    m = pd.read_csv(os.path.join(OUT, 'ST003255_DKD_vs_other_kidney_disease.tsv'), sep='\t')
    sig = m[m['q'] < 0.05].reindex(m[m['q'] < 0.05]['hedges_g'].abs()
                                   .sort_values(ascending=False).index)
    log('  q < 0.05 통과 %d개 중 상위 12개:' % len(sig))
    for _, r_ in sig.head(12).iterrows():
        log('    %-42s g = %+.2f  q = %.1e' % (r_['metabolite'][:40], r_['hedges_g'], r_['q']))
    sig.to_csv(os.path.join(OUT, 'experiments_survivors.tsv'), sep='\t', index=False)

    log('')
    log('=== 요약 ===')
    t = pd.DataFrame(rows)
    for _, r_ in t.iterrows():
        verdict = '무작위보다 유의하게 좋음' if r_['beats_random'] else '무작위와 구분 안 됨'
        log('  %-22s AUROC %.3f  (무작위 95%% 상한 %.3f) -> %s'
            % (r_['contrast'], r_['auroc_selected'], r_['auroc_random_p95'], verdict))
    log('')
    log('  전사체에서와 같은 교훈: 원시 AUROC가 아니라 무작위 분포 대비로 읽어야 합니다.')


if __name__ == '__main__':
    main()
