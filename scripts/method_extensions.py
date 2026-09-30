#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""교란을 선택 단계에 넣는 세 가지 방법을 같은 fold 에서 비교한다.

지금까지의 보정은 두 가지였다. 잔차화(handling score 축을 각 gene 에서 뺀 뒤 선택)와
penalty(RBS 점수에 1-|corr| 를 곱함). 둘 다 손잡이(lambda)를 손으로 정했고, 교란을
선택 모형 **안에** 넣지는 않았다. 세 가지를 더 해 본다.

    resid_learned   lambda 를 fold 안에서 고른다. 기준은 학습 cohort 만 쓴다 —
                    IEG module 이 상위 50에 1개 이하로 남는 lambda 중에서
                    inner-LODO AUROC 가 가장 높은 것. hold-out 도 specificity gate 도
                    보지 않는다. gate 는 ERCB 가 학습에 들어 있어 쓰면 누수다.
    joint_score     잔차화는 순차적이다 — handling 을 먼저 빼고 그다음 선택. 조건부로
                    바른 것은 handling 을 nuisance 로 둔 efficient score test 다.
                    logit(y) ~ gene + handling 에서 gene 에 대한 score 통계량을 쓴다.
    net_propagate   co-expression graph 위에서 기존 점수를 전파한다. 이것은 graph-Laplacian
                    정규화 선택(Li & Li 2008)이 **아니고**, 새 네트워크 방법의 성능 근거로
                    쓸 수 있는 것도 아니다. 이웃한 gene 이 서로를 끌어올리면 어떻게 되는지
                    보는 탐색적 ranking 이다. IEG module 이 강하게 뭉친 덩어리이므로 교란이
                    유지되거나 증폭될 것으로 예상한다.

비교 기준선은 같은 스크립트 안에 둔다(none, resid_fixed). 다른 표의 숫자와 섞어 읽지
않기 위해서다. 선택은 pooled |Hedges g| 랭킹으로 통일했다 — 여기서 묻는 것은 selector
가 아니라 교란을 다루는 방식이기 때문이다.

이 결과를 무엇으로 읽을 것인가
------------------------------
LODO fold 가 4개뿐이고 세 구현 모두 탐색적이다. 어느 것도 방법 기여로 주장할 단계가
아니며, 원고에 넣지 않는다. 지금 이 자료가 주는 가장 정확한 메시지는 하나다.

    네트워크 평활화는 AUROC 와 재현성 수치를 높일 수 있지만, 동시에 IEG·조달 관련
    신호를 유지하거나 강화할 수 있다. 단일 성능 지표의 개선은 교란 제거를 뜻하지 않는다.

이 문장도 정확한 조건부 모형과 독립 cohort 로 확인한 뒤에야 본문 주장이 될 수 있다.
"""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM
from dkd_deconfound import IEG, handling_score, residualise

OUT = 'results/method_extensions'
K = 50
LAMBDAS = (0.0, 0.25, 0.5, 0.75, 1.0)
ARMS = ('none', 'resid_fixed', 'resid_learned', 'joint_score', 'net_propagate')


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(X, y):
    a, b = X[y == 1], X[y == 0]
    if len(a) < 2 or len(b) < 2:
        return np.zeros(X.shape[1])
    sp = np.sqrt(((len(a) - 1) * a.var(0, ddof=1) + (len(b) - 1) * b.var(0, ddof=1))
                 / (len(a) + len(b) - 2))
    sp[sp == 0] = np.nan
    return np.nan_to_num((a.mean(0) - b.mean(0)) / sp)


def pooled_score(blocks, lam=None, ieg_cols=None):
    """cohort 별 |g| 의 평균. lam 이 주어지면 먼저 잔차화한다."""
    gs = []
    for Z, y in blocks:
        if lam:
            Z = residualise(Z, handling_score(Z, ieg_cols), y=y, lam=lam)
        gs.append(np.abs(hedges_g(Z, y)))
    return np.mean(gs, axis=0)


def joint_score(blocks, ieg_cols):
    """logit(y) ~ gene + handling 에서 gene 에 대한 efficient score test.

    처음 쓴 판은 귀무 모형(handling 만)의 잔차와 gene 의 가중 중심화 상관이었다. 그것은
    score test 가 아니다. nuisance 공변량에 대한 **가중 사영을 빼지 않아서** 분산이
    틀리고, 그 오차의 크기가 gene 마다 corr(gene, handling) 에 따라 달라진다. 즉 순위가
    바로 그 상관에 의해 체계적으로 왜곡된다 — 여기서 알고 싶은 것이 정확히 그 축인데.

    제대로 된 형태는 nuisance 설계 S = [1, handling] 에 대해

        r      = y - p                       귀무 모형 잔차
        W      = diag(p(1-p))                관측 가중
        X_perp = X - S (S'WS)^-1 S'W X       가중 사영을 뺀 나머지
        z      = X'r / sqrt(X_perp' W X_perp)

    이렇게 하면 handling 과 겹치는 부분이 분자와 분모 양쪽에서 제거되어, 남는 것이
    handling 으로 설명되지 않는 gene 의 기여가 된다.
    """
    out = []
    for Z, y in blocks:
        s = handling_score(Z, ieg_cols)
        S = np.column_stack([np.ones_like(s), s])
        m = LogisticRegression(max_iter=1000, C=1e6)
        m.fit(s.reshape(-1, 1), y)
        p = np.clip(m.predict_proba(s.reshape(-1, 1))[:, 1], 1e-6, 1 - 1e-6)
        r = y - p
        w = p * (1 - p)
        SW = S * w[:, None]
        A = np.linalg.pinv(S.T @ SW)            # (S'WS)^-1
        Zperp = Z - S @ (A @ (SW.T @ Z))        # 가중 사영을 뺀 나머지
        num = Z.T @ r
        den = np.sqrt((w[:, None] * Zperp ** 2).sum(0)) + 1e-12
        out.append(np.abs(num / den))
    return np.mean(out, axis=0)


def net_propagate_score(blocks, base, pool=2000, alpha=0.5, iters=10):
    """co-expression graph 위에서 기존 점수를 전파한다.

    graph-Laplacian 정규화 선택이 아니다. 그쪽은 벌점항을 목적함수에 넣고 계수를 함께
    추정하는 방법이고, 여기서는 이미 계산된 |g| 점수를 그래프에서 퍼뜨릴 뿐이다.
    탐색용이며, 네트워크 방법의 성능 근거로 인용할 수 없다.
    """
    idx = np.argsort(-base)[:pool]
    A = np.zeros((len(idx), len(idx)))
    for Z, _ in blocks:
        Y = Z[:, idx]
        Y = (Y - Y.mean(0)) / (Y.std(0, ddof=1) + 1e-12)
        A += np.abs((Y.T @ Y) / (len(Y) - 1))
    A /= len(blocks)
    np.fill_diagonal(A, 0.0)
    thr = np.percentile(A, 99.0)        # 상위 1% 만 간선으로
    A = np.where(A >= thr, A, 0.0)
    d = A.sum(1) + 1e-12
    W = A / np.sqrt(np.outer(d, d))
    s0 = base[idx].copy()
    s = s0.copy()
    for _ in range(iters):
        s = (1 - alpha) * s0 + alpha * (W @ s)
    full = np.zeros_like(base)
    full[idx] = s
    return full


def inner_lodo_auc(blocks, names, lam, ieg_cols, k=K):
    """학습 cohort 안에서만 도는 LODO. lambda 를 고르는 데 쓴다."""
    if len(blocks) < 2:
        return np.nan
    aucs = []
    for i in range(len(blocks)):
        tr = [blocks[j] for j in range(len(blocks)) if j != i]
        Zte, yte = blocks[i]
        if len(np.unique(yte)) < 2:
            continue
        sc = pooled_score(tr, lam=lam, ieg_cols=ieg_cols)
        cols = np.argsort(-sc)[:k]
        Xtr = np.vstack([b[0] for b in tr])
        ytr = np.concatenate([b[1] for b in tr])
        m = LogisticRegression(max_iter=5000, C=1.0, class_weight='balanced')
        m.fit(Xtr[:, cols], ytr)
        aucs.append(roc_auc_score(yte, m.predict_proba(Zte[:, cols])[:, 1]))
    return float(np.mean(aucs)) if aucs else np.nan


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)
    cohorts = GLOM + ['GSE142025']
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])
    gpos = {g: i for i, g in enumerate(genes)}
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = {s: g for g, s in zip(gs['entrez_id'], gs['symbol']) if isinstance(s, str)}
    ieg_cols = [gpos[sym2id[s]] for s in IEG if s in sym2id and sym2id[s] in gpos]
    ieg_set = set(ieg_cols)
    log('IEG module %d개, gene %d개, cohort %d개' % (len(ieg_cols), len(genes), len(cohorts)))

    spec = pd.read_csv('results/final_candidates.tsv', sep='\t')
    spec = spec.rename(columns={spec.columns[0]: 'entrez_id'})
    spec['entrez_id'] = spec['entrez_id'].astype(str)
    passes = dict(zip(spec['entrez_id'], spec['passes_specificity'].astype(bool)))
    proc = dict(zip(spec['entrez_id'], spec['procurement_flag'].astype(bool)))

    rows, orders = [], {a: {} for a in ARMS}
    for held in cohorts:
        tr = [c for c in cohorts if c != held]
        blocks = [(data[c][0].values, data[c][1]) for c in tr]
        Xtr = np.vstack([b[0] for b in blocks])
        ytr = np.concatenate([b[1] for b in blocks])
        Xte, yte = data[held][0].values, data[held][1]
        log('')
        log('== hold out %s ==' % held)

        base_none = pooled_score(blocks)
        scores = {'none': base_none,
                  'resid_fixed': pooled_score(blocks, lam=1.0, ieg_cols=ieg_cols)}

        # ---- lambda 학습. 학습 cohort 만 본다.
        cand = []
        for lam in LAMBDAS:
            sc = pooled_score(blocks, lam=lam or None, ieg_cols=ieg_cols)
            n_ieg = sum(1 for j in np.argsort(-sc)[:K] if j in ieg_set)
            cand.append((lam, n_ieg, inner_lodo_auc(blocks, tr, lam or None, ieg_cols)))
        ok = [c for c in cand if c[1] <= 1]
        pick = max(ok or cand, key=lambda c: (-c[1] if not ok else 0,
                                              c[2] if np.isfinite(c[2]) else -1))
        log('   lambda 후보 (lam, IEG수, inner AUROC): %s'
            % ', '.join('(%.2f, %d, %.3f)' % c for c in cand))
        log('   -> 고른 lambda = %.2f' % pick[0])
        scores['resid_learned'] = pooled_score(blocks, lam=pick[0] or None,
                                               ieg_cols=ieg_cols)
        scores['joint_score'] = joint_score(blocks, ieg_cols)
        scores['net_propagate'] = net_propagate_score(blocks, base_none)

        for arm in ARMS:
            order = np.argsort(-scores[arm])
            orders[arm][held] = order
            cols = order[:K]
            gid = genes[cols]
            m = LogisticRegression(max_iter=5000, C=1.0, class_weight='balanced')
            m.fit(Xtr[:, cols], ytr)
            auc = (roc_auc_score(yte, m.predict_proba(Xte[:, cols])[:, 1])
                   if len(np.unique(yte)) > 1 else np.nan)
            rows.append(dict(arm=arm, held_out=held, K=K, external_auroc=auc,
                             lam=pick[0] if arm == 'resid_learned' else np.nan,
                             n_ieg=sum(1 for j in cols if j in ieg_set),
                             procurement_frac=sum(1 for g in gid if proc.get(g, False)) / K,
                             spec_pass_frac=sum(1 for g in gid if passes.get(g, False)) / K))
            log('   %-14s AUROC %.3f  IEG %2d  procurement %2.0f%%  specificity %2.0f%%'
                % (arm, auc, rows[-1]['n_ieg'], 100 * rows[-1]['procurement_frac'],
                   100 * rows[-1]['spec_pass_frac']))

    res = pd.DataFrame(rows)
    for arm in ARMS:
        sets = [set(orders[arm][h][:K]) for h in cohorts]
        js = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
        res.loc[res['arm'] == arm, 'cross_fold_jaccard'] = float(np.mean(js))
    res.to_csv(os.path.join(OUT, 'extensions.tsv'), sep='\t', index=False)

    agg = (res.groupby('arm').agg(auroc=('external_auroc', 'mean'),
                                  ieg=('n_ieg', 'mean'),
                                  proc=('procurement_frac', 'mean'),
                                  spec=('spec_pass_frac', 'mean'),
                                  xfold=('cross_fold_jaccard', 'mean'))
           .reindex(ARMS))
    agg.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\t')
    log('')
    log('=' * 72)
    log('LODO 평균 (top %d)' % K)
    log('=' * 72)
    log('  %-14s %8s %6s %12s %12s %9s'
        % ('arm', 'AUROC', 'IEG', 'procurement', 'specificity', 'x-fold'))
    for a, r in agg.iterrows():
        log('  %-14s %8.3f %6.1f %11.0f%% %11.0f%% %9.3f'
            % (a, r['auroc'], r['ieg'], 100 * r['proc'], 100 * r['spec'], r['xfold']))
    lam = res[res['arm'] == 'resid_learned']['lam'].dropna()
    if len(lam):
        log('')
        log('  고른 lambda: %s' % ', '.join('%.2f' % v for v in lam))
    log('')
    log('  %s/ 에 저장' % OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
