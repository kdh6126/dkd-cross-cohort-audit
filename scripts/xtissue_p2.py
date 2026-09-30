#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""예측 P2: 즉시초기 모듈의 사례-대조 차이는 채취 비대칭을 따라가는가.

코호트마다 handling score(즉시초기 19유전자의 코호트 내 z 평균)를 표본별로 계산하고, 사례 대
대조의 Hedges' g 를 구한다. 신장 코호트도 같은 코드로 다시 계산해 세 조직을 한 표에 둔다.

g 하나만으로는 부족하다. 대장염처럼 전사체 전체가 크게 움직이는 질환에서는 아무 19유전자
묶음이나 큰 g 를 낼 수 있다. 그래서 같은 코호트에서 무작위 19유전자 묶음 2,000개의 |g| 분포를
만들고, 즉시초기 모듈의 |g| 가 그 분포의 몇 번째 백분위인지 함께 적는다. 모듈이 채취를
따라간다면 비대칭 코호트에서는 g 가 음수이면서 무작위 묶음보다 극단적이고, 대칭 코호트에서는
|g| 가 작아야 한다.

예측은 docs/xtissue/PREDICTIONS.md 에 데이터를 보기 전에 고정했다.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_cohort                     # noqa: E402
from dkd_deconfound import IEG                       # noqa: E402
from xtissue_config import COHORTS, arm_match        # noqa: E402

OUT = 'results/xtissue'
N_RANDOM = 2000

KIDNEY = {   # 신장 코호트는 원고의 표 그대로다. 사례는 모두 생검, 대조는 모두 비생검.
    'GSE30528': ('needle_biopsy', 'surgical_or_donor'),
    'GSE96804': ('needle_biopsy', 'surgical_or_donor'),
    'GSE104948': ('needle_biopsy', 'surgical_or_donor'),
    'GSE142025': ('needle_biopsy', 'surgical_or_donor'),
    'GSE294519': ('needle_biopsy', 'surgical_or_donor'),
    'GSE30529': ('needle_biopsy', 'surgical_or_donor'),
    'GSE104954': ('needle_biopsy', 'surgical_or_donor'),
}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def ieg_entrez():
    sym = {}
    with gzip.open('data/raw/annotation/Homo_sapiens.gene_info.gz', 'rt', encoding='utf-8',
                   errors='replace') as f:
        hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip('\n').split('\t')
            sym[p[ix['Symbol']]] = p[ix['GeneID']]
    return [sym[s] for s in IEG if s in sym]


def hedges(a, b):
    n1, n0 = len(a), len(b)
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / (n1 + n0 - 2))
    return (a.mean() - b.mean()) / (sp + 1e-12) * (1 - 3 / (4 * (n1 + n0) - 9))


def score(Z, cols):
    return Z[:, cols].mean(axis=1)


def one(tissue, cohort, root, proc, status, ieg_ids, rng):
    Z, y, _ = load_cohort(cohort, root=root)
    genes = [str(g) for g in Z.columns]
    pos = {g: i for i, g in enumerate(genes)}
    cols = [pos[g] for g in ieg_ids if g in pos]
    X = Z.values
    s = score(X, cols)
    g = hedges(s[y == 1], s[y == 0])
    boot = []
    for _ in range(1000):
        i1 = rng.choice(np.flatnonzero(y == 1), (y == 1).sum())
        i0 = rng.choice(np.flatnonzero(y == 0), (y == 0).sum())
        boot.append(hedges(s[i1], s[i0]))
    rnd = []
    for _ in range(N_RANDOM):
        c = rng.choice(X.shape[1], len(cols), replace=False)
        r = score(X, c)
        rnd.append(abs(hedges(r[y == 1], r[y == 0])))
    rnd = np.array(rnd)
    return dict(tissue=tissue, cohort=cohort, status=status,
                procurement_case=proc[0], procurement_control=proc[1], arms=arm_match(proc),
                n_case=int((y == 1).sum()), n_control=int((y == 0).sum()), n_ieg=len(cols),
                g_ieg=float(g), ci_lo=float(np.percentile(boot, 2.5)),
                ci_hi=float(np.percentile(boot, 97.5)),
                abs_g_random_median=float(np.median(rnd)),
                abs_g_percentile_vs_random=float(100 * (rnd < abs(g)).mean()))


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(20260914)
    ieg_ids = ieg_entrez()
    rows = []
    for c, proc in KIDNEY.items():
        rows.append(one('kidney_dkd', c, 'data/processed/harmonized', proc, 'manuscript',
                        ieg_ids, rng))
        log('%-10s %-11s g=%+.2f' % ('kidney', c, rows[-1]['g_ieg']))
    for tissue, cohorts in COHORTS.items():
        tab = pd.read_csv('results/xtissue/%s/cohorts.tsv' % tissue, sep='\t')
        ok = set(tab.loc[tab['n_case'].fillna(0) > 0, 'gse'])
        for c, cfg in cohorts.items():
            if c not in ok or cfg['label'] is None:
                continue
            if min(tab.loc[tab['gse'] == c, ['n_case', 'n_control']].values[0]) < 5:
                continue
            r = one(tissue, c, 'data/processed/xtissue/%s' % tissue, cfg['procurement'][:2],
                    cfg['status'].split(':')[0], ieg_ids, rng)
            rows.append(r)
            log('%-10s %-11s %-11s g=%+.2f [%+.2f, %+.2f]  |g| pct vs random=%.0f'
                % (tissue[:10], c, r['arms'], r['g_ieg'], r['ci_lo'], r['ci_hi'],
                   r['abs_g_percentile_vs_random']))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'p2_ieg_by_cohort.tsv'), sep='\t', index=False)

    # 예측 P2 판정
    log('\n=== P2 ===')
    verdict = []
    for tissue, g in df[df['tissue'] != 'kidney_dkd'].groupby('tissue'):
        m = g[g['arms'] == 'matched']
        a = g[g['arms'] == 'asymmetric']
        m_ok = bool((m['g_ieg'].abs() < 0.5).all()) if len(m) else None
        a_ok = (bool(((a['g_ieg'] < 0) & (a['g_ieg'].abs() > m['g_ieg'].abs().max())).all())
                if len(a) and len(m) else None)
        verdict.append(dict(tissue=tissue, n_matched=len(m), n_asymmetric=len(a),
                            matched_all_abs_g_below_0_5=m_ok,
                            asymmetric_negative_and_larger_than_every_matched=a_ok))
        log('%s: matched %d (all |g|<0.5: %s), asymmetric %d (negative and larger: %s)'
            % (tissue, len(m), m_ok, len(a), a_ok))
    pd.DataFrame(verdict).to_csv(os.path.join(OUT, 'p2_verdict.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
