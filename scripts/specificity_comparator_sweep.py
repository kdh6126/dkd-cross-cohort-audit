#!/usr/bin/env python
"""우리 논문의 특이성 게이트도 대조 질병을 묶은 것이 결과를 만든 것 아닌가.

ST003255 대사체에서, "다른 신장병"을 하나로 묶으면 29개가 통과하지만 질병별로 따로 보면
IgA신병증 16개 · 막성신병증 88개 · 고혈압신병증 2개로 답이 완전히 달라졌습니다. 통합이
결과를 만든 것입니다.

우리 논문의 특이성 게이트도 정확히 같은 구조입니다. ERCB에서 DN을 제외한 생검 진단 9개를
하나의 "other CKD"로 묶어 DKD와 비교했습니다. 같은 문제가 있다면 논문의 게이트도 흔들립니다.

    실험 1  진단별로 따로 게이트를 통과시켜 본다. 후보 30개가 9개 중 몇 개에서 통과하는가.
    실험 2  대조 진단 부분집합을 전수 조사해, 통과 유전자 수가 얼마나 흔들리는지 본다.
            코호트 구성 민감도 분석(3.2절)의 대조군 판(版)이다.
    실험 3  통합 게이트가 진단별 결과의 평균에 가까운가, 아니면 특정 진단이 끌고 가는가.

논문의 주장을 우리가 먼저 공격하는 것이 목적입니다.
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_cohort   # noqa: E402

ERCB = ['GSE104948', 'GSE104954']
OUT = 'results/specificity_sweep'
MIN_G = 0.5           # final_candidates.py와 동일한 문턱
MIN_N = 10            # 이보다 작은 진단군은 단독 비교에서 제외 (검정력 부족)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def hedges_g(a, b):
    n1, n0 = len(a), len(b)
    if n1 < 2 or n0 < 2:
        return np.full(a.shape[1], np.nan)
    sp = np.sqrt(((n1 - 1) * a.var(0, ddof=1) + (n0 - 1) * b.var(0, ddof=1)) / (n1 + n0 - 2))
    sp = np.where(sp == 0, np.nan, sp)
    return (a.mean(0) - b.mean(0)) / sp * (1 - 3 / (4 * (n1 + n0) - 9))


def norm(s):
    """두 코호트가 같은 진단을 다르게 표기한다. 소문자·공백 정규화로 맞춘다.

    밑줄도 공백으로 바꾼다. 이것을 빠뜨려 'living_donor'가 제외 목록의 'living donor'와
    매칭되지 않았고, 비생검 대조군이 '다른 신장병'으로 섞여 들어간 적이 있다.
    """
    return ' '.join(str(s).lower().replace('/', ' ').replace('-', ' ')
                    .replace('_', ' ').split())


def load():
    """ERCB 두 구획을 진단 라벨과 함께 읽는다."""
    frames = []
    for c in ERCB:
        X, y, p = load_cohort(c, labelled_only=False)
        X.columns = X.columns.astype(str)
        p = p.reset_index(drop=True)
        col = 'group' if 'group' in p.columns else 'stage'
        frames.append((c, X, p, col))
    return frames


def main():
    os.makedirs(OUT, exist_ok=True)
    frames = load()
    genes = frames[0][1].columns

    # 진단명 정규화 후, 두 구획에 공통으로 있는 비-DN 생검 진단만 쓴다
    per_cohort_dx = []
    for c, X, p, col in frames:
        dx = {norm(v): v for v in p[col].dropna().unique()}
        per_cohort_dx.append(dx)
    common = set(per_cohort_dx[0]) & set(per_cohort_dx[1])
    drop = {'dn', 'living donor', 'tumor nephrectomy', 'donor', 'nephrectomy', 'control'}
    comparators = sorted(common - drop)
    # 이름에 의존하지 않는 안전장치: 생검으로 얻은 군만 대조 질병이 될 수 있다
    biopsy = set()
    for (c, X, p, col), dxmap in zip(frames, per_cohort_dx):
        if 'control_type' in p.columns:
            for dxn, raw in dxmap.items():
                sel = p[col] == raw
                if sel.any() and p.loc[sel, 'control_type'].isna().all():
                    biopsy.add(dxn)
    if biopsy:
        removed = [d for d in comparators if d not in biopsy]
        if removed:
            log('  (비생검 군 제외: %s)' % ', '.join(removed))
        comparators = [d for d in comparators if d in biopsy]
    log('두 ERCB 구획에 공통인 비-DN 생검 진단 %d개' % len(comparators))

    # ---------------------------------------------------- 진단별 g
    gmat, sizes = {}, {}
    for dxn in comparators:
        gs, ns = [], 0
        for (c, X, p, col), dxmap in zip(frames, per_cohort_dx):
            M = X.values
            dkd = M[(p['label'] == 1).values]
            oth = M[(p[col] == dxmap[dxn]).values]
            ns += len(oth)
            if len(oth) >= 2:
                gs.append(hedges_g(dkd, oth))
        if gs:
            gmat[dxn] = np.nanmean(np.vstack(gs), axis=0)
            sizes[dxn] = ns
    used = [d for d in comparators if sizes.get(d, 0) >= MIN_N]
    log('그중 n>=%d 이라 단독 비교가 가능한 것 %d개' % (MIN_N, len(used)))
    log('')

    G = pd.DataFrame(gmat, index=genes)
    G.to_csv(os.path.join(OUT, 'g_by_diagnosis.tsv'), sep='\t')

    # ---------------------------------------------------- 실험 1
    log('=' * 76)
    log('실험 1 — 진단별로 따로 보면 몇 개가 통과하는가')
    log('=' * 76)
    log('  %-46s %5s %8s' % ('대조 진단', 'n', '통과 유전자'))
    rows = []
    for d in used:
        n_pass = int((G[d].abs() >= MIN_G).sum())
        rows.append(dict(diagnosis=d, n=sizes[d], n_pass=n_pass))
        log('  %-46s %5d %8d' % (d[:44], sizes[d], n_pass))
    pooled = pd.read_csv('results/final_candidates.tsv', sep='\t')
    pooled_pass = int(pooled['passes_specificity'].sum())
    log('  %-46s %5s %8d' % ('(통합: 9개 진단을 하나로)', '-', pooled_pass))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'per_diagnosis.tsv'), sep='\t', index=False)

    # ---------------------------------------------------- 실험 2
    log('')
    log('=' * 76)
    log('실험 2 — 대조 진단을 어떻게 고르느냐에 따라 얼마나 흔들리는가')
    log('=' * 76)
    sub = []
    for k in range(1, len(used) + 1):
        counts = []
        combos = list(itertools.combinations(used, k))
        if len(combos) > 60:
            rng = np.random.default_rng(0)
            idx = rng.choice(len(combos), 60, replace=False)
            combos = [combos[i] for i in idx]
        for cb in combos:
            gm = G[list(cb)].mean(axis=1)
            counts.append(int((gm.abs() >= MIN_G).sum()))
        sub.append(dict(k=k, n_subsets=len(combos), mean=np.mean(counts),
                        lo=int(np.min(counts)), hi=int(np.max(counts))))
        log('  대조 진단 %d개  조합 %3d개  통과 유전자 %5.0f  (범위 %d ~ %d)'
            % (k, len(combos), np.mean(counts), np.min(counts), np.max(counts)))
    pd.DataFrame(sub).to_csv(os.path.join(OUT, 'subset_sweep.tsv'), sep='\t', index=False)

    # ---------------------------------------------------- 실험 3
    log('')
    log('=' * 76)
    log('실험 3 — 후보 30개는 진단별로 얼마나 일관되게 통과하는가')
    log('=' * 76)
    dos = pd.read_csv('results/candidates_v2/master_candidate_table.tsv', sep='\t')
    gs_tab = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = {s: g for g, s in zip(gs_tab['entrez_id'], gs_tab['symbol'])}
    rows = []
    for sym in dos[dos.columns[0]]:
        eid = sym2id.get(sym)
        if eid is None or eid not in G.index:
            continue
        v = G.loc[eid, used].abs()
        rows.append(dict(gene=sym, n_pass=int((v >= MIN_G).sum()), n_dx=len(used),
                         median_g=float(G.loc[eid, used].median())))
    t = pd.DataFrame(rows).sort_values('n_pass', ascending=False)
    t.to_csv(os.path.join(OUT, 'candidate_consistency.tsv'), sep='\t', index=False)
    log('  후보 %d개 중 진단 %d개 전부에서 통과: %d개'
        % (len(t), len(used), int((t['n_pass'] == len(used)).sum())))
    log('  과반(%d개 이상)에서 통과            : %d개'
        % (len(used) // 2 + 1, int((t['n_pass'] >= len(used) // 2 + 1).sum())))
    log('  하나도 통과 못함                    : %d개' % int((t['n_pass'] == 0).sum()))
    log('')
    log('  상위 10개:')
    for _, r in t.head(10).iterrows():
        log('    %-10s %d/%d 진단에서 통과  (중앙 g = %+.2f)'
            % (r['gene'], r['n_pass'], r['n_dx'], r['median_g']))

    # ---------------------------------------------------- 판정
    log('')
    log('=' * 76)
    log('판정')
    log('=' * 76)
    counts = [int((G[d].abs() >= MIN_G).sum()) for d in used]
    spread, lo = max(counts), min(counts)
    log('  진단 하나만 대조로 쓰면 통과 유전자가 %d ~ %d개로 흔들립니다 (%.1f배).'
        % (lo, spread, spread / max(lo, 1)))
    frac_all = (t['n_pass'] == len(used)).mean()
    frac_half = (t['n_pass'] >= len(used) // 2 + 1).mean()
    log('  후보 30개 중 %.0f%%가 진단 전부에서, %.0f%%가 과반에서 통과합니다.'
        % (100 * frac_all, 100 * frac_half))
    log('')
    if frac_half >= 0.6:
        log('  -> 대사체와 달리, 전사체 게이트는 대조 진단 선택에 상당히 견고합니다.')
        log('     통합이 결과를 만들었다고 보기 어렵습니다.')
    else:
        log('  -> 전사체 게이트도 대조 진단 선택에 흔들립니다. 논문의 특이성 주장을')
        log('     "9개 진단을 통합했을 때"로 한정해 서술해야 합니다.')


if __name__ == '__main__':
    main()
