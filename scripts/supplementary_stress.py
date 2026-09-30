#!/usr/bin/env python
"""보조 스트레스 테스트 — 표본이 작아 본문급이 될 수 없는 데이터들.

미리 정한 기준: 이 둘은 보조자료 후보로 시작합니다. 결과가 아주 강하면 올리고,
아니면 그대로 보조로 둡니다. 좋은 결과만 고르지 않기 위해 기준을 먼저 적어둡니다.

    GSE1009      사구체 6샘플. 샘플명이 Control 1a/1b/2, Diabetes 1a/1b/2 입니다.
                 'a'와 'b'는 같은 환자의 기술 반복으로 보이므로, 실질 환자 수는 군당 2명
                 입니다. 이것을 6개로 세면 유사반복(pseudo-replication)이 됩니다.
                 그래서 두 가지로 계산해 둘 다 보고합니다.

    GSE131882    단일세포 3 vs 3. 원시 zUMIs .rds 파일만 등록되어 있어 파이썬에서 읽히지
                 않고, 이 환경에 R이 없습니다. 게다가 같은 질문에 대해 KPMP가 훨씬 큰
                 자료(DKD 27 vs 다른 신장병 10)를 이미 제공합니다. 처리 비용 대비 얻을
                 것이 없어 수행하지 않고, 그 판단을 기록으로 남깁니다.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from candidate_external_validation import read_series, collapse, hedges_g, enrich_vs_random

OUT = 'results/external_validation'
SEED = 0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)
    cand = pd.read_csv('results/candidates_v2/master_candidate_table.tsv', sep='\t')
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = dict(zip(gs['symbol'], gs['entrez_id']))
    cand_ids = {sym2id[s] for s in cand[cand.columns[0]] if s in sym2id}

    rows = []
    log('=' * 76)
    log('보조 1 — GSE1009 (사구체 6샘플, 오래된 플랫폼)')
    log('=' * 76)
    X, meta = read_series('GSE1009')
    Y = collapse(X, 'GPL8300')
    if Y is None:
        log('  GPL8300 매핑 실패. 건너뜁니다.')
    else:
        title = pd.Series(meta['!Sample_title'][0], index=X.columns)
        dm = title.str.contains('Diabetes', case=False)
        ct = title.str.contains('Control', case=False)
        log('  당뇨 %d · 대조 %d · 유전자 %d' % (dm.sum(), ct.sum(), len(Y)))
        log('  샘플명: %s' % ', '.join(title.values))

        # (1) 6개를 독립 표본으로 취급 — 유사반복. 참고용으로만.
        g6 = pd.Series({gid: hedges_g(Y.loc[gid, dm.values].values,
                                      Y.loc[gid, ct.values].values)
                        for gid in Y.index}).dropna()
        # (2) 기술 반복을 환자별로 평균 — 군당 2명. 이쪽이 정직한 계산.
        pid = title.str.extract(r'(\d+)')[0].values
        grp = np.where(dm.values, 'DM', 'CT')
        key = pd.Series([g + p for g, p in zip(grp, pid)], index=X.columns)
        Yp = Y.T.groupby(key.values).mean().T
        dm2 = pd.Series(Yp.columns).str.startswith('DM').values
        log('  기술 반복 병합 후: 당뇨 %d명 · 대조 %d명' % (dm2.sum(), (~dm2).sum()))
        g2 = pd.Series({gid: hedges_g(Yp.loc[gid, dm2].values, Yp.loc[gid, ~dm2].values)
                        for gid in Yp.index}).dropna()

        for label, gv in (('샘플 6개 (유사반복)', g6), ('환자 2+2 (정직한 계산)', g2)):
            cid = [i for i, gid in enumerate(gv.index) if gid in cand_ids]
            if len(cid) < 5:
                log('  %s: 매핑된 후보 %d개뿐이라 검정 불가' % (label, len(cid)))
                continue
            obs, nul, p = enrich_vs_random(gv.values, cid, rng, n_draw=2000)
            log('  %-22s 후보 |g| 중앙 %.3f  무작위 %.3f  순열 p = %.3f'
                % (label, obs, nul, p))
            rows.append(dict(dataset='GSE1009', variant=label, n_cand=len(cid),
                             median_g=obs, null_g=nul, p=p))
        log('')
        log('  주의 — 환자 2+2에서는 효과크기 자체가 거의 의미가 없습니다. 여기서 무엇이')
        log('  나오든 근거로 쓰지 않고, 유사반복이 결과를 얼마나 부풀리는지 보는 용도로만')
        log('  씁니다.')

    log('')
    log('=' * 76)
    log('보조 2 — GSE131882 (단일세포 3 vs 3)')
    log('=' * 76)
    log('  수행하지 않았습니다. 이유 셋:')
    log('    1. 원시 zUMIs .rds 파일만 등록되어 있고, pyreadr가 이 형식을 읽지 못합니다.')
    log('       ("The file contains an unrecognized object")')
    log('    2. 이 환경에 R이 없어 변환할 수단이 없습니다.')
    log('    3. 같은 질문에 KPMP가 훨씬 큰 자료를 제공합니다 — DKD 27명 vs 다른 신장병')
    log('       10명 vs 대조 18명, 이미 파이프라인에 들어가 있습니다. 3 vs 3을 추가로')
    log('       처리해 얻을 것이 없습니다.')
    log('  판단: 처리 비용 대비 정보 이득이 없어 제외. 데이터 자체는 보관합니다.')
    rows.append(dict(dataset='GSE131882', variant='수행 안 함', n_cand=0,
                     median_g=np.nan, null_g=np.nan, p=np.nan))

    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(OUT, 'supplementary.tsv'),
                                  sep='\t', index=False)
        log('')
        log('  wrote %s/supplementary.tsv' % OUT)


if __name__ == '__main__':
    main()
