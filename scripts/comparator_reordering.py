#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""대조군의 조달 방식을 바꾸면 후보 순위가 얼마나 재배열되는가.

3.6절의 표는 대조를 맞췄을 때 효과크기가 얼마나 줄어드는지를 중앙 |g| 비율로 보여줍니다.
줄어든다는 것과 **다른 유전자가 올라온다**는 것은 다른 말입니다. 신호가 고르게 줄면
순위는 그대로일 수 있고, 그러면 후보 목록은 바뀌지 않습니다.

그래서 여기서는 순위를 봅니다. 같은 사례 집합을 두 대조와 각각 비교하고, 각 비교에서
상위 K개를 뽑아 겹침을 셉니다.

    g_naive     사례 대 비생검 대조(신절제·공여자).  조달 방식이 사례와 다르다.
    g_matched   사례 대 다른 진단의 생검.            조달 방식이 사례와 같다.

한 가지를 분명히 해 둡니다. 이것은 **같은 사례 집합을 두 대조군과 각각 비교한 것**이지,
같은 환자에게서 두 가지 방식으로 조직을 채취해 짝지은 비교가 아닙니다. 공개 자료에
그런 짝은 없습니다.

전체 순위상관도 함께 냅니다. 상관이 높은데 상위 겹침이 낮으면, 대량의 유전자는 같은
방향으로 움직이되 **꼭대기만 갈린다**는 뜻입니다. 후보를 고르는 것은 꼭대기이므로
그 구분이 중요합니다.
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

SRC = 'results/ckd_generalization'
GSPACE = 'data/processed/gene_space_all6.tsv'
CAND_V1 = 'results/master_candidate_table.tsv'
CAND_V2 = 'results/candidates_v2/master_candidate_table.tsv'
OUT = 'results/comparator_reordering'
KS = (20, 50, 100)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)
    gs = pd.read_csv(GSPACE, sep='\t', dtype=str)
    id2sym = dict(zip(gs['entrez_id'], gs['symbol']))

    rows = []
    detail = {}
    for f in sorted(glob.glob(os.path.join(SRC, 'g_*.tsv'))):
        name = os.path.basename(f)[2:-4].replace('_', ' ')
        d = pd.read_csv(f, sep='\t', dtype={'entrez_id': str})
        ok = d['g_naive'].notna() & d['g_matched'].notna()
        d = d[ok]
        r = dict(diagnosis=name, n_genes=len(d),
                 spearman=float(stats.spearmanr(d['g_naive'], d['g_matched']).statistic))
        for K in KS:
            a = set(d.reindex(d['g_naive'].abs().sort_values(ascending=False).index)
                    .head(K)['entrez_id'])
            b = set(d.reindex(d['g_matched'].abs().sort_values(ascending=False).index)
                    .head(K)['entrez_id'])
            r['shared_top%d' % K] = len(a & b)
            r['jaccard_top%d' % K] = len(a & b) / len(a | b)
            if K == 50:
                detail[name] = (a, b)
        rows.append(r)

    t = pd.DataFrame(rows).sort_values('shared_top50')
    t.to_csv(os.path.join(OUT, 'top_overlap.tsv'), sep='\t', index=False)

    log('=' * 78)
    log('대조군 조달 방식을 바꿨을 때 상위 유전자가 얼마나 갈리는가')
    log('=' * 78)
    log('  같은 사례 집합을 두 대조군과 각각 비교한 것입니다. 같은 환자에서 두 방식으로')
    log('  채취해 짝지은 비교가 아닙니다.')
    log('')
    log('  %-32s %8s %8s %9s %9s'
        % ('진단', '상위20', '상위50', '상위100', '순위상관'))
    for _, r in t.iterrows():
        log('  %-32s %8d %8d %9d %9.2f'
            % (r['diagnosis'][:32], r['shared_top20'], r['shared_top50'],
               r['shared_top100'], r['spearman']))

    v1 = set(pd.read_csv(CAND_V1, sep='\t')['gene'])
    v2 = set(pd.read_csv(CAND_V2, sep='\t')['gene'])
    dn = t[t['diagnosis'] == 'dn'].iloc[0]
    a, b = detail['dn']
    log('')
    log('  당뇨병성 신장질환: 상위 50개 중 %d개만 겹칩니다.' % dn['shared_top50'])
    log('    비생검 대조에서만 : %s'
        % ', '.join(sorted(id2sym.get(x, x) for x in a - b)[:10]))
    log('    생검 대조에서만   : %s'
        % ', '.join(sorted(id2sym.get(x, x) for x in b - a)[:10]))
    log('    양쪽 공통         : %s'
        % (', '.join(sorted(id2sym.get(x, x) for x in a & b)) or '(없음)'))

    # 알파벳 순으로 잘라 보여주면 어느 쪽이 오르고 내리는지가 보이지 않는다. 후보
    # 유전자의 순위가 대조를 바꾸면 어디로 가는지를 그대로 적는 편이 훨씬 분명하다.
    dd = pd.read_csv(os.path.join(SRC, 'g_dn.tsv'), sep='\t',
                     dtype={'entrez_id': str}).dropna()
    dd['gene'] = dd['entrez_id'].map(id2sym)
    dd['rank_naive'] = dd['g_naive'].abs().rank(ascending=False)
    dd['rank_matched'] = dd['g_matched'].abs().rank(ascending=False)
    watch = dd[dd['gene'].isin(v1 | v2)].copy()
    watch['moved'] = watch['rank_naive'] - watch['rank_matched']
    watch = watch.sort_values('moved', ascending=False)
    watch[['gene', 'g_naive', 'rank_naive', 'g_matched', 'rank_matched', 'moved']].to_csv(
        os.path.join(OUT, 'candidate_rank_shift.tsv'), sep='\t', index=False)
    log('')
    log('  후보 유전자의 순위가 대조를 바꾸면 어디로 가는가 (당뇨병성 신장질환)')
    log('    %-9s %11s %10s' % ('유전자', '비생검대조', '생검대조'))
    for _, r in pd.concat([watch.head(5), watch.tail(5)]).iterrows():
        log('    %-9s %11.0f %10.0f   %s'
            % (r['gene'], r['rank_naive'], r['rank_matched'],
               '올라감' if r['moved'] > 0 else '내려감'))
    log('    오르는 것과 내려가는 것이 갈립니다. 신호가 고르게 줄어드는 것이 아닙니다.')

    log('')
    log('  최종 후보 목록: 보정 전 %d개, 보정 후 %d개, 공통 %d개'
        % (len(v1), len(v2), len(v1 & v2)))
    log('    보정 후에만 나타남: %s' % ', '.join(sorted(v2 - v1)))

    pd.DataFrame([dict(
        n_diagnoses=len(t), dn_shared_top50=int(dn['shared_top50']),
        dn_jaccard_top50=float(dn['jaccard_top50']),
        dn_spearman=float(dn['spearman']),
        min_shared_top50=int(t['shared_top50'].min()),
        max_shared_top50=int(t['shared_top50'].max()),
        median_shared_top50=float(t['shared_top50'].median()),
        candidates_v1=len(v1), candidates_v2=len(v2), candidates_shared=len(v1 & v2))]
    ).to_csv(os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)
    log('')
    log('  순위상관은 높은데 상위 겹침이 낮으면, 대량의 유전자는 같이 움직이되')
    log('  꼭대기만 갈린다는 뜻입니다. 후보를 고르는 것은 꼭대기입니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
