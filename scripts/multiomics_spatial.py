#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""진짜 멀티오믹스로 할 수 있는 한 가지 — 같은 참여자, 같은 부위, 두 층.

이 프로젝트는 "같은 환자에서 두 층을 재고 통합한다"는 뜻의 멀티오믹스를 DKD 에서 하지
못했습니다. 공개 데이터에 없기 때문입니다. 그런데 KPMP 를 다시 뒤져보니 하나가 남아
있습니다.

    지역 전사체   사구체 vs 세뇨관간질,  전체 참여자 풀,  유전자 26,485개
    지역 단백체   사구체 vs 세뇨관간질,  전체 참여자 풀,  단백질

두 층이 같은 조직 부위를, 같은 참여자 집단에서, 같은 대조로 쟀습니다. 이것은 정의상
멀티오믹스입니다. 다만 대조가 '질병 대 정상' 이 아니라 '부위 대 부위' 입니다.
그래서 이것으로 질병 바이오마커를 새로 만들 수는 없고, 대신 이렇게 씁니다.

    후보가 지목한 부위가 단백질 수준에서도 같은 부위인가.

전사체에서 사구체에 많다고 나온 유전자가 단백질에서도 사구체에 많으면, 그 후보의 부위
귀속은 두 층에서 일치하는 것입니다. 어긋나면 전사체만의 현상일 수 있습니다.

한계를 먼저 적습니다. 참여자 단위로 잇는 것이 아니라 집단 요약끼리 잇습니다. KPMP 공개
API 가 참여자별 값을 주지 않기 때문입니다. 그리고 단백체에는 질병별 층화가 없어
'DKD 에서 오르는가' 는 이 축으로 물을 수 없습니다.
"""
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kpmp_client import regional_proteomics   # noqa: E402
from kpmp_client import regional_transcriptomics_by_structure as rt_by_structure  # noqa: E402

OUT = 'results/multiomics_spatial'
CAND = 'results/candidates_v2/master_candidate_table.tsv'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)

    # ---------------------------------------------------------- 전사체 (한 번에)
    log('지역 전사체(사구체) 내려받는 중 ...')
    rt = pd.DataFrame(rt_by_structure('Glomerulus'))
    rt = rt[['geneSymbol', 'foldChange', 'adjPVal', 'sampleCount']].copy()
    rt.columns = ['gene', 'rt_fc', 'rt_q', 'rt_n']
    rt = rt.dropna(subset=['rt_fc']).drop_duplicates('gene')
    log('  유전자 %d개, 표본 %s' % (len(rt), rt['rt_n'].dropna().unique()[:3]))

    # ---------------------------------------------------------- 단백체 (후보만)
    cand = pd.read_csv(CAND, sep='\t')
    genes = list(cand[cand.columns[0]])
    log('후보 %d개의 단백체 조회 ...' % len(genes))
    rows = []
    for k, g in enumerate(genes, 1):
        try:
            r = regional_proteomics(g) or []
        except Exception:
            r = []
        glom = [x for x in r if str(x.get('region', '')).lower().startswith('glom')]
        if glom:
            x = glom[0]
            rows.append(dict(gene=g, rp_fc=x.get('foldChange'), rp_q=x.get('adjPVal'),
                             rp_n=x.get('sampleCount'), rp_peptides=x.get('numUniquePeptides')))
        if k % 10 == 0:
            log('  %d / %d' % (k, len(genes)))
        time.sleep(0.15)
    rp = pd.DataFrame(rows)
    log('  단백질이 검출된 후보 %d / %d' % (len(rp), len(genes)))

    # ---------------------------------------------------------- 합치기
    d = rp.merge(rt, on='gene', how='left')
    d['rt_dir'] = np.sign(d['rt_fc'])
    d['rp_dir'] = np.sign(d['rp_fc'])
    d['both_measured'] = d['rt_fc'].notna() & d['rp_fc'].notna()
    d['concordant'] = d['both_measured'] & (d['rt_dir'] == d['rp_dir'])
    d['both_significant'] = (d['rt_q'] < 0.05) & (d['rp_q'] < 0.05)
    d = d.merge(cand[[cand.columns[0], 'tier']].rename(columns={cand.columns[0]: 'gene'}),
                on='gene', how='left')
    d = d.sort_values(['concordant', 'rp_q'], ascending=[False, True])
    d.to_csv(os.path.join(OUT, 'candidates_two_layer.tsv'), sep='\t', index=False)

    # ---------------------------------------------------------- 귀무 비교
    # 후보가 특별히 일치하는 것인지, 아무 유전자나 그런지 본다.
    both = d[d['both_measured']]
    obs = float(both['concordant'].mean()) if len(both) else float('nan')
    log('')
    log('=' * 74)
    log('두 층에서 모두 측정된 후보 %d개 · 방향 일치 %d개 (%.0f%%)'
        % (len(both), int(both['concordant'].sum()), 100 * obs))

    # 배경: 전사체 전체에서 무작위로 같은 수를 뽑아 부호가 일치할 확률은 0.5 다.
    # 단백체를 전체 유전자에 대해 받을 수 없으므로, 이항검정으로 대신한다.
    from scipy import stats
    if len(both):
        p = stats.binomtest(int(both['concordant'].sum()), len(both), 0.5,
                            alternative='greater').pvalue
        log('  우연(50%%) 대비 이항검정 p = %.4g' % p)
    log('')
    sig = d[d['both_significant'] & d['concordant']]
    log('두 층 모두 유의(q<0.05)하고 방향이 같은 후보: %d개' % len(sig))
    if len(sig):
        log('  %-10s %-26s %9s %9s' % ('유전자', '등급', '전사체 FC', '단백체 FC'))
        for _, r in sig.iterrows():
            log('  %-10s %-26s %9.2f %9.2f'
                % (r['gene'], str(r['tier'])[:24], r['rt_fc'], r['rp_fc']))
    log('')
    log('  %s 에 저장했습니다.' % os.path.join(OUT, 'candidates_two_layer.tsv'))
    log('')
    log('  주의 — 이 축의 대조는 "사구체 대 세뇨관간질" 이지 "질병 대 정상" 이 아닙니다.')
    log('  따라서 여기서 일치한다는 것은 부위 귀속이 두 층에서 같다는 뜻이고,')
    log('  그 유전자가 DKD 바이오마커라는 뜻이 아닙니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
