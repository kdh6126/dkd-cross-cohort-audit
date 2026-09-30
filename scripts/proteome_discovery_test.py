#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""단백체 층 혼자서 후보를 발굴할 수 있는가 — 음성 결과.

proteome_meta.py 는 전사체에서 만든 후보 목록이 단백질 층에서 배경보다 자주 유의해진다는
것을 보였습니다(승산비 3.3). 그러면 자연히 다음 질문이 옵니다.

    단백체 자료만으로 후보 목록을 처음부터 만들 수는 없는가?

같은 두 코호트를 쓰고 전사체 정보를 전혀 넣지 않으면 됩니다. 두 코호트가 모두 잰
유전자 중에서 양쪽 다 q<0.05 이고 방향이 같은 것을 고르는 것이 가장 단순한 형태입니다.

결과는 음성입니다. 그리고 그 이유가 분명합니다.

    전사체 층    코호트 7개 · 환자 258명 · 유전자 9,900개
    단백체 층    코호트 2개 · 환자  45명 · 공통 유전자 170개

코호트가 둘이면 leave-one-dataset-out 이 "하나로 뽑아 하나로 채점" 이 되어 이 논문이
쓰는 절차가 성립하지 않습니다. 공통 유전자 170개는 교차 코호트 선택을 할 만한 공간이
아닙니다.

이 음성 결과를 남기는 이유는 두 가지입니다. 첫째, 보강과 발굴은 필요한 자료 규모가
다르다는 것을 수치로 보여 줍니다. 보강은 목록이 이미 있으므로 소수 유전자만 재도
되지만, 발굴은 목록을 만들어야 하므로 넓은 공간과 여러 코호트가 필요합니다. 둘째,
같은 자료로 두 가지를 다 했다고 주장하지 않는다는 것을 분명히 합니다.

검정은 주변분포를 지키는 순열로 합니다. 한쪽 코호트의 유전자 라벨만 섞으면 각 코호트의
유의 개수와 방향 분포는 그대로 두고 두 코호트 사이의 대응만 무너뜨릴 수 있습니다.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

SOMA = 'results/proteome_validation/protein_dkd_vs_control.tsv'
MS = 'results/proteome_ms/protein_ms_dkd_vs_control.tsv'
OUT = 'results/proteome_meta'
NPERM = 200000
SEED = 0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)

    s = pd.read_csv(SOMA, sep='\t')
    s = s[s['on_panel_and_measured']]
    m = pd.read_csv(MS, sep='\t')
    j = s[['gene', 'prot_g', 'prot_q', 'is_candidate']].merge(
        m[['gene', 'prot_g', 'prot_q']], on='gene', suffixes=('_a', '_b'))
    n = len(j)

    sa = (j['prot_q_a'] < 0.05).values
    sb = (j['prot_q_b'] < 0.05).values
    ga = j['prot_g_a'].values
    gb = j['prot_g_b'].values
    cc = np.sign(ga) == np.sign(gb)
    hit = sa & sb & cc
    obs = int(hit.sum())

    log('=' * 78)
    log('단백체 층 혼자서 발굴이 되는가')
    log('=' * 78)
    log('  두 코호트가 모두 잰 유전자      %d개' % n)
    log('  SOMAscan 에서 q<0.05            %d개' % sa.sum())
    log('  질량분석에서 q<0.05             %d개' % sb.sum())
    log('  방향이 일치                     %d개' % cc.sum())
    log('')
    log('  양쪽 유의 + 방향 일치           %d개' % obs)
    for _, r in j[hit].sort_values('prot_q_b').iterrows():
        log('    %-8s  g %+5.2f / %+5.2f   q %.3f / %.3f   %s'
            % (r['gene'], r['prot_g_a'], r['prot_g_b'], r['prot_q_a'], r['prot_q_b'],
               '전사체 후보' if r['is_candidate'] else '단백체 고유'))

    # 주변분포를 지키는 순열. 한쪽의 유전자 순서만 섞는다.
    idx = rng.permuted(np.tile(np.arange(n), (NPERM, 1)), axis=1)
    null = (sa & sb[idx] & (np.sign(ga)[None, :] == np.sign(gb)[idx])).sum(1)
    p = float((null >= obs).mean())
    exp_ind = n * sa.mean() * sb.mean() * 0.5
    p_binom = stats.binomtest(obs, n, sa.mean() * sb.mean() * 0.5,
                              alternative='greater').pvalue

    log('')
    log('  독립 가정 기대                  %.2f개' % exp_ind)
    log('  이항검정                        p = %.4f' % p_binom)
    log('  주변분포 보존 순열 %s회    p = %.4f  (귀무 평균 %.2f개)'
        % ('{:,}'.format(NPERM), p, null.mean()))
    log('')
    verdict = 'negative' if p >= 0.05 else 'positive'
    if verdict == 'negative':
        log('  음성입니다. 관측 %d개는 우연 기대와 구별되지 않습니다.' % obs)
        log('  단백체 층은 이 규모에서 후보를 발굴할 수 없습니다. 보강만 가능합니다.')
    else:
        log('  양성입니다. 서술을 다시 써야 합니다.')
    log('')
    log('  왜 그런가 — 두 층의 규모 차이')
    log('    전사체   코호트 7 · 환자 258 · 유전자 9,900')
    log('    단백체   코호트 2 · 환자  45 · 공통 유전자 %d' % n)

    j.assign(both_sig_concordant=hit).to_csv(
        os.path.join(OUT, 'discovery_attempt.tsv'), sep='\t', index=False)
    pd.DataFrame([dict(n_shared_genes=n, n_sig_somascan=int(sa.sum()),
                       n_sig_massspec=int(sb.sum()), n_concordant=int(cc.sum()),
                       n_both_sig_concordant=obs, expected_independent=exp_ind,
                       binom_p=p_binom, perm_p=p, perm_null_mean=float(null.mean()),
                       n_perm=NPERM, verdict=verdict)]).to_csv(
        os.path.join(OUT, 'discovery_summary.tsv'), sep='\t', index=False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
