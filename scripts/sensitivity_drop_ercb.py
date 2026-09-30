#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""GSE104948 을 발굴에서 빼면 특이성 결론이 유지되는가.

왜 필요한가
-----------
특이성 게이트는 ERCB 의 두 블록(GSE104948 사구체, GSE104954 세뇨관간질)에서 DKD 를
비당뇨 CKD 와 견주어 만듭니다. 그런데 GSE104948 은 발굴 코호트이기도 합니다. 같은 환자의
같은 잡음이 후보를 뽑는 데도, 그 후보가 게이트를 통과하는 데도 쓰일 수 있으므로 게이트는
완전히 독립인 검증이 아닙니다. 치명적 누수는 아니지만 부분적 중첩이고, 적어 두지 않으면
심사에서 먼저 지적당합니다.

그래서 한 번 떼어 봅니다. GSE104948 을 **발굴에서만** 빼고 나머지 코호트로 후보를 다시
만들면, 게이트는 발굴에 쓰이지 않은 자료로만 계산된 것이 됩니다. 그 상태에서

    게이트 통과율        게이트 적용 **전** 상위 N 개 중 몇 개가 통과하는가.
                         게이트를 통과한 표에서 재면 정의상 100% 라 아무 정보가 없다.
                         기준선은 무작위 유전자의 8.3% 다.
    GSE20602 외부 검사   ERCB 밖에서도 후보가 신경화증보다 DKD 쪽으로 더 움직이는가

이 둘이 유지되면 중첩이 결론을 만들어 낸 것이 아니라는 뜻입니다.

후보 목록이 그대로일 필요는 없습니다. 코호트를 하나 빼면 순위는 당연히 흔들립니다.
확인하는 것은 핵심 결론의 방향과 주요 후보군이 남는가입니다.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import candidate_external_validation as EV

TAB = {'main': 'results/candidates_v2/%s.tsv',
       'sens': 'results/sensitivity_no104948/%s.tsv'}
OUT = 'results/sensitivity_no104948'
TOPN = (50, 100, 200)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load(arm, which):
    d = pd.read_csv(TAB[arm] % which, sep='\t', dtype={'entrez_id': str})
    return d.set_index(d['entrez_id'].astype(str))


def gse20602(cand_ids, g_dkd, label):
    """후보가 신경화증보다 DKD 쪽으로 더 크게 움직이는가 (ERCB 밖)."""
    X, meta = EV.read_series('GSE20602')
    Y = EV.collapse(X, 'GPL96')
    src = EV.sample_annot(meta, X.columns, key='source_name')
    nsc = src.str.contains('NSC', case=False)
    ctl = src.str.contains('Nephrectomy', case=False)
    g_nsc = pd.Series({gid: EV.hedges_g(Y.loc[gid, nsc.values].values,
                                        Y.loc[gid, ctl.values].values)
                       for gid in Y.index}).dropna()
    common = [g for g in g_nsc.index if g in g_dkd.index]
    cid = [i for i, g in enumerate(common) if g in cand_ids]
    ratio = np.abs(g_dkd.loc[common].values) - np.abs(g_nsc.loc[common].values)
    tt = stats.mannwhitneyu(ratio[cid], np.delete(ratio, cid), alternative='greater')
    log('    %-24s 공통 %d개 중 후보 %d개   후보 중앙값 %+.2f   나머지 %+.2f   MWU p = %.1e'
        % (label, len(common), len(cid), np.nanmedian(ratio[cid]),
           np.nanmedian(np.delete(ratio, cid)), tt.pvalue))
    return dict(n_common=len(common), n_cand=len(cid),
                median_cand=float(np.nanmedian(ratio[cid])),
                median_rest=float(np.nanmedian(np.delete(ratio, cid))),
                mwu_p=float(tt.pvalue))


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    if not os.path.exists(TAB['sens'] % 'candidates_gated'):
        log('민감도 후보표가 없습니다. rederive_candidates.py 를 축소 코호트로 먼저 도세요.')
        return 1
    os.makedirs(OUT, exist_ok=True)

    gated = {a: load(a, 'candidates_gated') for a in TAB}
    allg = {a: load(a, 'candidates_all') for a in TAB}
    ms = set(gated['main']['symbol'])
    ss = set(gated['sens']['symbol'])

    log('=' * 78)
    log('GSE104948 을 발굴에서 뺀 민감도 분석')
    log('=' * 78)
    log('  최종 후보  본 분석 %d개, 민감도 %d개, 겹침 %d개 (Jaccard %.2f)'
        % (len(ms), len(ss), len(ms & ss), len(ms & ss) / len(ms | ss)))
    log('  두 목록에 모두 있는 유전자 (%d개):' % len(ms & ss))
    log('    %s' % ', '.join(sorted(ms & ss)))
    log('  민감도에서만: %s' % ', '.join(sorted(ss - ms)))
    log('')
    log('  게이트 통과율 — 게이트를 적용하기 **전** 상위 N 개 기준 (무작위 배경 8.3%)')
    for n in TOPN:
        log('    상위 %-4d   본 분석 %3.0f%%   민감도 %3.0f%%'
            % (n, 100 * allg['main'].head(n)['passes_specificity'].mean(),
               100 * allg['sens'].head(n)['passes_specificity'].mean()))
    log('  민감도 쪽 게이트는 발굴에 쓰이지 않은 자료로만 계산된 것입니다.')
    log('')
    log('  최종 후보의 |g| (DKD 대 비당뇨 CKD) 중앙값   본 분석 %.2f   민감도 %.2f'
        % (gated['main']['g_DKD_vs_otherCKD'].abs().median(),
           gated['sens']['g_DKD_vs_otherCKD'].abs().median()))
    log('')
    log('  ERCB 밖 독립 검증 (GSE20602, 신경화증 13 대 종양신절제 5)')

    rows = []
    for arm, label in (('main', 'main (GSE104948 포함)'), ('sens', 'sensitivity (제외)')):
        r = gse20602(set(gated[arm]['entrez_id'].astype(str)),
                     allg[arm]['g_DKD_vs_control'], label)
        r['arm'] = label
        for n in TOPN:
            r['gate_pass_top%d' % n] = float(allg[arm].head(n)['passes_specificity'].mean())
        r['n_final'] = len(gated[arm])
        r['n_overlap'] = len(ms & ss)
        rows.append(r)
    cols = ['arm', 'n_final', 'n_overlap'] + ['gate_pass_top%d' % n for n in TOPN] + \
           ['n_common', 'n_cand', 'median_cand', 'median_rest', 'mwu_p']
    pd.DataFrame(rows)[cols].to_csv(os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)
    log('')
    log('  %s 에 저장했습니다.' % os.path.join(OUT, 'summary.tsv'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
