#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""문헌 근거 장부 — 무엇을 언제 어떻게 물어서 그 등급이 나왔는가.

CrossRef 검증은 참고문헌이 실재하는지만 봅니다. "이 후보는 DKD 문헌에 희소하다" 는
주장 자체는 검증하지 않습니다. 그 주장은 PubMed 질의 결과에서 나오므로, 질의와 답을
그대로 남겨야 확인할 수 있습니다.

남기는 것은 넷입니다.

    질의       유전자마다 실제로 보낸 문자열 네 개를 그대로
    검색일     동결본에는 날짜가 없습니다. 그래서 지금 다시 물어 오늘 날짜와 함께 적습니다
    답         건수와 PMID
    제외 사유  기호가 다른 뜻과 충돌하면 건수를 믿을 수 없다고 적어 둔 이유

다시 묻는 이유가 하나 더 있습니다. PubMed 는 시간이 지나면 늘어납니다. 원고가 보고하는
것은 동결본이므로, 지금 값과 다르면 그 차이를 드러내야 합니다. 등급이 뒤집히는 유전자가
있으면 그것은 원고가 고쳐야 할 사실이지 숨길 것이 아닙니다.

질의 문자열은 literature_review.py 에서 가져옵니다. 여기서 다시 적으면 두 곳이 어긋납니다.
"""
import datetime as dt
import json
import os
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import literature_review as LR                      # noqa: E402

FROZEN = 'results/candidates_v2/literature_review.tsv'
OUT = 'results/proteome_meta'
EUTILS = LR.EUTILS


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def verdict_of(g, n_kid, n_dkd, n_bm):
    """literature_review.py 의 판정 분기를 그대로 옮긴 것."""
    if g in LR.AMBIGUOUS:
        return 'AMBIGUOUS SYMBOL - counts unreliable'
    if n_bm > 0:
        return 'ALREADY PROPOSED as DKD marker'
    if n_dkd >= 3:
        return 'established in DKD'
    if n_dkd > 0:
        return 'reported in DKD (few papers)'
    if n_kid > 0:
        return 'kidney-studied, NOT in DKD'
    return 'NOVEL - no kidney literature'


def count(term, retmax=8):
    url = EUTILS + '?' + urllib.parse.urlencode(
        {'db': 'pubmed', 'term': term, 'retmax': retmax, 'retmode': 'json'})
    for _ in range(3):
        try:
            d = json.loads(urllib.request.urlopen(url, timeout=60).read().decode())
            es = d.get('esearchresult', {})
            return int(es.get('count', 0)), es.get('idlist', [])
        except Exception:
            time.sleep(2)
    return -1, []


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)
    today = dt.date.today().isoformat()

    fz = pd.read_csv(FROZEN, sep='\t')
    log('=' * 78)
    log('문헌 근거 장부 — 후보 %d개, 조회일 %s' % (len(fz), today))
    log('=' * 78)
    log('  질의 틀 (literature_review.py 에서 그대로 가져옴)')
    log('    ANY       "<GENE>"[Title/Abstract]')
    log('    KIDNEY    ANY AND %s' % LR.KIDNEY_TERMS[:60] + ' ...')
    log('    DKD       ANY AND %s' % LR.DKD_TERMS[:60] + ' ...')
    log('    BIOMARKER ANY AND DKD AND %s' % LR.BIOMARKER_TERMS)
    log('')

    rows = []
    flips = []
    for _, r in fz.iterrows():
        g = str(r['gene'])
        gq = '"%s"[Title/Abstract]' % g
        q_kid = '%s AND %s' % (gq, LR.KIDNEY_TERMS)
        q_dkd = '%s AND %s' % (gq, LR.DKD_TERMS)
        q_bm = '%s AND %s AND %s' % (gq, LR.DKD_TERMS, LR.BIOMARKER_TERMS)
        n_kid, _ = count(q_kid)
        time.sleep(0.35)
        n_dkd, ids = count(q_dkd)
        time.sleep(0.35)
        n_bm, bm_ids = count(q_bm)
        time.sleep(0.35)

        # 판정 규칙은 literature_review.py 의 것을 그대로 옮긴다. 처음에 KIDNEY 질의를
        # 빼고 다시 적었더니 'kidney-studied, NOT in DKD' 분기가 없어져서, 건수가 하나도
        # 바뀌지 않은 유전자 4개가 등급이 뒤집힌 것처럼 나왔다. 규칙을 다시 쓰면 안 된다.
        now = verdict_of(g, n_kid, n_dkd, n_bm)

        same_class = (now == str(r['verdict']))
        if not same_class:
            flips.append((g, r['verdict'], now, int(r['n_dkd']), n_dkd))

        rows.append(dict(
            gene=g, queried_on=today,
            query_kidney=q_kid, query_dkd=q_dkd, query_biomarker=q_bm,
            n_kidney_today=n_kid,
            n_dkd_frozen=int(r['n_dkd']), n_dkd_today=n_dkd,
            n_dkd_biomarker_frozen=int(r['n_dkd_biomarker']), n_dkd_biomarker_today=n_bm,
            verdict_frozen=r['verdict'], verdict_today=now,
            class_unchanged=same_class,
            ambiguity_reason=LR.AMBIGUOUS.get(g, ''),
            pmids_dkd_today=';'.join(ids),
            pmids_biomarker_today=';'.join(bm_ids),
            pmids_dkd_frozen=str(r.get('dkd_pmids', '') or ''),
        ))
        log('  %-9s DKD %3d -> %3d   marker %2d -> %2d   %s'
            % (g, r['n_dkd'], n_dkd, r['n_dkd_biomarker'], n_bm,
               '' if same_class else '*** 등급 변동'))

    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'literature_ledger.tsv'), sep='\t', index=False)

    log('')
    log('  기호가 모호해 건수를 믿을 수 없다고 표시한 유전자 %d개' % len(LR.AMBIGUOUS))
    for g, why in LR.AMBIGUOUS.items():
        if g in set(fz['gene']):
            log('    %-9s %s' % (g, why))
    log('')
    if flips:
        log('  *** 동결본 이후 등급이 바뀐 유전자 %d개' % len(flips))
        for g, a, b, na, nb in flips:
            log('    %-9s %s  ->  %s   (DKD %d -> %d)' % (g, a, b, na, nb))
        log('    원고가 이 유전자를 신규로 부르고 있다면 고쳐야 합니다.')
    else:
        log('  동결본 이후 등급이 바뀐 유전자는 없습니다.')

    pd.DataFrame([dict(queried_on=today, n_genes=len(t), n_class_changed=len(flips),
                       genes_changed='|'.join(f[0] for f in flips))]).to_csv(
        os.path.join(OUT, 'literature_ledger_summary.tsv'), sep='\t', index=False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
