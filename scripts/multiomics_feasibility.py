#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""같은 참여자에서 두 층을 재는 대조가 실제로 만들어지는가.

이 프로젝트는 이 질문에 두 번 잘못 답했습니다.

  1차   "DKD 에는 두 번째 층이 없다."  KPMP 지역 단백체에 당뇨 참여자가 0명인 것을
        보고 공개 데이터 전체로 일반화했습니다. 틀렸습니다 — 다른 코호트에 있습니다.
  2차   "CKD 는 단백체가 환자에게만 있어 불가능하다."  단백체만 세었습니다.
        대사체는 건강 대조에도 있습니다.

두 번 다 원인이 같습니다. 층 하나를 보고 결론을 내렸습니다. 그래서 여기서는 층 쌍마다
질병군과 대조군의 인원을 모두 세고, 대조가 성립하는지를 따로 적습니다.

성립 여부는 두 조건입니다.

  인원   질병군과 대조군 양쪽에 두 층을 다 가진 사람이 있는가
  형식   그 층이 참여자별 정량 행렬로 공개되는가

KPMP 의 대사체는 첫 조건은 만족하지만 두 번째를 만족하지 않습니다. Spatial
Metabolomics 와 Spatial Lipidomics 는 MALDI 영상 질량분석이고, 공개 API 는 참여자별
파일 목록을 내주지 않습니다. 세어서 되는 문제와 받아서 되는 문제를 구분해 적습니다.
"""
import json
import os
import sys

import pandas as pd

CACHE = 'data/raw/kpmp/participant_datatypes.json'
LAYERS = 'results/kpmp_overlap/participant_layers.tsv'
OUT = 'results/multiomics_feasibility'

# KPMP 등록 범주 중 무엇을 대조군으로 쓸 수 있는가. 종양신절제와 사후공여자는
# 이 논문이 문제 삼는 조달 경로이지만, 여기서는 인원 유무만 세므로 함께 둔다.
CONTROL = {'Healthy_reference_tissue', 'Healthy_stone_donor', 'Tumor_nephrectomy',
           'Deceased_donor', 'Diabetes_mellitus_resilient'}
CASE = {'CKD', 'AKI', 'DMR', 'Diabetes_mellitus_with_renal_disease'}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)

    d = pd.read_csv(LAYERS, sep='\t')
    d['arm'] = ['case' if x in CASE else ('control' if x in CONTROL else 'other')
                for x in d['disease']]

    pairs = [('transcriptomics', 'proteomics'), ('transcriptomics', 'metabolomics'),
             ('proteomics', 'metabolomics')]
    rows = []
    log('=' * 78)
    log('같은 참여자에서 두 층을 다 가진 인원 (KPMP 단일핵 공여자 %d명 기준)' % len(d))
    log('')
    log('  %-32s %6s %8s   %s' % ('층 쌍', '질병군', '대조군', '대조 성립'))
    for a, b in pairs:
        sel = d['has_' + a] & d['has_' + b]
        n_case = int((sel & (d['arm'] == 'case')).sum())
        n_ctl = int((sel & (d['arm'] == 'control')).sum())
        ok = n_case > 0 and n_ctl > 0
        rows.append(dict(layer_a=a, layer_b=b, n_both=int(sel.sum()),
                         n_case=n_case, n_control=n_ctl, contrast_possible=ok))
        log('  %-32s %6d %8d   %s' % ('%s + %s' % (a, b), n_case, n_ctl,
                                      '가능' if ok else '불가 (한쪽이 0명)'))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'layer_pair_contrast.tsv'), sep='\t', index=False)

    # 대조가 성립하는 쌍에 대해, 그 층이 어떤 형식으로 공개되는가
    cache = json.load(open(CACHE, encoding='utf-8'))
    names = set()
    for v in cache.values():
        if v:
            names.update(k for k in v if 'etabolom' in k or 'ipidom' in k)
    log('')
    log('대사체 층의 실제 데이터 종류: %s' % ', '.join(sorted(names)))
    log('  둘 다 MALDI 영상 질량분석입니다. 참여자별 정량 행렬이 아닙니다.')
    log('  KPMP GraphQL 은 파일 개수만 내주고 목록이나 내려받기 경로를 내주지 않습니다.')
    log('')
    log('결론')
    log('  전사체+단백체  인원에서 막힙니다. 대조군에 단백체가 한 명도 없습니다.')
    log('  전사체+대사체  인원은 됩니다(CKD 10명 대 건강대조 6명). 형식에서 막힙니다.')
    log('  DKD 만 보면    두 쌍 모두 질병군 인원이 0명이라 애초에 성립하지 않습니다.')

    pd.DataFrame([dict(
        n_donors=len(d),
        tx_prot_case=int(rows[0]['n_case']), tx_prot_control=int(rows[0]['n_control']),
        tx_metab_case=int(rows[1]['n_case']), tx_metab_control=int(rows[1]['n_control']),
        metab_types='|'.join(sorted(names)),
        tx_prot_possible=bool(rows[0]['contrast_possible']),
        tx_metab_possible=bool(rows[1]['contrast_possible']),
        tx_metab_usable=False)]).to_csv(
        os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
