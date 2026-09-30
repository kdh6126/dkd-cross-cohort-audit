#!/usr/bin/env python
"""KPMP는 진짜 통합 다중 오믹스인가, 아니면 층마다 다른 사람인가.

질환별 집계로는 CKD에 전사체 72명과 단백체 14명이 있다는 것까지만 알 수 있습니다. 그런데
그 14명이 72명 안에 들어 있는지, 아니면 완전히 다른 사람들인지가 결정적입니다.

    겹치면    같은 환자에서 여러 층을 잰 것이므로 통합 다중 오믹스입니다. 샘플 단위로
              조인할 수 있고, 이 프로젝트가 계속 "공개 데이터에는 없다"고 적어온 바로
              그것이 존재한다는 뜻입니다.
    안 겹치면 층마다 다른 사람이므로 교차 오믹스입니다. 지금까지의 서술이 맞습니다.

참여자 단위로 확인합니다. snRNA h5ad에서 참여자 ID를 얻고, 각 ID에 대해 KPMP가 보유한
데이터 타입을 조회해 교차표를 만듭니다.
"""
import collections
import json
import os
import sys
import time
import urllib.request

import pandas as pd

OUT = 'results/kpmp_overlap'
DONORS = 'data/raw/kpmp/sn_donors.json'
CACHE = 'data/raw/kpmp/participant_datatypes.json'
ENDPOINT = 'https://atlas.kpmp.org/graphql'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def gql(query):
    r = urllib.request.Request(ENDPOINT, data=json.dumps({'query': query}).encode(),
                               headers={'Content-Type': 'application/json'})
    return json.load(urllib.request.urlopen(r, timeout=120))


def datatypes_for(rid):
    q = ('{getExperimentalStrategyCountsByParticipant(redcapId:"%s")'
         '{dataType count}}' % rid)
    try:
        d = gql(q)
        rows = (d.get('data') or {}).get('getExperimentalStrategyCountsByParticipant') or []
        # API는 참여자가 가진 것만 주는 것이 아니라 전체 카탈로그를 count와 함께 준다.
        # count가 0인 항목까지 '보유'로 세면 모든 참여자가 모든 층을 가진 것처럼 보인다.
        return {r['dataType']: r['count'] for r in rows if (r.get('count') or 0) > 0}
    except Exception:
        return None


def main():
    os.makedirs(OUT, exist_ok=True)
    donors = json.load(open(DONORS, encoding='utf-8'))
    cache = json.load(open(CACHE, encoding='utf-8')) if os.path.exists(CACHE) else {}

    todo = [d for d in donors if d not in cache]
    if todo:
        log('참여자 %d명 조회 (캐시 %d명)' % (len(todo), len(cache)))
        for k, rid in enumerate(todo, 1):
            cache[rid] = datatypes_for(rid)
            if k % 20 == 0:
                log('  %d / %d' % (k, len(todo)))
                json.dump(cache, open(CACHE, 'w', encoding='utf-8'), ensure_ascii=False)
            time.sleep(0.15)
        json.dump(cache, open(CACHE, 'w', encoding='utf-8'), ensure_ascii=False)

    ok = {k: v for k, v in cache.items() if v}
    log('')
    log('조회 성공 %d / %d명' % (len(ok), len(donors)))
    if not ok:
        log('  어떤 참여자도 조회되지 않았습니다. redcapId 형식이 다를 수 있습니다.')
        log('  h5ad의 donor_id가 KPMP redcapId와 같은 체계인지 확인이 필요합니다.')
        return 1

    # 어떤 데이터 타입들이 등장하는가
    all_types = collections.Counter()
    for v in ok.values():
        all_types.update(v.keys())
    log('')
    log('=' * 76)
    log('참여자가 보유한 데이터 타입')
    log('=' * 76)
    for t, n in all_types.most_common():
        log('  %-46s %3d명' % (t[:44], n))

    # 오믹스 층으로 묶는다
    def layer(t):
        s = t.lower()
        if 'proteom' in s:
            return 'proteomics'
        if 'metabolom' in s or 'lipidom' in s:
            return 'metabolomics'
        if 'single-nucleus' in s or 'snrna' in s or 'single-cell' in s or 'scrna' in s \
                or 'transcriptom' in s or 'rnaseq' in s or 'rna-seq' in s:
            return 'transcriptomics'
        if 'imag' in s or 'micro' in s or 'histolog' in s or 'ct ' in s:
            return 'imaging'
        return 'other'

    rows = []
    for rid, dt in ok.items():
        layers = {layer(t) for t in dt}
        rows.append(dict(participant=rid, disease=donors.get(rid),
                         n_datatypes=len(dt), layers='|'.join(sorted(layers)),
                         has_transcriptomics='transcriptomics' in layers,
                         has_proteomics='proteomics' in layers,
                         has_metabolomics='metabolomics' in layers))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(OUT, 'participant_layers.tsv'), sep='\t', index=False)

    # 참여자 한 명이 한 행인 표는 이 분석에 필요하지만 공개할 이유가 없다. 논문이 쓰는
    # 것은 "질환별로 몇 명이 어떤 층 조합을 가지는가"라는 집계뿐이다. 공개본에는 이쪽만
    # 넣는다 (make_release.py 가 participant_layers.tsv 를 제외한다).
    agg = (d.groupby(['disease', 'layers'])
           .agg(n_participants=('participant', 'size'),
                n_transcriptomics=('has_transcriptomics', 'sum'),
                n_proteomics=('has_proteomics', 'sum'),
                n_metabolomics=('has_metabolomics', 'sum'))
           .reset_index()
           .sort_values(['disease', 'n_participants'], ascending=[True, False]))
    agg.to_csv(os.path.join(OUT, 'layer_summary.tsv'), sep='\t', index=False)

    log('')
    log('=' * 76)
    log('층 조합별 참여자 수')
    log('=' * 76)
    for combo, n in d['layers'].value_counts().items():
        log('  %-48s %3d명' % (combo[:46], n))

    both = d[d['has_transcriptomics'] & d['has_proteomics']]
    log('')
    log('=' * 76)
    log('핵심 — 전사체와 단백체를 모두 가진 참여자')
    log('=' * 76)
    log('  전체: %d / %d명' % (len(both), len(d)))
    if len(both):
        for k, n in both['disease'].value_counts().items():
            log('    %-30s %3d명' % (k, n))
    log('')
    if len(both) >= 5:
        log('  같은 참여자에서 두 층이 측정됐습니다. 이것은 통합 다중 오믹스입니다.')
        log('  샘플 단위 조인이 가능하므로, 지금까지 "공개 데이터에는 없다"고 적어온')
        log('  서술을 CKD에 한해 고쳐야 합니다.')
    elif len(both):
        log('  겹치는 참여자가 %d명뿐입니다. 존재하지만 분석에 쓰기는 어렵습니다.' % len(both))
    else:
        log('  겹치는 참여자가 없습니다. 층마다 다른 사람이므로 교차 오믹스가 맞고,')
        log('  지금까지의 서술이 정확합니다.')

    log('')
    log('  한계 — snRNA h5ad에 등장하는 참여자만 조회했습니다. 단백체만 가진 참여자는')
    log('  이 목록에 애초에 없으므로, 위 숫자는 "전사체 보유자 중 단백체도 가진 비율"')
    log('  입니다. 전체 KPMP 참여자 기준 비율이 아닙니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
