#!/usr/bin/env python
"""데이터셋이 어디서 어떤 이름으로 공개됐고 어느 논문에 딸린 것인가.

설명자료는 지금 데이터셋을 우리가 붙인 약칭으로만 부르고 있습니다. 받는 사람이
GSE104948 을 직접 찾아보려면 공개된 원래 제목과 그 데이터를 낸 논문이 필요합니다.
그리고 GSE 와 ST 가 서로 다른 저장소라는 것도 어디에도 적혀 있지 않습니다.

두 저장소에 직접 물어 받아옵니다.

    GSE  NCBI Gene Expression Omnibus. E-utilities 로 공식 제목·플랫폼·표본 수·PMID.
    ST   Metabolomics Workbench. REST API 로 연구 제목·기관·표본 수.

받은 것은 results/dataset_provenance.tsv 에 저장하고 설명자료가 그것을 읽습니다.
손으로 옮겨 적으면 언젠가 어긋납니다.
"""
import json
import os
import subprocess
import urllib.parse
import sys
import time

import pandas as pd

OUT = 'results/dataset_provenance.tsv'
CACHE = 'data/raw/provenance_cache.json'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url):
    p = subprocess.run(['curl', '-sS', '-m', '40', url], capture_output=True)
    return p.stdout.decode('utf-8', 'replace')


def geo(acc):
    """GEO 시리즈 하나의 공식 제목과 딸린 논문."""
    e = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    # 대괄호를 그대로 두면 NCBI 가 질의를 받지 못하고 빈 결과를 준다. 인코딩해야 한다.
    term = urllib.parse.quote('%s[ACCN] AND gse[ETYP]' % acc)
    s = get(e + 'esearch.fcgi?db=gds&term=%s&retmode=json' % term)
    try:
        uid = json.loads(s)['esearchresult']['idlist'][0]
    except Exception:
        return None
    time.sleep(0.4)
    s = get(e + 'esummary.fcgi?db=gds&id=%s&retmode=json' % uid)
    try:
        d = json.loads(s)['result'][uid]
    except Exception:
        return None
    return dict(title=' '.join((d.get('title') or '').split()),
                platform='/'.join(d.get('gpl', '').split(';')),
                n_samples=d.get('n_samples'),
                pmid=';'.join(str(x) for x in (d.get('pubmedids') or [])),
                taxon=d.get('taxon', ''))


def workbench(acc):
    """Metabolomics Workbench 연구 하나."""
    s = get('https://www.metabolomicsworkbench.org/rest/study/study_id/%s/summary' % acc)
    try:
        d = json.loads(s)
    except Exception:
        return None
    if 'study_title' not in d:               # 여러 건이면 첫 건
        d = list(d.values())[0] if d else {}
    return dict(title=' '.join((d.get('study_title') or '').split()),
                platform=d.get('analysis_type', ''),
                n_samples=d.get('subject_species', ''),
                pmid='', taxon=d.get('subject_species', ''),
                institute=' '.join((d.get('institute') or '').split()))


def pub(pmid):
    """PMID 하나의 서지사항."""
    e = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
    s = get(e + '?db=pubmed&id=%s&retmode=json' % pmid)
    try:
        d = json.loads(s)['result'][pmid]
    except Exception:
        return ''
    au = (d.get('authors') or [{}])
    first = au[0].get('name', '') if au else ''
    return '%s et al., %s (%s)' % (first, d.get('source', ''),
                                   (d.get('pubdate', '') or '')[:4])


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    import sqlite3
    c = sqlite3.connect('db/dkd.sqlite')
    accs = [r[0] for r in c.execute('select accession from study order by study_id')]
    c.close()

    cache = json.load(open(CACHE, encoding='utf-8')) if os.path.exists(CACHE) else {}
    rows = []
    for a in accs:
        if a not in cache:
            log('  조회 %s ...' % a)
            cache[a] = (workbench(a) if a.startswith('ST') else geo(a)) or {}
            time.sleep(0.5)
            os.makedirs(os.path.dirname(CACHE), exist_ok=True)
            json.dump(cache, open(CACHE, 'w', encoding='utf-8'), ensure_ascii=False)
        d = dict(cache[a])
        d['accession'] = a
        d['repository'] = ('Metabolomics Workbench' if a.startswith('ST')
                           else 'NCBI Gene Expression Omnibus')
        if d.get('pmid'):
            first = d['pmid'].split(';')[0]
            if 'citation' not in d:
                d['citation'] = pub(first)
                cache[a] = d
                json.dump(cache, open(CACHE, 'w', encoding='utf-8'), ensure_ascii=False)
                time.sleep(0.4)
        rows.append(d)

    t = pd.DataFrame(rows)
    keep = [x for x in ('accession', 'repository', 'title', 'platform', 'n_samples',
                        'taxon', 'pmid', 'citation', 'institute') if x in t.columns]
    t = t[keep]
    t.to_csv(OUT, sep='\t', index=False)

    log('')
    log('%-11s %-9s %s' % ('접근번호', '표본', '공개된 제목'))
    for _, r in t.iterrows():
        log('%-11s %-9s %s' % (r['accession'], str(r.get('n_samples', ''))[:8],
                               str(r.get('title', ''))[:70]))
    log('')
    log('  %s 에 저장했습니다.' % OUT)
    n_pub = int(t['pmid'].astype(str).str.len().gt(0).sum()) if 'pmid' in t else 0
    log('  논문이 딸린 데이터셋 %d / %d' % (n_pub, len(t)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
