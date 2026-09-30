#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""조직을 넓히기 위한 후보 코호트의 표본 메타데이터를 모은다.

DKD 에서 논문의 두 주장을 가르려면 채취 방식이 다른 질환이 필요하다. 여기서는 두 질환을
본다.

    liver_masld  간 MASLD/MASH. 코호트마다 대조군 채취가 다르다(수술 절제·생체 공여 vs
                 같은 비만수술 중 쐐기 생검). 신장의 GSE162830 처럼 코호트 단위로 예측을
                 시험할 수 있다.
    colon_uc     궤양성 대장염. 사례와 대조가 모두 내시경 생검이라 채취 교란이 없어야 한다.

계열마다 series matrix 의 머리말만 스트리밍으로 읽고 표 본문에서 멈춘다. 행렬 전체를 받지
않으므로 빠르다. 결과는 사람이 읽고 채택 여부와 사례·대조 규칙, 채취 경로를 정하는 데 쓴다.
채택 결정은 이 스크립트가 하지 않는다. DKD 에서 제외한 세 계열은 전부 메타데이터를
사람이 읽어서 찾았다.
"""
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.request

HDR = {'User-Agent': 'dkd-xtissue/1.0 (mailto:etri.jhdh@gmail.com)'}
OUT = 'results/xtissue/candidates'

CANDIDATES = {
    'liver_masld': ['GSE89632', 'GSE48452', 'GSE61260', 'GSE63067', 'GSE164760', 'GSE83452',
                    'GSE66676', 'GSE24807', 'GSE17470', 'GSE37031', 'GSE126848', 'GSE130970',
                    'GSE162694', 'GSE213621', 'GSE163211'],
    'colon_uc': ['GSE75214', 'GSE59071', 'GSE87466', 'GSE16879', 'GSE38713', 'GSE36807',
                 'GSE9452', 'GSE47908', 'GSE53306', 'GSE13367', 'GSE48958', 'GSE107499',
                 'GSE179285', 'GSE22619', 'GSE11223', 'GSE10616'],
}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url, stream=False):
    for k in range(3):
        try:
            r = urllib.request.urlopen(urllib.request.Request(url, headers=HDR), timeout=90)
            return r if stream else r.read()
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(2 + 3 * k)
    raise err


def matrix_files(gse):
    stub = gse[:-3] + 'nnn'
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/' % (stub, gse)
    html = get(url).decode('utf-8', 'replace')
    files = sorted(set(re.findall(r'href="(%s[^"]*_series_matrix\.txt\.gz)"' % gse, html)))
    return [url + f for f in files]


def header(url):
    """gzip 을 스트리밍으로 풀면서 표 본문 직전까지만 읽는다."""
    r = get(url, stream=True)
    hdr = {}
    with gzip.GzipFile(fileobj=r) as gz:
        for raw in io.TextIOWrapper(gz, encoding='utf-8', errors='replace'):
            if raw.startswith('!series_matrix_table_begin'):
                break
            if not raw.startswith('!'):
                continue
            parts = raw.rstrip('\n').split('\t')
            hdr.setdefault(parts[0][1:], []).append([p.strip('"') for p in parts[1:]])
    r.close()
    return hdr


def summarise(gse, hdr):
    one = lambda k: ' '.join(v[0] for v in hdr.get(k, []) if v)
    samples = hdr.get('Sample_geo_accession', [[]])[0]
    n = len(samples)
    chars = {}
    for row in hdr.get('Sample_characteristics_ch1', []):
        for v in row:
            if ':' in v:
                k, val = v.split(':', 1)
                chars.setdefault(k.strip().lower(), []).append(val.strip())
    src = hdr.get('Sample_source_name_ch1', [[]])[0]
    titles = hdr.get('Sample_title', [[]])[0]
    return dict(
        gse=gse, n=n,
        platform=sorted(set(hdr.get('Sample_platform_id', [[]])[0])),
        title=one('Series_title'), summary=one('Series_summary')[:1500],
        design=one('Series_overall_design')[:1500],
        pubmed=sorted(set(sum(hdr.get('Series_pubmed_id', []), []))),
        source_counts=_counts(src),
        characteristics={k: _counts(v) for k, v in chars.items()},
        title_examples=titles[:6],
        samples=samples, sample_titles=titles, sample_source=src,
        sample_chars=hdr.get('Sample_characteristics_ch1', []),
    )


def _counts(vals):
    out = {}
    for v in vals:
        out[v] = out.get(v, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: -kv[1])[:25])


def main():
    os.makedirs(OUT, exist_ok=True)
    for domain, gses in CANDIDATES.items():
        allrec = []
        for gse in gses:
            try:
                recs = []
                for url in matrix_files(gse):
                    recs.append(summarise(gse, header(url)))
                    time.sleep(0.5)
                if not recs:
                    log('%-10s 행렬 파일 없음' % gse)
                    continue
                for r in recs:
                    allrec.append(r)
                    log('%-10s n=%-4d %-12s %s' % (gse, r['n'], ','.join(r['platform']),
                                                    r['title'][:70]))
            except Exception as e:  # noqa: BLE001
                log('%-10s 실패: %s' % (gse, e))
        json.dump(allrec, open(os.path.join(OUT, '%s.json' % domain), 'w', encoding='utf-8'),
                  ensure_ascii=False, indent=1)
        log('== %s: %d records' % (domain, len(allrec)))


if __name__ == '__main__':
    main()
