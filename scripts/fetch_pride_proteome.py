#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""PRIDE PXD041884 를 받아서 UniProt 기호까지 붙인다.

이 자료를 쓰는 이유는 proteome_ms_validation.py 의 머리말에 있습니다. 여기서는
받는 방법만 다룹니다.

PRIDE 는 v2 와 v3 두 가지 API 를 동시에 냅니다. 프로젝트 메타는 v2 가, 파일 목록은
v3 가 답합니다. v2 의 files/byProject 는 빈 응답을 줍니다. 둘을 섞어 씁니다.

UniProt 조회는 한 번에 40개씩 끊습니다. 200개를 한 질의에 넣으면 URL 이 너무 길어
400 을 돌려주는데, 오류 메시지가 그 이유를 말해 주지 않아 빈 결과처럼 보입니다.
"""
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request

ACC = 'PXD041884'
FILE = ('Secocnd_Extraction_GPQ_Norm_TS_All_Peptides_FFPE_'
        'Protein_Rollup_Workbook.xlsx')
DEST = 'data/raw/proteomics/pride_PXD041884'
BATCH = 40


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url, tries=3, timeout=120):
    for k in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=timeout).read()
        except Exception as e:
            if k == tries - 1:
                raise
            log('  재시도 %d (%s)' % (k + 1, e))
            time.sleep(3)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(DEST, exist_ok=True)
    import pandas as pd

    xls = os.path.join(DEST, FILE)
    if os.path.exists(xls):
        log('이미 받아 두었습니다: %s' % xls)
    else:
        j = json.loads(get('https://www.ebi.ac.uk/pride/ws/archive/v3/projects/'
                           '%s/files?pageSize=100' % ACC).decode())
        items = j if isinstance(j, list) else (j.get('_embedded', {}) or {}).get('files', [])
        url = None
        for f in items:
            if f.get('fileName') != FILE:
                continue
            for loc in (f.get('publicFileLocations') or []):
                v = str(loc.get('value') or '')
                if v.startswith('ftp://'):
                    url = 'https://' + v[len('ftp://'):]
        if not url:
            log('파일을 찾지 못했습니다. PRIDE 가 경로를 바꾼 것일 수 있습니다.')
            return 1
        log('내려받는 중: %s' % url)
        open(xls, 'wb').write(get(url, timeout=600))
        log('  %.1f MB' % (os.path.getsize(xls) / 1e6))

    out = os.path.join(DEST, 'uniprot_to_gene.tsv')
    if os.path.exists(out):
        log('기호 매핑이 이미 있습니다: %s' % out)
        return 0

    d = pd.read_excel(xls, sheet_name='All Protien Rolled Up')
    acc = sorted(set(d['Accession #'].astype(str).str.strip()))
    log('UniProt 기호를 %d개 조회합니다 (%d개씩)' % (len(acc), BATCH))
    m = {}
    for i in range(0, len(acc), BATCH):
        q = ' OR '.join('accession:%s' % a for a in acc[i:i + BATCH])
        u = ('https://rest.uniprot.org/uniprotkb/search?'
             + urllib.parse.urlencode({'query': q, 'fields': 'accession,gene_primary',
                                       'format': 'tsv', 'size': '500'}))
        t = pd.read_csv(io.StringIO(get(u).decode()), sep='\t')
        for _, r in t.iterrows():
            g = str(r.get('Gene Names (primary)') or '').strip()
            if g and g.lower() != 'nan':
                m[str(r['Entry']).strip()] = g.split()[0]
    pd.DataFrame(sorted(m.items()), columns=['uniprot', 'gene']).to_csv(
        out, sep='\t', index=False)
    log('  %d개 매핑, %s 에 저장' % (len(m), out))
    return 0


if __name__ == '__main__':
    sys.exit(main())
