#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""후보 코호트의 메타데이터를 사람이 읽을 수 있게 편다.

채택 판단에 필요한 것은 넷이다. 사례·대조 집단이 무엇으로 표시되는지, 대조군이 어떤
조직에서 왔는지(요약·설계 문장), 표본 수, 그리고 같은 연구실의 다른 계열과 표본을
공유하는지. 넷째는 표본 제목·환자 식별자를 계열끼리 맞대어 본다.
"""
import json
import re
import sys

dom = sys.argv[1]
recs = json.load(open('results/xtissue/candidates/%s.json' % dom, encoding='utf-8'))

CUE = re.compile(r'(biops\w*|resect\w*|surg\w*|donor\w*|bariatric|wedge|explant\w*|'
                 r'transplant\w*|endoscop\w*|colonoscop\w*|autops\w*|normal liver|'
                 r'healthy|control\w*|tumou?r|cancer|adjacent|intraoperat\w*|needle|'
                 r'percutaneous|cholecystectom\w*|hepatectom\w*|colectom\w*)', re.I)

for r in recs:
    print('=' * 110)
    print('%s  n=%d  %s  PMID %s' % (r['gse'], r['n'], ','.join(r['platform']),
                                     ','.join(r['pubmed'])))
    print('  TITLE  %s' % r['title'])
    text = r['summary'] + ' || ' + r['design']
    hits = sorted(set(m.group(0).lower() for m in CUE.finditer(text)))
    print('  CUES   %s' % ', '.join(hits))
    for sent in re.split(r'(?<=[.;])\s+', text):
        if CUE.search(sent):
            print('   > %s' % sent.strip()[:230])
    print('  SOURCE %s' % json.dumps(r['source_counts'], ensure_ascii=False)[:300])
    for k, v in r['characteristics'].items():
        print('  CHAR   %-28s %s' % (k[:28], json.dumps(v, ensure_ascii=False)[:230]))
    print('  TITLES %s' % ' | '.join(r['title_examples'])[:230])
