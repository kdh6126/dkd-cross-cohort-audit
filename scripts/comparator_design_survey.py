#!/usr/bin/env python
"""신장 대사체 연구는 대조군을 어떻게 잡는가 — 전수 조사.

논문의 핵심 권고는 "건강인 말고 같은 질병 안에서 비교하라"입니다. 그 권고가 얼마나 지켜지지
않고 있는지를 주장이 아니라 숫자로 대야 합니다. 지금까지는 "공개 DKD 데이터 11건 중 2건만
질병 대조를 가진다"는 전사체 쪽 수치뿐이었습니다.

이 스크립트는 Metabolomics Workbench의 신장 관련 사람 연구를 전수로 긁어, 각 연구의 비교군
구조를 네 가지로 분류합니다.

    healthy_control   환자 대 건강인. 질병 유무가 유일한 차이. 가장 흔하고 가장 약하다.
    within_disease    같은 질병 안에서 진행/중증도/아형을 비교. 논문이 권고하는 설계.
    disease_contrast  서로 다른 질병끼리 비교. ST003255가 여기.
    other             치료 전후, 시간 경과, 품질관리 등.

분류는 factors 필드의 문자열을 규칙으로 읽습니다. 완벽하지 않으므로 규칙을 코드에 그대로
남기고, 분류 결과를 파일로 내보내 사람이 검토할 수 있게 합니다.

동물 실험과 암 연구는 제외합니다. 신장이 등장하지만 이 질문과 무관합니다.
"""
import collections
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

OUT = 'results/design_survey'
TERMS = ['chronic kidney', 'kidney disease', 'renal', 'nephropathy', 'dialysis', 'eGFR']

EXCLUDE_TITLE = re.compile(
    r'carcinoma|cancer|tumou?r|lung|hepat|mouse|mice|rat\b|MoTrPAC|exercise|'
    r'streptozotocin|porcine|zebrafish', re.I)

HEALTHY = re.compile(r'health|normal|control(?!led)|non[- ]?disease|hc\b', re.I)
WITHIN = re.compile(r'progress|stage|early|late|severity|grade|mild|moderate|severe|'
                    r'baseline|follow.?up|responder|decline|rapid|slow', re.I)
# 질병명은 약어로 적히는 경우가 많다. 전체 이름만 찾으면 가장 중요한 연구를 놓친다 —
# ST003255의 군 이름은 IgAN / MN / HN / DKD / Normal 이라, 약어를 넣지 않으면
# 'IgAN' 하나만 잡혀 질병 대조가 아닌 것으로 분류된다.
# 질병명은 약어로 적히는 경우가 많다. 전체 이름만 찾으면 가장 중요한 연구를 놓친다 —
# ST003255의 군 이름은 IgAN / MN / HN / DKD / Normal 이라, 약어를 넣지 않으면
# 'IgAN' 하나만 잡혀 질병 대조가 아닌 것으로 분류된다.
# 단어 경계는 패턴 문자열에 직접 쓰지 않고 여기서 조립한다. 소스에 넣으면 편집 과정에서
# 제어문자로 바뀌어 조용히 아무것도 매치하지 않는 정규식이 된 적이 있다.
_ABBR = ['IgAN', 'MN', 'HN', 'DKD', 'DN', 'FSGS', 'MCD', 'SLE', 'ADPKD', 'AKI', 'ANCA']
_WORDS = ['nephropathy', 'glomerul', 'sclerosis', 'diabet', 'lupus', 'vasculitis',
          'membranous', 'minimal change', 'amyloid', 'nephritis']
_WB = chr(92) + 'b'
DISEASE = re.compile('|'.join([_WB + a + _WB for a in _ABBR] + _WORDS), re.I)
OTHER = re.compile(r'pre.?treatment|post.?treatment|placebo|QC|pooled|blank|time|visit|'
                   r'dose|supplement', re.I)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def fetch(url):
    return json.load(urllib.request.urlopen(url, timeout=120))


def scan():
    seen = {}
    for t in TERMS:
        u = ('https://www.metabolomicsworkbench.org/rest/study/study_title/%s/summary'
             % urllib.parse.quote(t))
        try:
            d = fetch(u)
        except Exception:
            continue
        if isinstance(d, dict) and 'study_id' in d:
            d = {'0': d}
        for r in (d.values() if isinstance(d, dict) else d):
            seen[r.get('study_id')] = r
        time.sleep(0.35)
    return seen


def classify(levels):
    """비교군 이름들을 보고 설계 유형을 정한다. 우선순위가 있다."""
    txt = ' | '.join(levels)
    has_h = bool(HEALTHY.search(txt))
    has_w = bool(WITHIN.search(txt))
    n_dis = len(set(m.group(0).lower() for m in DISEASE.finditer(txt)))
    if n_dis >= 2:
        return 'disease_contrast'
    if has_w and not has_h:
        return 'within_disease'
    if has_h and has_w:
        return 'within_disease+healthy'
    if has_h:
        return 'healthy_control'
    if OTHER.search(txt):
        return 'other'
    return 'unclear'


def main():
    os.makedirs(OUT, exist_ok=True)
    seen = scan()
    log('신장 관련 대사체 연구 %d건 수집' % len(seen))

    rows = []
    for sid, s in sorted(seen.items()):
        title = s.get('study_title') or ''
        if EXCLUDE_TITLE.search(title):
            rows.append(dict(study_id=sid, n=s.get('number_of_samples'), title=title[:70],
                             design='excluded', levels='', note='암/동물/무관'))
            continue
        try:
            f = fetch('https://www.metabolomicsworkbench.org/rest/study/study_id/%s/factors'
                      % sid)
        except Exception:
            rows.append(dict(study_id=sid, n=s.get('number_of_samples'), title=title[:70],
                             design='no_factors', levels='', note='factors 조회 실패'))
            continue
        recs = list(f.values()) if isinstance(f, dict) else f
        levels = sorted({str(r.get('factors')) for r in recs if isinstance(r, dict)})
        rows.append(dict(study_id=sid, n=s.get('number_of_samples'), title=title[:70],
                         design=classify(levels), levels=' ;; '.join(levels)[:220], note=''))
        time.sleep(0.3)

    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(OUT, 'kidney_metabolomics_designs.tsv'), sep='\t', index=False)

    human = d[~d['design'].isin(['excluded', 'no_factors'])]
    log('')
    log('=' * 76)
    log('설계 유형 분포 (암·동물 제외 %d건)' % len(human))
    log('=' * 76)
    c = human['design'].value_counts()
    for k, v in c.items():
        log('  %-24s %3d건  (%.0f%%)' % (k, v, 100 * v / len(human)))

    matched = c.get('within_disease', 0) + c.get('disease_contrast', 0) + \
        c.get('within_disease+healthy', 0)
    log('')
    log('  질병을 고정한 비교를 포함하는 연구: %d / %d (%.0f%%)'
        % (matched, len(human), 100 * matched / max(len(human), 1)))
    log('  건강인 대조만 쓰는 연구            : %d / %d (%.0f%%)'
        % (c.get('healthy_control', 0), len(human),
           100 * c.get('healthy_control', 0) / max(len(human), 1)))

    log('')
    log('=' * 76)
    log('질병을 고정한 설계를 가진 연구')
    log('=' * 76)
    good = human[human['design'].isin(['within_disease', 'disease_contrast',
                                       'within_disease+healthy'])]
    for _, r in good.sort_values('n', ascending=False).iterrows():
        log('  %-10s n=%-5s %-22s %s' % (r['study_id'], r['n'], r['design'], r['title'][:44]))

    unc = human[human['design'] == 'unclear']
    if len(unc):
        log('')
        log('=' * 76)
        log('분류 실패 — 사람이 확인해야 하는 연구')
        log('=' * 76)
        for _, r in unc.iterrows():
            log('  %-10s n=%-5s %s' % (r['study_id'], r['n'], r['title'][:52]))
            log('             군: %s' % (r['levels'][:96] or '(없음)'))

    log('')
    log('=' * 76)
    log('해석')
    log('=' * 76)
    log('  전사체 쪽에서는 공개 DKD 데이터 11건 중 2건만 질병 대조를 가졌습니다.')
    log('  대사체 쪽도 같은 방향인지 이 조사로 처음 확인했습니다.')
    log('')
    log('  한계 — 분류는 factors 문자열 규칙에 기댑니다. 등록자가 군 이름을 어떻게 적었느냐에')
    log('  좌우되고, 논문 본문을 읽지 않았습니다. 그래서 개별 판정보다 분포의 크기만 봅니다.')
    log('  분류 결과 전체를 %s 에 남겼으니 사람이 검토할 수 있습니다.' % OUT)


if __name__ == '__main__':
    main()
