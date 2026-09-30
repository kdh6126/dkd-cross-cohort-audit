#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""사전 등록 파일이 기록 이후 바뀌지 않았는지 확인한다.

예측(docs/xtissue/PREDICTIONS.md)과 사후 분석 계획(docs/xtissue/POSTHOC_PLAN.md)은 실행 전에
해시를 적어 두었다. 파일이 한 글자라도 바뀌었으면 여기서 멈춘다.

**이 검사가 보증하는 범위.** 여기서 비교하는 것은 로컬 파일과 로컬에 적어 둔 해시뿐이다.
따라서 "기록 이후 바뀌지 않았다"까지만 증명하고, 그 파일이 언제 쓰였는지나 데이터를 보기
전이었다는 사실은 증명하지 못한다. 제3자 등록·타임스탬프가 없으므로 원고에서도 "사전에
명시하고 해시로 고정한 예측"이라고만 적는다. 외부 등록을 하게 되면 그 식별자를 여기에
함께 기록하고 검사에 넣는다. 다만 이미 끝난 분석을 나중에 등록해 놓고 사전등록이었다고
적을 수는 없다. 그 경우에도 원고 표현은 "분석 전에 작성한 로컬 계획과 그 한계"로 남긴다.

바꿔야 할 이유가 생기면 파일을 고치지 말고 새 파일을 만들어 날짜와 이유를 적는다.
"""
import hashlib
import sys

PAIRS = [('docs/xtissue/PREDICTIONS.md', 'results/xtissue/predictions.sha256'),
         ('docs/xtissue/POSTHOC_PLAN.md', 'results/xtissue/posthoc_plan.sha256')]


def main():
    bad = 0
    for doc, rec in PAIRS:
        now = hashlib.sha256(open(doc, 'rb').read()).hexdigest()
        lines = open(rec, encoding='utf-8').read().split()
        was, when = lines[0], lines[-1]
        ok = now == was
        print('%-32s %s  (recorded %s)' % (doc, 'unchanged' if ok else '*** CHANGED', when))
        bad += 0 if ok else 1
    if bad:
        print('사전 등록 파일이 기록 이후 바뀌었습니다. 판정을 믿을 수 없습니다.')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
