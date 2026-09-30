#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""조판된 쪽마다 넘침과 그림-캡션 분리를 찾는다.

그림 폭을 바꾸면 한 쪽만 고쳐도 다른 쪽이 밀린다. 눈으로 30쪽을 넘기기 전에 기계가
먼저 세 가지를 센다.

    1  본문 영역을 넘어간 글자        오른쪽·왼쪽 여백으로 삐져나온 글
    2  그림 없는 그림 캡션            캡션만 다음 쪽으로 넘어간 경우
    3  글자 없는 쪽                   그림만 남고 본문이 밀려난 쪽

여백은 문서마다 다르므로 쪽 크기에서 추정하지 않고 인자로 받는다. 기본값은 이
저장소의 두 판에 맞춰 두었다.

    python scripts/page_layout_audit.py

이 검사도 한계가 있다. 겹침은 글자 상자가 서로 닿는지로만 보고, 선과 점의 겹침은
보지 못한다. 눈으로 보는 일을 대신하지 못한다.
"""
import argparse
import os
import re
import sys

import pymupdf

DOCS = [
    # (경로, 좌우 여백 pt, 위아래 여백 pt)
    ('submission/dkd-manuscript.pdf', 60, 50),
    ('submission_bib/bib-manuscript.pdf', 30, 40),
    ('submission_bib/bib-supplement.pdf', 50, 50),
]
# 진짜 캡션만: 줄 맨 앞의 "Fig. 3", "Figure S1:" 꼴. 본문 속 "Figure 3 of the main
# text" 같은 교차참조는 캡션이 아니다.
CAPTION = re.compile(r'(?m)^\s*(?:Fig\.?|Figure)\s*(S?\d+)\s*[.:]?\s+[A-Z(]')


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def audit(path, mx, my, slack=2.0):
    bad = []
    doc = pymupdf.open(path)
    for i, page in enumerate(doc, 1):
        r = page.rect
        box = pymupdf.Rect(mx - slack, my - slack, r.width - mx + slack, r.height - my + slack)
        spans = [sp for blk in page.get_text('dict')['blocks']
                 for line in blk.get('lines', []) for sp in line.get('spans', [])
                 if sp['text'].strip()]
        # 좌우로 삐져나온 글만 본다. 쪽번호와 머리글은 아래위 여백에 있는 것이 정상이다.
        def overshoot(sp):
            return max(sp['bbox'][2] - box.x1, box.x0 - sp['bbox'][0])
        out = [sp for sp in spans if overshoot(sp) > 0]
        if out:
            worst = max(out, key=overshoot)
            bad.append((i, 'outside text area',
                        '%.1fpt: %s' % (overshoot(worst), worst['text'][:40])))
        text = page.get_text()
        caps = CAPTION.findall(text)
        drawings = len(page.get_drawings())
        if caps and drawings < 20:
            bad.append((i, 'caption without a figure', 'caption %s, %d drawings'
                        % (', '.join(sorted(set(caps))[:3]), drawings)))
        if not spans and drawings:
            bad.append((i, 'page with no text', '%d drawings' % drawings))
    return len(doc), bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slack', type=float, default=2.0, help='여백을 넘어도 넘긴 것으로 보지 않을 pt')
    a = ap.parse_args()
    total = 0
    for path, mx, my in DOCS:
        if not os.path.exists(path):
            log('%-42s 파일 없음' % path)
            continue
        n, bad = audit(path, mx, my, a.slack)
        log('')
        log('=== %s (%d쪽) ===' % (path, n))
        if not bad:
            log('  넘침·분리 없음')
        for page, kind, detail in bad:
            log('  %3d쪽  %-26s %s' % (page, kind, detail))
        total += len(bad)
    log('')
    log('의심 항목 %d건. 글자 상자 기준이므로 선·점의 겹침은 눈으로 확인해야 한다.' % total)
    return 1 if total else 0


if __name__ == '__main__':
    sys.exit(main())
