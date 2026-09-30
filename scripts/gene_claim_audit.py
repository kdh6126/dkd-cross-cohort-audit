#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""본문이 유전자 이름에 붙인 주장이 후보표와 맞는가.

기존 검사들은 숫자만 봤습니다. "20 of 30 already proposed as DKD markers" 의 20 은
맞는지 확인하면서, 그 뒤에 이름을 댄 다섯 유전자가 정말 그 20개 안에 있는지는 보지
않았습니다. 실제로 MOXD1 이 그 목록에 잘못 들어가 있었습니다 — 마커로 제안한 논문이
0건인 유전자입니다. 원고를 눈으로 읽다가 걸렸고, 통과하던 검사는 60개가 넘었습니다.

그래서 이름을 댄 주장을 표와 대조합니다. 두 가지 형태를 봅니다.

    등급    "\\textit{A} and \\textit{B} are tier 2"  형태
    판정    "already proposed as DKD markers ... including \\textit{A}, \\textit{B}"  형태

문장을 전부 이해하려는 것이 아닙니다. 기계가 확인할 수 있는 형태만 골라서 봅니다.
확인하지 못한 문장이 남는 것은 감수하고, 확인한 것은 확실하게 합니다.
"""
import os
import re
import sys

import pandas as pd

TEX = 'submission/dkd-manuscript.tex'
CAND = 'results/candidates_v2/master_candidate_table.tsv'
E = chr(92)
RX = E + E          # 정규식 안에서 백슬래시 하나를 뜻한다


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def genes_in(seg):
    return re.findall(RX + r'textit\{([A-Z0-9]+)\}', seg)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    s = open(TEX, encoding='utf-8').read()
    d = pd.read_csv(CAND, sep='\t').set_index('gene')
    bad = 0
    checked = 0

    log('=' * 78)
    log('본문이 이름을 댄 유전자 주장과 후보표 대조')
    log('=' * 78)

    # ---- 1. 등급 주장
    log('')
    log('  등급 주장')
    for m in re.finditer(r'((?:' + RX + r'textit\{[A-Z0-9]+\}[,\s]*(?:and\s+)?){1,6})'
                         r'(?:are\s+)?tier[~\s]*(\d)', s):
        for g in genes_in(m.group(1)):
            if g not in d.index:
                log('    *** %-9s 후보표에 없습니다' % g)
                bad += 1
                continue
            want = m.group(2)
            got = str(d.loc[g, 'tier'])[0]
            checked += 1
            ok = got == want
            bad += 0 if ok else 1
            log('    %-9s 본문 %s등급 · 표 %s등급  %s'
                % (g, want, got, '' if ok else '*** 불일치'))

    # ---- 2. '이미 마커로 제안됨' 목록
    log('')
    log('  이미 DKD 마커로 제안되었다는 목록')
    for m in re.finditer(r'already proposed as DKD markers(.{0,320}?)(?:' + RX
                         + r'cite|\.\s)', s, re.S):
        for g in genes_in(m.group(1)):
            if g not in d.index:
                continue
            v = str(d.loc[g, 'verdict'])
            checked += 1
            ok = 'ALREADY PROPOSED' in v
            bad += 0 if ok else 1
            log('    %-9s %s  %s' % (g, v[:44], '' if ok else '*** 이 판정이 아닙니다'))

    log('')
    log('  확인한 주장 %d건, 어긋난 것 %d건' % (checked, bad))
    if bad:
        log('  *** 본문이 표와 다른 말을 합니다.')
    else:
        log('  이름을 댄 주장은 전부 표와 일치합니다.')
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
