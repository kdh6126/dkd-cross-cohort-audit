#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""초록 단어 수를 센다. BMC Bioinformatics 는 350단어가 한도다.

손으로 세면 틀리고, LaTeX 명령과 수식을 단어로 세면 더 틀린다. 명령과 중괄호를 걷어낸
뒤 센다. 편집할 때마다 이 값을 확인하려고 스크립트로 둔다.
"""
import io
import re
import sys

TEX = 'submission/dkd-manuscript.tex'
LIMIT = 350
B = chr(92)


def main():
    s = io.open(TEX, encoding='utf-8').read()
    i = s.find(B + 'abstract{')
    if i < 0:
        print('초록을 찾지 못했습니다.', file=sys.stderr)
        return 1
    depth, k = 1, i + len(B + 'abstract{')
    while k < len(s) and depth:
        if s[k] == '{':
            depth += 1
        elif s[k] == '}':
            depth -= 1
        k += 1
    t = s[i + len(B + 'abstract{'):k - 1]
    t = re.sub(re.escape(B) + r'[a-zA-Z]+', ' ', t)      # 명령 제거
    t = re.sub(r'[{}$]', ' ', t)                          # 중괄호와 수식 기호
    t = re.sub(r'\s+', ' ', t).strip()
    n = len(t.split())
    print('초록 %d단어 (한도 %d, 여유 %d)' % (n, LIMIT, LIMIT - n))
    for part in ('Background:', 'Results:', 'Conclusions:'):
        j = t.find(part.replace(':', ''))
        if j >= 0:
            print('  %-14s 시작 위치 %d단어째' % (part, len(t[:j].split())))
    return 0 if n <= LIMIT else 1


if __name__ == '__main__':
    sys.exit(main())
