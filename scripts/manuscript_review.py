#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""투고 직전 원고 전반 점검.

사람이 눈으로 훑으면 놓치는 것들을 기계가 세어 준다. 다섯 가지를 본다.

    1  줄표(---)  어디에 몇 개 쓰였는가. 본문과 캡션을 나눠 센다.
    2  그림 인용  label 마다 ref 가 있는가. 표도 같이 본다.
    3  절 길이    본문 단어 수. 너무 짧은 절은 서술이 빠진 것일 수 있다.
    4  수식       수식 환경이 몇 개이고 어디에 있는가.
    5  참고문헌   bib 에 있는데 인용 안 된 것, 인용됐는데 bib 에 없는 것.

수정은 하지 않는다. 무엇이 걸리는지만 보고한다.
"""
import os
import re
import sys
from collections import Counter

TEX = 'submission/dkd-manuscript.tex'
BIB = 'submission/dkd-references.bib'
B = chr(92)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load():
    return open(TEX, encoding='utf-8').read()


def strip_env(s, name):
    """caption / footnotetext 같은 중괄호 환경의 내용을 뽑아내고 본문에서 지운다."""
    out, rest, i = [], [], 0
    key = B + name + '{'
    while True:
        j = s.find(key, i)
        if j < 0:
            rest.append(s[i:])
            break
        rest.append(s[i:j])
        k, depth = j + len(key), 1
        while k < len(s) and depth:
            if s[k] == '{':
                depth += 1
            elif s[k] == '}':
                depth -= 1
            k += 1
        out.append(s[j + len(key):k - 1])
        i = k
    return ''.join(out), ''.join(rest)


def sections(body):
    """(번호, 제목, 본문) 목록. 표·그림 환경은 빼고 센다."""
    parts = re.split(r'%s(?:sub)?section\{' % re.escape(B), body)
    out = []
    for p in parts[1:]:
        depth, k = 1, 0
        while k < len(p) and depth:
            if p[k] == '{':
                depth += 1
            elif p[k] == '}':
                depth -= 1
            k += 1
        title = re.sub(r'\s+', ' ', p[:k - 1])
        title = re.sub(re.escape(B) + r'label\{[^}]*\}', '', title).strip()
        text = p[k:]
        for env in ('figure', 'table', 'tabular'):
            text = re.sub(r'%sbegin\{%s\*?\}.*?%send\{%s\*?\}'
                          % (re.escape(B), env, re.escape(B), env), ' ', text, flags=re.S)
        text = re.sub(re.escape(B) + r'[a-zA-Z]+', ' ', text)
        out.append((title, len(text.split())))
    return out


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    s = load()
    caps, rest1 = strip_env(s, 'caption')
    foots, body = strip_env(rest1, 'footnotetext')
    bad = 0

    log('=' * 76)
    log('1. 줄표(---)')
    log('=' * 76)
    log('  본문        %d개' % body.count('---'))
    log('  캡션        %d개' % caps.count('---'))
    log('  각주        %d개' % foots.count('---'))
    log('  합계        %d개' % s.count('---'))
    if body.count('---'):
        log('')
        for m in re.finditer(r'.{58}---.{58}', body, re.S):
            log('    ' + re.sub(r'\s+', ' ', m.group(0)))
        bad += body.count('---')

    log('')
    log('=' * 76)
    log('2. 그림과 표의 인용')
    log('=' * 76)
    for kind in ('fig', 'tab'):
        labs = re.findall(re.escape(B) + r'label\{(%s:[^}]+)\}' % kind, s)
        for l in labs:
            n = len(re.findall(re.escape(B) + r'ref\{' + re.escape(l) + r'\}', s))
            mark = '' if n else '   *** 인용 없음'
            log('  %-26s %d회%s' % (l, n, mark))
            if not n:
                bad += 1
    inc = re.findall(re.escape(B) + r'includegraphics\[[^]]*\]\{([^}]+)\}', s)
    log('  포함된 그림 파일: %s' % ', '.join(inc))

    log('')
    log('=' * 76)
    log('3. 절 길이 (본문 단어 수, 표·그림 제외)')
    log('=' * 76)
    for title, n in sections(body):
        mark = '   *** 짧음' if n < 120 else ''
        log('  %-62s %4d%s' % (title[:60], n, mark))
        if n < 120:
            bad += 1

    log('')
    log('=' * 76)
    log('4. 수식')
    log('=' * 76)
    envs = Counter(re.findall(re.escape(B) + r'begin\{(equation\*?|align\*?|gather\*?)\}', s))
    log('  수식 환경: %s' % (dict(envs) or '없음'))
    log('  번호 붙은 수식 label: %s'
        % (re.findall(re.escape(B) + r'label\{(eq:[^}]+)\}', s) or '없음'))
    inline = len(re.findall(r'(?<!\$)\$(?!\$)[^$]+\$', s))
    log('  인라인 수학 조각 약 %d개' % inline)

    log('')
    log('=' * 76)
    log('5. 참고문헌')
    log('=' * 76)
    keys = set(re.findall(r'@[a-zA-Z]+\{([^,]+),', open(BIB, encoding='utf-8').read()))
    cited = set()
    for m in re.findall(re.escape(B) + r'cite[a-z]*\{([^}]+)\}', s):
        cited |= {x.strip() for x in m.split(',')}
    log('  bib 항목 %d개, 본문 인용 %d개' % (len(keys), len(cited)))
    miss = sorted(cited - keys)
    unused = sorted(keys - cited)
    if miss:
        log('  *** 인용됐으나 bib 에 없음: %s' % ', '.join(miss))
        bad += len(miss)
    if unused:
        log('  *** bib 에 있으나 인용 안 됨: %s' % ', '.join(unused))
        bad += len(unused)
    if not miss and not unused:
        log('  모든 항목이 짝을 이룹니다.')

    log('')
    log('=' * 76)
    log('걸리는 항목 %d건' % bad)
    return 0


if __name__ == '__main__':
    sys.exit(main())
