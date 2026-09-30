#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""조판된 원고를 문장 수준에서 훑는다.

기존 검사들은 숫자와 참조만 봤습니다. 글이 맞는지는 아무도 보지 않았고, 실제로 두 곳이
틀린 채 남아 있었습니다. 마침표 뒤가 소문자로 시작하는 문장이 둘이었는데, 하나는 본문
한계 절에 있었고 하나는 그림 1 안의 글이라 원고 소스를 아무리 봐도 나오지 않았습니다.
그래서 소스가 아니라 **만들어진 PDF** 를 읽습니다. 그림 안의 글도 그래야 보입니다.

찾는 것

    마침표 뒤 소문자로 시작하는 문장   참고문헌 목록은 형식상 그런 곳이 많아 제외한다
    같은 단어가 두 번 이어진 곳        표에서 오는 것은 걸러낸다
    한 문장에 줄표가 세 개 이상        읽기 어렵다는 신호

기계가 확신할 수 있는 것만 봅니다. 문체를 판정하려는 것이 아닙니다.
"""
import os
import re
import sys

DOC = 'submission/dkd-manuscript.pdf'
# 참고문헌은 "Nat Rev Genet. 2022;23(3):169-181." 처럼 마침표 뒤 소문자가 정상이다.
REF_MARK = re.compile(r'(?:PMID|doi\.org|https?:|\d{4};\d|@|E-mail)')
SKIP_WORDS = {'et', 'al', 'cf', 'vs', 'e', 'i'}


def source_text():
    """원고 소스와 그림 스크립트의 글을, 명령을 걷어내고 한 줄로 잇는다.

    PDF 에서 뽑은 글은 그림 글자와 본문이 이어 붙습니다. "...is a marker. control obtained
    differently..." 처럼요. 앞은 본문 문장이고 뒤는 그림 축 이름인데, 뽑아낸 글만 보면
    "마침표 뒤 소문자" 로 보입니다.

    진짜 오타라면 그 두 낱말이 **소스에도 이어져 있습니다**. 조판이 만들어낸 인접이라면
    소스 어디에도 그 순서가 없습니다. 그래서 소스에서 확인합니다.
    """
    buf = []
    for f in ['submission/dkd-manuscript.tex'] + [
            os.path.join('scripts', x) for x in os.listdir('scripts')
            if x.endswith('.py') and ('figure' in x or x == 'make_figures.py')]:
        if os.path.exists(f):
            buf.append(open(f, encoding='utf-8', errors='replace').read())
    t = ' '.join(buf)
    E = chr(92)
    t = re.sub(re.escape(E) + r'[a-zA-Z]+\*?', ' ', t)          # LaTeX 명령
    t = re.sub(r'[{}$~' + re.escape(E) + r']', ' ', t)
    return re.sub(r'\s+', ' ', t)


def figure_labels():
    """그림 안의 축 이름. 추출된 글에서 본문 문장 바로 뒤에 붙어 나온다.

    그림 글자도 검사 대상이어야 합니다 — 실제로 그림 1의 오타를 이 검사가 잡았습니다.
    다만 그림의 축 이름이 본문 문장 끝에 이어 붙으면 "마침표 뒤 소문자" 로 보입니다.
    이름을 candidate_figure.py 에서 읽어 두 경우를 구분합니다. 그 파일에서 축 이름을
    바꾸면 여기도 따라 바뀝니다.
    """
    f = os.path.join('scripts', 'candidate_figure.py')
    if not os.path.exists(f):
        return set()
    src = open(f, encoding='utf-8', errors='replace').read()
    return set(re.findall(r"e\['([a-z][^']+)'\]", src))


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    try:
        import pypdf
    except ImportError:
        log('pypdf 가 없어 건너뜁니다.')
        return 0
    if not os.path.exists(DOC):
        log('%s 가 없습니다. 먼저 조판하세요.' % DOC)
        return 1

    r = pypdf.PdfReader(DOC)
    t = re.sub(r'\s+', ' ', ' '.join(p.extract_text() for p in r.pages))
    # 참고문헌 목록 뒤는 형식이 달라 같은 규칙을 적용할 수 없다
    cut = t.find('References')
    body = t[:cut] if cut > 0 else t
    bad = 0

    log('=' * 78)
    log('조판된 원고의 문장 검사')
    log('=' * 78)

    labels = figure_labels()
    # 소스와 견줄 때 공백을 지운다. 양쪽정렬된 줄은 추출하면 공백이 사라져
    # "Shrinkingisnotthewholeofit.asignal" 처럼 나오기 때문이다. 처음에 공백을 요구했다가
    # 본문 오타를 통째로 놓쳤다.
    src_ns = re.sub(r'\s+', '', source_text())
    log('')
    log('  마침표 뒤 소문자  (그림 축 이름 %d개는 제외)' % len(labels))
    n = 0
    for m in re.finditer(r'[a-z0-9]\.\s*([a-z][a-z]{2,})', body):
        seg = body[max(0, m.start() - 70):m.end() + 20]
        after = body[m.end() - len(m.group(1)):m.end() + 40]
        pair = re.sub(r'\s+', '', body[m.start():m.end()])       # "x.word"
        if (m.group(1) in SKIP_WORDS or REF_MARK.search(seg)
                or any(after.startswith(lab) for lab in labels)
                or pair not in src_ns):
            continue
        log('    *** ...%s...' % seg)
        n += 1
    bad += n
    log('    %d건' % n)

    log('')
    log('  같은 단어 반복')
    n = 0
    # 대소문자를 구분한다. "by diagnosis Diagnosis Patients" 처럼 캡션 끝과 표 머리글이
    # 붙어서 생기는 것은 실제 반복이 아니고, 그런 경우 대소문자가 다르다.
    for m in re.finditer(r'\b([A-Za-z]{4,})\s+\1\b', body):
        seg = body[max(0, m.start() - 50):m.end() + 25]
        # 표에서 온 것은 숫자와 슬래시가 섞여 있다
        if re.search(r'\d\s*/\s*\d', seg):
            continue
        log('    *** ...%s...' % seg)
        n += 1
    bad += n
    log('    %d건' % n)

    log('')
    log('  한 문장에 줄표 3개 이상')
    # 표·그림이 문단 사이에 들어가면 추출된 문장 안으로 캡션이 끼어든다. 그러면 원문에는
    # 줄표가 두 개뿐인데 추출본에서는 세 개로 세어진다. 실제로 쪽 나눔이 바뀌었을 때
    # 그런 거짓 양성이 났다. 그래서 이 스크립트가 다른 검사에서 쓰는 방식과 같이,
    # 문장이 원문에 통째로 있는지 먼저 확인한다. 없으면 추출 과정의 합성이다.
    src_ns = re.sub(r'\s+', '', source_text()).replace('---', chr(8212))
    n = 0
    for sent in re.split(r'(?<=[.!?]) ', body):
        if sent.count(chr(8212)) < 3:
            continue
        if re.sub(r'\s+', '', sent) not in src_ns:
            continue
        log('    *** %s' % sent[:140])
        n += 1
    bad += n
    log('    %d건' % n)

    log('')
    if bad:
        log('  *** %d건. 조판된 문서를 눈으로 확인하세요.' % bad)
    else:
        log('  걸리는 곳이 없습니다.')
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
