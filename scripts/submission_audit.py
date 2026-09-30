#!/usr/bin/env python
"""투고 위생 감사 - 번호와 참조와 부록 연결이 실제로 맞는가.

내용 검증(verify_provenance.py)과는 다른 축입니다. 저기는 "숫자가 결과 파일과 같은가"를
묻고, 여기는 "심사자가 원고를 열었을 때 가리키는 것이 실제로 거기 있는가"를 묻습니다.
그림 한 장을 넣거나 빼면 뒤의 번호가 전부 밀리는데, LaTeX는 \\ref 를 알아서 고쳐주기 때문에
빌드가 조용히 성공합니다. 조용히 어긋나는 것은 다음 셋입니다.

    1  본문이 그림을 번호 순서대로 인용하지 않는 경우. BMC는 순서 인용을 요구합니다.
    2  "Figure 3" 처럼 번호를 손으로 적어둔 곳. \\ref 가 아니므로 재번호를 따라가지 않습니다.
    3  Additional file N 을 언급하는데 그 파일이 제출 폴더에 없는 경우.

전부 원고 파일과 디스크에서 직접 읽습니다. 종료 코드 1이면 아직 보낼 수 없습니다.
"""
import collections
import os
import re
import sys

TEX = 'submission/dkd-manuscript.tex'
BS = '\\'
E = re.escape(BS)
NL = chr(10)


def say(*a):
    # 출력이 파이프로 잡히면 콘솔 인코딩이 cp949 가 되어 일부 문자를 못 찍고 죽는다.
    # 감사 결과가 인코딩 때문에 실패로 보이면 안 되므로 찍을 수 없는 글자는 버린다.
    enc = sys.stdout.encoding or 'utf-8'
    print(*[str(x).encode(enc, 'replace').decode(enc) for x in a])


# ---------------------------------------------------------------- BMC 형식


def braced(s, start):
    """{ 부터 짝이 맞는 } 까지. 중첩 명령이 들어 있어 정규식으로는 못 자른다."""
    depth = 0
    for k in range(start, len(s)):
        if s[k] == '{':
            depth += 1
        elif s[k] == '}':
            depth -= 1
            if depth == 0:
                return s[start + 1:k]
    return ''


# BMC 연구논문의 Declarations 는 이 일곱 개를 이 순서로 요구한다.
# (Authors' information 은 선택이라 넣지 않는다.)
DECLARATIONS = ['Ethics approval and consent to participate',
                'Consent for publication',
                'Availability of data and materials',
                'Competing interests',
                'Funding',
                "Authors' contributions",
                'Acknowledgements']
BODY_ORDER = ['Background', 'Methods', 'Results', 'Discussion', 'Conclusions']
ABSTRACT_HEADS = ['Background:', 'Results:', 'Conclusions:']


def broken_commands(s):
    r"""백슬래시가 먹혀 반쪽만 남은 명령을 찾는다.

    실제로 있었던 일입니다. 스크립트가 Section~\ref{sec:external} 을 쓰면서 \r 이
    캐리지리턴으로 해석되어 'Section~' + 줄바꿈 + 'ef{sec:external}' 이 됐습니다.
    \ref 자체가 사라졌으므로 LaTeX 은 미해결 참조로 보지 않고, 'efsec:external' 을
    본문에 그대로 찍은 채 빌드가 성공했습니다. 눈으로만 보면 넘어갑니다.

    같은 방식으로 먹히는 것은 \r \t \b \f \n \a \v 입니다. 남는 꼬리는 아래와 같습니다.
    """
    bad = 0
    say('=== 먹힌 명령 ===')
    tails = ('ef', 'extbf', 'extit', 'exttt', 'extrm', 'extsuperscript',
             'egin', 'ibliography', 'race', 'ho', 'ightarrow', 'imes', 'rac',
             'eft', 'ight', 'ootnote', 'ewcommand', 'ame', 'erb')
    hits = []
    for m in re.finditer(r'(?m)^(' + '|'.join(tails) + r')(?=[{~\s])', s):
        hits.append((s[:m.start()].count(NL) + 1, m.group(1)))
    # 줄 끝의 물결표는 그 자리에 있던 명령이 잘려 나갔다는 뜻이다
    for m in re.finditer(r'(?m)~$', s):
        hits.append((s[:m.start()].count(NL) + 1, '~ (줄 끝)'))
    # 제어문자가 그대로 박힌 경우
    for m in re.finditer(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', s):
        hits.append((s[:m.start()].count(NL) + 1, '제어문자 %r' % m.group(0)))
    if hits:
        for line, what in hits:
            say('  *** %d행  %s' % (line, what))
        bad += len(hits)
    else:
        say('  없음')
    say('')
    return bad


def ascii_only(path):
    """투고 폼에 붙여넣는 문서는 ASCII 로만 쓴다.

    커버레터에 유니코드 마이너스와 위첨자가 들어 있었습니다. UTF-8 로는 정상이지만
    투고 시스템이나 메일을 거치면 깨져서, 정제된 원고와 나란히 놓였을 때 인상이
    나빠집니다. 원고(.tex)는 LaTeX 이 조판하므로 여기서 검사하지 않습니다.
    """
    import unicodedata
    if not os.path.exists(path):
        return 0
    t = open(path, encoding='utf-8').read()
    bad = sorted({c for c in t if ord(c) > 126})
    say('=== %s 의 ASCII 여부 ===' % os.path.basename(path))
    if bad:
        for c in bad:
            say('  *** U+%04X %s  %r' % (ord(c), unicodedata.name(c, '?'), c))
        say('')
        return len(bad)
    say('  ASCII 전용 OK')
    say('')
    return 0


def journal_format(s, body):
    """저널 가이드라인 준수 검사.

    구조가 맞는지는 눈으로 한 번 보면 알 수 있지만, 손을 댈 때마다 어긋납니다. 실제로
    약어 목록에 SVA 를 적어놓고 본문에서는 계속 'surrogate-variable adjustment' 로 풀어
    써서, 정의만 있고 쓰이지 않는 약어가 남아 있었습니다.
    """
    bad = 0
    say('=== BMC 형식 ===')

    m = re.search(E + r'title\[([^]]*)\]', s)
    if m:
        n = len(m.group(1))
        i = s.index(BS + "title[")
        full = ' '.join(braced(s, s.index('{', i + len(BS + 'title[') + n)).split())
        ok = len(full) <= 150
        say('  제목 %d자 (BMC 상한 150) %s' % (len(full), 'OK' if ok else '*** 초과'))
        bad += 0 if ok else 1
    else:
        say('  *** 짧은 제목이 없습니다. 러닝 헤드가 비게 됩니다.')
        bad += 1

    i = s.find(BS + 'abstract{')
    if i < 0:
        say('  *** 초록이 없습니다')
        bad += 1
    else:
        ab = braced(s, s.index('{', i))
        heads = re.findall(E + r'textbf\{([^}]*)\}', ab)
        ok = heads[:3] == ABSTRACT_HEADS
        say('  초록 소제목 %s %s' % (heads[:3], 'OK' if ok else '*** BMC는 %s' % ABSTRACT_HEADS))
        bad += 0 if ok else 1
        words = len(re.sub(r'[{}$]', '', re.sub(E + r'[a-zA-Z]+\*?', '', ab)).split())
        say('  초록 %d단어 (상한 350) %s' % (words, 'OK' if words <= 350 else '*** 초과'))
        bad += 0 if words <= 350 else 1
        if 'cite' in ab:
            say('  *** 초록에 인용이 있습니다. BMC 초록은 인용을 넣지 않습니다.')
            bad += 1

    m = re.search(E + r'keywords\{([^}]*)\}', s, re.S)
    n_kw = len([x for x in m.group(1).split(',') if x.strip()]) if m else 0
    ok = 3 <= n_kw <= 10
    say('  키워드 %d개 (BMC 3~10) %s' % (n_kw, 'OK' if ok else '*** 범위 밖'))
    bad += 0 if ok else 1

    got = [x for x in re.findall(E + r'section\{([^}]*)\}', s) if x in BODY_ORDER]
    ok = got == BODY_ORDER
    say('  본문 절 순서 %s %s' % (got, 'OK' if ok else '*** BMC는 %s' % BODY_ORDER))
    bad += 0 if ok else 1

    numbered = [x for x in re.findall(E + r'section\{([^}]*)\}', s)
                if x not in BODY_ORDER]
    if numbered:
        say('  *** 후미 절이 번호를 받습니다 (section* 이어야 함): %s' % numbered)
        bad += 1
    else:
        say('  후미 절은 전부 번호 없음 OK')

    heads = re.findall(E + r'bmhead\{([^}]*)\}', s[s.find('Declarations'):])
    decl = [h for h in heads if h in DECLARATIONS]
    ok = decl == DECLARATIONS
    say('  Declarations %d/%d개 %s' % (len(decl), len(DECLARATIONS),
                                       'OK' if ok else '*** 누락/순서'))
    if not ok:
        for h in DECLARATIONS:
            if h not in decl:
                say('    빠짐: %s' % h)
        bad += 1

    # 약어 목록에 적어놓고 본문에서 쓰지 않으면 정의만 남는다.
    # body 는 주석을 지운 문자열이라 s 의 인덱스를 그대로 쓰면 자르는 위치가 밀린다.
    # 처음에 그렇게 짜서, 목록 자체가 '본문'에 포함되는 바람에 미사용 약어를 못 잡았다.
    i = s.find('List of abbreviations')
    ib = body.find('List of abbreviations')
    if i > 0 and ib > 0:
        j = s.index(BS + 'backmatter', i)
        listed = re.findall(r'([A-Z][A-Za-z0-9-]{1,9}):', s[i:j])
        head = body[:ib]
        unused = [a for a in listed
                  if not re.search(r'(?<![A-Za-z])' + re.escape(a) + r'(?![A-Za-z])', head)]
        say('  약어 %d개, 본문 미사용 %d개 %s'
            % (len(listed), len(unused), unused if unused else 'OK'))
        bad += 1 if unused else 0

    ok = 'sn-vancouver-num' in s
    say('  서지 양식 sn-vancouver-num %s' % ('OK' if ok else '*** BMC는 Vancouver 번호식'))
    bad += 0 if ok else 1
    say('')
    return bad


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    s = open(TEX, encoding='utf-8').read()
    here = os.path.dirname(TEX)
    bad = 0

    # 주석 줄은 본문이 아니다. \% 는 주석이 아니므로 살린다.
    body = re.sub(r'(?m)(?<!' + E + r')%.*$', '', s)

    for kind, pref in (('figure', 'fig:'), ('table', 'tab:')):
        env = re.findall('(?s)' + E + r'begin\{' + kind + r'\}.*?'
                         + E + r'end\{' + kind + r'\}', s)
        declared = []
        for k, e in enumerate(env, 1):
            m = re.search(E + r'label\{(' + pref + r'[^}]+)\}', e)
            if m:
                declared.append(m.group(1))
            else:
                say('  *** %s %d 번에 label 이 없습니다' % (kind, k))
                bad += 1

        cited = []
        for m in re.finditer(E + r'ref\{(' + pref + r'[^}]+)\}', body):
            if m.group(1) not in cited:
                cited.append(m.group(1))

        say('=== %s %d개 ===' % (kind, len(env)))
        say('  선언 순서  %s' % ', '.join(x.split(':', 1)[1] for x in declared))
        say('  인용 순서  %s' % ', '.join(x.split(':', 1)[1] for x in cited))

        if [x for x in cited if x in declared] != [x for x in declared if x in cited]:
            say('  *** 인용 순서가 번호 순서와 다릅니다. BMC는 순서대로 인용을 요구합니다.')
            bad += 1
        else:
            say('  순서 OK')
        for l in declared:
            if l not in cited:
                say('  *** %s 는 실려 있는데 본문에서 한 번도 인용되지 않습니다' % l)
                bad += 1
        for l in cited:
            if l not in declared:
                say('  *** %s 를 인용하는데 그런 %s 가 없습니다' % (l, kind))
                bad += 1
        say('')

    say('=== 손으로 적은 그림/표 번호 ===')
    hard = re.findall(r'(?i)(?<!' + E + r'ref\{)\b(?:Figure|Fig\.|Table)~?\s*(\d+)', body)
    if hard:
        say('  *** %s - \\ref 가 아니므로 재번호를 따라가지 않습니다' % collections.Counter(hard))
        bad += 1
    else:
        say('  없음')
    say('')

    say('=== 제출 그림이 최신 원본과 같은가 ===')
    # 제출 폴더의 Fig*.pdf 는 results/figures 의 작업 그림을 복사한 것이다. 원본을
    # 다시 그리고 스테이징을 다시 돌리지 않으면 제출본만 낡은 채로 남는데, 어떤 검사도
    # 그것을 보지 않았다. 실제로 그림 8이 그렇게 낡아 있었다.
    src_dir = os.path.join(os.path.dirname(here), 'results', 'figures')
    if os.path.isdir(src_dir):
        import hashlib

        def h(path):
            return hashlib.sha256(open(path, 'rb').read()).hexdigest()

        srcs = {}
        for f in os.listdir(src_dir):
            if f.endswith('.pdf'):
                srcs.setdefault(h(os.path.join(src_dir, f)), []).append(f)
        stale = 0
        for f in sorted(os.listdir(here)):
            if not re.match(r'^Fig\d+\.pdf$', f):
                continue
            hit = srcs.get(h(os.path.join(here, f)))
            say('  %-10s %s' % (f, (', '.join(hit) if hit else
                                    '*** results/figures 의 어떤 그림과도 일치하지 않습니다')))
            if not hit:
                stale += 1
        if stale:
            say('  *** %d개가 낡았습니다. 원본을 다시 그렸다면 스테이징을 다시 도세요.'
                % stale)
            bad += stale
    say('')

    say('=== 그림 파일 ===')
    for g in re.findall(E + r'includegraphics\[[^]]*\]\{([^}]+)\}', s):
        hit = [x for x in (g + '.pdf', g + '.png', g)
               if os.path.exists(os.path.join(here, x))]
        say('  %-28s %s' % (g, hit[0] if hit else '*** 제출 폴더에 없습니다'))
        if not hit:
            bad += 1
    say('')

    say('=== Additional file ===')
    ment = collections.Counter(re.findall(r'(?i)Additional file~?s?\s*(\d+)', body))
    disk = sorted(f for f in os.listdir(here) if f.lower().startswith('additional'))
    say('  본문 언급  %s' % (dict(ment) or '없음'))
    say('  실제 파일  %s' % (disk or '*** 없음'))
    for n in sorted(ment):
        if not any(n in f for f in disk):
            say('  *** Additional file %s 를 언급하는데 대응 파일이 없습니다' % n)
            bad += 1
    for f in disk:
        m = re.search(r'(\d+)', f)
        if m and m.group(1) not in ment:
            say('  *** %s 를 제출하는데 본문이 언급하지 않습니다' % f)
            bad += 1
    say('')

    bad += broken_commands(s)
    bad += ascii_only(os.path.join(here, 'COVER_LETTER.md'))
    bad += journal_format(s, body)

    say('=== 절 참조 ===')
    lab = set(re.findall(E + r'label\{(sec:[^}]+)\}', s))
    ref = set(re.findall(E + r'ref\{(sec:[^}]+)\}', body))
    for r in sorted(ref - lab):
        say('  *** %s 를 인용하는데 그런 절이 없습니다' % r)
        bad += 1
    say('  인용 %d개, 대상 %d개, 미해결 %d개' % (len(ref), len(lab), len(ref - lab)))
    say('')

    say('=' * 64)
    if bad:
        say('감사 실패 %d건. 이대로 보내면 심사자가 없는 것을 가리키게 됩니다.' % bad)
        return 1
    say('감사 통과 - 원고가 가리키는 것이 전부 실제로 있습니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
