#!/usr/bin/env python
"""참고문헌이 실제로 존재하는가 — DOI 를 CrossRef 에 조회해 대조한다.

지어낸 참고문헌은 게재 후에 발견되면 철회 사유가 됩니다. 사람이 26개를 눈으로 확인하는
것은 신뢰할 수 없으므로, DOI 를 등록기관에 직접 물어 제목·저자·저널·연도를 맞춰봅니다.

맞춰보는 것은 넷입니다.

    존재      DOI 가 CrossRef 에 있는가. 없으면 그 문헌은 없는 것이다.
    제목      .bib 의 제목과 등록된 제목이 같은가. 단어 겹침 비율로 본다.
    저자      첫 저자 성이 같은가.
    연도      같은가. 온라인 선공개 때문에 1년 차이는 경고로만 둔다.

DOI 가 없는 항목(학회지 미등록, arXiv)은 여기서 판정할 수 없으므로 따로 표시하고
사람이 확인해야 합니다. 종료 코드 1이면 그대로 투고해서는 안 됩니다.
"""
import json
import os
import re
import subprocess
import sys
import time

BIB = 'submission/dkd-references.bib'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def entries(path):
    s = open(path, encoding='utf-8').read()
    out = []
    for m in re.finditer(r'@(\w+)\{([^,]+),(.*?)\n\}', s, re.S):
        body = m.group(3)

        def f(k):
            # 끝의 줄바꿈을 필수로 두면 항목의 '마지막' 필드를 못 읽는다. 엔트리 정규식이
            # 그 줄바꿈을 이미 먹었기 때문이다. 그래서 doi 가 맨 끝에 적힌 문헌은
            # 'DOI 없음' 으로 보고되고 조회 자체를 건너뛰었다. 지어낸 문헌이 하필 그
            # 자리에 있으면 검사가 통째로 지나친다.
            mm = re.search(k + r'\s*=\s*[{"](.*?)[}"],?\s*(?:\n|$)', body, re.S)
            return ' '.join(mm.group(1).split()) if mm else ''
        out.append(dict(key=m.group(2).strip(), type=m.group(1),
                        title=f('title').strip('{}'), author=f('author'),
                        journal=f('journal'), year=f('year'), doi=f('doi')))
    return out


def crossref(doi):
    """curl 로 조회한다. 실패하면 None."""
    url = 'https://api.crossref.org/works/' + doi
    try:
        p = subprocess.run(['curl', '-sS', '-m', '30', '-H',
                            'User-Agent: dkd-reference-check (mailto:etri.jhdh@gmail.com)',
                            url], capture_output=True)
        d = json.loads(p.stdout.decode('utf-8', 'replace'))
        return d.get('message') if d.get('status') == 'ok' else None
    except Exception:
        return None


def words(s):
    return set(re.findall(r'[a-z0-9]+', s.lower())) - {
        'a', 'an', 'the', 'of', 'in', 'for', 'and', 'to', 'on', 'with', 'by'}


def check_pdf_names(refs):
    """조판된 PDF 의 텍스트 층에 저자명이 온전히 남는가.

    화면으로 맞아 보이는 것과 텍스트 층이 맞는 것은 다릅니다. fontenc 없이 조판하면
    LaTeX 이 u 위에 움라우트를 겹쳐 그리고, PDF 에는 'uller' 만 남습니다. 복사도, 색인도,
    출판사의 저자명 자동 추출도 전부 틀어집니다. 그림과 달리 이건 눈으로는 절대 안 보입니다.
    """
    pdf = 'submission/dkd-manuscript.pdf'
    if not os.path.exists(pdf):
        return 0
    try:
        import pypdf
    except ImportError:
        log('  pypdf 가 없어 건너뜁니다')
        return 0
    t = ''.join((p.extract_text() or '')
                for p in pypdf.PdfReader(pdf).pages)
    # bib 의 악센트 문자를 담은 성을 뽑아 PDF 에서 그대로 찾는다
    want = set()
    # bst 는 저자 여섯 명까지 찍고 et al. 로 줄인다(sn-vancouver-num.bst 의 #6).
    # 일곱째 이후는 PDF 에 없는 것이 맞으므로 찾지 않는다.
    for r in refs:
        typeset = ' and '.join(r['author'].split(' and ')[:6])
        for name in re.split(r' and |,', typeset):
            name = name.strip()
            if name and any(ord(c) > 126 for c in name):
                want.add(name)
        for m in re.finditer(r'([A-Z][a-z]*\{' + E_BS + r'["\'`^~][a-z]\}[a-z]*)', typeset):
            want.add(m.group(1))
    log('=== PDF 텍스트 층의 악센트 저자명 ===')
    bad = 0
    # LaTeX 이스케이프는 악센트를 '떼어낸' 형태가 아니라 '붙인' 형태로 바꿔야 한다.
    # 처음에 떼어낸 형태(Buhlmann)로 찾았더니, PDF 에 Bühlmann 이 제대로 있는데도
    # 없다고 보고했다. 검사기가 틀리면 통과도 실패도 믿을 수 없다.
    ACC = {'"': '̈', "'": '́', '`': '̀', '^': '̂', '~': '̃'}

    def unescape(nm):
        import unicodedata as ud

        def one(m):
            return ud.normalize('NFC', m.group(2) + ACC[m.group(1)])
        return re.sub(r'\{' + E_BS + r'(["\'`^~])([a-zA-Z])\}', one, nm)

    for name in sorted(want):
        shown = unescape(name)
        if shown in t:
            log('  OK   %s' % shown)
        else:
            log('  ***  %s  — PDF 텍스트 층에 없습니다 (악센트가 분해되었을 수 있음)' % shown)
            bad += 1
    if not want:
        log('  악센트가 있는 저자명이 없습니다')
    log('')
    return bad


E_BS = re.escape(chr(92))


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    refs = entries(BIB)
    log('참고문헌 %d개' % len(refs))
    log('')
    bad_pdf = check_pdf_names(refs)

    bad, nodoi = bad_pdf, []
    log('%-17s %-6s %-7s %-7s %-6s %s' % ('key', '존재', '제목', '첫저자', '연도', '등록된 제목'))
    log('-' * 100)
    for r in refs:
        if not r['doi']:
            nodoi.append(r)
            log('%-17s %-6s  DOI 없음 — 사람이 확인해야 합니다' % (r['key'], '?'))
            continue
        m = crossref(r['doi'])
        time.sleep(0.3)
        if m is None:
            log('%-17s %-6s' % (r['key'], '없음') + '  *** DOI 가 CrossRef 에 없습니다: ' + r['doi'])
            bad += 1
            continue

        got_title = ' '.join((m.get('title') or [''])[0].split())
        a = words(r['title'])
        b = words(got_title)
        overlap = len(a & b) / max(1, len(a))
        t_ok = overlap >= 0.7

        first = (m.get('author') or [{}])[0].get('family', '')
        want = re.split(r',| and ', r['author'])[0].strip()
        # bib 은 'Family, Given' 이므로 첫 토큰이 성이다
        # 등록에 저자가 없는 옛 문헌(예: Mantel & Haenszel 1959)은 대조할 수 없다.
        # 못 본 것과 틀린 것은 다르므로 '~' 로 표시하고 실패로 세지 않는다.
        a_ok = (first.split()[-1].lower() in want.lower()) if first else None

        got_year = None
        for k in ('published-print', 'published-online', 'issued'):
            if m.get(k, {}).get('date-parts'):
                got_year = m[k]['date-parts'][0][0]
                break
        y_ok = str(got_year) == r['year']
        y_near = got_year is not None and abs(int(got_year) - int(r['year'] or 0)) <= 1

        mark = lambda ok, near=False: 'OK' if ok else ('~' if near else '***')
        log('%-17s %-6s %-7s %-7s %-6s %s'
            % (r['key'], 'OK', mark(t_ok), mark(bool(a_ok), a_ok is None), mark(y_ok, y_near), got_title[:44]))
        if not t_ok:
            log('    제목 불일치  bib: %s' % r['title'][:80])
            log('              등록: %s' % got_title[:80])
        if a_ok is False:
            log('    첫 저자 불일치  bib: %s  /  등록: %s' % (want, first))
        if not y_ok:
            log('    연도  bib: %s  /  등록: %s' % (r['year'], got_year))
        bad += (0 if t_ok else 1) + (0 if a_ok is not False else 1) + (0 if y_ok or y_near else 1)

    log('')
    if nodoi:
        log('DOI 가 없어 자동 확인 불가 %d건:' % len(nodoi))
        for r in nodoi:
            log('  %-17s %s (%s, %s)' % (r['key'], r['title'][:56], r['journal'], r['year']))
    log('')
    log('=' * 78)
    if bad:
        log('불일치 %d건. 지어낸 참고문헌이거나 서지사항이 틀렸습니다.' % bad)
        return 1
    log('DOI 가 있는 %d건은 전부 CrossRef 의 등록 내용과 일치합니다.' % (len(refs) - len(nodoi)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
