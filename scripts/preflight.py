#!/usr/bin/env python
"""Refuse to call the submission package ready while anything is unfinished.

Two classes of defect have already occurred once each and are cheap to catch mechanically:

    placeholders      Bracketed fields nobody filled. Harmless in a working draft, but they read
                      as unfinished boilerplate to an editor.
    internal markers  Notes written for ourselves that were never meant to leave the building.
                      A "delete before sending" heading only works if somebody deletes it.

Exit status is 1 when anything is outstanding, so this can gate a release.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 제출 폴더가 둘이다. BMC 판과 BiB 판을 따로 검사한다. 한쪽만 통과한 것을 양쪽 통과로
# 읽어서는 안 된다(BMC 만 검사하고 BiB 폴더를 검사했다고 착각한 적이 있다).
TARGETS = {
    'bmc': dict(
        dir='submission',
        outbound=['dkd-manuscript.tex', 'dkd-references.bib', 'COVER_LETTER.md'],
        expected=(r'^(dkd-manuscript' + chr(92) + r'.(tex|pdf|bbl)|dkd-references' + chr(92) + r'.bib|'
                  r'COVER_LETTER' + chr(92) + r'.md|Fig' + chr(92) + r'd+' + chr(92) + r'.pdf|table_[a-z_]+' + chr(92) + r'.tex|'
                  r'Additional_file_' + chr(92) + r'd+' + chr(92) + r'.(pdf|xlsx)|sn-[a-z-]+' + chr(92) + r'.(cls|bst)|'
                  r'bst|sn-article-template)$'),
    ),
    'bib': dict(
        dir='submission_bib',
        outbound=['bib-manuscript.tex', 'bib-supplement.tex', 'dkd-references.bib',
                  'COVER_LETTER_BiB.md'],
        expected=(r'^(bib-(manuscript|supplement)' + chr(92) + r'.(tex|pdf|bbl)|dkd-references' + chr(92) + r'.bib|'
                  r'COVER_LETTER_BiB' + chr(92) + r'.md|figures|'
                  r'Supplementary_File_' + chr(92) + r'd+' + chr(92) + r'.(pdf|xlsx)|bib-vancouver' + chr(92) + r'.bst)$'),
    ),
}
SUB = os.path.join(ROOT, 'submission')
OUTBOUND = TARGETS['bmc']['outbound']

# files that live in submission/ but are ours, not the journal's
INTERNAL = {'BUILD.md', 'SUBMISSION_CHECKLIST.md', 'TITLE_OPTIONS.md'}

PLACEHOLDER = re.compile(r'\[[A-Z][A-Z0-9 ()/,\'`.-]{2,}\]')

# phrases that must never appear in an outbound file
BANNED = [
    'delete before sending',
    'TO BE COMPLETED',
    'TO BE ADDED',
    'ON ACCEPTANCE',
    'TODO',
    'FIXME',
    'XXX',
    'Notes on tone',
    'lorem ipsum',
]

# sn-jnl.cls sample values and config tokens. Matched case-sensitively, because folding these
# into BANNED made 'FAMILY' match the word 'family' in ordinary prose.
UNEDITED = [
    r'\sur{Author}',
    r'\fnm{First}',
    'iauthor@gmail.com',
    'corresponding@example.org',
    r'\orgdiv{Division}',
    'FAMILY',
    'ORG/REPO',
]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def scan(path):
    """Return (placeholders, banned) as lists of (line_no, text)."""
    ph, bn = [], []
    with open(path, encoding='utf-8') as fh:
        for n, line in enumerate(fh, 1):
            if line.lstrip().startswith('%%'):      # LaTeX comment: ours, not the journal's
                continue
            for m in PLACEHOLDER.finditer(line):
                ph.append((n, m.group(0)))
            low = line.lower()
            for b in BANNED:
                if b.lower() in low:
                    bn.append((n, b))
            for u in UNEDITED:
                if u in line:
                    bn.append((n, u))
    return ph, bn


def non_ascii_comments(path):
    """나가는 tex 안의 우리말 주석. 편집자가 원본을 열면 보인다."""
    if not path.endswith('.tex'):
        return []
    out = []
    for n, line in enumerate(open(path, encoding='utf-8'), 1):
        t = line.strip()
        if t.startswith('%') and any(ord(c) > 127 for c in t):
            out.append((n, t[:60]))
    return out


def main():
    global SUB, OUTBOUND
    target = 'bmc'
    for a in sys.argv[1:]:
        if a.startswith('--target'):
            target = a.split('=', 1)[1] if '=' in a else sys.argv[sys.argv.index(a) + 1]
    cfg = TARGETS[target]
    SUB = os.path.join(ROOT, cfg['dir'])
    OUTBOUND = cfg['outbound']
    problems = 0
    log('=== outbound files (%s: %s) ===' % (target, cfg['dir']))
    for name in OUTBOUND:
        path = os.path.join(SUB, name)
        if not os.path.exists(path):
            log('  MISSING  %s' % name)
            problems += 1
            continue
        ph, bn = scan(path)
        ko = non_ascii_comments(path)
        if not ph and not bn and not ko:
            log('  clean    %s' % name)
            continue
        log('  %-24s %d placeholder(s), %d internal marker(s), %d non-English comment(s)'
            % (name, len(ph), len(bn), len(ko)))
        for n, t in ph:
            log('      L%-5d placeholder  %s' % (n, t))
        for n, t in bn:
            log('      L%-5d INTERNAL     "%s"' % (n, t))
        for n, t in ko:
            log('      L%-5d COMMENT      %s' % (n, t))
        problems += len(ph) + len(bn) + len(ko)

    # An internal file inside submission/ is a hazard only if it could be attached by mistake;
    # flag its presence so the list stays deliberate rather than accidental.
    log('')
    log('=== files in %s/ that must NOT be uploaded ===' % cfg['dir'])
    # 제출 폴더에 무엇이 있어야 하는지를 정해 두고, 그 밖의 것은 전부 보고한다.
    # 이름이 한글인 낡은 사본이 한 달 넘게 이 폴더에 있었는데 아무도 보지 않았다.
    # 사람이 이 폴더를 열어 올리므로, 낯선 파일은 그 자체로 사고 위험이다.
    expected = re.compile(cfg['expected'])
    build = re.compile(r'^P\d+_[A-Za-z0-9_]+\.pdf$')
    stray = 0
    for name in sorted(os.listdir(SUB)):
        if name in INTERNAL:
            log('  internal %s' % name)
        elif build.match(name):
            log('  build     %s  (figure source an Additional file embeds)' % name)
        elif not expected.match(name):
            log('  *** unexpected  %s' % name)
            stray += 1
    if stray:
        log('  %d file(s) in submission/ are neither outbound nor known internal notes.'
            % stray)
        problems += stray

    log('')
    if problems:
        log('NOT READY: %d item(s) outstanding across the outbound files.' % problems)
        log('Fill them, or move the text out of submission/ if it was never meant to go.')
        return 1
    log('READY: no placeholders or internal markers in any outbound file.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
