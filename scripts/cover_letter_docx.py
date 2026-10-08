#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Fill the journal cover-letter .docx from COVER_LETTER_BiB.md.

The markdown letter is the audited text (preflight reads it; run_all --submit hashes it). The
journal system wants a document, so this copies the markdown paragraphs into the formatted
template in submission_bib/CoverLetter.docx, paragraph by paragraph, keeping the template's
styles. Re-run after editing the markdown; never edit the .docx by hand or the two drift.

    python scripts/cover_letter_docx.py            # write submission_bib/CoverLetter.docx
    python scripts/cover_letter_docx.py --check    # print what would go in, write nothing
"""
import argparse
import copy
import os
import re
import sys

import docx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
sys.path.insert(0, ROOT)
import submission_config as cfg  # noqa: E402

# --target 로 고른다. GigaScience 판(2026-10-08)은 같은 서식 파일을 submission_giga/ 에 두고 채운다.
TARGETS = {
    'bib': ('submission_bib/COVER_LETTER_BiB.md', 'submission_bib/CoverLetter.docx', 'Briefings in Bioinformatics'),
    'giga': ('submission_giga/COVER_LETTER_Giga.md', 'submission_giga/CoverLetter.docx', 'GigaScience'),
}
MD, DOCX, JOURNAL = TARGETS['bib']
# Contact details the markdown letter does not carry. The phone and fax are the corresponding
# author's office numbers as given on the author's previous cover letter.
PHONE = '+82-31-739-7273'
FAX = '+82-31-739-7264'


def md_paragraphs():
    s = open(MD, encoding='utf-8').read()
    # markers sit on their own lines inside a paragraph; drop the line break with them
    s = re.sub(r'\n?<!--.*?-->', '', s, flags=re.S)
    s = re.sub(r'(?m)^# .*$', '', s)
    s = s.replace('**', '')
    paras = [' '.join(p.split()) for p in re.split(r'\n\s*\n', s) if p.strip()]
    # bullets: one paragraph each, numbered like the template's points
    out = []
    for p in paras:
        if p.startswith('- '):
            for k, item in enumerate(re.split(r'(?:^|\s)- ', p)):
                if item.strip():
                    out.append('%d. %s' % (len([x for x in out if re.match(r'\d+\. ', x)]) + 1,
                                           item.strip()))
        else:
            out.append(p)
    return out


def split_letter(paras):
    """date, addressee lines, greeting, body paragraphs, closing, signature lines."""
    i_greet = next(i for i, p in enumerate(paras) if p.startswith('Dear '))
    i_close = next(i for i, p in enumerate(paras) if p.startswith('Yours sincerely'))
    body = paras[i_greet + 1:i_close]
    sig = paras[i_close + 1:]
    return paras[i_greet], body, paras[i_close], sig


def title_from(body):
    m = re.search(r'"([^"]+)"', ' '.join(body))
    if not m:
        raise SystemExit('cover letter carries no quoted title')
    return m.group(1).rstrip('.')


def set_text(p, text):
    """Replace a paragraph's text, keeping the first run's formatting."""
    runs = p.runs
    if not runs:
        p.add_run(text)
        return
    runs[0].text = text
    for r in runs[1:]:
        r._r.getparent().remove(r._r)


def set_lines(p, lines):
    """Multi-line paragraph (signature block) with soft breaks, formatting from the first run."""
    runs = p.runs
    rpr = copy.deepcopy(runs[0]._r.rPr) if runs and runs[0]._r.rPr is not None else None
    for r in runs:
        r._r.getparent().remove(r._r)
    for k, line in enumerate(lines):
        r = p.add_run()
        if rpr is not None:
            r._r.insert(0, copy.deepcopy(rpr))
        if k:
            r.add_break()
        r.add_text(line)


def clone_after(p, text):
    new = copy.deepcopy(p._p)
    p._p.addnext(new)
    q = docx.text.paragraph.Paragraph(new, p._parent)
    set_text(q, text)
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--target', choices=sorted(TARGETS), default='bib')
    a = ap.parse_args()
    global MD, DOCX, JOURNAL
    MD, DOCX, JOURNAL = TARGETS[a.target]

    greeting, body, closing, sig = split_letter(md_paragraphs())
    title = title_from(body)
    corr = next(x for x in cfg.AUTHORS if x.get('corresponding'))
    af = cfg.AFFILIATIONS[corr['affiliation']]
    authors = ['%s %s' % (x['given'], x['family']) for x in cfg.AUTHORS]
    coauthors = [n for n, x in zip(authors, cfg.AUTHORS) if x is not corr]
    author_line = ', '.join(authors[:-1]) + ' and ' + authors[-1]
    org = af['organisation']
    org_short = org + ' (ETRI)' if 'Electronics and Telecommunications' in org else org
    signature = ['%s %s (corresponding author)' % (corr['given'], corr['family']),
                 af['division'], org_short, '%s, %s' % (af['city'], af['country']), corr['email']]
    behalf = 'on behalf of the co-authors %s and %s' % (', '.join(coauthors[:-1]), coauthors[-1])
    table = {
        'Title:': title,
        'Submitted to:': JOURNAL,
        'Authors:': author_line,
        'Affiliation:': '%s, %s' % (org_short, af['country']),
        'Contact Address:': '%s, %s, %s, %s' % (af['division'], org_short, af['city'], af['country']),
        'Contact Email:': corr['email'],
        'Phone:': PHONE,
        'Fax:': FAX,
    }

    if a.check:
        print(greeting)
        for p in body:
            print('  -', p[:90])
        print(closing)
        print('  |', ' / '.join(signature))
        print('  |', behalf)
        for k, v in table.items():
            print('  %-17s %s' % (k, v))
        return 0

    d = docx.Document(DOCX)
    P = d.paragraphs
    i_greet = next(i for i, p in enumerate(P) if p.text.startswith('Dear '))
    i_thanks = next(i for i, p in enumerate(P) if p.text.startswith('Thank you'))
    i_close = i_thanks + 1
    i_sig, i_behalf = i_close + 1, i_close + 2
    # 'Sincerely' in the untouched template, 'Yours sincerely' once this script has run
    if not (P[i_close].text.lower().endswith('sincerely')
            and P[i_behalf].text.startswith('on behalf')):
        raise SystemExit('template paragraph order changed; refuse to guess')

    set_text(P[i_greet], greeting)
    # body: drop the template's paragraphs, clone the first one for each of ours
    template = P[i_greet + 1]
    for p in P[i_greet + 2:i_thanks]:
        p._p.getparent().remove(p._p)
    set_text(template, body[0])
    prev = template
    for text in body[1:]:
        prev = clone_after(prev, text)
    set_text(P[i_close], closing.rstrip(','))
    set_lines(P[i_sig], signature)
    set_text(P[i_behalf], behalf)

    for row in d.tables[0].rows:
        key = row.cells[0].text.strip()
        if key not in table:
            raise SystemExit('unknown table row %r' % key)
        cell = row.cells[1]
        set_text(cell.paragraphs[0], table[key])
        for extra in cell.paragraphs[1:]:
            extra._p.getparent().remove(extra._p)

    d.save(DOCX)
    print('wrote %s: %d body paragraphs, %d table rows' % (DOCX, len(body), len(table)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
