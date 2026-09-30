#!/usr/bin/env python
"""Write submission_config.py into the manuscript and the cover letter.

The placeholders sit in four files: the BMC manuscript and cover letter, and the BiB manuscript and
cover letter. Filling them by hand is how one gets left behind, so this does it from a single
config and reports every substitution it made. It is idempotent: it replaces the generated regions wholesale, so it can be re-run after
correcting a value.

    python scripts/fill_submission.py --check     # report what would change, write nothing
    python scripts/fill_submission.py             # write

Regions are delimited by markers in the source files. Text outside them is never touched.
"""
import argparse
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import submission_config as cfg   # noqa: E402

TEX = os.path.join(ROOT, 'submission', 'dkd-manuscript.tex')
LETTER = os.path.join(ROOT, 'submission', 'COVER_LETTER.md')
# BiB 판도 같은 값을 쓴다. 여기를 빼 두면 BMC 만 채워 놓고 다 됐다고 착각하게 된다.
BIB_TEX = os.path.join(ROOT, 'submission_bib', 'bib-manuscript.tex')
BIB_LETTER = os.path.join(ROOT, 'submission_bib', 'COVER_LETTER_BiB.md')


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ------------------------------------------------------------------ generators

def author_block():
    lines = []
    for a in cfg.AUTHORS:
        star = '*' if a.get('corresponding') else ''
        mail = r'\email{%s}' % a['email'] if a.get('email') else ''
        lines.append(r'\author%s[%d]{\fnm{%s} \sur{%s}}%s'
                     % (star, a['affiliation'], a['given'], a['family'], mail))
    for num, af in sorted(cfg.AFFILIATIONS.items()):
        star = '*' if any(a.get('corresponding') and a['affiliation'] == num
                          for a in cfg.AUTHORS) else ''
        lines.append(r'\affil%s[%d]{\orgdiv{%s}, \orgname{%s}, \orgaddress{\city{%s}, '
                     r'\country{%s}}}'
                     % (star, num, af['division'], af['organisation'], af['city'],
                        af['country']))
    return '\n'.join(lines)


def availability_block():
    r = cfg.REPOSITORY
    return (
        r'\textbf{Processed data.} The harmonised per-cohort expression matrices, the frozen '
        '9,900-gene\nEntrez space, the probe-to-gene mapping, and the sample-level phenotype '
        'tables including the\nprocurement annotation are deposited at \\url{%s} under %s. The '
        'cohort audit is\ndistributed with them as a relational database whose views encode the '
        'exclusion and overlap\nrules stated in Section~\\ref{sec:cohorts}, so that the audit is '
        'enforced by the schema rather\nthan by prose.\n\n'
        r'\textbf{Code.} All analysis code is available at \url{%s} under %s, archived at '
        '\\url{%s}.' % (r['doi'], r['data_licence'], r['url'], r['code_licence'], r['doi'])
    )


def funding_block():
    """Join the funder statements. Each entry supplies its own sentence-fragment wording,
    because Korean funding bodies mandate specific phrasing that differs between them."""
    parts = ['%s (No. %s)' % (f['funder'], f['grant']) for f in cfg.FUNDING]
    if len(parts) == 1:
        joined, noun = parts[0], 'funder'
    else:
        joined = '; by '.join(parts)
        noun = 'funders'
    return ('This research was supported by %s. The %s had no role in the design of the study, '
            'in the collection, analysis and interpretation of data, or in writing the '
            'manuscript.' % (joined, noun))


# CRediT roles rendered as the verb phrase BMC's Authors' contributions section expects.
CREDIT_PHRASE = {
    'conceptualization':   'conceived the study',
    'methodology':         'designed the methodology',
    'software':            'implemented the analysis pipeline',
    'validation':          'designed and ran the validation and control analyses',
    'formal analysis':     'performed the formal analysis',
    'investigation':       'carried out the analyses',
    'resources':          'provided study resources',
    'data curation':       'curated and harmonised the data',
    'writing original':    'wrote the original draft',
    'writing review':      'reviewed and edited the manuscript',
    'visualization':       'produced the figures',
    'supervision':         'supervised the work',
    'project administration': 'administered the project',
    'funding acquisition': 'acquired funding',
}


def _phrase(contribution):
    """Accept either free text or a list of CRediT role keys."""
    if isinstance(contribution, str):
        return contribution
    unknown = [r for r in contribution if r not in CREDIT_PHRASE]
    if unknown:
        raise SystemExit('unknown CRediT role(s): %s\nvalid keys: %s'
                         % (', '.join(unknown), ', '.join(sorted(CREDIT_PHRASE))))
    parts = [CREDIT_PHRASE[r] for r in contribution]
    if len(parts) == 1:
        return parts[0]
    return ', '.join(parts[:-1]) + ' and ' + parts[-1]


def contributions_block():
    """Group authors who share an identical role set.

    Six authors spelled out one by one repeats the same clause list four times and runs to a
    paragraph; BMC's own examples group co-equal contributors instead.
    """
    groups = []
    for a in cfg.AUTHORS:
        roles = tuple(a['contribution']) if not isinstance(a['contribution'], str)             else (a['contribution'],)
        for g in groups:
            if g['roles'] == roles:
                g['who'].append(a['initials'])
                break
        else:
            groups.append(dict(roles=roles, who=[a['initials']]))

    out = []
    for g in groups:
        who = g['who']
        if len(who) == 1:
            subject = who[0]
        elif len(who) == 2:
            subject = '%s and %s' % (who[0], who[1])
        else:
            subject = '%s and %s' % (', '.join(who[:-1]), who[-1])
        out.append('%s %s.' % (subject, _phrase(list(g['roles']))))
    return ' '.join(out) + ' All authors read and approved the final manuscript.'


def acknowledgements_block():
    """KPMP's mandated attribution, plus the people the data actually came from."""
    k = cfg.KPMP
    kpmp = (
        'The results reported here are in whole or in part based upon data generated by the '
        r'Kidney Precision Medicine Project. Accessed %s. \url{%s}. The KPMP study is '
        'supported by the National Institute of Diabetes and Digestive and Kidney Diseases '
        '(NIDDK). We acknowledge the essential contribution of the patient participants who '
        'consented to kidney biopsy for research.' % (k['accessed'], k['url']))
    others = (
        'We also thank the European Renal cDNA Bank consortium and the investigators who '
        'deposited the public datasets analysed here, without which this study would not have '
        'been possible.')
    return kpmp + '\n\n' + others


def letter_signature():
    a = next((x for x in cfg.AUTHORS if x.get('corresponding')), cfg.AUTHORS[0])
    af = cfg.AFFILIATIONS[a['affiliation']]
    orcid = ' · ORCID %s' % a['orcid'] if a.get('orcid') else ''
    return ('%s %s\n%s, %s\n%s, %s\n%s%s\non behalf of all authors'
            % (a['given'], a['family'], af['division'], af['organisation'], af['city'],
               af['country'], a['email'], orcid))


def letter_reviewers():
    revs = cfg.COVER.get('suggested_reviewers') or []
    if not revs:
        return 'We do not request the exclusion of any reviewer.'
    items = '; '.join('%s (%s, %s) - %s' % (r['name'], r['affiliation'], r['email'],
                                            r['expertise']) for r in revs)
    return ('We suggest the following as reviewers with relevant expertise and no conflict with '
            'the authors: %s. We request the exclusion of no reviewer.' % items)


# ------------------------------------------------------------------ region replacement

def replace_region(text, name, new, path):
    """Replace between `%% BEGIN name` / `%% END name`, inserting markers on first run."""
    if path.endswith('.tex'):
        begin, end = '%%%% BEGIN %s' % name, '%%%% END %s' % name
    else:
        begin, end = '<!-- BEGIN %s -->' % name, '<!-- END %s -->' % name
    pat = re.compile(re.escape(begin) + r'.*?' + re.escape(end), re.S)
    block = '%s\n%s\n%s' % (begin, new, end)
    if pat.search(text):
        return pat.sub(lambda _: block, text, count=1), True
    return text, False


REGIONS_TEX = [
    ('AUTHORS', author_block),
    ('AVAILABILITY', availability_block),
    ('FUNDING', funding_block),
    ('CONTRIBUTIONS', contributions_block),
    ('ACKNOWLEDGEMENTS', acknowledgements_block),
]
def bib_availability_block():
    r = cfg.REPOSITORY
    return ('phenotype tables with procurement annotation and the cohort audit database are '
            'deposited at\n' + chr(92) + 'url{%s}, and all code is available at\n'
            % r['doi'] + chr(92) + 'url{%s}.' % r['url'])


REGIONS_BIB_TEX = [
    ('BIBAVAIL', bib_availability_block),
]
REGIONS_BIB_LETTER = [
    ('DATE', lambda: cfg.COVER['date']),
    ('REPO', lambda: 'audit database and code regenerating every figure and table are released '
                     'at %s\n(%s).' % (cfg.REPOSITORY['url'], cfg.REPOSITORY['doi'])),
]
REGIONS_LETTER = [
    ('TITLE', lambda: '**"%s."**' % cfg.PROJECT['title']),
    ('DATE', lambda: cfg.COVER['date']),
    ('REPO', lambda: 'the harmonised matrices, the audit database, and all code are deposited '
                     'at %s (%s)' % (cfg.REPOSITORY['url'], cfg.REPOSITORY['doi'])),
    ('REVIEWERS', letter_reviewers),
    ('SIGNATURE', letter_signature),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true', help='report only, write nothing')
    ap.add_argument('--allow-partial', action='store_true',
                    help='fill what is known and write the remaining placeholders through, '
                         'so the near-final document can be read before the DOI exists')
    args = ap.parse_args()

    left = cfg.placeholders()
    if left and not args.allow_partial:
        log('submission_config.py still has unfilled values: %s' % ', '.join(left))
        log('Fill them, or pass --allow-partial to write the rest through anyway.')
        return 1
    if left:
        log('PARTIAL: these are still placeholders and will appear verbatim in the output:')
        log('    %s' % ', '.join(left))
        log('preflight.py will keep failing until they are filled. This is for reading the')
        log('near-final document, not for submitting it.')
        log('')

    total = 0
    for path, regions in ((TEX, REGIONS_TEX), (LETTER, REGIONS_LETTER),
                          (BIB_TEX, REGIONS_BIB_TEX), (BIB_LETTER, REGIONS_BIB_LETTER)):
        if not os.path.exists(path):
            log('  MISSING %s' % path)
            return 1
        text = open(path, encoding='utf-8').read()
        original = text
        for name, gen in regions:
            text, found = replace_region(text, name, gen(), path)
            if found:
                log('  filled   %-14s in %s' % (name, os.path.basename(path)))
                total += 1
            else:
                log('  NO MARKER %-13s in %s -- add %s ... %s around the block'
                    % (name, os.path.basename(path), 'BEGIN ' + name, 'END ' + name))
        if not args.check and text != original:
            open(path, 'w', encoding='utf-8').write(text)

    log('')
    log('%d region(s) %s' % (total, 'would be filled' if args.check else 'filled'))
    if not args.check:
        log('Now run:  python scripts/preflight.py')
    return 0


if __name__ == '__main__':
    sys.exit(main())
