"""Single source of truth for everything the submission needs and the analysis cannot supply.

Fill this in once. Then:

    python scripts/fill_submission.py     # writes it into the manuscript and cover letter
    python scripts/preflight.py           # confirms nothing was missed
    python scripts/make_release.py        # regenerates CITATION.cff and .zenodo.json

Editing the manuscript by hand instead means touching thirteen separate places, which is how a
stray placeholder survives to submission. Everything below propagates from here.

This file lives outside submission/ so it can never be uploaded to the journal by accident.
"""

# ---------------------------------------------------------------------------- authors
# Order is the author order on the paper. Exactly one entry must have corresponding=True.
# `initials` is what BMC's Authors' contributions section uses; they must be distinct.
#
# NAMES ARE ROMANISED IN REVISED ROMANIZATION AND SHOULD BE CONFIRMED BY EACH AUTHOR.
# Whatever is published becomes that author's permanent PubMed index entry, and several of these
# names have a common alternative spelling (Dohun/Do-Hun, Jinho/Jin-Ho, Seonghee/Sung-Hee,
# Sangyeop/Sang-Yeop). Correcting one is a one-line change here followed by fill_submission.py.
#
# `email` is optional for non-corresponding authors: the typeset paper carries only the
# corresponding author's address, and BMC collects the rest through the submission system.
# `contribution` accepts EITHER a free-text phrase OR a list of CRediT role keys, which
# fill_submission.py renders into BMC's initials-based prose. The fourteen valid keys, with
# what each one corresponds to in this particular project:
#
#   'conceptualization'      asking whether cross-cohort reproducibility measures the methods
#   'methodology'            the LODO design, the stability selector, the handling-score
#                            correction, the block-size simulation
#   'software'               the 57-script pipeline and its single entry point
#   'data curation'          the cohort audit, harmonisation, the frozen gene space, the
#                            hybrid database
#   'formal analysis'        the benchmark, bootstrap intervals, cohort-composition sweep
#   'investigation'          running the analyses
#   'validation'             null controls, the KPMP confounder-matched control, GSE175759 and
#                            GSE162830 replication, the 32-variant WGCNA check
#   'visualization'          figures P1-P6
#   'resources'              clinical and domain input, nephrology interpretation
#   'supervision'            overseeing the work
#   'project administration' coordinating it
#   'funding acquisition'    the three grants
#   'writing original'       the first draft
#   'writing review'         revising it
#
# Example:  contribution=['methodology', 'software', 'writing original'],
#
# ICMJE, which BMC follows, requires every listed author to meet all four criteria: a
# substantial contribution; drafting or critically revising; final approval; and accountability.
# Every author should therefore carry at least one writing role.
AUTHORS = [
    dict(
        given='Dohun', family='Kim', initials='DK',
        email='dohun@etri.re.kr',
        orcid=None,
        corresponding=True,
        title='senior researcher',
        affiliation=1,
        contribution=['conceptualization', 'methodology', 'software', 'data curation',
                      'formal analysis', 'investigation', 'validation', 'visualization',
                      'resources', 'writing original', 'writing review',],
    ),
    dict(
        given='Minho', family='Bae', initials='MB',
        email='minkkang@etri.re.kr', orcid=None, corresponding=False, affiliation=1,
        title='senior researcher',
        contribution=['methodology', 'data curation', 'writing review',],
    ),
    dict(
        given='Jinho', family='Park', initials='JP',
        email='jinho.park@etri.re.kr', orcid=None, corresponding=False, affiliation=1,
        title='senior researcher',
        contribution=['methodology', 'data curation', 'writing review',],
    ),
    dict(
        given='Hayoung', family='Lee', initials='HL',
        email='underzero11@etri.re.kr', orcid=None, corresponding=False, affiliation=1,
        title='researcher',
        contribution=['methodology', 'data curation', 'writing review',],
    ),
    # 이니셜이 이하영(HL)과 겹치지 않도록 HJL 로 둔다. BMC 기여 문단이 이니셜로 쓰인다.
    dict(
        given='Hyeokjin', family='Lim', initials='HJL',
        email='hyeokjin.lim@etri.re.kr', orcid=None, corresponding=False, affiliation=1,
        title='senior researcher',
        contribution=['methodology', 'data curation', 'writing review',],
    ),
    dict(
        given='Seonghee', family='Lee', initials='SL',
        email='slee0003@etri.re.kr', orcid=None, corresponding=False, affiliation=1,
        title='principal researcher',
        contribution=['funding acquisition', 'project administration', 'supervision',
                      'writing review',],
    ),
]

# Referenced by the `affiliation` number in each author entry.
AFFILIATIONS = {
    1: dict(
        division='Autonomous Intelligence DX Research Section',
        organisation='Electronics and Telecommunications Research Institute',
        city='Daejeon',
        country='Republic of Korea',
    ),
}

# ---------------------------------------------------------------------------- project
PROJECT = dict(
    title='Cross-cohort feature-selection benchmarks in diabetic kidney disease are decided by '
          'cohort composition and a shared procurement confounder',
    # Zenodo/CITATION 이 쓰는 제목. BiB 판 제목이다. 위 title 은 BMC 판·커버레터용.
    release_title='Auditing cross-cohort feature-selection benchmarks: cohort composition and '
                  'tissue-context-dependent confounding across kidney, liver and colon transcriptomes',
    running_head='Cohort composition and a shared procurement confounder decide cross-cohort benchmarks',
    version='1.0.0',
    year='2026',
    keywords=['feature selection', 'benchmarking', 'reproducibility', 'batch effects',
              'pre-analytical variables', 'diabetic kidney disease', 'transcriptomics'],
)

# ---------------------------------------------------------------------------- repository
# Fill after creating the repository and minting the DOI. The manuscript's Availability of
# data and materials statement is built from these.
REPOSITORY = dict(
    url='https://github.com/kdh6126/dkd-cross-cohort-audit',
    doi='https://doi.org/10.5281/zenodo.23051697',
    code_licence='the MIT licence',
    data_licence='CC BY 4.0',
)

# ---------------------------------------------------------------------------- funding
# Each entry supplies its own wording because the funding bodies mandate different phrasing.
# fill_submission.funding_block() joins them into one statement.
FUNDING = [
    dict(
        # NRF's template names the funded programme; without it, this is the standard
        # programme-less form the foundation also accepts.
        funder='the National Research Foundation of Korea (NRF) grant funded by the Korean '
               'government (MSIT)',
        grant='RS-2026-25524613',
    ),
    dict(
        funder='an internal research grant of the Electronics and Telecommunications Research '
               'Institute (ETRI)',
        grant='25YT1400',
    ),
    dict(
        funder='the Culture, Sports and Tourism R\\&D Program through the Korea Creative '
               'Content Agency grant funded by the Ministry of Culture, Sports and Tourism in '
               '2024',
        grant='RS-2024-00398310',
    ),
]

# ---------------------------------------------------------------------------- attribution
# KPMP mandates a specific acknowledgement in any publication using its data, in AMA form, with
# the access date. Omitting it is a data-use violation, not a stylistic choice.
KPMP = dict(
    accessed='August 24, 2026',      # when the atlas and API were retrieved
    url='https://www.kpmp.org',
)

# ---------------------------------------------------------------------------- cover letter
COVER = dict(
    # Bracketed so that, if it is ever written through unfilled, preflight's placeholder
    # pattern still catches it in the outbound file.
    date='[SUBMISSION DATE]',
    # Optional. BMC asks for suggested reviewers; leave the list empty to drop the paragraph.
    suggested_reviewers=[
        # dict(name='', affiliation='', email='', expertise='benchmarking methodology'),
    ],
)


def placeholders():
    """Return the tokens that are still unfilled, so callers can refuse to proceed."""
    import json
    blob = json.dumps([AUTHORS, AFFILIATIONS, PROJECT, REPOSITORY, FUNDING, COVER], default=str)
    tokens = ['GIVEN', 'FAMILY', 'DIVISION', 'DEPARTMENT', 'CITY', 'ORG/REPO', 'XXXXXXX',
              'GRANT NUMBER', 'ROLES TO BE ASSIGNED', 'SUBMISSION DATE']
    return [t for t in tokens if t.strip('"\'') in blob]
