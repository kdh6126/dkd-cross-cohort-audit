#!/usr/bin/env python
"""Stage the public release tree under release/.

Nothing here publishes anything. It assembles the directory that would be pushed to a
repository and archived to Zenodo, so the contents can be inspected before that decision.

Three rules decide what ships:

    redistributable   Anything we derived ourselves: harmonised matrices, phenotype tables,
                      result tables, figures, code, schema, docs.
    not redistributed Third-party primary data with its own terms or its own canonical home:
                      the KPMP .h5ad atlases (4.6 GB), the GEO series matrices (1.6 GB), the
                      GWAS catalog (664 MB), platform annotations (209 MB). run_all.py fetches
                      all of these, so excluding them costs a reader nothing but bandwidth.
    regenerable       Caches and derived stores that any run rebuilds: expr_cache.npz, the
                      SQLite file, the Parquet mirror. Shipping them invites drift against the
                      code that builds them.

Fill AUTHORS and PROJECT below once; LICENSE, CITATION.cff and .zenodo.json are generated from
them, so the author list cannot disagree with itself across files.
"""
import json
import os
import re
import shutil
import sys

# Author and project metadata come from submission_config.py, the same file the manuscript and
# cover letter are filled from, so CITATION.cff and .zenodo.json cannot disagree with the paper.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import submission_config as cfg   # noqa: E402

AUTHORS = [dict(family=a['family'], given=a['given'], orcid=a.get('orcid'),
                affiliation='%s, %s, %s, %s' % (cfg.AFFILIATIONS[a['affiliation']]['division'],
                                                cfg.AFFILIATIONS[a['affiliation']]['organisation'],
                                                cfg.AFFILIATIONS[a['affiliation']]['city'],
                                                cfg.AFFILIATIONS[a['affiliation']]['country']))
           for a in cfg.AUTHORS]
PROJECT = dict(
    # 릴리스 제목은 목표 저널(BiB) 판을 따른다. PROJECT['title'] 은 BMC 커버레터가 쓴다.
    title=cfg.PROJECT.get('release_title', cfg.PROJECT['title']),
    repo_url=cfg.REPOSITORY['url'],
    doi=cfg.REPOSITORY['doi'],
    year=cfg.PROJECT['year'],
    version=cfg.PROJECT['version'],
    keywords=cfg.PROJECT['keywords'],
    grant='; '.join('%s (No. %s)' % (f['funder'], f['grant']) for f in cfg.FUNDING),
)


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def additional_files():
    """원고가 선언한 Additional file 을 그대로 담는다.

    손으로 적은 목록은 부록이 늘면 반드시 뒤처집니다. 실제로 부록 6이 생겼는데 공개본에는
    1~5만 들어가 있었습니다. 그래서 원고에서 번호를 읽고, 그 번호의 파일을 찾습니다.
    선언은 있는데 파일이 없으면 그것은 알려야 할 결함이므로 그대로 둡니다.
    """
    tex = os.path.join(ROOT, 'submission', 'dkd-manuscript.tex')
    if not os.path.exists(tex):
        return []
    body = open(tex, encoding='utf-8').read()
    nums = sorted({int(n) for n in
                   re.findall(r'textbf\{Additional file (\d+)\.\}', body)})
    out = []
    for n in nums:
        for ext in ('.pdf', '.xlsx'):
            rel = 'submission/Additional_file_%d%s' % (n, ext)
            if os.path.exists(os.path.join(ROOT, rel)):
                out.append((rel, rel))
                break
        else:
            print('원고가 Additional file %d 을 선언했는데 파일이 없습니다.' % n,
                  file=sys.stderr)
    return out
OUT = os.path.join(ROOT, 'release')

# (source, destination) — directories are copied whole, subject to EXCLUDE
INCLUDE = [
    ('run_all.py', 'run_all.py'),
    ('requirements.txt', 'requirements.txt'),
    # scripts/fill_submission.py imports this; shipping the script without it makes the
    # release crash on import rather than merely lack a feature.
    ('submission_config.py', 'submission_config.py'),
    ('README.md', 'README.md'),
    ('REPRODUCE.md', 'REPRODUCE.md'),
    ('scripts', 'scripts'),
    ('docs', 'docs'),
    ('db/schema.sql', 'db/schema.sql'),
    # the whole processed tree: the frozen gene space and probe map live at its top level and
    # were missed when this listed only harmonized/, which left `genespace` unreproducible from
    # a fresh clone and omitted the very file the paper offers as a resource.
    ('data/processed', 'data/processed'),
    ('data/metadata', 'data/metadata'),
    ('results', 'results'),
    # 투고본은 BiB 판 하나다. BMC 형식 확장판(submission/)은 공개본에 넣지 않는다 -
    # 투고하지 않는 원고가 공개되면 혼동을 낳는다. 반려되어 BMC 로 갈 때 다시 넣는다.
    # 커버레터는 아래 EXCLUDE 로 뺀다.
    ('submission_bib', 'submission_bib'),
    # 부록은 아래 additional_files() 가 원고에서 읽어 붙인다. 손으로 적어 두었더니
    # 여섯 번째 부록이 생겼을 때 목록이 그대로 남아, 공개본에만 빠져 있었다.
]

# matched against the path relative to ROOT, with forward slashes
EXCLUDE_SUFFIX = ('.pyc', '.npz', '.pptx')
EXCLUDE_CONTAINS = ('__pycache__', '/figures/.ipynb', 'db/columnar', 'dkd.sqlite',
                    # 커버레터는 편집자에게만 가는 문서다(심사자 추천 등). 공개하지 않는다.
                    'COVER_LETTER', '.aux', '.blg', '.out', '.log',
                    # 옛 마크다운 초안. 투고하지 않는 또 하나의 원고라 공개본에서 뺀다.
                    'docs/manuscript/',
                    # 투고 진행 메모(체크리스트·제목 후보·빌드 메모)는 내부 문서다.
                    'docs/submission-notes/',
                    'GSE142025_expr_matrix.tsv',
                    # KPMP 참여자 한 명이 한 행인 표. 공개 Atlas 에서 얻은 것이라 위법은
                    # 아니지만, 논문이 쓰는 것은 질환별 집계뿐이므로 공개할 이유가 없다.
                    # 같은 폴더의 layer_summary.tsv 가 필요한 것을 전부 담는다.
                    'kpmp_overlap/participant_layers.tsv',
                    # 투고 폴더에서 치운 옛 원고 사본을 모아 둔 자리다. 기록으로는
                    # 남기되 공개 릴리스에 옛 원고가 같이 실릴 이유는 없다.
                    'submission-notes/old/')

MIT = """MIT License

Copyright (c) {year} {holders}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

DATA_LICENCE = """# Licence for the data in this repository

## What this covers

The derived data files distributed here -- the harmonised expression matrices under
`data/processed/harmonized/`, the sample-level phenotype and metadata tables under
`data/metadata/`, and the result tables under `results/` -- are released under the
**Creative Commons Attribution 4.0 International licence (CC BY 4.0)**.

Full text: https://creativecommons.org/licenses/by/4.0/legalcode
Summary:   https://creativecommons.org/licenses/by/4.0/

Attribute by citing the associated publication (see `CITATION.cff`).

## What this does not cover

These files are *derived* from primary data that we did not generate and do not redistribute.
The primary data remain under the terms of their original sources, and those terms govern any
use of them:

- **Gene Expression Omnibus series** (GSE30528, GSE30529, GSE96804, GSE104948, GSE104954,
  GSE142025, GSE294519, GSE162830, GSE175759, GSE142153, GSE1009, GSE111154, GSE20602,
  GSE5406). https://www.ncbi.nlm.nih.gov/geo/
- **Kidney Precision Medicine Project** single-nucleus atlas and Atlas API.
  https://atlas.kpmp.org/ -- subject to KPMP's own data use terms. Only the open-access,
  de-identified portion was used; no controlled-access data was requested or held. The `.h5ad`
  atlas files are **not** included in this repository; `run_all.py --stage download` retrieves
  them from KPMP. Results derived from KPMP are distributed as aggregates: the modality table
  is summarised by diagnosis and participant identifiers are not reproduced.
- **Metabolomics Workbench** studies ST003255, ST001411, ST000691, ST004483, ST004442 and
  ST002145. https://www.metabolomicsworkbench.org/ -- the measurement tables are **not**
  included here; only per-metabolite contrast results computed from them.
- **Mendeley Data 83k89shdx5** -- SOMAscan proteomics of human kidney cortex biopsies
  (23 diabetic kidney disease, 10 healthy), CC BY 4.0.
  https://data.mendeley.com/datasets/83k89shdx5 -- the measurement workbook is **not**
  included here; only the per-protein contrast computed from it.
- **PRIDE / ProteomeXchange PXD041884** -- LC-MS/MS of archival formalin-fixed human kidney
  cortex (5 diabetic, 7 non-diabetic). https://www.ebi.ac.uk/pride/archive/projects/PXD041884
  -- the protein rollup workbook is **not** included here; `run_all.py --stage pride_fetch`
  retrieves it. Only the per-protein contrast is distributed.
- **NHGRI-EBI GWAS Catalog**. https://www.ebi.ac.uk/gwas/
- **Platform annotation files** from GEO and NCBI `gene_info`.

The code in this repository is licensed separately, under the MIT licence in `LICENSE`.
"""

GITIGNORE = """# raw third-party data: fetched by `python run_all.py --stage download`
data/raw/

# regenerable stores and caches
db/dkd.sqlite
db/columnar/
**/expr_cache.npz
__pycache__/
*.pyc

# stage logs
logs/

# LaTeX build products
submission/*.aux
submission/*.bbl
submission/*.blg
submission/*.log
submission/*.out
"""


RELEASE_NOTE = """

---

## What is and is not in this repository

Everything derived here is included: the harmonised expression matrices, the frozen 9,900-gene
Entrez space and probe map, the sample-level phenotype tables with the procurement annotation,
every result table the paper quotes, the figures, the database schema, and all code.

Three categories are deliberately absent.

**Primary third-party data** (`data/raw/`, about 7 GB) is not redistributed. The KPMP
single-nucleus atlases carry their own data use terms, and the GEO series, GWAS catalog and
platform annotations all have canonical homes that should remain the source of truth. Fetch them
with:

```bash
python run_all.py --stage download
```

**Regenerable stores** are omitted so they cannot drift from the code that builds them: the
SQLite database and its Parquet mirror (`python run_all.py --stage database`, 3 minutes) and the
single-nucleus expression caches, which their own stages rebuild.

**The status deck** is an internal progress artefact and is not part of the published work.

As a result, a fresh clone reports five stages as outstanding: `download`, `database`, `deck`,
and the two exploratory stages that declare no output files. Every stage backing a number in the
manuscript is already satisfied by the files shipped here, so
`python run_all.py --paper --dry-run` will show them all as present.
"""

def log(*a):
    print(*a, file=sys.stderr, flush=True)


def excluded(rel):
    rel = rel.replace(os.sep, '/')
    return rel.endswith(EXCLUDE_SUFFIX) or any(c in rel for c in EXCLUDE_CONTAINS)


def copy_tree(src, dst):
    n = bytes_ = 0
    for dirpath, dirnames, filenames in os.walk(src):
        dirnames[:] = [d for d in dirnames if not excluded(os.path.join(dirpath, d))]
        for fn in filenames:
            s = os.path.join(dirpath, fn)
            rel = os.path.relpath(s, src)
            if excluded(s):
                continue
            d = os.path.join(dst, rel)
            os.makedirs(os.path.dirname(d), exist_ok=True)
            shutil.copy2(s, d)
            n += 1
            bytes_ += os.path.getsize(s)
    return n, bytes_


def holders():
    return ', '.join('%s %s' % (a['given'], a['family']) for a in AUTHORS)


def write_citation(path):
    lines = ['cff-version: 1.2.0',
             'message: "If you use this software or data, please cite the associated paper."',
             'title: "%s"' % PROJECT['title'],
             'version: "%s"' % PROJECT['version'],
             'license: MIT',
             'repository-code: "%s"' % PROJECT['repo_url'],
             'doi: "%s"' % PROJECT['doi'].replace('https://doi.org/', ''),
             'keywords:'] + ['  - "%s"' % k for k in PROJECT['keywords']] + ['authors:']
    for a in AUTHORS:
        lines.append('  - family-names: "%s"' % a['family'])
        lines.append('    given-names: "%s"' % a['given'])
        lines.append('    affiliation: "%s"' % a['affiliation'])
        if a.get('orcid'):
            lines.append('    orcid: "https://orcid.org/%s"' % a['orcid'])
    open(path, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')


def write_zenodo(path):
    creators = []
    for a in AUTHORS:
        c = dict(name='%s, %s' % (a['family'], a['given']), affiliation=a['affiliation'])
        if a.get('orcid'):
            c['orcid'] = a['orcid']
        creators.append(c)
    meta = dict(
        title=PROJECT['title'],
        upload_type='software',
        creators=creators,
        keywords=PROJECT['keywords'],
        license='mit',
        version=PROJECT['version'],
        description=(
            'Code and derived data for a study of how much a cross-cohort feature-selection '
            'benchmark depends on which cohorts are in it, using diabetic kidney disease as the '
            'case. Includes harmonised expression matrices on a frozen 9,900-gene Entrez space, '
            'a cohort audit identifying three public series as patient-level duplicates of '
            'others, and a single entry point that regenerates every figure and table in the '
            'associated paper. Primary GEO and KPMP data are not redistributed; the pipeline '
            'fetches them.'),
        grants=[],
        # DOI 와 저장소가 서로를 가리켜야 같은 tagged release 임을 확인할 수 있다.
        related_identifiers=[dict(identifier=PROJECT['repo_url'],
                                  relation='isSupplementTo',
                                  scheme='url')],
        notes='Funded by %s.' % PROJECT['grant'],
    )
    json.dump(meta, open(path, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)


def hygiene(root):
    """공개해서는 안 될 것이 섞여 들어갔는가. 사람 눈으로 세 번 놓친 것들이다.

    셋 다 실제로 있었던 일입니다.

        원자료      남의 데이터를 그대로 재배포하면 GEO/KPMP/Workbench 약관 문제가 된다.
                    파생물만 나가야 한다.
        개인 단위   KPMP 참여자 한 명이 한 행인 표가 저장소와 부록에 둘 다 들어가 있었다.
                    법 문제 이전에, 논문이 쓰지 않는 것을 공개할 이유가 없다.
        비ASCII     한국어 파일명이 결과 폴더에 있었다. GitHub zip 과 일부 도구에서 깨진다.

    출처를 적어두는 것과 실제로 안 넣는 것은 다르므로, 적어둔 대로인지 파일로 확인한다.
    """
    import re as _re
    raw = _re.compile(r'datatable|series_matrix|\.h5ad$|(^|/)raw/', _re.I)
    private = _re.compile(r'participant_layers')
    bad = {'raw data': [], 'participant-level': [], 'non-ASCII path': []}
    for dp, _, fs in os.walk(root):
        for f in fs:
            rel = os.path.relpath(os.path.join(dp, f), root).replace(os.sep, '/')
            if raw.search(rel) and not rel.endswith('.py'):
                bad['raw data'].append(rel)
            if private.search(rel):
                bad['participant-level'].append(rel)
            if any(ord(c) > 127 for c in rel):
                bad['non-ASCII path'].append(rel)
    hits = sum(len(v) for v in bad.values())
    if not hits:
        log('  release hygiene OK - no raw third-party data, no participant-level tables, '
            'no non-ASCII paths.')
        return 0
    log('')
    log('  RELEASE HYGIENE FAILED - %d file(s) must not be published:' % hits)
    for kind, paths in bad.items():
        for p in paths[:6]:
            log('    %-18s %s' % (kind, p))
        if len(paths) > 6:
            log('    %-18s ... and %d more' % (kind, len(paths) - 6))
    log('  Add them to EXCLUDE_CONTAINS, or stop producing them.')
    return 1


def main():
    # Clear the contents rather than the directory itself: on Windows a shell or editor sitting
    # in release/ holds a handle on it, and rmtree(OUT) then fails after having already emptied
    # it, leaving a half-staged tree.
    os.makedirs(OUT, exist_ok=True)
    for name in os.listdir(OUT):
        # release/ is itself the git repository that is pushed to GitHub. Wiping .git here once
        # destroyed the local history and made a later `git add` run in the project root.
        if name == '.git':
            continue
        path = os.path.join(OUT, name)
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path, ignore_errors=True)
        else:
            try:
                os.remove(path)
            except OSError:
                pass

    total_n = total_b = 0
    for src_rel, dst_rel in INCLUDE:   # BMC 부록(additional_files())은 공개본에 넣지 않는다
        src = os.path.join(ROOT, src_rel)
        dst = os.path.join(OUT, dst_rel)
        if not os.path.exists(src):
            log('  (missing, skipped) %s' % src_rel)
            continue
        if os.path.isdir(src):
            n, b = copy_tree(src, dst)
        else:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)
            n, b = 1, os.path.getsize(src)
        total_n += n
        total_b += b
        log('  %-34s %5d files  %7.1f MB' % (dst_rel, n, b / 1048576))

    # the README is copied verbatim from the project root, so the release-specific note is
    # appended here rather than maintained as a second copy that would drift
    rp = os.path.join(OUT, 'README.md')
    if os.path.exists(rp):
        with open(rp, 'a', encoding='utf-8') as fh:
            fh.write(RELEASE_NOTE)

    open(os.path.join(OUT, 'LICENSE'), 'w', encoding='utf-8').write(
        MIT.format(year=PROJECT['year'], holders=holders()))
    open(os.path.join(OUT, 'LICENSE-DATA.md'), 'w', encoding='utf-8').write(DATA_LICENCE)
    open(os.path.join(OUT, '.gitignore'), 'w', encoding='utf-8').write(GITIGNORE)
    write_citation(os.path.join(OUT, 'CITATION.cff'))
    write_zenodo(os.path.join(OUT, '.zenodo.json'))

    log('')
    log('  staged %d files, %.1f MB into release/' % (total_n, total_b / 1048576))

    if hygiene(OUT):
        return 1

    placeholders = cfg.placeholders()
    if placeholders:
        log('')
        log('  NOT READY TO PUBLISH - unfilled placeholders: %s' % ', '.join(placeholders))
        log('  Edit submission_config.py and re-run.')


if __name__ == '__main__':
    sys.exit(main() or 0)
