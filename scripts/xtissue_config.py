# -*- coding: utf-8 -*-
"""조직 확장 코호트의 정의. 조화·감사·분석 스크립트가 모두 이 한 파일을 읽는다.

각 코호트에 적는 것은 다섯이다.

    platform       GEO 플랫폼. 주석 해석 경로를 정한다.
    label          표본 하나(dict: title, source, chars)를 받아 1(사례), 0(대조), None(제외)을
                   돌려주는 규칙. 메타데이터를 사람이 읽고 정했다(results/xtissue/candidates).
    subject        같은 사람의 여러 표본을 하나로 줄일 때 쓰는 식별자 규칙. 없으면 None.
    procurement    사례 팔과 대조 팔을 각각 어떻게 얻었는가, 그리고 그 근거 문장.
    status         core(구성 스윕에 들어감), p2(코호트 단위 교란 검정에만), excluded(이유와 함께).

채취 코드는 넷 중 하나다.
    needle_biopsy          경피 바늘 생검
    endoscopic_biopsy      내시경 집게 생검
    intraoperative         다른 수술(비만수술 등) 도중의 쐐기·바늘 생검
    surgical_or_donor      절제 표본, 생체 공여자 수술 조직, 적출 조직, 상업용 조직 RNA
두 팔의 코드가 같으면 matched, 다르면 asymmetric 이다.
"""
import re


def src(s):
    return s['source'].lower()


def ch(s, key):
    for k, v in s['chars'].items():
        if k.lower() == key.lower():
            return v.lower()
    return ''


COHORTS = {
    # ------------------------------------------------------------------ colon, UC
    'colon_uc': {
        'GSE75214': dict(
            platform='GPL6244',
            label=lambda s: (1 if 'inflamed colonic mucosa of active uc' in src(s) else
                             0 if 'normal colonic mucosa of control' in src(s) else None),
            subject=None,
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy',
                         'Mucosal biopsies were obtained at endoscopy from the colon (series text).'),
            status='core'),
        'GSE59071': dict(
            platform='GPL6244', label=None, subject=None,
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy', 'series text'),
            status='excluded: colon subset of GSE75214 (identical sample titles and group sizes); '
                   'confirmed or refuted by the value-vector audit'),
        'GSE87466': dict(
            platform='GPL13158',
            label=lambda s: (1 if 'ulcerative colitis' in ch(s, 'disease') else
                             0 if ch(s, 'disease') == 'normal' else None),
            subject=None,
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy',
                         'Mucosal biopsy samples (87 UC and 21 normal) were collected (series text).'),
            status='core'),
        'GSE179285': dict(
            platform='GPL6480',
            label=lambda s: (1 if src(s) == 'uc inflamed sigmoid colon' else
                             0 if src(s).startswith('controls') and 'sigmoid' in src(s) else None),
            subject=lambda s: s['title'].split('_')[0],
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy',
                         'Patients ... or normal healthy controls underwent ileocolonoscopy; biopsies '
                         'placed in RNAlater at the clinical site (series text).'),
            status='core'),
        'GSE11223': dict(
            platform='GPL1708',
            label=lambda s: (1 if src(s) == 'uc inflamed sigmoid colon' else
                             0 if src(s) == 'normal uninflamed sigmoid colon' else None),
            subject=lambda s: ch(s, 'patient'),
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy',
                         'Transcriptional profiling of colon epithelial biopsies from UC patients and '
                         'healthy control donors; indication for procedure recorded (series text).'),
            status='core'),
        'GSE47908': dict(
            platform='GPL570',
            label=lambda s: (1 if ch(s, 'disease state') in ('left-sided coltis', 'pancolitis') else
                             0 if ch(s, 'disease state') == 'control' else None),
            subject=None,
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy',
                         'Colonic mucosal biopsy from patient / control individual (sample source).'),
            status='core'),
        'GSE38713': dict(
            platform='GPL570',
            label=lambda s: (1 if 'active disease (involved mucosa)' in src(s) else
                             0 if 'non-inflammatory control' in src(s) else None),
            subject=None,
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy',
                         'Colonic biopsies from patients ... as well as non-inflammatory controls '
                         '(series text).'),
            status='p2'),
        'GSE53306': dict(
            platform='GPL14951',
            label=lambda s: (1 if ch(s, 'disease state') == 'ulcerative colitis active' else
                             0 if ch(s, 'disease state') == 'healthy control' else None),
            subject=None,
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy',
                         'Single endoscopic pinch biopsies (series text).'),
            status='p2'),
        'GSE16879': dict(
            platform='GPL570',
            # 'uc' in source 는 'mucosal' 에도 걸려 크론병까지 사례로 들어갔다. 원문을 통째로 맞춘다.
            # 대조도 disease == control 로 두면 회장 대조 6개가 섞인다. 대장 대조만 쓴다.
            label=lambda s: (1 if (src(s).startswith('colonic mucosal biopsy from uc ')
                                   and 'before first infliximab' in src(s)) else
                             0 if src(s) == 'colonic mucosal biopsy from control individual'
                             else None),
            subject=None,
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy',
                         'Mucosal biopsies were obtained at endoscopy (series text).'),
            status='p2'),
        'GSE48958': dict(
            platform='GPL6244',
            # 'active disease' 는 'inactive disease' 에도 걸린다. 특성 필드를 정확히 비교한다.
            label=lambda s: (1 if ch(s, 'disease activity') == 'active' else
                             0 if 'control individual' in src(s) else None),
            subject=None,
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy',
                         'Colonic mucosal biopsies were obtained during endoscopy (series text).'),
            status='p2'),
        'GSE13367': dict(
            platform='GPL570',
            label=lambda s: (None if 'colonocytes' in src(s) else
                             1 if 'uc inflamed' in s['title'].lower() else
                             0 if s['title'].lower().endswith('control') else None),
            subject=None,
            procurement=('endoscopic_biopsy', 'endoscopic_biopsy',
                         'Adjacent mucosal colonic biopsies were obtained endoscopically (series text).'),
            status='p2'),
    },

    # ------------------------------------------------------------------ liver, MASH
    'liver_masld': {
        'GSE89632': dict(
            platform='GPL14951',
            label=lambda s: (1 if ch(s, 'diagnosis') == 'nash' else
                             0 if ch(s, 'diagnosis') == 'hc' else None),
            subject=None,
            procurement=('needle_biopsy', 'surgical_or_donor',
                         'patients with NAFLD ... and 24 healthy living liver donors as HC '
                         '(series text).'),
            status='core'),
        'GSE48452': dict(
            platform='GPL11532',
            label=lambda s: (None if 'after surgery' in src(s) else
                             1 if ch(s, 'group') == 'nash' else
                             0 if ch(s, 'group') == 'control' else None),
            subject=None,
            procurement=('needle_or_intraoperative_biopsy', 'surgical_or_donor',
                         'Ahrens et al. 2013 is not open access. Horvath et al. 2014 (GSE61260, same '
                         'Kiel group and platform): "Normal control samples were recruited from samples '
                         'obtained for exclusion of liver malignancy during major oncological surgery"; '
                         'NAFLD samples "obtained percutaneously ... or intraoperatively".'),
            status='core'),
        'GSE61260': dict(
            platform='GPL11532',
            label=lambda s: (1 if ch(s, 'diseasestatus') == 'nash' else
                             0 if ch(s, 'diseasestatus') == 'normal control' else None),
            subject=None,
            procurement=('needle_or_intraoperative_biopsy', 'surgical_or_donor',
                         'Horvath et al. 2014, PMC4217403: "Normal control samples were recruited from '
                         'samples obtained for exclusion of liver malignancy during major oncological '
                         'surgery"; "Liver samples were obtained percutaneously ... or intraoperatively".'),
            status='excluded: the series matrix table is empty and the deposited normalised matrix '
                   'labels columns Sample1-Sample134 with no stated mapping to GSM accessions, so '
                   'case and control labels cannot be attached without guessing; GSE48452 from the '
                   'same group is used instead, and the two are audited for shared samples'),
        'GSE162694': dict(
            platform='GPL21290', counts='GSE162694_raw_counts.csv.gz',
            # NAS 5 이상을 확정 NASH 로 둔다. 섬유화 4기는 이식 적출 간 8명이 섞여 있어 뺀다.
            label=lambda s: (None if ch(s, 'fibrosis stage') == '4' else
                             0 if ch(s, 'fibrosis stage') == 'normal liver histology' else
                             1 if ch(s, 'nas score').isdigit() and int(ch(s, 'nas score')) >= 5
                             else None),
            count_column=lambda s: s['title'].split(' ')[-1],
            subject=None,
            procurement=('intraoperative', 'intraoperative',
                         'PMC8433177: "The majority of subjects (N = 133) underwent bariatric surgery '
                         'and had standard of care wedge liver biopsies performed intra-operatively, 8 '
                         'subjects had NAFLD cirrhosis and underwent liver transplantation" (fibrosis '
                         'stage 4 excluded for that reason).'),
            status='core'),
        'GSE126848': dict(
            platform='GPL18573', counts='GSE126848_Gene_counts_raw.txt.gz',
            label=lambda s: (1 if ch(s, 'disease') == 'nash' else
                             0 if ch(s, 'disease') == 'healthy' else None),
            count_column=None,     # 열 이름이 표본 번호라 제목과 직접 맞지 않는다; 아래 규칙으로 맞춘다
            subject=None,
            procurement=('needle_biopsy', 'needle_biopsy',
                         'Sample source "Liver needle biopsy" for all 57 samples, including healthy '
                         'normal-weight volunteers (sample source field and series text).'),
            status='core'),
        'GSE163211': dict(
            platform='GPL29503', label=None, subject=None,
            procurement=('intraoperative', 'intraoperative', 'PMC8710788'),
            status='excluded: NanoString panel of 795 target genes, not genome-wide'),
        'GSE130970': dict(
            platform='GPL16791', label=None, subject=None,
            procurement=('needle_biopsy', 'needle_biopsy', 'series text'),
            status='excluded: no histologically normal control arm (NAS 0 in 4 samples only)'),
        'GSE83452': dict(
            platform='GPL16686',
            label=lambda s: (1 if src(s) == 'nash liver baseline' else
                             0 if src(s) == 'no nash liver baseline' else None),
            subject=None,
            procurement=('intraoperative', 'intraoperative',
                         'bariatric cohort; liver biopsies at surgery in both arms (series text). '
                         'Controls are obese without NASH, not normal liver: flagged.'),
            status='core'),
        'GSE66676': dict(
            platform='GPL6244',
            label=lambda s: (1 if ch(s, 'histology') in ('borderline nash', 'definite nash') else
                             0 if ch(s, 'histology') == 'no nafld' else None),
            subject=None,
            procurement=('intraoperative', 'intraoperative',
                         'liver biopsies were obtained intra-operatively; Liver wedge biopsy '
                         '(series text and sample type).'),
            status='p2'),
        'GSE63067': dict(
            platform='GPL570',
            label=lambda s: (1 if ch(s, 'disease status') == 'non-alcoholic steatohepatitis' else
                             0 if ch(s, 'disease status') == 'healthy' else None),
            subject=None,
            procurement=('unresolved', 'unresolved', 'not stated in the series text'),
            status='p2'),
        'GSE24807': dict(
            platform='GPL2895',
            label=lambda s: (1 if 'nash' in ch(s, 'disease state') else
                             0 if 'normal' in ch(s, 'disease state') else None),
            subject=None,
            procurement=('needle_biopsy', 'surgical_or_donor',
                         'Patient liver biopsy vs Total RNA purchased from ADMET Technologies, '
                         'extracted from healthy liver tissue (sample source).'),
            status='p2'),
        'GSE17470': dict(
            platform='GPL2895', label=None, subject=None,
            procurement=('needle_biopsy', 'surgical_or_donor', 'sample source'),
            status='excluded: patient titles (P53, P55, P59, ...) repeat those of GSE24807'),
        'GSE37031': dict(
            platform='GPL14877',
            label=lambda s: (1 if 'nash' in ch(s, 'disease state') else
                             0 if ch(s, 'disease state') == 'control' else None),
            subject=None,
            procurement=('unresolved', 'unresolved', 'Liver biopsy sample for both arms (source); '
                                                      'route not stated'),
            status='p2'),
    },
}


def arm_match(proc):
    a, b = proc[0], proc[1]
    if 'unresolved' in (a, b):
        return 'unresolved'
    return 'matched' if a == b else 'asymmetric'
