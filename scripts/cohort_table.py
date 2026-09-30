#!/usr/bin/env python
"""본문에 들어갈 "쓴 데이터" 표를 조화 산출물에서 직접 만든다.

원고에는 '제외한 계열' 표는 있는데 '실제로 쓴 코호트' 표가 없었습니다. 순서가 뒤집혀
있습니다. 심사자가 가장 먼저 확인하려는 것은 무엇을 뺐는지가 아니라 무엇을 썼는지입니다.

읽는 곳은 data/processed/harmonized/ 입니다. 요약 파일(harmonized_summary.tsv)을 읽지
않는 이유가 있습니다. 그 파일은 harmonize.py 가 만드는데 GSE294519 는 나중에
add_gse294519.py 로 따로 붙어서 요약에 들어가지 않았습니다. 다섯 발견 코호트 중 하나가
빠진 요약을 표의 근거로 삼으면 표도 같이 틀립니다. 조화된 파일 자체가 유일하게 믿을 수
있는 기록입니다.

이 표가 논증을 대신하기도 합니다. 채취 방식 열을 보면, 모든 발견 코호트에서 환자는
생검이고 대조군은 신절제 또는 생체 공여자입니다. 결과를 읽기 전에 그 사실이 보입니다.
"""
import glob
import os
import sys

import pandas as pd

H = 'data/processed/harmonized'
OUT = 'results/cohort_table.tsv'
TEX = 'submission/table_cohorts.tex'
BS = chr(92)

# 어느 코호트가 어떤 역할인지. 이 목록은 실행 코드에서 가져온다.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_sensitivity import ALL5   # noqa: E402

NICE_PROC = {'biopsy': 'biopsy', 'tumor_nephrectomy': 'nephrec.',
             'living_donor': 'donor', '': 'not stated'}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def overlap_table():
    """남긴 코호트 사이에 같은 환자가 몇 명이나 겹치는가.

    본문은 "GSE104948 과 GSE104954 는 DN 환자 5명을 공유한다"고만 적고 있었습니다.
    그건 DKD 환자에 한한 수이고, 두 계열은 실제로 55명을 공유합니다 — 환자 5, 대조 1,
    그리고 다른 진단 생검 49명입니다. 이 논문이 바로 그 다른 진단 생검을 대조군으로
    쓰므로, 49명을 빼고 말하면 남은 중복을 실제보다 작게 보이게 합니다.
    """
    import itertools
    import sqlite3
    c = sqlite3.connect('db/dkd.sqlite')
    d = pd.read_sql(
        'select st.accession s, sm.subject_id, sm.label, sm.compartment '
        'from sample sm join study st on sm.study_id = st.study_id '
        'join study_audit a on a.study_id = st.study_id '
        'where a.usable = 1 and sm.subject_id is not null', c)
    c.close()
    g = d.groupby('s')['subject_id'].apply(set)
    rows = []
    for a, b in itertools.combinations(sorted(g.index), 2):
        shared = g[a] & g[b]
        if not shared:
            continue
        # label 은 1=DKD, 0=비생검 대조, -1=다른 신장병 생검이다. 한 참여자를 두 구획에서
        # 재면 두 행의 label 이 같으므로 처음 값을 쓴다. max() 를 쓰면 -1 과 0 이 섞였을 때
        # 0 으로 올라가 분류가 뒤집힌다.
        lab = (d[d['subject_id'].isin(shared)].groupby('subject_id')['label']
               .agg(lambda x: x.dropna().iloc[0] if x.notna().any() else float('nan')))
        rows.append(dict(series_a=a, series_b=b, shared=len(shared),
                         dkd_cases=int((lab == 1).sum()),
                         non_biopsy_controls=int((lab == 0).sum()),
                         other_biopsy=int((lab == -1).sum()),
                         unlabelled=int(lab.isna().sum())))
    o = pd.DataFrame(rows)
    o.to_csv('results/subject_overlap.tsv', sep='\t', index=False)

    L = [BS + 'begin{table}[h]',
         BS + 'caption{Subjects appearing in more than one retained series. These are the same '
         'patients profiled in two compartments, not duplicate cohorts.}'
         + BS + 'label{tab:overlap}',
         BS + 'begin{tabular}{@{}llrrrr@{}}', BS + 'toprule',
         'Series A & Series B & Shared & DKD & Non-biopsy & Other ' + BS * 2,
         '        &          & subjects & cases & controls & biopsy ' + BS * 2,
         BS + 'midrule']
    for _, r in o.iterrows():
        L.append(' & '.join([r['series_a'], r['series_b'], str(r['shared']),
                             str(r['dkd_cases']), str(r['non_biopsy_controls']),
                             str(r['other_biopsy'])]) + ' ' + BS * 2)
    L += [BS + 'botrule', BS + 'end{tabular}',
          BS + 'footnotetext{Other biopsy are the non-diabetic kidney-disease biopsies that '
          'Section~' + BS + 'ref{sec:procurement} uses as the procurement-matched comparator. '
          'They are counted here because the leave-one-dataset-out folds are not independent '
          'for those subjects either.}',
          BS + 'end{table}']
    open('submission/table_overlap.tex', 'w', encoding='utf-8').write(
        (chr(10)).join(L) + chr(10))
    log('')
    log('  남은 중복:')
    for _, r in o.iterrows():
        log('    %-10s %-10s %3d명 (환자 %d · 비생검 대조 %d · 다른 진단 생검 %d)'
            % (r['series_a'], r['series_b'], r['shared'], r['dkd_cases'],
               r['non_biopsy_controls'], r['other_biopsy']))


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    rows = []
    for f in sorted(glob.glob(os.path.join(H, '*_pheno.tsv'))):
        name = os.path.basename(f)[:-len('_pheno.tsv')]
        p = pd.read_csv(f, sep='\t', dtype=str)
        p['label'] = pd.to_numeric(p['label'], errors='coerce')
        expr = os.path.join(H, name + '_expr.tsv')
        n_genes = sum(1 for _ in open(expr, encoding='utf-8')) - 1

        case = p[p['label'] == 1]
        ctrl = p[p['label'] == 0]
        # 대조군의 채취 방식은 control_type 에, 환자의 채취 방식은 그것이 비어 있다는
        # 사실에 들어 있다. 환자는 전원 생검이다.
        ct = sorted({str(v) for v in ctrl['control_type'].fillna('') if str(v).strip()})
        rows.append(dict(
            cohort=name,
            role='discovery' if name in ALL5 else 'compartment pair',
            compartment=p['compartment'].dropna().iloc[0] if p['compartment'].notna().any()
            else '',
            platform='/'.join(sorted(set(p['platform'].dropna()))),
            n_case=len(case), n_control=len(ctrl),
            n_other=int((p['label'].isna()).sum()),
            genes=n_genes,
            case_procurement='biopsy',
            control_procurement=', '.join(NICE_PROC.get(x, x) for x in ct) or 'not stated'))

    d = pd.DataFrame(rows).sort_values(['role', 'cohort'], ascending=[False, True])
    d.to_csv(OUT, sep='\t', index=False)

    # ---------------------------------------------------------------- LaTeX
    # 전체 폭이 본문 단을 39mm 넘겼다. 글자를 줄이는 것만으로는 모자라서 구획과 역할을
    # 줄임말로 쓰고 각주에서 푼다. 표에서 잘리는 것보다 각주 한 줄이 낫다.
    SHORT_C = {'glomerulus': 'glom.', 'tubulointerstitium': 'tubulo.',
               'whole_cortex': 'cortex'}
    SHORT_R = {'discovery': 'discovery', 'compartment pair': 'compartment'}
    L = [BS + 'begin{table}[h]',
         BS + 'caption{Transcriptome cohorts used. Cases are kidney biopsies in every cohort; '
         'controls are not. This is the design feature Section~' + BS + 'ref{sec:procurement} '
         'tests.}' + BS + 'label{tab:cohorts}',
         BS + 'footnotesize',
         BS + 'setlength{' + BS + 'tabcolsep}{4pt}',
         BS + 'begin{tabular}{@{}llllrrl@{}}',
         BS + 'toprule',
         'Cohort & Role & Compart. & Platform & Cases & Controls & Control tissue ' + BS * 2,
         BS + 'midrule']
    for _, r in d.iterrows():
        L.append(' & '.join([
            r['cohort'],
            SHORT_R.get(r['role'], r['role']),
            SHORT_C.get(r['compartment'], r['compartment'].replace('_', ' ')),
            r['platform'].replace('_', ' '),
            str(r['n_case']),
            str(r['n_control']),
            r['control_procurement']]) + ' ' + BS * 2)
    L += [BS + 'botrule', BS + 'end{tabular}',
          BS + 'footnotetext{Compartments: glom., glomerulus; tubulo., tubulointerstitium; '
          'cortex, whole renal cortex. Control tissue: nephrec., tumour nephrectomy; donor, '
          'living kidney donor. All cohorts are expressed on the same frozen 9,900-gene '
          'Entrez space. Cohorts marked \\emph{compartment} are the second compartment of donors '
          'already represented by an ERCB or Woroniecka series; they enter the procurement and '
          'comparator analyses but not the cohort-composition sweep, which requires independent '
          'cohorts.}',
          BS + 'end{table}']
    open(TEX, 'w', encoding='utf-8').write((chr(10)).join(L) + chr(10))

    log('%-12s %-18s %-20s %6s %8s  %s'
        % ('코호트', '역할', '구획', '환자', '대조군', '대조군 조직'))
    for _, r in d.iterrows():
        log('%-12s %-18s %-20s %6d %8d  %s'
            % (r['cohort'], r['role'], r['compartment'], r['n_case'], r['n_control'],
               r['control_procurement']))
    overlap_table()

    log('')
    log('  코호트 %d개, 전부 9,900 유전자 공간. %s 에 표를 썼습니다.' % (len(d), TEX))
    disc = d[d['role'] == 'discovery']
    log('  발견 코호트 %d개 — 스윕이 요구하는 %d개와 %s'
        % (len(disc), len(ALL5), '일치' if len(disc) == len(ALL5) else '*** 불일치'))
    return 0 if len(disc) == len(ALL5) else 1


if __name__ == '__main__':
    sys.exit(main())
