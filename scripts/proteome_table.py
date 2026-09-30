#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""본문 표 — 두 단백체 코호트에서 잰 후보.

값은 results/proteome_validation/ 과 results/proteome_ms/ 에서 읽습니다. 손으로
옮겨 적지 않습니다. 코호트마다 재는 후보가 다르므로 한 표 안에 블록 두 개로 냅니다.
"""
import os
import sys

import numpy as np
import pandas as pd

BS = chr(92)
NL = chr(10)
OUT = 'submission/table_proteome.tex'
SRC = [('SOMAscan aptamer panel; 23 DKD, 10 healthy (Mendeley 83k89shdx5)',
        'results/proteome_validation/protein_dkd_vs_control.tsv', 'is_candidate'),
       ('LC-MS/MS of FFPE cortex; 5 DKD, 7 non-diabetic (PRIDE PXD041884)',
        'results/proteome_ms/protein_ms_dkd_vs_control.tsv', 'is_cand')]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    meta = pd.read_csv('results/proteome_meta/summary.tsv', sep='\t').iloc[0]

    L = [BS + 'begin{table}[h]',
         BS + 'caption{Candidates measured at the protein level in two patient-independent '
         'datasets. Neither shares patients with the other or with any transcriptome cohort '
         'used here, and the two use different measurement principles. Transcriptome effect is the '
         'mean across the harmonised cohorts; protein effect and $q$ are from the proteomic '
         'cohort named in each block.}' + BS + 'label{tab:proteome}',
         BS + 'begin{tabular}{@{}lrrrl@{}}', BS + 'toprule',
         'Gene & Transcriptome $g$ & Protein $g$ & Protein $q$ & Direction ' + BS * 2]

    for title, path, flag in SRC:
        d = pd.read_csv(path, sep='\t')
        c = d[d[flag].astype(bool)].sort_values('prot_q')
        L += [BS + 'midrule',
              BS + 'multicolumn{5}{@{}l}{' + BS + 'textit{' + title + '}} ' + BS * 2,
              BS + 'midrule']
        for _, r in c.iterrows():
            agree = 'agrees' if (r['prot_g'] > 0) == (r['rna_g'] > 0) else 'differs'
            q = r['prot_q']
            qs = ('$<0.001$' if q < 1e-3 else '%.3f' % q)
            L.append('%s & %+.2f & %+.2f & %s & %s %s'
                     % (r['gene'], r['rna_g'], r['prot_g'], qs, agree, BS * 2))

    L += [BS + 'botrule', BS + 'end{tabular}',
          BS + 'footnotetext{Neither cohort is genome-wide. The aptamer panel carries 841 of '
          'the 9,900 genes in the frozen space and 7 of the 30 candidates; the mass-spectrometry '
          'cohort quantifies 1,142 and 7. Twelve distinct candidates are measured in total, two '
          '(\\textit{LUM}, \\textit{MMP7}) in both. Absence means the protein was not measured, '
          'not that it failed. Pooled across the two cohorts, candidates reach $q<0.05$ more '
          'often than the proteins that are measured but not candidates (Cochran--Mantel--'
          'Haenszel $p = %.3f$; permutation without replacement $p = %.3f$; permutation with '
          'the background matched on abundance and peptide count $p = %.3f$).}'
          % (meta['cmh_p'], meta['perm_p_sig'], meta['perm_p_matched']),
          BS + 'end{table}']
    open(OUT, 'w', encoding='utf-8').write(NL.join(L) + NL)

    for title, path, flag in SRC:
        d = pd.read_csv(path, sep='\t')
        c = d[d[flag].astype(bool)]
        log('  %-16s 후보 %d개 · q<0.05 %d개 · 방향일치 %d개'
            % (title.split(';')[0][:16], len(c), int((c['prot_q'] < 0.05).sum()),
               int(((c['prot_g'] > 0) == (c['rna_g'] > 0)).sum())))
    log('  %s 에 표를 썼습니다.' % OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
