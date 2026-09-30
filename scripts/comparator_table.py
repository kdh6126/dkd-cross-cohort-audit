#!/usr/bin/env python
"""비교군이 정본 WGCNA 와 무엇을 공유하고 어디서 갈라지는지 — 한 표로.

이 표는 원래 부록에만 있었습니다. 그런데 심사자가 가장 먼저 묻는 것이 "당신들이 WGCNA
라고 부르는 것이 WGCNA 인가" 입니다. 그 답을 부록에 두면, 답이 있는데도 없는 것처럼
읽힙니다. 본문으로 올립니다.

부록에도 같은 표가 들어갑니다. 두 곳에 손으로 적으면 언젠가 어긋나므로 여기 한 곳에서
정의하고, 원고용 LaTeX 과 부록용 행을 같이 내보냅니다.

soft-threshold 값은 손으로 적지 않고 conventional_pipeline.wgcna_hub 의 기본 인자에서
읽습니다. 부록에서 격자의 최빈값을 beta 로 착각해 본문과 어긋나는 값을 적은 적이 있습니다.
"""
import inspect
import os
import sys

BS = chr(92)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def used_beta():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from conventional_pipeline import wgcna_hub
    return int(inspect.signature(wgcna_hub).parameters['power'].default)


def rows():
    b = used_beta()
    return [
        ('Similarity measure', 'Pearson correlation', 'Pearson correlation', 'same'),
        ('Adjacency', 'signed or unsigned power', 'unsigned power', 'restricted'),
        ('Soft-threshold power', 'chosen by scale-free fit', 'fixed at $\\beta = %d$' % b,
         'differs'),
        ('Module detection', 'TOM distance, dynamic tree cut',
         'clustering on correlation distance', 'differs'),
        ('Hub score', 'kME, eigengene correlation', 'intramodular connectivity kIN', 'differs'),
        ('Module eigengene', 'first principal component', 'not computed', 'absent'),
    ]


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    L = [BS + 'begin{table}[h]',
         BS + 'caption{What the correlation-module comparator shares with the canonical WGCNA '
         'pipeline, and where it departs. Additional file~1 tests whether any departure changes '
         'a conclusion.}' + BS + 'label{tab:comparator}',
         BS + 'footnotesize',
         BS + 'setlength{' + BS + 'tabcolsep}{4pt}',
         BS + 'begin{tabular}{@{}p{24mm}p{41mm}p{41mm}p{15mm}@{}}',
         BS + 'toprule',
         'Element & Canonical WGCNA & Used here & Status ' + BS * 2,
         BS + 'midrule']
    for r in rows():
        L.append(' & '.join(r) + ' ' + BS * 2)
    L += [BS + 'botrule', BS + 'end{tabular}',
          BS + 'footnotetext{The departures are why the manuscript calls this a '
          'correlation-module comparator rather than WGCNA. We do not claim it reproduces the '
          'reference package.}',
          BS + 'end{table}']
    open('submission/table_comparator.tex', 'w', encoding='utf-8').write(
        (chr(10)).join(L) + chr(10))
    log('  비교군 대응표 %d행, beta = %d' % (len(rows()), used_beta()))
    log('  submission/table_comparator.tex 에 썼습니다. 부록 1 도 같은 rows() 를 씁니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
