# -*- coding: utf-8 -*-
"""발표 자료가 쓰는 현황 수치를 결과 파일에서 읽는다.

덱과 Word 문서가 같은 값을 쓰게 하려고 따로 뒀습니다. 손으로 적으면 분석이 바뀔 때
두 자료가 서로 다른 숫자를 말하게 되고, 그 어긋남은 발표장에서 드러납니다.
"""
import os
import sys

import pandas as pd


def facts():
    """발표에 쓸 숫자를 결과 파일과 코드에서 직접 읽는다."""
    f = {}
    g = pd.read_csv('results/gate_characteristics.tsv', sep='\t')
    f['bg'] = 100 * g['background_rate'].iloc[0]
    f['fdr'] = 100 * g['gate_fdr_estimate'].iloc[0]
    s = pd.read_csv('results/sensitivity_no104948/summary.tsv', sep='\t')
    v = s[~s['arm'].str.startswith('main')].iloc[0]
    m = s[s['arm'].str.startswith('main')].iloc[0]
    f['sens_keep'] = int(v['n_overlap'])
    f['sens_p'] = float(v['mwu_p'])
    f['main_p'] = float(m['mwu_p'])
    for key, path in (('g3', 'results/comparison_glom3.tsv'),
                      ('a4', 'results/comparison_all4.tsv')):
        d = pd.read_csv(path, sep='\t')
        d = d[d['K'] == 100]
        f[key + '_rbs'] = float(d[d['method'] == 'rbs']['cross_fold_jaccard'].iloc[0])
        f[key + '_next'] = float(d[d['method'] != 'rbs']['cross_fold_jaccard'].max())
    sys.path.insert(0, os.getcwd())
    import run_all
    f['n_stage'] = len(run_all.S)
    f['n_paper'] = len(run_all.PAPER)
    try:
        import pypdf
        f['pages'] = len(pypdf.PdfReader('submission/dkd-manuscript.pdf').pages)
    except Exception:
        f['pages'] = 0
    return f
