#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""KPMP 공여자 한 명의 메타데이터를 h5ad 에서 뽑는다.

두 가지를 고치려고 만듭니다.

  1  등록 범주에는 DKD 가 없습니다. 그래서 카탈로그가 "KPMP 에 당뇨병성 신장질환
     참여자가 없다" 고 표시했는데, 그것은 틀렸습니다. h5ad 에는 diabetes_history 가
     있고, CKD 로 등록된 공여자 중 상당수가 당뇨 병력이 있습니다. 등록 범주와 실제
     병력은 다른 것입니다.
  2  나이와 성별이 공개 자료에 없다고 적었는데, 그것도 KPMP 에 한해서는 틀렸습니다.
     obs 에 sex 와 development_stage 가 있습니다. 교란 축을 그 둘까지 넓힐 수 있습니다.

h5ad 는 3 GB 입니다. 행렬은 건드리지 않고 obs 의 범주형 열 몇 개만 읽습니다. 세포 하나가
한 행이므로 공여자 단위로 접습니다 — 한 공여자의 모든 세포가 같은 값을 가지는지도
함께 확인해서, 아니면 그 사실을 보고합니다.
"""
import os
import sys

import numpy as np
import pandas as pd

H5 = 'data/raw/kpmp/KPMP_v1.5_snRNA_human_kidney.h5ad'
OUT = 'results/kpmp_overlap/donor_meta.tsv'
COLS = ['donor_id', 'disease_category', 'diabetes_history', 'disease',
        'sex_ontology_term_id', 'development_stage_ontology_term_id']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def read_col(obs, name):
    """범주형이면 codes 를 풀고, 아니면 그대로 읽는다."""
    import h5py
    g = obs[name]
    if isinstance(g, h5py.Group):
        cat = np.array(g['categories']).astype(str)
        return pd.Series(cat[np.array(g['codes'])])
    return pd.Series(np.array(g).astype(str))


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    if not os.path.exists(H5):
        log('%s 가 없습니다. run_all.py --stage download 를 먼저 도세요.' % H5)
        return 1
    import h5py

    with h5py.File(H5, 'r') as f:
        obs = f['obs']
        have = [c for c in COLS if c in obs]
        d = pd.DataFrame({c: read_col(obs, c) for c in have})
    log('세포 %d개에서 열 %d개를 읽었습니다.' % (len(d), len(have)))

    # 한 공여자의 세포가 서로 다른 값을 가지면 접는 것이 거짓이 된다. 확인하고 적는다.
    incons = []
    for c in have:
        if c == 'donor_id':
            continue
        k = d.groupby('donor_id')[c].nunique()
        if (k > 1).any():
            incons.append('%s (%d명)' % (c, int((k > 1).sum())))
    if incons:
        log('  *** 공여자 안에서 값이 갈리는 열: %s' % ', '.join(incons))

    g = d.groupby('donor_id').agg(lambda x: x.iloc[0]).reset_index()
    g['n_cells'] = d.groupby('donor_id').size().values
    g = g.rename(columns={'sex_ontology_term_id': 'sex',
                          'development_stage_ontology_term_id': 'age_stage'})
    # 온톨로지 식별자는 사람이 못 읽는다. 알려진 것만 옮기고 나머지는 그대로 둔다.
    g['sex'] = g['sex'].replace({'PATO:0000384': 'male', 'PATO:0000383': 'female'})
    g['age_stage'] = g['age_stage'].str.replace('HsapDv:', '', regex=False)

    # 카탈로그가 쓸 최종 판정. 등록 범주가 아니라 병력으로 당뇨 여부를 정한다.
    g['is_diabetic'] = g.get('diabetes_history', pd.Series(dtype=str)) \
        .astype(str).str.strip().str.lower().isin(['yes', 'true', '1'])
    g['dkd'] = g['is_diabetic'] & g.get(
        'disease_category', pd.Series(dtype=str)).astype(str).str.upper().eq('CKD')

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    g.to_csv(OUT, sep='\t', index=False)

    log('')
    log('공여자 %d명' % len(g))
    for c in ('disease_category', 'diabetes_history', 'sex'):
        if c in g:
            log('  %-18s %s' % (c, dict(g[c].value_counts())))
    log('')
    log('  당뇨 병력 있음        %d명' % int(g['is_diabetic'].sum()))
    log('  그중 CKD 로 등록      %d명  <- 이들이 실질적인 DKD 다' % int(g['dkd'].sum()))
    log('')
    log('  등록 범주에는 DKD 항목이 없습니다. 범주만 보고 "DKD 참여자가 없다" 고')
    log('  말하면 틀립니다. 병력을 봐야 합니다.')
    log('  %s 에 저장했습니다.' % OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
