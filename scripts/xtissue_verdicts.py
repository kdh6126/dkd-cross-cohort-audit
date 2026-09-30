#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""사전 등록한 세 예측의 판정을 한 표로 모은다.

판정은 docs/xtissue/PREDICTIONS.md 에 적힌 문장 그대로 한다. 결과를 보고 기준을 바꾸지 않는다.
기준이 너무 강했다고 판단되는 곳은 판정은 그대로 두고 비고에 적는다.
"""
import os

import numpy as np
import pandas as pd

OUT = 'results/xtissue'


def main():
    rows = []
    # ---------------------------------------------------------------- P1
    for tissue in ('colon_uc', 'liver_masld'):
        p = os.path.join(OUT, tissue, 'sweep', 'by_size.tsv')
        if not os.path.exists(p):
            rows.append(dict(prediction='P1', tissue=tissue, verdict='not run'))
            continue
        b = pd.read_csv(p, sep='\t')
        r3 = b[b['size'] == 3].iloc[0]
        spans = bool(r3['jaccard_spans_zero'])
        aur_fewer = int(min(r3['auroc_n_negative'], r3['n'] - r3['auroc_n_negative'])) < \
            int(min(r3['jaccard_n_negative'], r3['n'] - r3['jaccard_n_negative']))
        rows.append(dict(
            prediction='P1', tissue=tissue,
            observed='3-cohort Jaccard diff %+.3f to %+.3f (%d of %d negative); AUROC diff %+.3f '
                     'to %+.3f (%d of %d negative)'
                     % (r3['jaccard_min'], r3['jaccard_max'], r3['jaccard_n_negative'], r3['n'],
                        r3['auroc_min'], r3['auroc_max'], r3['auroc_n_negative'], r3['n']),
            verdict='held' if (spans and aur_fewer) else
                    ('partly held: Jaccard spans zero, AUROC not more stable' if spans else
                     'failed: Jaccard difference does not span zero')))
    # ---------------------------------------------------------------- P2
    v = pd.read_csv(os.path.join(OUT, 'p2_verdict.tsv'), sep='\t')
    for _, r in v.iterrows():
        ok_m = r['matched_all_abs_g_below_0_5']
        ok_a = r['asymmetric_negative_and_larger_than_every_matched']
        rows.append(dict(prediction='P2', tissue=r['tissue'],
                         observed='matched cohorts %d, asymmetric %d' % (r['n_matched'],
                                                                         r['n_asymmetric']),
                         verdict='held' if (ok_m is True and ok_a in (True, None)) else 'failed'))
    # ---------------------------------------------------------------- P3
    counts = {}
    for tissue in ('kidney_dkd', 'liver_masld', 'colon_uc'):
        p = os.path.join(OUT, tissue, 'sweep', 'subsets.tsv')
        if os.path.exists(p):
            s = pd.read_csv(p, sep='\t')
            full = s[s['size'] == s['size'].max()]
            counts[tissue] = dict(full[['method', 'n_ieg_topK']].values)
    if len(counts) == 3:
        k, l, c = (counts[t].get('RBS', np.nan) for t in ('kidney_dkd', 'liver_masld', 'colon_uc'))
        order_ok = k > l > c
        colon_ok = c <= 1
        corr = {t: counts[t].get('RBS', np.nan) - counts[t].get('RBS_orth', np.nan) for t in counts}
        rows.append(dict(prediction='P3', tissue='all',
                         observed='uncorrected RBS top-50 IEG: kidney %.1f, liver %.1f, colon %.1f; '
                                  'reduction by correction: kidney %.1f, liver %.1f, colon %.1f'
                                  % (k, l, c, corr['kidney_dkd'], corr['liver_masld'],
                                     corr['colon_uc']),
                         verdict='held' if (order_ok and colon_ok) else 'failed'))
    else:
        rows.append(dict(prediction='P3', tissue='all', verdict='not run (missing sweeps: %s)'
                         % sorted(set(['kidney_dkd', 'liver_masld', 'colon_uc']) - set(counts))))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'verdicts.tsv'), sep='\t', index=False)
    print(df.to_string(index=False))


if __name__ == '__main__':
    main()
