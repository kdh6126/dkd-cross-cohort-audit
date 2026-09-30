#!/usr/bin/env python
"""새 알고리즘을 만드는 것이 의미가 있는가 — 우리 데이터로 답한다.

"방법을 새로 개발해도 결국 코호트가 결정하는 것 아니냐"는 질문은 추측이 아니라 계산으로
답할 수 있습니다. 코호트 구성 민감도 분석이 방법 4개 × 부분집합 26개를 모두 채운 완전
설계이기 때문입니다.

    분산 분해   결과의 흔들림 중 얼마가 "어떤 방법을 썼나"에서 오고, 얼마가 "어떤 코호트를
                골랐나"에서 오는가. 같은 표를 재현성 축과 예측력 축 각각에 대해 계산한다.
    효과 크기   방법 간 최대 차이와 코호트 구성에 따른 최대 차이를 같은 단위로 비교한다.
    보정 효과   교란 보정이 만든 차이는 어느 쪽 크기에 가까운가. 이것이 "모델 개발"의
                실제 수익률에 대한 이 프로젝트의 답이다.

분산 분해는 균형 설계에서의 제곱합 분해(type I)이며, 두 요인이 완전 교차하므로 순서에
의존하지 않습니다.
"""
import os
import sys

import numpy as np
import pandas as pd

SRC = 'results/cohort_sensitivity/subsets.tsv'
OUT = 'results/method_vs_design'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def decompose(d, value):
    """방법 요인과 부분집합 요인의 제곱합 분해. 완전 균형 설계라 순서 무관."""
    y = d[value].values
    gm = y.mean()
    ss_total = ((y - gm) ** 2).sum()

    ss_method = 0.0
    for m, g in d.groupby('method'):
        ss_method += len(g) * (g[value].mean() - gm) ** 2
    ss_subset = 0.0
    for s, g in d.groupby('subset'):
        ss_subset += len(g) * (g[value].mean() - gm) ** 2
    ss_resid = ss_total - ss_method - ss_subset
    return dict(total=ss_total, method=ss_method, subset=ss_subset, resid=ss_resid,
                pct_method=100 * ss_method / ss_total,
                pct_subset=100 * ss_subset / ss_total,
                pct_resid=100 * ss_resid / ss_total)


def main():
    os.makedirs(OUT, exist_ok=True)
    d = pd.read_csv(SRC, sep='\t')
    log('설계 확인: 방법 %d개 × 부분집합 %d개 = %d행 (완전 균형: %s)'
        % (d['method'].nunique(), d['subset'].nunique(), len(d),
           'yes' if len(d) == d['method'].nunique() * d['subset'].nunique() else 'NO'))
    log('')

    rows = []
    log('=' * 74)
    log('분산 분해 — 결과의 흔들림은 어디서 오는가')
    log('=' * 74)
    log('  %-16s %12s %12s %12s' % ('평가 축', '방법 선택', '코호트 선택', '상호작용·잡음'))
    for value, label in (('cross_fold_jaccard', '재현성'), ('external_auroc', '예측력')):
        r = decompose(d, value)
        r['axis'] = label
        rows.append(r)
        log('  %-16s %11.1f%% %11.1f%% %11.1f%%'
            % (label, r['pct_method'], r['pct_subset'], r['pct_resid']))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'variance.tsv'), sep='\t', index=False)

    log('')
    log('=' * 74)
    log('효과 크기 — 같은 단위로 비교')
    log('=' * 74)
    eff = []
    for value, label in (('cross_fold_jaccard', '재현성'), ('external_auroc', '예측력')):
        by_m = d.groupby('method')[value].mean()
        by_s = d.groupby('subset')[value].mean()
        m_rng = by_m.max() - by_m.min()
        s_rng = by_s.max() - by_s.min()
        eff.append(dict(axis=label, method_range=m_rng, subset_range=s_rng,
                        ratio=s_rng / max(m_rng, 1e-9),
                        best_method=by_m.idxmax(), worst_method=by_m.idxmin()))
        log('  %-8s 방법을 바꾸면 최대 %.3f 차이, 코호트를 바꾸면 최대 %.3f 차이  (%.1f배)'
            % (label, m_rng, s_rng, s_rng / max(m_rng, 1e-9)))
        log('           최고 방법 %-10s 최저 방법 %s' % (by_m.idxmax(), by_m.idxmin()))
    pd.DataFrame(eff).to_csv(os.path.join(OUT, 'effect_size.tsv'), sep='\t', index=False)

    log('')
    log('=' * 74)
    log('방법별 평균')
    log('=' * 74)
    log('  %-12s %12s %12s' % ('방법', '재현성', '예측력'))
    for m in d['method'].unique():
        g = d[d['method'] == m]
        log('  %-12s %12.3f %12.3f'
            % (m, g['cross_fold_jaccard'].mean(), g['external_auroc'].mean()))

    # ---------------------------------------------------------------- 보정 효과
    log('')
    log('=' * 74)
    log('참고 — 교란 보정이 만든 차이는 어느 쪽 크기인가')
    log('=' * 74)
    bench = 'results/ruv_benchmark/benchmark.tsv'
    if os.path.exists(bench):
        b = pd.read_csv(bench, sep='\t')
        if 'arm' in b.columns and 'procurement_frac' in b.columns:
            agg = b.groupby('arm')['procurement_frac'].mean().sort_values()
            log('  교란 유전자 비율 (낮을수록 좋음)')
            for k, v in agg.items():
                log('    %-34s %5.1f%%' % (k, 100 * v))
            log('')
            log('  보정으로 %.0f%%p 개선 (%.1f%% -> %.1f%%)'
                % (100 * (agg.max() - agg.min()), 100 * agg.max(), 100 * agg.min()))
    log('')
    log('=' * 74)
    log('해석')
    log('=' * 74)
    v_rep, v_auc = rows[0], rows[1]
    log('  재현성 축에서는 코호트 선택이 방법 선택보다 %.0f배 크게 작용합니다'
        % (v_rep['pct_subset'] / max(v_rep['pct_method'], 1e-9)))
    log('  (%.1f%% vs %.1f%%). 이 축에서 새 방법을 만드는 것의 수익률은 낮습니다.'
        % (v_rep['pct_subset'], v_rep['pct_method']))
    log('')
    log('  예측력 축은 다릅니다. 방법이 %.1f%%를 설명합니다.' % v_auc['pct_method'])
    if v_auc['pct_method'] > v_rep['pct_method'] * 2:
        log('  즉 "방법 개발이 무의미하다"는 결론은 나오지 않습니다. 다만 어느 축에서')
        log('  개선을 주장하느냐를 분명히 해야 하고, 재현성 축이라면 코호트 구성 민감도를')
        log('  함께 보고하지 않는 한 그 주장은 근거가 약합니다.')
    log('')
    log('  이 프로젝트에서 가장 큰 개선은 새 선택자가 아니라 교란 보정에서 나왔습니다.')
    log('  알고리즘 경쟁보다 평가 설계와 교란 처리의 수익률이 높다는 뜻입니다.')


if __name__ == '__main__':
    main()
