#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""선택 알고리즘 8종 각각에 교란 보정 6종을 걸어 전후를 비교한다.

기존 보정 벤치마크는 선택기 하나(relieff)에 대해서만 돌렸다. 그러면 "보정이 듣는다" 는
말이 그 선택기에만 해당하는지 알 수 없다. 여덟 개 전부에 같은 보정을 걸어 보면, 보정
효과가 선택기에 의존하는지 아닌지가 드러난다.

부트스트랩 수는 낮춘다. 여기서 묻는 것은 개별 수치의 정밀도가 아니라 보정 전후의
방향이고, 여덟 번을 돌려야 하기 때문이다. B 는 결과 파일에 함께 적는다.
"""
import os
import subprocess
import sys

import pandas as pd

SELECTORS = ['relieff', 'univariate', 'lasso', 'elastic_net', 'mrmr', 'rf', 'boruta']
OUT = 'results/correction_sweep'
B = '40'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for sel in SELECTORS:
        d = os.path.join(OUT, sel)
        f = os.path.join(d, 'benchmark.tsv')
        if not os.path.exists(f):
            log('== %s 실행' % sel)
            # capture_output=True 로 파이프를 잡았더니 자식이 192분 동안 CPU 2초만 쓰고
            # 멈춰 있었다. 자식의 출력은 파일로 보낸다 — 파이프를 아예 만들지 않으면
            # 버퍼가 찰 일이 없고, 로그도 나중에 그대로 읽을 수 있다.
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, 'run.log'), 'w', encoding='utf-8') as fh:
                r = subprocess.run([sys.executable, 'scripts/ruv_benchmark.py',
                                    '--base', sel, '--B', B, '--out', d],
                                   cwd=root, stdout=fh, stderr=subprocess.STDOUT)
            if r.returncode != 0:
                log('   실패 (자세한 내용은 %s/run.log)' % d)
                continue
        t = pd.read_csv(f, sep='\t')
        t = t[t['K'] == 50]
        g = t.groupby('arm').agg(auroc=('external_auroc', 'mean'),
                                 ieg=('n_ieg', 'mean'),
                                 proc=('procurement_frac', 'mean'),
                                 spec=('spec_pass_frac', 'mean'),
                                 xfold=('cross_fold_jaccard', 'mean')).reset_index()
        g.insert(0, 'selector', sel)
        rows.append(g)
        log('   %s 완료' % sel)
    if not rows:
        return 1
    df = pd.concat(rows, ignore_index=True)
    df['B'] = int(B)
    df.to_csv(os.path.join(OUT, 'sweep.tsv'), sep='\t', index=False)

    log('')
    log('=== 선택기별 보정 전후 (채취 관련 피처 비율) ===')
    p = df.pivot_table(index='selector', columns='arm', values='proc')
    log(p.round(3).to_string())
    log('')
    log('%s/sweep.tsv 에 저장' % OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
