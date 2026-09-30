#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""후보 30개 목록이 단백체 결과와 문헌 등급을 보기 전에 확정되었는가.

본문은 "분석 계획을 첫 코호트에서 정하고 바꾸지 않았다" 고 적습니다. 문장으로는 누구나
그렇게 쓸 수 있으므로 확인 가능한 형태로 남깁니다.

파일 시각은 증거가 되지 못합니다. 다시 돌릴 때마다 바뀝니다. 대신 **의존 방향**을 봅니다.
후보 목록을 만드는 단계가 단백체 결과를 입력으로 받지 않는다면, 그 결과를 보고 목록을
고치는 것은 이 파이프라인에서 불가능합니다.

이 감사의 첫 판은 틀렸습니다. 경로 문자열 주위 240자를 훑어 "쓰기" 를 판정했더니,
옆줄의 다른 to_csv 에 걸려 읽기만 하는 스크립트를 생산자로 지목했습니다. 그리고 진짜
생산자는 찾지 못했습니다. 실제로 그때 드러난 것은 다른 문제였습니다 — 이 파일을 만드는
단계가 run_all.py 에 아예 없었습니다. 손으로 한 번 만들어 두고 쓰고 있었던 것입니다.
그 틈은 master_table_v2 단계를 넣어 막았고, 여기서는 판정을 run_all.py 의 선언에서
직접 읽습니다. 문자열을 훑지 않습니다.
"""
import hashlib
import os
import sys

import pandas as pd

CAND = 'results/candidates_v2/master_candidate_table.tsv'
GATED = 'results/candidates_v2/candidates_gated.tsv'
DOWNSTREAM_PREFIX = ('proteome', 'multiomics')
OUT = 'results/proteome_meta'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 16), b''):
            h.update(b)
    return h.hexdigest()


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, root)
    import run_all

    stages = list(run_all.S)
    names = [st.name for st in stages]
    outs = {st.name: [p.replace(chr(92), '/') for p in st.outputs] for st in stages}
    cmds = {st.name: ' '.join(str(x) for x in st.cmd) for st in stages}

    producer = next((n for n in names if CAND in outs[n]), None)
    log('=' * 78)
    log('후보 목록 동결 감사')
    log('=' * 78)
    n = len(pd.read_csv(CAND, sep='\t'))
    log('  대상      %s' % CAND)
    log('  후보 수   %d' % n)
    log('  SHA-256   %s' % sha(CAND))
    log('  입력      %s  (SHA-256 %s)' % (GATED, sha(GATED)[:32]))
    log('')

    if producer is None:
        log('  *** 이 파일을 산출물로 선언한 단계가 없습니다.')
        log('      저장소를 새로 복제하면 만들어지지 않습니다.')
        verdict = 'undeclared'
        ok_order = False
        reads_down = []
    else:
        log('  만드는 단계   %s' % producer)
        log('    %s' % cmds[producer])
        # 그 단계의 명령줄이 하류 결과 경로를 입력으로 받는가
        reads_down = [p for p in cmds[producer].split()
                      if p.startswith('results/')
                      and any(k in p for k in ('proteome', 'multiomics'))]
        down = [n_ for n_ in names if n_.startswith(DOWNSTREAM_PREFIX)]
        pi = names.index(producer)
        log('')
        log('  하류 단계가 뒤에 오는가')
        ok_order = True
        for d in down:
            after = names.index(d) > pi
            ok_order &= after
            log('    %-18s 위치 %-3d  %s' % (d, names.index(d), '뒤' if after else '앞 ***'))
        verdict = ('circular' if reads_down else
                   'frozen' if ok_order else 'order-violation')

    log('')
    if verdict == 'frozen':
        log('  통과. 후보 목록을 만드는 단계는 단백체 결과를 입력으로 받지 않고,')
        log('  단백체 단계는 전부 그 뒤에 옵니다. 결과를 보고 목록을 고치는 것은')
        log('  이 의존 구조에서 불가능합니다.')
    else:
        log('  *** %s' % verdict)

    pd.DataFrame([dict(target=CAND, n_candidates=n, sha256=sha(CAND),
                       input_gated=GATED, input_sha256=sha(GATED),
                       producer_stage=producer or '',
                       producer_cmd=cmds.get(producer, ''),
                       producer_reads_downstream='|'.join(reads_down),
                       downstream_all_after=ok_order, verdict=verdict)]).to_csv(
        os.path.join(OUT, 'freeze_summary.tsv'), sep='\t', index=False)
    return 0 if verdict == 'frozen' else 1


if __name__ == '__main__':
    sys.exit(main())
