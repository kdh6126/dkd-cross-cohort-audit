#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""--paper 가 원고를 정말로 처음부터 다시 만드는가.

원고는 "python run_all.py --paper 가 이 원고의 모든 그림과 표를 다시 만든다" 고 적습니다.
그 문장이 참이려면 PAPER 목록이 **원고 산출물의 의존성 폐쇄집합**이어야 합니다. 손으로
유지하는 목록은 단계가 늘 때마다 조용히 어긋나고, 실제로 어긋나 있었습니다.

앞선 감사가 이것을 놓친 이유가 중요합니다. 그 감사는 "쓰는 스크립트가 하나도 없는 파일"
만 실패로 봤습니다. 그런데 원래 결함은 달랐습니다. master_candidate_table.py 라는
**쓰는 스크립트는 있었고**, 그것을 그 경로로 실행하는 **단계가 없었을** 뿐입니다.
스크립트가 아니라 단계가 생산자입니다.

그래서 여기서는 생산자를 Stage.outputs 에서만 읽습니다. 선언이 곧 생산 책임입니다.

    1  파일 -> 그 파일을 outputs 로 선언한 단계.  정확히 하나여야 한다.
    2  단계 -> 그 단계의 스크립트가 읽는 results/ 파일.  정적 분석.
    3  씨앗 = 제출물(그림, 표, 부록, 원고)을 만드는 단계.
    4  씨앗에서 읽기를 따라 올라가며 닫는다.
    5  그 폐쇄집합이 PAPER 안에 전부 들어 있어야 한다.
    6  거꾸로, PAPER 에 있으나 폐쇄집합 밖인 단계도 --paper 가 실제로 실행한다.
       그 단계의 입력에도 생산자가 있어야 한다.

6번이 처음에 빠져 있었고, 그 틈으로 결함이 하나 들어왔습니다. figures 단계는
results/figures/F*.pdf 만 쓰므로 씨앗이 아니고, 그 산출을 읽는 단계도 없어서 폐쇄집합에
닿지 않습니다. 그런데 PAPER 목록에는 들어 있으므로 --paper 가 실행합니다. 그 단계가
읽던 results/comparison_glom3.tsv 는 손으로 만든 파일이어서 생산 단계가 없었고, 새로
복제한 저장소에서는 거기서 죽었을 것입니다. 폐쇄집합만 검사하면 이것이 보이지 않습니다.

디렉터리 단위로 봐주지 않습니다. 한 단계가 디렉터리에 파일 다섯 개를 쓰면서 하나만
선언하면, 나머지 넷은 생산자가 없는 것으로 봅니다. 실제로 그 틈으로 결함이 들어왔습니다.
"""
import os
import re
import sys

SKIP = {'provenance_gap_audit.py', 'paper_closure_audit.py',
        'make_release.py', 'preflight.py', 'fill_submission.py',
        # BiB 판은 BMC 원고에서 파생된다. BMC 원고가 그것에 기대지 않으므로 씨앗이 아니다.
        'bib_build.py',
        # 그림 글자 크기 검사는 제출물을 읽기만 하고 아무것도 만들지 않는다.
        'figure_typography.py'}
# verify_provenance.py 는 원고에 적힌 수치와 결과 파일을 잇는 다리다. 그것이 읽는
# 파일은 곧 원고가 의존하는 파일이므로 반드시 폐쇄집합의 씨앗이어야 한다. 처음에
# 이것을 제외 목록에 넣었더니, 원고의 GWAS 음성 결과를 뒷받침하는 gwas_layer.tsv 의
# 생산 단계가 PAPER 에 없는데도 감사가 통과했다.
BRIDGE = 'verify_provenance.py'
ART = re.compile(r'submission/(Fig\d+\.pdf|table_[a-z_]+\.tex|Additional_file_\d|'
                 r'dkd-manuscript\.tex)')


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def script_of(stage):
    for x in stage.cmd:
        if str(x).endswith('.py'):
            return os.path.basename(str(x))
    return ''


def reads_of(script, cache={}):
    """그 스크립트가 여는 results/ 경로. 쓰기로 여는 것은 뺀다."""
    if script in cache:
        return cache[script]
    path = os.path.join('scripts', script)
    if not os.path.exists(path):
        cache[script] = set()
        return cache[script]
    s = open(path, encoding='utf-8', errors='replace').read()
    W = r"(?:to_csv|to_excel|savefig|to_parquet|np\.savez[a-z_]*)\s*\("
    wtext = ' '.join(s[m.end():m.end() + 200] for m in re.finditer(W, s))
    wtext += ' ' + ' '.join(s[m.start():m.start() + 160]
                            for m in re.finditer(r"open\s*\([^)]*['\"][wa]", s))
    consts = dict(re.findall(r"^([A-Z][A-Z0-9_]*)\s*=\s*'(results/[^']+)'", s, re.M))
    # argparse 의 --out 기본값은 읽는 경로가 아니라 쓰는 경로다.
    wtext += ' ' + ' '.join(
        m.group(1) for m in re.finditer(
            r"add_argument\(\s*'--out[a-z-]*'[^)]*default\s*=\s*'(results/[^']+)'", s))
    out = set()
    for m in re.finditer(r"['\"](results/[A-Za-z0-9_./-]+\.(?:tsv|csv|xlsx|json|npz|pdf))['\"]",
                         s):
        p = m.group(1)
        names = [k for k, v in consts.items() if v == p]
        if p in wtext or any(k in wtext for k in names):
            continue
        out.add(p)
    cache[script] = out
    return out


def inputs_of(stage):
    """그 단계가 읽는 results/ 파일.

    스크립트 본문만 보면 안 된다. 입력을 명령줄로 받는 단계가 있고, 그 경우 스크립트
    안의 기본값은 다른 경로를 가리킨다. 실제로 master_table_v2 가 그렇다 — 스크립트
    기본값은 v1 경로이고, v2 경로는 --lit/--kpmp/--pattern 인자로만 들어온다. 명령줄을
    빼고 세면 폐쇄집합이 엉뚱한 단계를 가리킨다.
    """
    out = set(reads_of(script_of(stage)))
    declared = {str(o).replace(chr(92), '/') for o in stage.outputs}
    # stage_figures 는 그림 경로를 os.path.join 으로 조립하므로 소스만 훑어서는
    # 무엇을 읽는지 보이지 않는다. 실제로 새 그림을 넣었을 때 그 생산 단계가 폐쇄집합에
    # 들어오지 않았다. 이 단계는 results/figures 아래의 선언된 그림을 전부 읽는 것으로
    # 본다. 구조상 참이다 — 원고가 가리키는 그림을 골라 복사하는 단계이기 때문이다.
    if script_of(stage) == 'stage_figures.py':
        out |= {p for p in FIGURE_OUTPUTS if p.endswith('.pdf')}
    # 명령줄에 나오는 results/ 토큰 중, 선언된 산출물이 그 아래에 있으면 그것은
    # 출력 디렉터리 인자(--out results/baselines_all4)이지 입력이 아니다.
    outdirs = {os.path.dirname(d) for d in declared}
    for x in stage.cmd:
        x = str(x).replace(chr(92), '/')
        if not x.startswith('results/') or x in declared:
            continue
        if x.rstrip('/') in outdirs:
            continue
        out.add(x)
    return out


FIGURE_OUTPUTS = set()      # main() 에서 채운다


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    sys.path.insert(0, root)
    import run_all

    stages = {st.name: st for st in run_all.S}
    FIGURE_OUTPUTS.update(p.replace(chr(92), '/') for st in run_all.S for p in st.outputs
                          if 'results/figures/' in p.replace(chr(92), '/'))
    paper = set(run_all.PAPER)

    # ---- 1. 파일마다 선언된 생산자 단계
    producer = {}
    dup = []
    for st in run_all.S:
        for p in st.outputs:
            p = p.replace(chr(92), '/')
            if p in producer:
                dup.append((p, producer[p], st.name))
            producer[p] = st.name

    log('=' * 78)
    log('--paper 가 원고를 처음부터 다시 만드는가')
    log('=' * 78)
    log('  단계 %d개 · PAPER %d개 · 선언된 산출물 %d개'
        % (len(run_all.S), len(paper), len(producer)))
    bad = 0
    if dup:
        log('')
        log('  *** 같은 파일을 두 단계가 선언합니다')
        for p, a, b in dup:
            log('    %-52s %s, %s' % (p, a, b))
        bad += len(dup)

    # 입력이 디렉터리일 때(--baselines results/baselines_all4), 그 안의 파일을 선언한
    # 단계가 생산자다. 확장자 있는 파일만 찾으면 이런 단계 간 연결이 끊긴다.
    def owner_of(f):
        o = producer.get(f)
        if o is not None:
            return o
        if not os.path.splitext(f)[1]:
            under = [producer[q] for q in producer if q.startswith(f.rstrip('/') + '/')]
            if len(set(under)) == 1:
                return under[0]
        return None

    # ---- 3. 씨앗: 제출물을 만드는 단계
    seeds = set()
    missing_bridge = {}
    for st in run_all.S:
        sc = script_of(st)
        if sc in SKIP or not sc:
            continue
        body = open(os.path.join('scripts', sc), encoding='utf-8',
                    errors='replace').read() if os.path.exists(
            os.path.join('scripts', sc)) else ''
        if ART.search(body) or any('submission' in str(o) for o in st.outputs):
            seeds.add(st.name)
    log('')
    log('  제출물을 만드는 단계 %d개: %s' % (len(seeds), ', '.join(sorted(seeds))))

    # 원고 수치가 읽는 파일의 생산자도 씨앗에 넣는다
    bridge_files = reads_of(BRIDGE)
    for f in sorted(bridge_files):
        owner = owner_of(f)
        if owner is None:
            missing_bridge.setdefault(f, set()).add(BRIDGE)
        else:
            seeds.add(owner)
    log('  원고 수치가 읽는 파일 %d개 -> 그 생산 단계도 씨앗에 넣습니다'
        % len(bridge_files))

    # ---- 4. 읽기를 따라 닫는다
    need = set(seeds)
    missing_producer = dict(missing_bridge)
    frontier = list(seeds)
    while frontier:
        nm = frontier.pop()
        st = stages.get(nm)
        if st is None:
            continue
        for f in inputs_of(st):
            owner = owner_of(f)
            if owner is None:
                missing_producer.setdefault(f, set()).add(nm)
                continue
            if owner not in need:
                need.add(owner)
                frontier.append(owner)

    log('')
    log('  원고 산출물의 의존성 폐쇄집합: %d개 단계' % len(need))

    # ---- 5b. PAPER 자체가 의존성에 대해 닫혀 있는가.
    # 씨앗에서 출발한 폐쇄집합만 보면 부족하다. --paper 는 PAPER 전체를 실행하므로,
    # 폐쇄집합 밖에 있는 PAPER 단계(figures 처럼 원고가 읽지 않는 산출을 만드는 단계)의
    # 입력도 생산자가 있어야 하고, 그 생산자 또한 PAPER 안에 있어야 한다. 앞의 조건만
    # 검사했더니 compare_glom3 이 생산자로 존재하되 PAPER 밖이라 --paper 가 여전히
    # 그 파일 없이 figures 에 닿는 상태가 남았다.
    outside = sorted(paper - need)
    log('  PAPER 에 있으나 폐쇄집합 밖인 단계 %d개 (--paper 가 실행하므로 함께 검사)'
        % len(outside))
    pneed = set(paper)
    pfront = list(paper)
    while pfront:
        nm = pfront.pop()
        st = stages.get(nm)
        if st is None:
            continue
        for f in inputs_of(st):
            owner = owner_of(f)
            if owner is None:
                missing_producer.setdefault(f, set()).add(nm)
                continue
            if owner not in pneed:
                pneed.add(owner)
                pfront.append(owner)
    unclosed = sorted(pneed - paper)
    if unclosed:
        log('')
        log('  *** PAPER 단계가 입력으로 기대는데 PAPER 에 없는 단계 %d개' % len(unclosed))
        for g in unclosed:
            log('    %-20s %s' % (g, stages[g].desc[:60]))
        bad += len(unclosed)

    gap = sorted(need - paper)
    log('')
    if gap:
        log('  *** PAPER 에 빠진 단계 %d개' % len(gap))
        for g in gap:
            log('    %-20s %s' % (g, stages[g].desc[:60]))
        bad += len(gap)
    else:
        log('  폐쇄집합이 전부 PAPER 안에 있습니다.')

    if missing_producer:
        log('')
        log('  *** 폐쇄집합이 읽는데 선언된 생산자 단계가 없는 파일 %d개'
            % len(missing_producer))
        for f, who in sorted(missing_producer.items()):
            log('    %-54s <- %s' % (f, ', '.join(sorted(who))))
        bad += len(missing_producer)

    log('')
    if bad:
        log('  *** %d건. --paper 가 원고를 처음부터 다시 만든다고 말할 수 없습니다.' % bad)
    else:
        log('  통과. --paper 는 원고 산출물의 의존성을 전부 포함합니다.')
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
