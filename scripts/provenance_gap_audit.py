#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""읽히지만 아무 단계도 만들지 않는 결과 파일을 찾는다.

이 감사가 생긴 이유가 있습니다. results/candidates_v2/master_candidate_table.tsv 는
원고의 후보 30개가 나오는 표이고 단백체 스크립트 전부가 읽는데, run_all.py 의 어떤
단계도 그것을 산출물로 선언하지 않았습니다. 손으로 한 번 만들어 두고 그대로 쓴 것입니다.

저장소를 새로 복제해 run_all.py 를 돌리면 그 파일이 없으므로 하류가 전부 실패합니다.
검증 스크립트들은 파일이 **있는지**와 값이 맞는지만 보았기 때문에 이것을 잡지 못했습니다.
있는 파일의 출처를 묻지 않았던 것입니다.

그래서 방향을 뒤집어 봅니다.

    스크립트가 여는 results/ 아래 경로를 전부 모은다
    run_all.py 가 산출물로 선언한 경로를 전부 모은다
    읽히는데 선언되지 않은 것이 틈이다

덮임 판정은 **파일 단위**로 합니다. 처음에는 "선언된 경로의 부모 디렉터리에 속하면
덮인 것" 으로 보았는데, 그것이 이 감사를 무력화했습니다. results/comparison_all4.tsv 가
선언되어 있으면 부모가 'results' 이므로, results/ 바로 아래 **모든** 파일이 자동으로
덮인 것이 되어 버립니다. 실제로 그 규칙 때문에 results/comparison_glom3.tsv 가 아무
단계도 만들지 않는데도 통과했습니다. 게다가 선언된 산출물 128개 중 디렉터리 선언은
0개였으므로, 그 예외는 지켜 주는 것 없이 거짓 음성만 만들고 있었습니다.

그래서 확장자가 없는 선언(=진짜 디렉터리 선언)만 접두사로 인정합니다.
"""
import os
import re
import sys

import pandas as pd

OUT = 'results/proteome_meta'
SKIP = ('verify_provenance.py', 'provenance_gap_audit.py', 'make_release.py',
        'preflight.py', 'submission_audit.py')


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def declared_outputs():
    """run_all.py 를 정규식으로 읽지 않고 그대로 불러온다.

    처음에는 Stage(...) 블록을 정규식으로 잘랐는데 0개를 찾았습니다. 선언 형식이
    한 줄짜리와 여러 줄짜리로 섞여 있어서입니다. 모듈을 불러오면 그런 문제가 없습니다.
    """
    sys.path.insert(0, os.getcwd())
    import run_all
    out = set()
    for st in run_all.S:
        for p in st.outputs:
            out.add(p.replace('\\', '/'))
    return out


def read_paths():
    """스크립트가 여는 results/ 경로를, 그 스크립트가 만드는 것인지와 함께 모은다.

    쓰기 판정은 창(window) 검색으로 하면 이웃한 to_csv 에 걸려 거짓 양성이 납니다.
    실제로 그렇게 한 번 틀렸습니다. 그래서 경로 문자열이 쓰기 호출의 인자 위치에
    있는지를 봅니다. 상수에 담아 두고 나중에 쓰는 형태(OUT = '...' 뒤에
    to_csv(OUT))도 흔하므로, 그 상수 이름이 쓰기 호출에 나타나는지도 함께 봅니다.
    """
    got = {}
    W = r"(?:to_csv|to_excel|savefig|to_parquet|np\.savez[a-z_]*)\s*\("
    for f in sorted(os.listdir('scripts')):
        if not f.endswith('.py') or f in SKIP:
            continue
        s = open(os.path.join('scripts', f), encoding='utf-8', errors='replace').read()
        # 이 파일에서 쓰기 호출의 인자로 등장하는 텍스트를 전부 모은다
        wtext = ' '.join(s[m.end():m.end() + 200] for m in re.finditer(W, s))
        wtext += ' ' + ' '.join(s[m.start():m.start() + 160]
                                for m in re.finditer(r"open\s*\([^)]*['\"][wa]", s))
        # argparse 의 --out 기본값은 읽는 경로가 아니라 쓰는 경로다. 이것을
        # 읽기로 세면 results/comparison.tsv 처럼 있지도 않은 파일이 틈으로 잡힌다.
        wtext += ' ' + ' '.join(
            m.group(1) for m in re.finditer(
                r"add_argument\(\s*'--out[a-z-]*'[^)]*default\s*=\s*'(results/[^']+)'", s))
        # 상수 이름 -> 값
        consts = dict(re.findall(r"^([A-Z][A-Z0-9_]*)\s*=\s*'(results/[^']+)'", s, re.M))
        for m in re.finditer(
                r"['\"](results/[A-Za-z0-9_./-]+\.(?:tsv|csv|xlsx|json|npz))['\"]", s):
            p = m.group(1)
            names = [k for k, v in consts.items() if v == p]
            written = (p in wtext) or any(k in wtext for k in names)
            got.setdefault(p, {'readers': set(), 'writers': set()})
            (got[p]['writers'] if written else got[p]['readers']).add(f)
        # 디렉터리 상수 + 파일명 결합형: os.path.join(OUT, 'x.tsv')
        for m in re.finditer(r"os\.path\.join\(\s*([A-Z][A-Z0-9_]*)\s*,\s*'([^']+\.tsv)'", s):
            d = dict(re.findall(r"^([A-Z][A-Z0-9_]*)\s*=\s*'(results/[^']*)'", s, re.M))
            if m.group(1) in d:
                p = d[m.group(1)].rstrip('/') + '/' + m.group(2)
                got.setdefault(p, {'readers': set(), 'writers': set()})
                got[p]['writers'].add(f)
    return got


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)

    dec = declared_outputs()
    # 확장자가 없는 선언만 디렉터리로 본다. 파일 선언의 부모를 넣으면
    # results/ 아래 전부가 덮인 것이 되어 감사가 아무것도 잡지 못한다.
    dec_dirs = {p.rstrip('/') for p in dec if not os.path.splitext(p)[1]}
    used = read_paths()
    # 재현 경로 밖(탐색용)의 스크립트만 읽는 파일은 알려 두되 실패로 세지 않는다.
    # 논문의 결과가 그 파일에 기대지 않기 때문이다.
    sys.path.insert(0, os.getcwd())
    import run_all
    paper = set(run_all.PAPER)
    paper_scripts = {os.path.basename(str(x)) for st in run_all.S
                     if st.name in paper for x in st.cmd if str(x).endswith('.py')}

    rows = []
    for p, v in sorted(used.items()):
        covered = p in dec or any(p.startswith(d + '/') for d in dec_dirs)
        rows.append(dict(path=p, declared=covered, exists=os.path.exists(p),
                         n_readers=len(v['readers']), n_writers=len(v['writers']),
                         readers='|'.join(sorted(v['readers'])),
                         writers='|'.join(sorted(v['writers']))))
    t = pd.DataFrame(rows)
    # 재현 경로 = PAPER 폐쇄집합의 스크립트. 그 밖에서만 읽히면 원고가 기대지
    # 않으므로 알려만 두고 실패로 세지 않는다.
    t['only_exploratory'] = [
        bool(v['readers']) and not (v['readers'] & paper_scripts)
        for _, v in sorted(used.items())]
    gap = t[(~t['declared']) & (t['n_writers'] == 0) & (t['n_readers'] > 0)]
    hard = gap[~gap['only_exploratory']]

    log('=' * 78)
    log('산출 출처 감사')
    log('=' * 78)
    log('  스크립트가 여는 results/ 경로 %d개' % len(t))
    log('  run_all.py 가 선언한 산출물 %d개 (디렉터리 %d개)' % (len(dec), len(dec_dirs)))
    log('')
    log('  읽히지만 어떤 스크립트도 쓰지 않고 단계 선언도 없는 파일: %d개' % len(gap))
    for _, r in gap.iterrows():
        log('    %-52s %s 읽는 곳 %d개%s'
            % (r['path'], '있음' if r['exists'] else '없음', r['n_readers'],
               '  (탐색용 단계만 읽음)' if r['only_exploratory'] else ''))
        log('        %s' % r['readers'][:140])
    if not len(gap):
        log('    없습니다. 읽히는 것은 전부 어딘가가 만듭니다.')
    log('')
    log('  그중 재현 경로에 영향을 주는 것: %d개' % len(hard))

    t.to_csv(os.path.join(OUT, 'provenance_gap.tsv'), sep='\t', index=False)
    log('')
    log('  %s 에 전체 표를 썼습니다.' % os.path.join(OUT, 'provenance_gap.tsv'))
    return 0 if not len(hard) else 1


if __name__ == '__main__':
    sys.exit(main())
