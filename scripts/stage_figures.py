#!/usr/bin/env python
"""그림을 인용된 자리로 옮기고, 번호를 인용 순서에 맞추고, BMC 이름으로 제출 폴더에 깐다.

세 가지를 한꺼번에 맞춥니다.

    자리  그림 환경을 원고 끝에 몰아 두면 LaTeX 이 전부 마지막 몇 쪽에 밀어 넣습니다.
          독자는 그림 하나를 보려고 여러 쪽을 넘겨야 합니다. 각 그림을 그것을 처음 인용한
          문단 바로 뒤로 옮기면 LaTeX 이 그 근처 쪽에 띄웁니다.
    순서  BMC 는 그림을 본문에 나오는 순서대로 번호 매길 것을 요구합니다. 작성 순서대로
          두면 절을 옮기거나 그림을 끼워 넣을 때 조용히 어긋납니다. 빌드는 성공하므로
          눈으로는 잡히지 않습니다.
    이름  제출 파일은 Fig1.pdf, Fig2.pdf 처럼 번호 이름이어야 합니다. 작업용 이름
          (P3_procurement)은 무엇을 그린 것인지 알려주므로 results/ 에서는 그대로 두고,
          제출 폴더로 옮길 때만 번호로 바꿉니다.

인용 순서와 인용 위치는 원고에서 직접 읽습니다. 목록을 여기 적어두면 그 목록이 또
어긋납니다. 원본 이름은 그림 환경 안에 주석으로 남겨, 두 번 돌려도 같은 결과가 나옵니다.
"""
import os
import re
import shutil
import sys

TEX = 'submission/dkd-manuscript.tex'
SRC = 'results/figures'
DST = 'submission'
BS = '\\'
E = re.escape(BS)
NL = chr(10)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def mask_comments(s):
    """주석을 같은 길이의 공백으로 덮는다.

    주석을 '지운' 문자열에서 위치를 찾아 원본에 쓰면 안 됩니다. 길이가 달라 위치가
    앞으로 밀립니다. 실제로 그렇게 짜서, Methods 에서 인용되는 표 셋이 Background 로
    옮겨졌습니다. 길이를 유지하면 두 문자열의 인덱스가 같습니다.
    """
    out = list(s)
    for m in re.finditer(r'(?m)(?<!' + E + r')%.*$', s):
        for k in range(m.start(), m.end()):
            out[k] = ' '
    return ''.join(out)


def mask_figures(s):
    """그림 환경을 같은 길이의 공백으로 덮는다.

    그림 캡션 안에서 표를 인용하는 일이 있습니다. 개요 그림의 캡션이 본문 표 하나를
    가리키자, 그 표가 Results 에서 Background 의 그림 뒤로 끌려왔습니다. 그림 번호는
    이미 그림 환경을 떼어낸 뒤에 세므로 같은 일이 없는데, 표 번호만 캡션 속 인용을
    본문 인용으로 착각했습니다. 주석과 같은 방식으로 덮어 길이를 유지합니다.
    """
    out = list(s)
    for m in re.finditer('(?s)' + E + r'begin\{figure\}.*?' + E + r'end\{figure\}', s):
        for k in range(m.start(), m.end()):
            out[k] = ' '
    return ''.join(out)


def paragraph_end(s, pos):
    """pos 뒤에서 문단이 끝나는 지점. 표 환경 안이면 그 뒤로 넘긴다.

    표 캡션 안에서 그림을 인용하는 경우가 있고, 그 자리에 그림 환경을 끼우면 표 안에
    그림이 박혀 빌드가 깨진다.
    """
    tables = [(m.start(), m.end()) for m in
              re.finditer('(?s)' + E + r'begin\{table\}.*?' + E + r'end\{table\}', s)]
    while True:
        j = s.find(NL + NL, pos)
        if j < 0:
            return len(s)
        j += 1
        for a, b in tables:
            if a < j < b:
                pos = b
                break
        else:
            return j


def reorder_tables(s):
    """인용 순서와 어긋난 표를 인용된 자리로 옮긴다.

    표는 본문에 직접 써 넣으므로 대개 인용 자리에 있습니다. 그런데 표를 새로 끼워 넣으면
    그 뒤 문단에서 먼저 인용되는 표가 생겨 순서가 어긋납니다. BMC 는 순서 인용을
    요구하고, 빌드는 조용히 성공하므로 눈으로는 잡히지 않습니다.

    그림과 달리 표는 전부 옮기지 않고, 어긋난 것만 옮깁니다. 제자리에 있는 표를 건드릴
    이유가 없습니다.
    """
    body = mask_comments(mask_figures(s))
    cited = []
    for m in re.finditer(E + r'ref\{(tab:[^}]+)\}', body):
        if m.group(1) not in cited:
            cited.append(m.group(1))
    envs = list(re.finditer('(?s)' + E + r'begin\{table\}.*?' + E + r'end\{table\}', s))

    # 순서가 맞아도 옮긴다. 순서만 보던 판에서는, Methods 에서 인용되는 표 셋이
    # Background 안에 놓여 있는데도 순서가 맞다는 이유로 그냥 두었습니다. 표가 어느
    # 절에 있는지는 순서와 별개이고, 인용 문단 뒤가 언제나 맞는 자리입니다.
    # 전부 떼어내고 각자의 인용 문단 뒤에 다시 넣는다
    blocks = {}

    def take(m):
        lab = re.search(E + r'label\{(tab:[^}]+)\}', m.group(0))
        if not lab:
            return m.group(0)
        blocks[lab.group(1)] = m.group(0)
        return ''

    stripped = re.sub('(?s)' + E + r'begin\{table\}.*?' + E + r'end\{table\}' + NL + '?',
                      take, s)
    stripped = re.sub(NL + '{3,}', NL + NL, stripped)
    sbody = mask_comments(mask_figures(stripped))
    moved = []
    for lab in cited:
        if lab not in blocks:
            continue
        r = re.search(E + r'ref\{' + re.escape(lab) + r'\}', sbody)
        if not r:
            continue
        moved.append((paragraph_end(stripped, r.end()), blocks.pop(lab)))
    # 인용되지 않은 표는 원래 자리를 알 수 없으므로 끝에 남긴다
    tail = ''.join(NL + b + NL for b in blocks.values())
    out = stripped
    for pos, blk in sorted(moved, key=lambda x: -x[0]):
        out = out[:pos] + NL + blk + NL + out[pos:]
    n = len(re.findall('(?s)' + E + r'begin\{table\}.*?' + E + r'end\{table\}', out + tail))
    if n != len(envs):
        log('표 재배치 결과가 %d개뿐입니다 (있어야 할 수 %d). 표는 손대지 않습니다.'
            % (n, len(envs)))
        return s, 0
    return out + tail, len(moved)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    s = open(TEX, encoding='utf-8').read()
    s, n_tab = reorder_tables(s)
    if n_tab:
        open(TEX, 'w', encoding='utf-8').write(s)
        log('표 %d개를 인용된 자리로 옮겼습니다.' % n_tab)

    # ---------------------------------------------------------- 1. 떼어낸다
    blocks = {}

    def take(m):
        lab = re.search(E + r'label\{(fig:[^}]+)\}', m.group(0))
        if not lab:
            raise SystemExit('label 없는 그림 환경이 있습니다. 먼저 고치세요.')
        blocks[lab.group(1)] = m.group(0)
        return ''

    stripped = re.sub('(?s)' + E + r'begin\{figure\}.*?' + E + r'end\{figure\}' + NL + '?',
                      take, s)
    stripped = re.sub(r'(?m)^%' + r' source: .*$' + NL, '', stripped)
    stripped = re.sub(NL + '{3,}', NL + NL, stripped)

    body = mask_comments(stripped)
    order = []
    for m in re.finditer(E + r'ref\{(fig:[^}]+)\}', body):
        if m.group(1) not in order:
            order.append(m.group(1))

    missing = [l for l in blocks if l not in order] + [l for l in order if l not in blocks]
    if missing:
        log('선언과 인용이 맞지 않습니다: %s' % ', '.join(missing))
        return 1

    # ---------------------------------------------------------- 2. 자리를 찾는다
    mapping, inserts = [], []
    for k, lab in enumerate(order, 1):
        blk = blocks[lab]
        cur = re.search(E + r'includegraphics\[[^]]*\]\{([^}]+)\}', blk).group(1)
        note = re.search(r'%' + r' source: ([A-Za-z0-9_]+)', blk)
        old = note.group(1) if note else cur
        new = 'Fig%d' % k
        mapping.append((k, lab, old, new))

        blk = re.sub(r'(' + E + r'includegraphics\[[^]]*\]\{)[^}]+\}',
                     lambda m: m.group(1) + new + '}', blk, count=1)
        blk = re.sub(r'(?m)^%' + r' source: .*$' + NL, '', blk)
        # 원본 이름은 환경 안에 남긴다. 밖에 두면 다음 실행에서 블록을 잘라낼 때
        # 딸려오지 않아 이름을 잃는다.
        blk = blk.replace(NL, NL + '%' + ' source: ' + old + NL, 1)
        # [h] 만 주면 그 자리에 못 넣을 때 LaTeX 이 그림을 문서 끝으로 던진다.
        # 치환은 반드시 함수로 준다. re.sub 의 치환 '문자열' 에서 backslash-b 는
        # 백스페이스 문자로 해석되어 begin{figure} 자체가 사라진다.
        blk = re.sub(E + r'begin\{figure\}(\[[^]]*\])?',
                     lambda m: BS + 'begin{figure}[htbp]', blk, count=1)

        r = re.search(E + r'ref\{' + re.escape(lab) + r'\}', body)
        inserts.append((paragraph_end(stripped, r.end()), blk))

    # ---------------------------------------------------------- 3. 끼워 넣는다
    # 뒤에서부터 넣어야 앞쪽 위치가 밀리지 않는다
    out = stripped
    for pos, blk in sorted(inserts, key=lambda x: -x[0]):
        out = out[:pos] + NL + blk + NL + out[pos:]

    # 쓰기 전에 세어 본다. 한 번은 치환에서 begin{figure} 가 통째로 사라진 채로
    # '성공' 을 찍고 원고를 망가뜨렸다. 성공 메시지는 증거가 아니다.
    n_env = len(re.findall('(?s)' + E + r'begin\{figure\}.*?' + E + r'end\{figure\}', out))
    if n_env != len(order):
        log('조립 결과에 그림 환경이 %d개뿐입니다 (있어야 할 수 %d). 원고를 쓰지 않습니다.'
            % (n_env, len(order)))
        return 1
    open(TEX, 'w', encoding='utf-8').write(out)

    # ---------------------------------------------------------- 4. 파일을 깐다
    for f in os.listdir(DST):
        if re.match(r'^(P\d_|Fig\d+\.pdf)', f) and f.endswith('.pdf'):
            os.remove(os.path.join(DST, f))
    staged = []
    for k, lab, old, new in mapping:
        src = os.path.join(SRC, old + '.pdf')
        if not os.path.exists(src):
            log('원본이 없습니다: %s' % src)
            return 1
        shutil.copy2(src, os.path.join(DST, new + '.pdf'))
        staged.append((new, old, lab))

    # 어떤 작업 그림이 어떤 제출 그림이 되었는지 남긴다. 이것이 없으면 나중에 원본을
    # 다시 그려도 제출본이 낡은 채로 남고, 어떤 검사도 알아채지 못한다.
    import csv
    with open(os.path.join(SRC, 'staged_map.tsv'), 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter=chr(9))
        w.writerow(['submission_name', 'source_name', 'label'])
        w.writerows(staged)

    total = out.count(NL) + 1
    log('%-5s %-13s %-24s %-11s %s' % ('번호', '라벨', '작업 이름', '제출 이름', '들어간 행'))
    for (k, lab, old, new), (pos, _) in zip(mapping, inserts):
        log('%-5d %-13s %-24s %-11s %d / %d'
            % (k, lab.split(':', 1)[1], old, new + '.pdf',
               stripped[:pos].count(NL) + 1, total))
    log('')
    log('  각 그림은 처음 인용한 문단 뒤로 갔습니다. 번호는 인용 순서에서 읽었습니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
