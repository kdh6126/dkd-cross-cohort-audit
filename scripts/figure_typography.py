#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""그림 글자가 지면에서 몇 pt 로 찍히는가 — 계산으로 확인한다.

코드에 fontsize=6 이라고 적혀 있어도 최종 크기는 삽입 배율에 달렸다. 그래서 두 값을
곱한다.

    그림 PDF 안의 글자 크기(원본 pt) x (지면에서의 폭 / 그림 원본 폭)

지면 폭은 각 문서의 클래스에 직접 물어본다(\\the\\textwidth). 삽입 폭은 tex 의
includegraphics 옵션에서 읽는다. 본문 두 판과 BiB 보충자료, BMC 부록 5번까지 본다.

**이 검사가 보증하지 않는 것.** 글자 크기만 계산한다. 겹침, 잘림, 범례 가독성, 선 굵기,
해상도는 보지 않으며, 조판된 쪽을 눈으로 확인하는 일을 대신하지 못한다. 기준 6pt 도
이 프로젝트가 정한 가독성 하한이고, 저널 그림 규정 전체를 대신하지 않는다.

    python scripts/figure_typography.py --min 6
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

import pymupdf

B = chr(92)
TARGETS = {
    'bmc': dict(tex='submission/dkd-manuscript.tex', dir='submission',
                cls=r'\documentclass[pdflatex,sn-vancouver-num]{sn-jnl}', cwd='submission'),
    'bib': dict(tex='submission_bib/bib-manuscript.tex', dir='submission_bib',
                cls=r'\documentclass[unnumsec,webpdf,contemporary,large,numbered]'
                    r'{oup-authoring-template}', cwd='submission_bib',
                extra=['submission_bib/figures']),
    # BiB 보충자료: 그림을 원본 이름으로 불러오므로 results/figures 에서도 찾는다
    'bib-supp': dict(tex='submission_bib/bib-supplement.tex', dir='submission_bib',
                     cls=(r'\documentclass[11pt]{article}'
                          r'\usepackage[a4paper,margin=22mm]{geometry}'),
                     cwd='submission_bib', extra=['submission_bib/figures']),
    # BMC 부록 5번(후보 30개 전체 그림)은 같은 article 조판에 25mm 여백이다
    'bmc-af5': dict(tex=None, dir='submission',
                    cls=(r'\documentclass[11pt]{article}'
                         r'\usepackage[a4paper,margin=25mm]{geometry}'),
                    cwd='submission', extra=['results/figures'],
                    manual={'P7_candidates_full': (1.0, 'textwidth', False)}),
}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def widths(cfg):
    """그 클래스에서의 textwidth 와 두 단 전체 폭(pt)을 실제로 물어본다."""
    src = (cfg['cls'] + B + 'begin{document}' + B + 'typeout{TW=' + B + 'the' + B + 'textwidth}'
           + B + 'typeout{LW=' + B + 'the' + B + 'linewidth}' + B + 'end{document}')
    d = tempfile.mkdtemp(prefix='tw_')
    open(os.path.join(d, 'probe.tex'), 'w', encoding='utf-8').write(src)
    p = subprocess.run(['pdflatex', '-interaction=nonstopmode', 'probe.tex'],
                       cwd=cfg['cwd'] if False else d, capture_output=True,
                       env=dict(os.environ, TEXINPUTS=os.path.abspath(cfg['dir']) + os.pathsep))
    out = p.stdout.decode('utf-8', 'replace')
    tw = re.search(r'TW=([\d.]+)pt', out)
    return float(tw.group(1)) if tw else None


def placements(tex):
    """라벨 대신 그림 파일 이름으로 삽입 폭을 읽는다. figure* 는 두 단 전체 폭이다."""
    s = open(tex, encoding='utf-8').read()
    out = {}
    for m in re.finditer(r'(?s)' + re.escape(B) + r'begin\{figure(\*?)\}(.*?)'
                         + re.escape(B) + r'end\{figure\*?\}', s):
        star, body = m.group(1), m.group(2)
        g = re.search(re.escape(B) + r'includegraphics\[([^]]*)\]\{([^}]+)\}', body)
        if not g:
            continue
        opt, name = g.group(1), g.group(2)
        f = re.search(r'width=([\d.]*)' + re.escape(B) + r'(textwidth|columnwidth|linewidth)', opt)
        frac = float(f.group(1)) if f and f.group(1) else 1.0
        base = f.group(2) if f else 'textwidth'
        out[name] = (frac, base, bool(star))
    return out


def min_font(pdf):
    sizes = []
    for page in pymupdf.open(pdf):
        for blk in page.get_text('dict')['blocks']:
            for line in blk.get('lines', []):
                for sp in line.get('spans', []):
                    if sp['text'].strip():
                        sizes.append(sp['size'])
    return (min(sizes), len(sizes)) if sizes else (None, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min', type=float, default=6.0, help='저널 최소 글자 크기(pt)')
    a = ap.parse_args()
    bad = 0
    for name, cfg in TARGETS.items():
        tw = widths(cfg)
        if tw is None:
            log('%s: textwidth 를 읽지 못했습니다' % name)
            bad += 1
            continue
        # 두 단 레이아웃의 figure* 는 textwidth 전체를 쓴다. 한 단은 columnwidth.
        colw = (tw - 20.0) / 2 if name == 'bib' else tw
        log('')
        log('=== %s (textwidth %.1fpt, columnwidth %.1fpt, 기준 %.1fpt) ==='
            % (name, tw, colw, a.min))
        items = cfg.get('manual') or placements(cfg['tex'])
        for fig, (frac, base, star) in sorted(items.items()):
            cands = [os.path.join(d, fig + '.pdf')
                     for d in [cfg['dir']] + list(cfg.get('extra', []))]
            p = next((c for c in cands if os.path.exists(c)), None)
            if p is None:
                log('  %-22s 파일 없음' % fig)
                bad += 1
                continue
            native_w = pymupdf.open(p)[0].rect.width
            target = frac * (tw if (star or base == 'textwidth') else colw)
            scale = target / native_w
            mn, n = min_font(p)
            eff = mn * scale if mn else None
            # 부동소수 비교로 6.0 이 미달로 잡히던 것을 막는다
            flag = '' if (eff is None or eff >= a.min - 0.05) else '  *** 기준 미달'
            log('  %-22s 원본 %.0fpt 폭, 지면 %.0fpt (배율 %.2f), 최소 글자 %.1f -> %.1fpt%s'
                % (fig, native_w, target, scale, mn or 0, eff or 0, flag))
            if eff is not None and eff < a.min - 0.05:
                bad += 1
    log('')
    log('기준 미달 %d건 (기준 %.1fpt 는 이 프로젝트가 정한 가독성 하한이다)' % (bad, a.min))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
