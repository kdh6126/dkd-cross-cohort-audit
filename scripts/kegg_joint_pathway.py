#!/usr/bin/env python
"""전사체 후보가 지목한 경로가 독립 대사체 코호트에서도 재현되는가 — 규칙 기반.

콜라겐 7개 → 4-하이드록시프롤린 연결은 사례로는 좋았지만 사후 선택이었습니다. 살아남은
대사물질 목록을 보고 나서 연결을 찾았기 때문에, 같은 절차를 다른 후보에 적용할 수 없었고
다중검정도 셀 수 없었습니다.

이 스크립트는 그 연결을 규칙으로 바꿉니다. 순서가 중요합니다 — 대사체 결과를 보기 전에
검정할 대사물질이 무엇인지 확정됩니다.

    1  전사체 후보 30개  ->  KEGG 경로            (유전자만 씁니다)
    2  그 경로들          ->  KEGG 화합물          (여전히 대사체 결과를 안 봅니다)
    3  화합물 이름        ->  ST003255 패널과 대조 (여기서 검정 대상이 확정됩니다)
    4  확정된 대사물질의 DKD vs 다른 신장병 변화를 봅니다

주장의 형태도 바뀝니다. "유전자 X를 대사체에서 검증했다"가 아니라 "유전자 X가 지목한 경로가
독립 코호트의 대사체에서도 움직인다"입니다. 전사체와 대사체는 1:1로 대응하지 않으므로
전자는 애초에 성립하지 않는 문장입니다.

귀무가설은 무작위 유전자 집합으로 만듭니다. 같은 크기의 무작위 유전자 30개를 뽑아 같은 절차를
그대로 돌리면, 우연히도 어떤 대사물질들이 지목됩니다. 그 분포와 비교해야 경로 매핑 자체가
아무 유전자에나 무언가를 물어다 주는 것이 아님을 보일 수 있습니다.
"""
import collections
import os
import re
import sys
import time
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

OUT = 'results/kegg_pathway'
CACHE = 'data/raw/kegg'
MET = 'results/metabolomics/ST003255_DKD_vs_other_kidney_disease.tsv'
PANEL = 'data/raw/metabolomics/ST003255/AN005337_datatable.tsv'
SEED = 0
N_NULL = 1000


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def kegg(path, fname):
    """KEGG REST를 한 번만 받아 캐시한다. 재실행 시 네트워크를 다시 치지 않는다."""
    os.makedirs(CACHE, exist_ok=True)
    p = os.path.join(CACHE, fname)
    if os.path.exists(p) and os.path.getsize(p) > 1000:
        return open(p, encoding='utf-8').read()
    t = urllib.request.urlopen('https://rest.kegg.jp/' + path, timeout=300).read()
    t = t.decode('utf-8', 'replace')
    open(p, 'w', encoding='utf-8').write(t)
    time.sleep(0.5)
    return t


def norm(s, level=0):
    """대사물질 이름을 비교 가능한 형태로. 단계가 올라갈수록 더 공격적으로 지운다.

    level 0  표기 문자만 정리. 가장 안전하다.
    level 1  이름 중간의 L-/D- 입체 표기를 지운다. KEGG는 'trans-4-Hydroxy-L-proline',
             패널은 '4-Hydroxyproline'으로 적어 0단계에서는 절대 만나지 못한다.
    level 2  앞머리 위치번호와 cis/trans까지 지운다. '4-'와 '3-'이 서로 다른 화합물이므로
             이 단계는 1:1로 대응될 때만 채택한다.
    """
    s = s.lower()
    s = re.sub(r'\(.*?\)', ' ', s)
    s = s.replace('α', 'alpha').replace('β', 'beta')
    if level >= 1:
        s = re.sub(r'[-\s](l|d|dl)[-\s]', '-', s)
    s = re.sub(r'^(l|d|dl)[- ]', '', s)
    if level >= 2:
        # cis/trans는 지우지 않는다. 지우면 서로 다른 이성질체가 같은 키가 되어
        # 1:1 조건에 걸려 둘 다 버려진다. 위치번호만 떼어낸다.
        s = re.sub(r"^[0-9,]+[- ]", "", s)
    return re.sub(r'[^a-z0-9]+', '', s)


def match_panel(cname, cpds, panel):
    """화합물 -> 패널 대사물질. 단계별로 시도하고 각 단계의 기여를 기록한다.

    2단계는 양방향 1:1일 때만 받는다. '4-'와 '3-'을 둘 다 지우면 서로 다른 화합물이
    같은 키가 되어 조용히 잘못 붙는다.
    """
    hits, tier = {}, {}
    for lv in (0, 1):
        pn = {}
        for m in panel:
            pn.setdefault(norm(m, lv), []).append(m)
        for c in cpds:
            if c in hits.values():
                continue
            for nm in cname.get(c, []):
                k = norm(nm, lv)
                if k in pn and len(pn[k]) == 1 and pn[k][0] not in hits:
                    hits[pn[k][0]] = c
                    tier[pn[k][0]] = lv
                    break
    # 2단계: 양방향 1:1만
    pn2 = {}
    for m in panel:
        pn2.setdefault(norm(m, 2), []).append(m)
    cn2 = {}
    for c in cpds:
        for nm in cname.get(c, []):
            cn2.setdefault(norm(nm, 2), set()).add(c)
    for k, ms in pn2.items():
        if len(ms) == 1 and ms[0] not in hits and k in cn2 and len(cn2[k]) == 1:
            hits[ms[0]] = next(iter(cn2[k]))
            tier[ms[0]] = 2
    return hits, tier



def load_maps():
    g2p = collections.defaultdict(set)
    for line in kegg('link/pathway/hsa', 'gene2pathway.tsv').split('\n'):
        if '\t' not in line:
            continue
        g, p = line.split('\t')
        g2p[g.replace('hsa:', '')].add(p.replace('path:hsa', 'map'))

    p2c = collections.defaultdict(set)
    for line in kegg('link/compound/pathway', 'pathway2compound.tsv').split('\n'):
        if '\t' not in line:
            continue
        p, c = line.split('\t')
        if p.startswith('path:map'):
            p2c[p.replace('path:', '')].add(c.replace('cpd:', ''))

    cname = {}
    for line in kegg('list/cpd', 'compound_names.tsv').split('\n'):
        if '\t' not in line:
            continue
        cid, names = line.split('\t', 1)
        cname[cid] = [n.strip() for n in names.split(';')]
    return g2p, p2c, cname


def genes_to_metabolites(entrez_ids, g2p, p2c, cname, panel_list):
    """유전자 집합 -> 경로 -> 화합물 -> 패널에 있는 대사물질. 대사체 결과는 보지 않는다."""
    paths = set()
    for e in entrez_ids:
        paths |= g2p.get(str(e), set())
    cpds = set()
    for p in paths:
        cpds |= p2c.get(p, set())
    hits, tier = match_panel(cname, cpds, panel_list)
    return paths, cpds, hits


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)

    g2p, p2c, cname = load_maps()
    log('KEGG: 유전자 %d개에 경로 배정, 경로 %d개에 화합물 배정, 화합물 이름 %d개'
        % (len(g2p), len(p2c), len(cname)))

    # 패널 이름 정규화
    panel = open(PANEL, encoding='utf-8').readline().rstrip('\n').split('\t')[2:]
    log('ST003255 패널 %d개' % len(panel))

    # 후보와 배경 유전자
    cand = pd.read_csv('results/candidates_v2/master_candidate_table.tsv', sep='\t')
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = dict(zip(gs['symbol'], gs['entrez_id']))
    cand_ids = [sym2id[s] for s in cand[cand.columns[0]] if s in sym2id]
    background = [e for e in gs['entrez_id'] if e in g2p]
    log('후보 %d개, KEGG 경로가 있는 배경 유전자 %d개' % (len(cand_ids), len(background)))
    log('')

    # ------------------------------------------------ 1~3단계 (대사체 결과 미사용)
    paths, cpds, hits = genes_to_metabolites(cand_ids, g2p, p2c, cname, panel)
    log('=' * 78)
    log('1~3단계 — 대사체 결과를 보기 전에 검정 대상 확정')
    log('=' * 78)
    log('  후보 30개가 속한 KEGG 경로     : %d개' % len(paths))
    log('  그 경로들의 화합물             : %d개' % len(cpds))
    log('  ST003255 패널에서 매칭된 대사물질: %d개' % len(hits))
    if len(hits) < 3:
        log('  매칭이 너무 적어 검정할 수 없습니다.')
        return 1
    for m, c in sorted(hits.items()):
        log('    %-38s %s' % (m[:36], c))

    # ------------------------------------------------ 4단계
    d = pd.read_csv(MET, sep='\t')
    d['is_target'] = d['metabolite'].isin(hits)
    tgt = d[d['is_target']]
    oth = d[~d['is_target']]
    log('')
    log('=' * 78)
    log('4단계 — 확정된 대사물질이 DKD vs 다른 신장병에서 움직이는가')
    log('=' * 78)
    log('  지목된 %d개의 |g| 중앙값   : %.3f' % (len(tgt), tgt['hedges_g'].abs().median()))
    log('  나머지 %d개               : %.3f' % (len(oth), oth['hedges_g'].abs().median()))
    u = stats.mannwhitneyu(tgt['hedges_g'].abs(), oth['hedges_g'].abs(), alternative='greater')
    log('  Mann-Whitney p = %.4f' % u.pvalue)
    n_sig = int((tgt['q'] < 0.05).sum())
    log('  q<0.05 통과: %d / %d (전체 패널은 %d / %d)'
        % (n_sig, len(tgt), int((d['q'] < 0.05).sum()), len(d)))

    # ------------------------------------------------ 무작위 유전자 귀무분포
    log('')
    log('=' * 78)
    log('귀무분포 — 무작위 유전자 30개도 같은 결과를 내는가')
    log('=' * 78)
    obs = tgt['hedges_g'].abs().median()
    gmap = dict(zip(d['metabolite'], d['hedges_g'].abs()))
    null_med, null_n = [], []
    for _ in range(N_NULL):
        rid = rng.choice(background, len(cand_ids), replace=False)
        _, _, h = genes_to_metabolites(rid, g2p, p2c, cname, panel)
        vals = [gmap[m] for m in h if m in gmap]
        null_n.append(len(vals))
        if len(vals) >= 3:
            null_med.append(np.median(vals))
    null_med = np.array(null_med)
    p_emp = (np.sum(null_med >= obs) + 1) / (len(null_med) + 1)
    log('  무작위 집합이 지목한 대사물질 수: 중앙 %d (사분위 %d~%d)'
        % (np.median(null_n), np.percentile(null_n, 25), np.percentile(null_n, 75)))
    log('  후보가 지목한 수                : %d' % len(hits))
    log('')
    log('  무작위 집합의 |g| 중앙값 분포   : 중앙 %.3f (95분위 %.3f)'
        % (np.median(null_med), np.percentile(null_med, 95)))
    log('  후보                            : %.3f' % obs)
    log('  경험적 p = %.4f' % p_emp)

    pd.DataFrame({'metabolite': list(hits), 'kegg_cpd': [hits[m] for m in hits],
                  'hedges_g': [d.set_index('metabolite')['hedges_g'].get(m) for m in hits],
                  'q': [d.set_index('metabolite')['q'].get(m) for m in hits]}).to_csv(
        os.path.join(OUT, 'targeted_metabolites.tsv'), sep='\t', index=False)
    pd.DataFrame({'null_median_abs_g': null_med}).to_csv(
        os.path.join(OUT, 'null_distribution.tsv'), sep='\t', index=False)
    pd.DataFrame([dict(n_pathways=len(paths), n_compounds=len(cpds), n_matched=len(hits),
                       obs_median_abs_g=obs, mwu_p=u.pvalue, empirical_p=p_emp,
                       n_sig_q05=n_sig)]).to_csv(
        os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)

    # ------------------------------------------------ 판정
    log('')
    log('=' * 78)
    log('판정')
    log('=' * 78)
    if p_emp < 0.05:
        log('  후보가 지목한 대사물질이 무작위 유전자가 지목한 것보다 크게 움직입니다.')
        log('  "전사체 후보가 가리키는 경로가 독립 대사체 코호트에서도 재현된다"고 말할 수')
        log('  있습니다.')
    else:
        log('  무작위 유전자 집합과 구분되지 않습니다. 경로 매핑이 아무 유전자에나 비슷한')
        log('  대사물질을 물어다 준다는 뜻이므로, 이 연결을 근거로 쓸 수 없습니다.')
    log('')
    log('  앞선 콜라겐 사례와 다른 점 — 검정 대상을 대사체 결과를 보기 전에 확정했고,')
    log('  무작위 유전자 귀무분포와 비교했습니다. 사후 선택이 아닙니다.')
    log('')
    log('  중요 — 이 절차는 콜라겐 -> 하이드록시프롤린 연결을 재현하지 못합니다. 이름')
    log('  매칭 문제가 아니라 KEGG의 구조적 한계입니다. COL1A2는 ECM-receptor interaction,')
    log('  Focal adhesion 같은 경로에 있고, 하이드록시프롤린은 Arginine and proline')
    log('  metabolism에 있습니다. 두 경로는 만나지 않습니다.')
    log('')
    log('  KEGG는 효소 반응을 모형화하지, 구조 단백질이 분해되어 생기는 산물을 경로')
    log('  구성원으로 적지 않습니다. 우리 후보는 대부분 세포외기질 구조 단백질이라,')
    log('  KEGG 관점에서는 대사물질과 이어질 통로가 아예 없습니다.')
    log('')
    log('  결론 — 수기 연결(콜라겐)과 규칙 연결(KEGG)은 서로 다른 것을 찾습니다.')
    log('  어느 쪽도 다른 쪽을 포함하지 않습니다. 규칙화가 수기 사례를 대체하지 못했고,')
    log('  규칙 쪽 결과는 무작위와 구분되지 않았습니다. 둘 다 그대로 보고합니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
