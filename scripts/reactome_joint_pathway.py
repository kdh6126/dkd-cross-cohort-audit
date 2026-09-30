#!/usr/bin/env python
"""KEGG가 못 이은 것을 Reactome은 잇는가 — 같은 절차, 다른 온톨로지.

KEGG 판(kegg_joint_pathway.py)은 두 가지를 보여줬습니다. 후보가 지목한 대사물질이 무작위
유전자가 지목한 것과 구분되지 않았고(경험적 p = 0.26), 수기로 찾았던 콜라겐 -> 하이드록시
프롤린 연결을 재현하지 못했습니다. 후자는 이름 매칭이 아니라 구조 문제였습니다. KEGG는 효소
반응을 모형화하므로 세포외기질 구조 단백질이 분해되어 생기는 산물을 경로 구성원으로 적지
않습니다.

Reactome은 그 빈틈을 메울 것으로 기대할 만합니다. 'Collagen degradation'(R-HSA-1442490)이
명시적 경로로 있고, COL1A2가 실제로 그 안에 있습니다. 그래서 같은 절차를 온톨로지만 바꿔
돌립니다.

    1  후보 30개  ->  Reactome 경로        (Entrez 기준, NCBI2Reactome)
    2  그 경로들  ->  ChEBI 화합물          (ChEBI2Reactome)
    3  ChEBI 이름 ->  ST003255 패널과 대조  (여기서 검정 대상 확정)
    4  확정된 대사물질의 변화를 봅니다

되면 콜라겐 사례가 수기에서 규칙으로 승격됩니다. 안 되면 "두 온톨로지 모두에서 이 연결은
표현되지 않는다"는 더 강한 음성이 되고, 그 자체가 보고할 값입니다.

이름은 EBI OLS API에서 표제어와 동의어를 함께 받아 캐시합니다.
"""
import collections
import json
import os
import re
import sys
import time
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kegg_joint_pathway import norm, match_panel   # 매칭 규칙을 공유해야 비교가 성립한다

OUT = 'results/reactome_pathway'
RAW = 'data/raw/reactome'
NAMES = os.path.join(RAW, 'chebi_names.json')
MET = 'results/metabolomics/ST003255_DKD_vs_other_kidney_disease.tsv'
PANEL = 'data/raw/metabolomics/ST003255/AN005337_datatable.tsv'
SEED = 0
N_NULL = 1000


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load_maps():
    g2p = collections.defaultdict(set)
    with open(os.path.join(RAW, 'NCBI2Reactome_All_Levels.txt'),
              encoding='utf-8', errors='replace') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 6 and f[5] == 'Homo sapiens':
                g2p[f[0]].add(f[1])
    p2c = collections.defaultdict(set)
    pname = {}
    with open(os.path.join(RAW, 'ChEBI2Reactome_All_Levels.txt'),
              encoding='utf-8', errors='replace') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 6 and f[5] == 'Homo sapiens':
                p2c[f[1]].add(f[0])
                pname[f[1]] = f[3]
    return g2p, p2c, pname


def chebi_names(ids):
    """ChEBI ID -> [표제어, 동의어...]. OLS에서 받아 캐시한다."""
    cache = json.load(open(NAMES, encoding='utf-8')) if os.path.exists(NAMES) else {}
    todo = [i for i in ids if i not in cache]
    if todo:
        log('  ChEBI 이름 %d개 조회 (캐시 %d개)' % (len(todo), len(cache)))
        for k, cid in enumerate(todo, 1):
            u = ('https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms?obo_id=CHEBI:%s'
                 % cid)
            try:
                d = json.load(urllib.request.urlopen(u, timeout=60))
                t = d['_embedded']['terms'][0]
                cache[cid] = [t.get('label') or ''] + list(t.get('synonyms') or [])
            except Exception:
                cache[cid] = []
            if k % 100 == 0:
                log('    %d / %d' % (k, len(todo)))
                json.dump(cache, open(NAMES, 'w', encoding='utf-8'), ensure_ascii=False)
            time.sleep(0.05)
        json.dump(cache, open(NAMES, 'w', encoding='utf-8'), ensure_ascii=False)
    return {i: cache.get(i, []) for i in ids}


def genes_to_metabolites(entrez_ids, g2p, p2c, cname, panel):
    paths = set()
    for e in entrez_ids:
        paths |= g2p.get(str(e), set())
    cpds = set()
    for p in paths:
        cpds |= p2c.get(p, set())
    hits, _ = match_panel(cname, cpds, panel)
    return paths, cpds, hits


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)
    g2p, p2c, pname = load_maps()
    log('Reactome: 유전자 %d개, 화합물을 가진 경로 %d개' % (len(g2p), len(p2c)))

    panel = open(PANEL, encoding='utf-8').readline().rstrip('\n').split('\t')[2:]
    cand = pd.read_csv('results/candidates_v2/master_candidate_table.tsv', sep='\t')
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = dict(zip(gs['symbol'], gs['entrez_id']))
    cand_syms = list(cand[cand.columns[0]])
    cand_ids = [sym2id[s] for s in cand_syms if s in sym2id]
    background = [e for e in gs['entrez_id'] if e in g2p]
    log('후보 %d개, Reactome 경로가 있는 배경 유전자 %d개' % (len(cand_ids), len(background)))

    # 후보와 무작위 집합이 건드리는 화합물 전체를 미리 모아 한 번에 이름을 받는다
    paths_c = set()
    for e in cand_ids:
        paths_c |= g2p.get(str(e), set())
    cpds_c = set()
    for p in paths_c:
        cpds_c |= p2c.get(p, set())
    pool = set(cpds_c)
    probe = [rng.choice(background, len(cand_ids), replace=False) for _ in range(40)]
    for rid in probe:
        for e in rid:
            for p in g2p.get(str(e), set()):
                pool |= p2c.get(p, set())
    log('  이름이 필요한 ChEBI 화합물 %d개' % len(pool))
    cname = chebi_names(sorted(pool))

    # ---------------------------------------------------- 1~3단계
    paths, cpds, hits = genes_to_metabolites(cand_ids, g2p, p2c, cname, panel)
    log('')
    log('=' * 78)
    log('1~3단계 — 대사체 결과를 보기 전에 검정 대상 확정')
    log('=' * 78)
    log('  후보가 속한 Reactome 경로 : %d개' % len(paths))
    log('  그 경로들의 ChEBI 화합물  : %d개' % len(cpds))
    log('  패널에서 매칭된 대사물질  : %d개' % len(hits))
    for m, c in sorted(hits.items()):
        log('    %-38s CHEBI:%s' % (m[:36], c))

    # 콜라겐 연결이 살아났는지 명시적으로 확인
    log('')
    coll = [p for p in paths if 'collagen' in pname.get(p, '').lower()]
    log('  후보가 속한 콜라겐 관련 경로 %d개: %s'
        % (len(coll), ', '.join(sorted(pname[p] for p in coll))[:100]))
    hyp = [m for m in hits if 'ydroxyprol' in m]
    log('  하이드록시프롤린이 지목되었는가: %s' % ('예 -> ' + ', '.join(hyp) if hyp else '아니오'))

    if len(hits) < 3:
        log('  매칭이 너무 적어 검정할 수 없습니다.')
        return 1

    # ---------------------------------------------------- 4단계
    d = pd.read_csv(MET, sep='\t')
    d['is_target'] = d['metabolite'].isin(hits)
    tgt, oth = d[d['is_target']], d[~d['is_target']]
    obs = tgt['hedges_g'].abs().median()
    log('')
    log('=' * 78)
    log('4단계 — 지목된 대사물질이 DKD vs 다른 신장병에서 움직이는가')
    log('=' * 78)
    log('  지목 %d개 |g| 중앙 %.3f  ·  나머지 %d개 %.3f'
        % (len(tgt), obs, len(oth), oth['hedges_g'].abs().median()))
    log('  q<0.05 통과 %d / %d' % (int((tgt['q'] < 0.05).sum()), len(tgt)))

    # ---------------------------------------------------- 귀무분포
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
    log('')
    log('=' * 78)
    log('귀무분포 — 무작위 유전자 %d개' % len(cand_ids))
    log('=' * 78)
    log('  무작위가 지목한 대사물질 수 중앙 %d · 후보 %d' % (np.median(null_n), len(hits)))
    log('  무작위 |g| 중앙값 분포: 중앙 %.3f (95분위 %.3f)'
        % (np.median(null_med), np.percentile(null_med, 95)))
    log('  후보: %.3f   경험적 p = %.4f' % (obs, p_emp))

    pd.DataFrame([dict(ontology='Reactome', n_pathways=len(paths), n_compounds=len(cpds),
                       n_matched=len(hits), obs_median_abs_g=obs, empirical_p=p_emp,
                       collagen_pathways=len(coll),
                       hydroxyproline_captured=bool(hyp))]).to_csv(
        os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)
    pd.DataFrame({'metabolite': list(hits), 'chebi': [hits[m] for m in hits]}).to_csv(
        os.path.join(OUT, 'targeted_metabolites.tsv'), sep='\t', index=False)

    log('')
    log('=' * 78)
    log('KEGG와 비교')
    log('=' * 78)
    kp = 'results/kegg_pathway/summary.tsv'
    if os.path.exists(kp):
        k = pd.read_csv(kp, sep='\t').iloc[0]
        log('  %-10s %8s %10s %12s %10s' % ('온톨로지', '경로', '지목 대사물질', '|g| 중앙', '경험적 p'))
        log('  %-10s %8d %10d %14.3f %10.4f'
            % ('KEGG', k['n_pathways'], k['n_matched'], k['obs_median_abs_g'],
               k['empirical_p']))
        log('  %-10s %8d %10d %14.3f %10.4f'
            % ('Reactome', len(paths), len(hits), obs, p_emp))
    log('')
    if hyp:
        log('  Reactome이 콜라겐 연결을 규칙으로 재현했습니다. 수기 사례가 승격됩니다.')
    else:
        log('  Reactome도 콜라겐 -> 하이드록시프롤린을 잇지 못했습니다. 두 온톨로지 모두에서')
        log('  이 연결이 표현되지 않으므로, 수기 사례는 규칙으로 승격될 수 없습니다.')
        log('  이유는 같습니다 — 두 온톨로지 다 콜라겐 분해를 단백질에서 펩타이드까지만')
        log('  적고, 개별 아미노산까지 내려가지 않습니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
