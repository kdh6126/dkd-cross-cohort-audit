#!/usr/bin/env python
"""후보 30개를, 지금까지 한 번도 쓰지 않은 공개 데이터 세 개에서 검증한다.

각 데이터가 서로 다른 질문에 답합니다. 같은 질문을 세 번 묻는 것은 검증이 아닙니다.

    GSE20602   신경화증(NSC) 사구체 13 vs 종양신절제 5.
               DKD 데이터가 아니지만, 대조군 유형이 우리 GSE30528과 같습니다. 그래서
               "DKD vs 종양신절제"와 "신경화증 vs 종양신절제"를 같은 축에서 비교할 수
               있습니다. 후보가 DKD 특이적이라면 DKD 쪽에서만 크게 움직여야 합니다.
               논문이 한계로 적은 "특이성 게이트가 ERCB에만 의존한다"를 직접 겨냥합니다.

    GSE142153  말초혈액 단핵구 40샘플 (정상 / DN / 말기신부전).
               조직이 아니라 혈액입니다. 후보가 신장 섬유화 유전자라면 혈액에서는 보이지
               않아야 합니다. 반대로 보인다면 임상적으로 훨씬 쓸모 있지만, 조직 특이성
               주장은 약해집니다. 어느 쪽이든 정보입니다.

    GSE111154  초기 당뇨병성 신증 4 vs 대조 4.
               검정력이 매우 낮습니다. "초기에도 이미 변해 있는가"를 묻되, 음성이 나와도
               부재의 증거로 읽지 않습니다.

세 검증 모두 사후에 문턱을 정하지 않습니다. 후보군 대 무작위 유전자군의 비교로 봅니다.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

OUT = 'results/external_validation'
SEED = 0
RAW = 'data/raw/geo'
ANN = 'data/raw/annotation'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ------------------------------------------------------------------ 읽기

def read_series(acc):
    """series matrix에서 발현 행렬과 샘플 특성을 읽는다."""
    path = os.path.join(RAW, acc, '%s_series_matrix.txt.gz' % acc)
    meta, rows, header = {}, [], None
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        in_tab = False
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('!series_matrix_table_begin'):
                in_tab = True
                continue
            if line.startswith('!series_matrix_table_end'):
                break
            if in_tab:
                p = [x.strip('"') for x in line.split('\t')]
                if header is None:
                    header = p
                else:
                    rows.append(p)
            elif line.startswith('!Sample_'):
                p = [x.strip('"') for x in line.split('\t')]
                meta.setdefault(p[0], []).append(p[1:])
    X = pd.DataFrame(rows).set_index(0)
    X.columns = header[1:]
    X = X.apply(pd.to_numeric, errors='coerce')
    return X, meta


def probe2entrez(platform):
    """GEO .annot 또는 platform 테이블에서 프로브->Entrez 사전."""
    p = os.path.join(ANN, '%s.annot.gz' % platform)
    if os.path.exists(p):
        m, hdr, started = {}, None, False
        i_id = i_gid = -1
        with gzip.open(p, 'rt', encoding='utf-8', errors='replace') as fh:
            # 헤더는 !platform_table_begin 다음 줄이다. 그 앞에 '^Annotation' 처럼
            # ! 로도 # 로도 시작하지 않는 줄이 있어, 첫 줄을 헤더로 잡으면 깨진다.
            for line in fh:
                if not started:
                    if line.startswith('!platform_table_begin'):
                        started = True
                    continue
                if line.startswith('!platform_table_end'):
                    break
                f = line.rstrip('\n').split('\t')
                if hdr is None:
                    hdr = f
                    if 'ID' not in hdr or 'Gene ID' not in hdr:
                        return {}
                    i_id, i_gid = hdr.index('ID'), hdr.index('Gene ID')
                    continue
                if len(f) > i_gid and f[i_gid] and '///' not in f[i_gid]:
                    m[f[i_id]] = f[i_gid]
        return m
    if platform == 'GPL17586':
        # 이 플랫폼은 Entrez 컬럼이 없고 gene_assignment 문자열 안에 들어 있다.
        # 본 파이프라인에 이미 파서가 있으므로 재구현하지 않고 재사용한다.
        try:
            from build_gene_space import map_gpl17586
            return {k: str(v) for k, v in map_gpl17586().items()}
        except Exception as e:
            log('  (GPL17586 매핑 재사용 실패: %s)' % str(e)[:60])
            return {}
    return {}


def collapse(X, platform):
    """max-mean 규칙으로 프로브를 유전자로 접는다. 본 파이프라인과 동일한 규칙."""
    m = probe2entrez(platform)
    if not m:
        return None
    g = pd.Series({p: m.get(p) for p in X.index})
    keep = g.notna()
    X, g = X[keep.values], g[keep]
    means = X.mean(axis=1)
    order = means.sort_values(ascending=False).index
    X = X.loc[order]
    g = g.loc[order]
    first = ~g.duplicated()
    Y = X[first.values]
    Y.index = g[first.values].astype(str)
    return Y


def sample_annot(meta, columns, key=None):
    """샘플별 주석을 한 문자열로 합친다.

    GEO는 특성을 여러 !Sample_characteristics_ch1 줄에 나눠 적는다. 첫 줄만 보면
    엉뚱한 필드를 읽는다 — GSE142153에서 첫 줄은 조직, 진단은 두 번째 줄이었다.
    """
    parts = []
    for k, blocks in meta.items():
        if key and key not in k:
            continue
        for b in blocks:
            if len(b) == len(columns):
                parts.append(pd.Series(b, index=columns))
    if not parts:
        return pd.Series([''] * len(columns), index=columns)
    return pd.concat(parts, axis=1).agg(' | '.join, axis=1)


def hedges_g(a, b):
    n1, n0 = len(a), len(b)
    if n1 < 2 or n0 < 2:
        return np.nan
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / (n1 + n0 - 2))
    if not np.isfinite(sp) or sp == 0:
        return np.nan
    return (a.mean() - b.mean()) / sp * (1 - 3 / (4 * (n1 + n0) - 9))


def enrich_vs_random(values, cand_idx, rng, n_draw=2000):
    """후보군의 |g| 중앙값이 같은 크기 무작위군보다 큰가. 사후 문턱을 피한다."""
    v = np.abs(values)
    obs = np.nanmedian(v[cand_idx])
    k = len(cand_idx)
    null = np.array([np.nanmedian(v[rng.choice(len(v), k, replace=False)])
                     for _ in range(n_draw)])
    p = (np.sum(null >= obs) + 1) / (n_draw + 1)
    return obs, float(np.nanmedian(null)), p


# ------------------------------------------------------------------ 검증

def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)

    cand = pd.read_csv('results/candidates_v2/master_candidate_table.tsv', sep='\t')
    gs = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = dict(zip(gs['symbol'], gs['entrez_id']))
    cand_ids = [sym2id[s] for s in cand[cand.columns[0]] if s in sym2id]
    log('후보 %d개 중 유전자 공간에 매핑된 것 %d개' % (len(cand), len(cand_ids)))
    log('')

    results = []

    # ================================================== 1. GSE20602
    log('=' * 76)
    log('검증 1 — 신경화증도 같은 방향으로 움직이는가 (ERCB 밖 독립 특이성)')
    log('=' * 76)
    X, meta = read_series('GSE20602')
    Y = collapse(X, 'GPL96')
    src = sample_annot(meta, X.columns, key='source_name')
    nsc = src.str.contains('NSC', case=False)
    ctl = src.str.contains('Nephrectomy', case=False)
    log('  신경화증 %d, 종양신절제 %d, 유전자 %d' % (nsc.sum(), ctl.sum(), len(Y)))
    g_nsc = pd.Series({gid: hedges_g(Y.loc[gid, nsc.values].values,
                                     Y.loc[gid, ctl.values].values)
                       for gid in Y.index if gid in Y.index})
    g_nsc = g_nsc.dropna()

    # 같은 대조군 유형(종양신절제)을 쓴 우리 DKD 코호트
    fc = pd.read_csv('results/final_candidates.tsv', sep='\t')
    fc = fc.rename(columns={fc.columns[0]: 'entrez_id'})
    fc['entrez_id'] = fc['entrez_id'].astype(str)
    g_dkd = fc.set_index('entrez_id')['g_DKD_vs_control']

    common = [g for g in g_nsc.index if g in g_dkd.index]
    cid = [i for i, g in enumerate(common) if g in cand_ids]
    log('  두 데이터에 공통인 유전자 %d개, 그중 후보 %d개' % (len(common), len(cid)))
    dv = g_dkd.loc[common].values
    nv = g_nsc.loc[common].values
    r = stats.pearsonr(dv, nv)
    log('')
    log('  전체 유전자에서 DKD 효과와 신경화증 효과의 상관: r = %.3f (p = %.1e)'
        % (r[0], r[1]))
    log('  -> 두 질병이 공유하는 "신장이 나빠서" 성분이 %.0f%%입니다.' % (100 * r[0] ** 2))

    ratio = np.abs(dv) - np.abs(nv)
    obs, nul, p = enrich_vs_random(ratio, cid, rng)
    log('')
    log('  후보군의 |g(DKD)| - |g(NSC)| 중앙값 : %+.3f' % np.nanmedian(ratio[cid]))
    log('  무작위군                             : %+.3f' % np.nanmedian(ratio))
    tt = stats.mannwhitneyu(ratio[cid], np.delete(ratio, cid), alternative='greater')
    log('  Mann-Whitney p = %.2e  -> %s' % (tt.pvalue,
        '후보가 DKD 쪽으로 더 치우침' if tt.pvalue < 0.05 else '무작위와 구분 안 됨'))
    pd.DataFrame({'entrez_id': common, 'g_DKD_vs_control': dv,
                  'g_NSC_vs_control': nv,
                  'is_candidate': [g in cand_ids for g in common]}).to_csv(
        os.path.join(OUT, 'gse20602_nsc.tsv'), sep='\t', index=False)
    results.append(dict(dataset='GSE20602', question='ERCB 밖 독립 특이성',
                        stat='MWU p', value=tt.pvalue,
                        verdict='DKD 특이' if tt.pvalue < 0.05 else '구분 안 됨'))

    # ================================================== 2. GSE142153
    log('')
    log('=' * 76)
    log('검증 2 — 혈액에서도 보이는가 (조직 특이성 / 임상 접근성)')
    log('=' * 76)
    X2, meta2 = read_series('GSE142153')
    Y2 = collapse(X2, 'GPL6480')
    ch = sample_annot(meta2, X2.columns, key='characteristics')
    dn = ch.str.contains('diabetic nephropathy', case=False)
    hc = ch.str.contains('healthy', case=False)
    log('  DN %d, 정상 %d, 유전자 %d' % (dn.sum(), hc.sum(), len(Y2)))
    if dn.sum() >= 2 and hc.sum() >= 2:
        g_bl = pd.Series({gid: hedges_g(Y2.loc[gid, dn.values].values,
                                        Y2.loc[gid, hc.values].values)
                          for gid in Y2.index}).dropna()
        cid2 = [i for i, g in enumerate(g_bl.index) if g in cand_ids]
        log('  혈액에서 측정된 후보 %d개' % len(cid2))
        obs, nul, p = enrich_vs_random(g_bl.values, cid2, rng)
        log('  후보군 |g| 중앙값 %.3f  vs  무작위 %.3f   순열 p = %.3f' % (obs, nul, p))
        log('  -> %s' % ('혈액에서도 변합니다. 조직 특이성이 완전하지 않다는 뜻이자, '
                         '임상적으로는 유리한 신호입니다.' if p < 0.05 else
                         '혈액에서는 무작위 유전자와 구분되지 않습니다. 조직 국한 신호로 '
                         '보입니다.'))
        pd.DataFrame({'entrez_id': g_bl.index, 'g_DN_vs_healthy_blood': g_bl.values,
                      'is_candidate': [g in cand_ids for g in g_bl.index]}).to_csv(
            os.path.join(OUT, 'gse142153_blood.tsv'), sep='\t', index=False)
        results.append(dict(dataset='GSE142153', question='혈액에서 보이는가',
                            stat='permutation p', value=p,
                            verdict='혈액에서도 변함' if p < 0.05 else '조직 국한'))

    # ================================================== 3. GSE111154
    log('')
    log('=' * 76)
    log('검증 3 — 초기 당뇨병성 신증에서도 이미 변해 있는가 (검정력 낮음)')
    log('=' * 76)
    X3, meta3 = read_series('GSE111154')
    Y3 = collapse(X3, 'GPL17586')
    if Y3 is None or len(Y3) == 0:
        log('  GPL17586 프로브->유전자 매핑 실패. 이 검증은 건너뜁니다.')
    else:
        src3 = sample_annot(meta3, X3.columns, key='source_name')
        edn = src3.str.contains('early', case=False)
        ct3 = src3.str.contains('Non-diabetic', case=False)
        log('  초기 DN %d, 대조 %d, 유전자 %d' % (edn.sum(), ct3.sum(), len(Y3)))
        g_e = pd.Series({gid: hedges_g(Y3.loc[gid, edn.values].values,
                                       Y3.loc[gid, ct3.values].values)
                         for gid in Y3.index}).dropna()
        cid3 = [i for i, g in enumerate(g_e.index) if g in cand_ids]
        log('  측정된 후보 %d개' % len(cid3))
        if len(cid3) >= 5:
            obs, nul, p = enrich_vs_random(g_e.values, cid3, rng)
            log('  후보군 |g| 중앙값 %.3f  vs  무작위 %.3f   순열 p = %.3f' % (obs, nul, p))
            log('  주의 — 각 군 4명입니다. 음성이 나와도 "초기에는 변하지 않는다"의 증거로')
            log('  읽으면 안 됩니다. 검정력이 부족해 아무것도 못 볼 수 있습니다.')
            results.append(dict(dataset='GSE111154', question='초기에도 변하는가',
                                stat='permutation p', value=p,
                                verdict='초기에도 변함' if p < 0.05 else '결론 보류(검정력)'))

    # ================================================== 요약
    log('')
    log('=' * 76)
    log('요약')
    log('=' * 76)
    t = pd.DataFrame(results)
    t.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)
    for _, r_ in t.iterrows():
        log('  %-12s %-24s %s = %.2e  ->  %s'
            % (r_['dataset'], r_['question'], r_['stat'], r_['value'], r_['verdict']))
    log('')
    log('  wrote %s/' % OUT)


if __name__ == '__main__':
    main()
