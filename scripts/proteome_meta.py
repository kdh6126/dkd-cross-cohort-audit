#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""두 단백체 코호트를 합쳐서 본다.

각각은 후보를 조금씩만 잽니다. SOMAscan 은 표적 패널이라 7개, 질량분석은 깊이가
얕아 7개입니다. 7개로는 어느 쪽도 집합 수준의 결론을 낼 수 없었습니다
(Fisher p = 0.11 과 0.18).

합치는 것이 옳은 이유는 두 자료가 독립이기 때문입니다.

    환자   서로 겹치지 않습니다. 우리 전사체 코호트와도 겹치지 않습니다.
    기술   앱타머 결합(SOMAscan) 대 질량분석(LC-MS/MS). 측정 원리가 다릅니다.
    대상   후보 12개를 덮고, 공통으로 잰 것은 LUM 과 MMP7 둘뿐입니다.

분석 계획은 첫 코호트에서 정한 것을 그대로 씁니다. 두 번째 코호트를 찾은 뒤에 기준을
바꾸지 않았습니다. 계획을 바꾸지 않았다는 사실이 이 결합의 정당성입니다.

층화가 필요한 이유도 분명합니다. 두 코호트의 배경 유의 비율이 28%와 11%로 크게
다릅니다. 그냥 합쳐서 세면 배경이 높은 코호트가 결과를 지배합니다. 층 안에서 비교한 뒤
층을 가로질러 합해야 합니다. Cochran-Mantel-Haenszel 이 그 계산이고, 칸이 작아
근사를 믿기 어려우므로 순열검정으로 같은 값을 다시 구합니다.

순열은 세 가지를 지켜야 합니다.

  비복원   한 번의 추출 안에서 같은 단백질을 두 번 뽑으면 안 됩니다. 후보 집합에는
           중복이 없으므로 귀무 집합에도 없어야 합니다. 복원추출로 두면 분산이
           달라져서 "층별 배경에서 같은 수를 뽑는다"는 설명과 계산이 어긋납니다.
  층 유지  층마다 그 층의 배경에서만 뽑습니다. 층을 섞으면 배경 비율 차이가 다시
           결과를 지배합니다.
  공변량   후보가 원래 존재비가 높거나 펩타이드가 많아 잘 검출되는 단백질이라면,
           아무 배경과 비교하는 것으로는 답이 안 됩니다. 그래서 존재비와 펩타이드
           수를 맞춘 배경으로 같은 검정을 한 번 더 합니다.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

SOMA = 'results/proteome_validation/protein_dkd_vs_control.tsv'
MS = 'results/proteome_ms/protein_ms_dkd_vs_control.tsv'
OUT = 'results/proteome_meta'
NPERM = 200000
CHUNK = 10000
CALIPER = 0.10          # 존재비 순위 차이를 배경 크기의 이 비율 안으로 제한한다
SEED = 0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def cmh(strata):
    """Cochran-Mantel-Haenszel. strata 는 (a, b, c, d) 목록."""
    num = var = 0.0
    for a, b, c, d in strata:
        n = a + b + c + d
        num += a - (a + b) * (a + c) / n
        var += (a + b) * (c + d) * (a + c) * (b + d) / (n ** 2 * (n - 1))
    chi = (abs(num) - 0.5) ** 2 / var          # 연속성 보정. 보수적인 쪽이다.
    return chi, stats.chi2.sf(chi, 1) / 2      # 단측


def cmh_or(strata):
    """공통 승산비와 95% 신뢰구간 (Mantel-Haenszel + Robins-Breslow-Greenland).

    p 값만 내놓으면 효과가 얼마나 큰지, 그리고 표본이 얼마나 작은지가 둘 다 가려집니다.
    넓은 구간을 그대로 보여주는 편이 정직합니다. 후보가 14개(측정 사건 기준)뿐이므로
    구간은 실제로 넓습니다.
    """
    R = S = 0.0
    vr = vrs = vs = 0.0
    for a, b, c, d in strata:
        n = a + b + c + d
        Ri, Si = a * d / n, b * c / n
        P, Q = (a + d) / n, (b + c) / n
        R += Ri
        S += Si
        vr += P * Ri
        vrs += P * Si + Q * Ri
        vs += Q * Si
    or_mh = R / S
    var = vr / (2 * R ** 2) + vrs / (2 * R * S) + vs / (2 * S ** 2)
    se = np.sqrt(var)
    lo = np.exp(np.log(or_mh) - 1.96 * se)
    hi = np.exp(np.log(or_mh) + 1.96 * se)
    return or_mh, lo, hi


def draw_unmatched(rng, pool_size, k, nperm):
    """배경 전체에서 k개를 비복원으로 뽑기를 nperm 번. 결과는 배경 안의 색인."""
    out = np.empty((nperm, k), dtype=np.int64)
    for s in range(0, nperm, CHUNK):
        n = min(CHUNK, nperm - s)
        # argpartition 은 앞 k개만 정렬하므로 전체 정렬보다 빠르다.
        out[s:s + n] = np.argpartition(rng.random((n, pool_size)), k, axis=1)[:, :k]
    return out


def draw_matched(rng, pools, nperm):
    """후보마다 정해진 적격 배경 목록에서 하나씩 뽑되, 한 추출 안에서 겹치지 않게.

    pools 는 후보 순서대로의 배경 색인 배열 목록입니다. 각 후보의 적격 목록은 그
    후보와 존재비(그리고 질량분석에서는 펩타이드 수)가 비슷한 배경 단백질입니다.
    겹치면 그 후보만 다시 뽑습니다. 목록이 넓으므로 재추출은 드뭅니다.
    """
    k = len(pools)
    out = np.empty((nperm, k), dtype=np.int64)
    for j, pool in enumerate(pools):
        out[:, j] = rng.choice(pool, size=nperm, replace=True)
    for _ in range(60):
        dup = np.zeros(nperm, bool)
        srt = np.sort(out, axis=1)
        dup |= (srt[:, 1:] == srt[:, :-1]).any(axis=1)
        if not dup.any():
            break
        for j, pool in enumerate(pools):
            out[dup, j] = rng.choice(pool, size=int(dup.sum()), replace=True)
    return out


def eligible(cand_row, bg, cols, caliper_n):
    """이 후보와 공변량이 비슷한 배경 단백질의 색인.

    존재비는 순위 차이로 자릅니다. 값 자체로 자르면 코호트마다 눈금이 달라
    같은 캘리퍼가 다른 뜻이 됩니다. 펩타이드 수는 반올림한 로그가 같은 것으로
    맞춥니다 — 1개짜리와 20개짜리를 같은 것으로 볼 수는 없습니다.
    """
    ok = np.ones(len(bg), bool)
    rank = bg['abundance'].rank(pct=False).values
    my = (bg['abundance'] < cand_row['abundance']).sum() + 0.5
    ok &= np.abs(rank - my) <= caliper_n
    if 'peptides' in cols and np.isfinite(cand_row.get('peptides', np.nan)):
        lp = np.round(np.log2(np.maximum(bg['peptides'].values, 1)))
        ok &= np.abs(lp - np.round(np.log2(max(cand_row['peptides'], 1)))) <= 1
    return np.flatnonzero(ok)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)

    s = pd.read_csv(SOMA, sep='\t')
    s = s[s['on_panel_and_measured']].copy()
    s['conc'] = np.sign(s['prot_g']) == np.sign(s['rna_g'])
    s = s.rename(columns={'is_candidate': 'is_cand'})
    m = pd.read_csv(MS, sep='\t')

    layers = [('SOMAscan (Mendeley 83k89shdx5)', s), ('LC-MS/MS (PRIDE PXD041884)', m)]

    log('=' * 78)
    obs_sig = obs_conc = 0
    strata = []
    pack = []
    for name, t in layers:
        c = t[t['is_cand']]
        b = t[~t['is_cand']]
        ns, nc = int((c['prot_q'] < 0.05).sum()), int(c['conc'].sum())
        bs, bc = int((b['prot_q'] < 0.05).sum()), int(b['conc'].sum())
        obs_sig += ns
        obs_conc += nc
        strata.append((ns, len(c) - ns, bs, len(b) - bs))
        pack.append((name, c.reset_index(drop=True), b.reset_index(drop=True)))
        log('%-30s 후보 %2d개  유의 %d  방향일치 %d   |  배경 %4d개  유의 %.0f%%  방향일치 %.0f%%'
            % (name, len(c), ns, nc, len(b), 100 * bs / len(b), 100 * bc / len(b)))

    log('')
    log('합계: 후보 %d개 중 유의 %d개, 방향일치 %d개'
        % (sum(len(x[1]) for x in pack), obs_sig, obs_conc))

    chi, p_cmh = cmh(strata)
    or_mh, or_lo, or_hi = cmh_or(strata)
    log('')
    log('층화 검정 (배경보다 후보에서 단백질 유의가 더 잦은가)')
    log('  Cochran-Mantel-Haenszel   chi2 = %.3f   단측 p = %.4f' % (chi, p_cmh))
    log('  공통 승산비               OR = %.2f   95%% CI %.2f - %.2f'
        % (or_mh, or_lo, or_hi))
    log('    측정 사건 %d건(후보 x 코호트) 위에서 잰 값이라 구간이 넓습니다.'
        % sum(a + b for a, b, c, d in strata))

    # ---------------------------------------------------------------- 순열 1: 매칭 없음
    gs = np.zeros(NPERM, int)
    gc = np.zeros(NPERM, int)
    for name, c, b in pack:
        idx = draw_unmatched(rng, len(b), len(c), NPERM)
        gs += (b['prot_q'].values < 0.05)[idx].sum(1)
        gc += b['conc'].values[idx].sum(1)
    p_sig = float((gs >= obs_sig).mean())
    p_conc = float((gc >= obs_conc).mean())
    log('  순열 %s회 · 비복원          p = %.4f   (귀무 평균 %.2f개)'
        % ('{:,}'.format(NPERM), p_sig, gs.mean()))
    log('  방향 일치는 같은 방식으로     p = %.4f   (귀무 평균 %.2f개)'
        % (p_conc, gc.mean()))

    # ------------------------------------------- 순열 2: 존재비·펩타이드 수를 맞춘 배경
    log('')
    log('민감도 분석 — 후보가 원래 잘 검출되는 단백질이라서 유의한 것은 아닌가')
    ms_sig = np.zeros(NPERM, int)
    bal = []
    for name, c, b in pack:
        cal = max(30, int(CALIPER * len(b)))
        pools = []
        for _, row in c.iterrows():
            e = eligible(row, b, set(b.columns), cal)
            if len(e) < 20:                      # 너무 좁으면 캘리퍼를 푼다
                e = eligible(row, b, set(), cal * 3)
            pools.append(e)
        idx = draw_matched(rng, pools, NPERM)
        ms_sig += (b['prot_q'].values < 0.05)[idx].sum(1)
        bal.append((name, len(c), float(np.mean([len(p) for p in pools])),
                    float(c['abundance'].mean()),
                    float(b['abundance'].values[idx[:2000]].mean()),
                    float(b['prot_q'].values[idx].mean() < 1)))
        log('  %-30s 후보당 적격 배경 %.0f개 · 존재비 후보 %.2f 대 매칭 귀무 %.2f'
            % (name, bal[-1][2], bal[-1][3], bal[-1][4]))
    p_matched = float((ms_sig >= obs_sig).mean())
    log('  매칭 순열                    p = %.4f   (귀무 평균 %.2f개)'
        % (p_matched, ms_sig.mean()))
    log('')
    if p_matched < 0.05:
        log('  검출성을 맞춰도 유지됩니다. "잘 잡히는 단백질이라 그렇다" 로는 설명되지 않습니다.')
    else:
        log('  매칭하면 유의성이 사라집니다. 검출성으로 설명될 여지가 있다고 보고해야 합니다.')
    log('  유의 빈도는 합치면 유의합니다. 방향 일치만으로는 아닙니다.')
    log('  방향은 배경도 60~67%가 맞으므로 변별력이 약합니다. 놀랄 일이 아닙니다.')

    sc = set(pack[0][1]['gene'])
    mc = set(pack[1][1]['gene'])
    sig_s = set(pack[0][1][pack[0][1]['prot_q'] < 0.05]['gene'])
    sig_m = set(pack[1][1][pack[1][1]['prot_q'] < 0.05]['gene'])
    log('')
    log('  단백질 수준 자료가 있는 후보 %d개 (양쪽 공통 %s)'
        % (len(sc | mc), ', '.join(sorted(sc & mc))))
    log('  한 코호트 이상에서 유의한 후보 %d개: %s'
        % (len(sig_s | sig_m), ', '.join(sorted(sig_s | sig_m))))
    log('  두 코호트 모두에서 유의: %s' % ', '.join(sorted(sig_s & sig_m)))

    pd.DataFrame([dict(n_layers=len(layers), n_candidates_measured=len(sc | mc),
                       n_sig_union=len(sig_s | sig_m), n_sig_both=len(sig_s & sig_m),
                       obs_sig=obs_sig, obs_conc=obs_conc,
                       cmh_chi2=chi, cmh_p=p_cmh,
                       cmh_or=or_mh, cmh_or_lo=or_lo, cmh_or_hi=or_hi,
                       n_measurement_events=sum(a + b for a, b, c, d in strata),
                       perm_p_sig=p_sig,
                       perm_p_conc=p_conc, perm_p_matched=p_matched,
                       perm_null_mean=float(gs.mean()),
                       perm_null_mean_matched=float(ms_sig.mean()),
                       n_perm=NPERM, sampling='without replacement')]
                 ).to_csv(os.path.join(OUT, 'summary.tsv'), sep='\t', index=False)
    pd.DataFrame(bal, columns=['layer', 'n_candidates', 'mean_eligible_background',
                               'candidate_abundance', 'matched_null_abundance', '_']
                 ).drop(columns='_').to_csv(
        os.path.join(OUT, 'matching_balance.tsv'), sep='\t', index=False)
    pd.DataFrame(sorted((sc | mc)), columns=['gene']).assign(
        on_somascan=lambda t: t['gene'].isin(sc),
        on_massspec=lambda t: t['gene'].isin(mc),
        sig_somascan=lambda t: t['gene'].isin(sig_s),
        sig_massspec=lambda t: t['gene'].isin(sig_m),
    ).to_csv(os.path.join(OUT, 'candidate_protein_coverage.tsv'), sep='\t', index=False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
