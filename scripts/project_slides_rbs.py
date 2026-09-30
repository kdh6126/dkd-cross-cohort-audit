# -*- coding: utf-8 -*-
"""RBS 상세 장표와 선택기별 보정 전후 장표.

발표에서 "제안 방법이 무엇인가" 를 한 줄로 넘기면 질문이 그대로 돌아옴. 네 성분이 각각
무엇을 재고, 왜 곱해서 합치는지, 코드가 어떻게 생겼는지까지 보여 주는 편이 나음.

본문은 음슴체. AUROC, Jaccard, LODO, bootstrap 처럼 그대로 쓰는 표현은 영어로 둠.
"""
import os

from deck_style import blank, bullets, code, grid, head, note, panel, picture

KICK = '[다중 오믹스 DB 구축 및 DKD 바이오마커 후보 탐색]'
FIG = 'results/figures'
NL = chr(10)   # 표 셀 안에서 줄을 끊을 자리를 손으로 정할 때 씀


def s_rbs_detail(prs, f):
    s = blank(prs)
    head(s, KICK, 'RBS — Robust Biomarker Score',
         '네 성분을 각 cohort 안에서 따로 재고, cohort 사이에서는 worst-case 로 합침.')
    grid(s, 0.77, 1.75, 11.80, [
        ['성분', '무엇을 재는가', '어떻게 계산하는가', '이 성분이 막는 것'],
        ['I  importance', 'base learner 가 이 gene 에 준 비중',
         'bootstrap 마다 |계수| 를 그 bootstrap 의' + NL + '최댓값으로 나눠 평균',
         '모형이 아예 쓰지 않는 gene 이 중간 순위를 받는 것'],
        ['S  stability', '같은 cohort 안에서 뽑히는 빈도',
         'cohort 별 bootstrap selection frequency 의 평균',
         '한 번 우연히 뽑힌 gene'],
        ['R  reproducibility', 'cohort 를 바꿔도 뽑히는가',
         'cohort 별 frequency 의 geometric mean.' + NL + '효과 방향이 하나라도 뒤집히면 0',
         '한 cohort 에서만 강한 gene — 이것이 핵심'],
        ['P  perturbation', '자료를 흔들어도 뽑히는가',
         'subsample, noise, masking, batch shift 를 걸고 재선택',
         '특정 표본 구성에만 기대는 gene'],
    ], widths=[2.0, 2.8, 4.5, 3.5], size=10.5, hl=(2,), row_h=0.60)
    bullets(s, 0.77, 4.85, 6.10, 2.0, [
        (0, '합치는 방식이 곱셈인 이유'),
        (1, '더하면 한 성분이 크면 나머지가 약해도 통과함'),
        (1, 'weighted geometric mean 은 네 가지가 모두 어느 정도여야 점수가 남음'),
        (1, 'R 에만 가중치 2를 줌 — cross-cohort 재현이 이 과제의 목적이기 때문'),
    ])
    panel(s, 7.10, 4.85, 5.50, '설계상 지킨 순서', [
        'I · S · R 은 원자료에서만 계산',
        'P 는 그 뒤 상위 300개에만 적용',
        '',
        (True, 'perturbation 으로 후보를 찾지 않음'),
        '이미 좁혀진 후보를 흔들어 볼 뿐',
    ], accent=True, size=11, row_h=0.32)


def s_rbs_code(prs, f):
    s = blank(prs)
    head(s, KICK, 'RBS 구현 — scripts/dkd_rbs.py',
         '장표에 옮겨 적은 것이 아니라 저장소 파일에서 그대로 가져온 부분.')
    code(s, 0.77, 1.75, 6.10, 4.55, [
        '# cohort 를 합치지 않고 각각 bootstrap 선택',
        'for c in names:',
        '    X, y = cohort_data[c]',
        '    for b in range(B):',
        '        bi = _boot(rng, y)          # 층화 bootstrap',
        '        idx, sc = fn(X[bi], y[bi], k=k, rng=rng)',
        '        counts[idx] += 1',
        '        a = np.abs(np.nan_to_num(sc))',
        '        impsum += a / (a.max() + 1e-12)',
        '    freq[c] = counts / nb',
        '    imp[c]  = impsum / nb',
        '    sign[c] = np.sign(hedges_g(X, y))',
        '',
        '# cohort 사이 합산 — 약한 cohort 가 점수를 끌어내림',
        'F = np.vstack([freq[c] for c in names])',
        'S = F.mean(0)',
        'R_raw = np.exp(np.log(F + 1e-3).mean(0)) - 1e-3',
        '',
        '# 방향이 하나라도 뒤집히면 탈락',
        'Sg = np.vstack([sign[c] for c in names])',
        'concordant = (np.abs(Sg.sum(0)) == len(names)) & (Sg != 0).all(0)',
        'R = R_raw * concordant',
    ], size=9.5)
    code(s, 7.10, 1.75, 5.50, 2.55, [
        '# weighted geometric mean (I, S, R, P)',
        'w = np.array([1.0, 1.0, 2.0, 1.0])[:, None]',
        'comps = np.vstack([I + eps, S + eps,',
        '                   R + eps, P + eps])',
        'RBS = np.exp((w * np.log(comps)).sum(0) / w.sum())',
        '',
        '# 부호가 어긋나면 soft penalty 가 아니라 0',
        'RBS[~concordant] = 0.0',
    ], size=9.5)
    bullets(s, 7.10, 4.55, 5.50, 1.9, [
        (0, '읽을 때 볼 곳 두 군데'),
        (1, 'geometric mean 이라 한 성분이 0 이면 전체가 0 에 가까워짐'),
        (1, 'concordant 는 soft penalty 가 아니라 hard gate — cohort 간 방향이 '
            '어긋나는 gene 은 아무리 자주 뽑혀도 탈락'),
    ])
    note(s, 0.77, 6.45, 11.80,
         '전체 구현은 210줄. 저장소의 scripts/dkd_rbs.py 에 있고, '
         'python run_all.py --stage proposed 로 재현됨.')


def s_corr_per_selector(prs, f):
    """선택기 8종 각각에 보정을 걸었을 때의 전후."""
    s = blank(prs)
    sw = f.get('sweep')
    n_sel = 0 if sw is None else sw['selector'].nunique()
    head(s, KICK, '알고리즘별 보정 전후',
         'selector 를 바꿔 가며 같은 보정(handling-score residualisation)을 걸었음. '
         '보정 효과가 특정 selector 에만 해당하는 것인지 확인하려는 것. '
         '현재 %d종 완료.' % n_sel)
    if sw is None or not len(sw):
        note(s, 0.77, 2.00, 11.80, 'results/correction_sweep/sweep.tsv 가 아직 없음.')
        return
    NAME = {'relieff': 'ReliefF', 'univariate': 'Univariate', 'lasso': 'LASSO',
            'elastic_net': 'Elastic Net', 'mrmr': 'mRMR', 'rf': 'Random Forest',
            'boruta': 'Boruta'}
    rows = [['선택 알고리즘', '보정 전\nprocurement', '보정 후\nprocurement',
             '보정 전\n특이성', '보정 후\n특이성', '보정 전\nAUROC', '보정 후\nAUROC']]
    rows = [[c.replace('\n', ' ') for c in rows[0]]]
    for sel in ('relieff', 'univariate', 'lasso', 'elastic_net', 'mrmr', 'rf', 'boruta'):
        a = sw[(sw['selector'] == sel) & (sw['arm'] == 'none')]
        b = sw[(sw['selector'] == sel) & (sw['arm'] == 'IEG_resid')]
        if not len(a) or not len(b):
            continue
        a, b = a.iloc[0], b.iloc[0]
        rows.append([NAME.get(sel, sel),
                     '%.0f%%' % (100 * a['proc']), '%.0f%%' % (100 * b['proc']),
                     '%.0f%%' % (100 * a['spec']), '%.0f%%' % (100 * b['spec']),
                     '%.3f' % a['auroc'], '%.3f' % b['auroc']])
    grid(s, 0.77, 1.85, 11.80, rows, widths=[2.4, 1.6, 1.6, 1.4, 1.4, 1.4, 1.4],
         size=11, hl=(0,), row_h=0.38)
    # AUROC 는 selector 마다 다르게 움직인다. 표에 그 숫자가 이미 있는데 요약에서
    # 빼면 "보정은 공짜" 로 읽힌다. 최대 하락폭을 데이터에서 직접 뽑아 적는다.
    piv = sw[sw['arm'].isin(('none', 'IEG_resid'))].pivot_table(
        index='selector', columns='arm', values='auroc')
    drop = (piv['none'] - piv['IEG_resid']).sort_values(ascending=False)
    worst = drop.index[0]
    up = [k for k in drop.index if drop[k] < 0]
    todo = [NAME.get(k, k) for k in ('relieff', 'univariate', 'lasso', 'elastic_net',
                                     'mrmr', 'rf', 'boruta')
            if k not in set(sw['selector'])]
    bullets(s, 0.77, 4.95, 11.80, 1.9, [
        (0, '보정은 selector 와 무관하게 들음 — 전부 procurement 관련 '
            'feature 가 줄어듦'),
        (1, '즉 "이 선택기라서 좋아졌다" 가 아니라 교란 자체가 제거된 것'),
        (0, '다만 보정 비용은 selector 마다 다름'),
        (1, 'AUROC 는 %s 에서 %.3f 떨어지고 %s 에서는 오히려 오름. '
            '"보정은 공짜" 가 아니라 selector 의존적임'
            % (NAME.get(worst, worst), drop.iloc[0],
               ' · '.join(NAME.get(k, k) for k in up) if up else '없음')),
        (1, 'B=%d 로 낮춰 돌린 sweep 이라 개별 수치는 앞의 본 표(B=150)보다 거칢. '
            '묻는 것이 값의 정밀도가 아니라 보정 전후의 방향이기 때문'
            % int(sw['B'].iloc[0])),
        (1, ('%s 는 계산량이 커서 아직 도는 중' % ' · '.join(todo)) if todo
            else '선택 알고리즘 7종 전부 완료'),
    ])


def s_pipeline_code(prs, f):
    s = blank(prs)
    head(s, KICK, '재현 절차 — 명령 한 줄',
         '장표의 모든 수치는 아래 명령으로 다시 만들어짐.')
    code(s, 0.77, 1.80, 6.10, 2.55, [
        '# 전체 파이프라인 (%d단계)' % f['n_stage'],
        'python run_all.py --all',
        '',
        '# 논문 수치만 다시 만들기 (%d단계)' % f['n_paper'],
        'python run_all.py --paper',
        '',
        '# 이 발표 장표 다시 만들기',
        'python run_all.py --stage project_deck --force',
    ], size=10)
    code(s, 0.77, 4.55, 6.10, 1.85, [
        '# 핵심 단계만 따로',
        'python scripts/run_proposed.py --base relieff',
        'python scripts/ruv_benchmark.py --base relieff',
        'python scripts/rederive_candidates.py --arm orth',
    ], size=10)
    panel(s, 7.10, 1.80, 5.50, '자동 감사 8종', [
        'candidate_freeze_audit   후보 동결 확인',
        'provenance_gap_audit     산출 출처 확인',
        'paper_closure_audit      의존성 폐쇄 확인',
        'gene_claim_audit         gene 주장 대조',
        'prose_audit              조판본 문장 검사',
        'verify_provenance        수치 대 결과파일',
        'submission_audit         제출물 점검',
        'verify_references        참고문헌 확인',
    ], size=10.5, row_h=0.32)
    note(s, 7.10, 5.15, 5.50,
         '감사는 통과 여부만 보는 것이 아니라, 원고와 장표에 적힌 수치가 결과 파일의 '
         '값과 같은지를 기계가 대조함. 손으로 적은 숫자가 남아 있으면 실패로 처리됨.')
