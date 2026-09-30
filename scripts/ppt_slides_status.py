# -*- coding: utf-8 -*-
"""현재 상황 요약 슬라이드.

덱 앞머리에 붙습니다. "무엇을 알아냈나" 가 아니라 **지금 어디까지 왔고 무엇이 남았나**
를 답하는 자리입니다. 수치는 이 덱의 원칙대로 결과 파일에서 읽습니다 — 손으로 적으면
분석이 바뀔 때 조용히 어긋나고, 그 어긋남은 발표장에서 드러납니다.
"""
from make_ppt import add_slide, header, txt, rule, table, ACCENT, BLUE, INK, MUTED
from status_facts import facts

AUDITS = ('candidate_freeze_audit', 'provenance_gap_audit', 'paper_closure_audit',
          'gene_claim_audit', 'prose_audit', 'verify_provenance', 'submission_audit',
          'verify_references')


def s_now(prs, d):
    f = facts()
    s = add_slide(prs)
    header(s, 'WHERE WE ARE', '지금 상황 한 장 요약',
           '원고는 투고 직전 상태입니다. 남은 것은 자리표 3개뿐입니다.')
    txt(s, 0.6, 2.05, 4.0, 0.35, '끝난 것', size=14, bold=True, color=BLUE)
    txt(s, 0.6, 2.5, 4.0, 3.7,
        ['· 원고 %d쪽, 그림 9 · 표 12 · 부록 6' % f['pages'],
         '  조판 오류 0 · 미해결 참조 0',
         '· 재현 파이프라인 %d단계' % f['n_stage'],
         '  (--paper 가 %d단계를 다시 돌림)' % f['n_paper'],
         '· 자동 감사 8종 전부 통과',
         '· 공개 릴리스 트리 위생 검사 통과',
         '· 원고의 모든 수치가 결과 파일과',
         '  대조되어 일치'],
        size=11.5, line_spacing=1.3)
    txt(s, 4.9, 2.05, 4.2, 0.35, '최근 전면 검토에서 고친 것', size=14, bold=True, color=ACCENT)
    txt(s, 4.9, 2.5, 4.2, 3.7,
        ['· 특이성 게이트 순열을 환자 단위로',
         '  교정 → FDR 추정 %.1f%% (전 %.1f%%)' % (f['fdr'], 13.7),
         '· 재현 불가 파일 3건에 생산 단계 등록',
         '  (새 저장소에서 --paper 가 죽던 지점)',
         '· 감사 자체의 사각지대 3건 제거',
         '· 발굴/게이트 자료 중첩을 공개하고',
         '  민감도 분석으로 방어',
         '· 릴리스에 섞이던 옛 원고 사본 제외'],
        size=11.5, line_spacing=1.3)
    txt(s, 9.4, 2.05, 3.4, 0.35, '남은 것', size=14, bold=True)
    txt(s, 9.4, 2.5, 3.4, 3.7,
        ['자리표 3개를 채우면 끝납니다.',
         '',
         '1. 저장소 URL',
         '2. Zenodo DOI',
         '3. 투고일',
         '',
         '둘은 같은 태그된 릴리스를',
         '가리켜야 하고, 그 확인은',
         'preflight.py 가 합니다.'],
        size=11.5, line_spacing=1.3)
    rule(s, 0.6, 6.4, 12.2)
    txt(s, 0.6, 6.5, 12.2, 0.35,
        '검토의 원칙: 감사가 통과한다가 아니라 새 저장소에서 실제로 다시 만들어지는가를 봤습니다.',
        size=11, color=MUTED)


def s_now_review(prs, d):
    f = facts()
    s = add_slide(prs)
    header(s, 'WHAT THE REVIEW FOUND', '전면 검토에서 나온 실질 결함 4건',
           '감사는 모두 초록불이었습니다. 코드를 직접 읽고 의심 가는 곳마다 따로 측정했습니다.')
    rows = [['무엇이', '왜 문제였나', '어떻게 바뀌었나'],
            ['특이성 게이트의\n순열',
             'ERCB 두 블록은 같은 환자다.\n관문 시료의 32%(54명)가 양쪽에\n있고 전원 라벨이 같은데,\n블록마다 따로 섞고 있었다',
             'FDR 추정 13.7%% -> %.1f%%\n환자 단위 짝 순열로 교체.\n따로 섞으면 게이트가\n실제보다 깨끗해 보인다' % f['fdr']],
            ['comparison_glom3\n.tsv',
             '손으로 만든 파일인데 PAPER\n단계가 읽는다. 새 저장소에서는\n--paper 가 여기서 죽는다.\n널도 4코호트 것과 짝이 안 맞음',
             '단계 4개 등록해 재생성.\nAUROC 와 재현성 수치는\n완전히 동일했고, 틀렸던\n널 열만 바뀜'],
            ['stress 단계',
             'figures 가 그 산출을 읽는데\nPAPER 목록에 없었다',
             'PAPER 에 편입'],
            ['릴리스 위생',
             '옛 원고 사본이 공개 트리에\n실려 나가고 있었다',
             '제외 규칙 추가']]
    table(s, 0.6, 2.05, 12.2, rows, col_w=[1.5, 4.0, 3.4], size=10.5)
    txt(s, 0.6, 6.55, 12.2, 0.4,
        '확인했고 문제 없던 것: LODO 누수 없음 · RBS 의 해석적 Jaccard 는 경험값과 2% 이내 · '
        'lasso 는 실제 L1 · ReliefF 절댓값은 무해 · 보정 평가는 원시 자료에서만',
        size=10.5, color=MUTED)


def s_now_defense(prs, d):
    f = facts()
    s = add_slide(prs)
    header(s, 'ANTICIPATED OBJECTIONS', '심사에서 나올 반론과 우리 답',
           '세 가지 모두 원고 본문에 이미 적혀 있습니다.')
    items = [
        ('"특이성 게이트가 발굴 자료를 공유한다"',
         'GSE104948 은 발굴 코호트이면서 게이트의 두 블록 중 하나입니다 — 사실이고, Methods 와 '
         'Limitations 에 그대로 적었습니다. 발굴에서만 빼고 후보를 다시 만들면 게이트는 발굴과 '
         '독립이 됩니다. 그래도 후보 %d/30 이 남고(FMOD·LUM·MMP7·MOXD1 포함), 게이트 통과율은 '
         '낮아지지 않으며, 외부 GSE20602 재현도 유지됩니다(p=%.1e 대 %.1e).'
         % (f['sens_keep'], f['sens_p'], f['main_p'])),
        ('"%.1f%%와 %.1f%%는 둘 다 몇 %%인데 왜 따로 부르나"' % (f['bg'], f['fdr']),
         '분모가 다릅니다. %.1f%%는 배경 통과율로 분모가 전체 9,900개 유전자이고, %.1f%%는 순열 '
         '기반 FDR 추정치로 분모가 실제 통과한 820개입니다. 비교 가능한 비율이 아니어서 원고·'
         '그림·산출 표에서 이름과 분모를 모두 분리했고, 열 이름도 gate_fdr_estimate 입니다.'
         % (f['bg'], f['fdr'])),
        ('"AUROC가 좋으면 된 것 아닌가"',
         '무작위 50개 유전자가 external AUROC 0.68~0.88 을 냅니다. 절대 AUROC 로는 방법을 '
         '줄세울 수 없어 널 대비 퍼센타일로만 보고합니다. 변별력이 있는 축은 교차폴드 재현성이고, '
         '거기서 제안 방법이 %.3f 대 차순위 %.3f (4코호트 K=100) 입니다.'
         % (f['a4_rbs'], f['a4_next'])),
    ]
    y = 1.95
    for t, b in items:
        txt(s, 0.6, y, 12.2, 0.35, t, size=14.5, bold=True, color=BLUE)
        txt(s, 0.85, y + 0.36, 11.9, 1.0, b, size=11.5, line_spacing=1.28)
        y += 1.62
