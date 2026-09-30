#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""두 단백체 코호트가 정확히 무엇인가 — 원논문과 저장소 메타데이터를 그대로 옮긴다.

이 표가 필요한 이유는 앞서 이 자료들을 두 번 부정확하게 서술했기 때문입니다.

  "same case/control contrast"   두 코호트의 대조군 정의가 같지 않습니다. 한쪽은
                                 건강 참여자, 다른 쪽은 비당뇨 병리 사례입니다.
  "procurement asymmetry is not  조달 경로가 같다는 것은 설계에서 생검 대 신절제
   built in"                     대비가 없다는 뜻이지, 허혈시간이나 고정 시간을
                                 실제로 쟀다는 뜻이 아닙니다.

그래서 여기서는 '무엇이 문헌에 적혀 있는가' 와 '무엇이 적혀 있지 않은가' 를 나눠
적습니다. 적혀 있지 않은 것을 추론으로 채우지 않습니다.

한 가지 더 밝혀야 할 것이 있습니다. SOMAscan 코호트의 원논문(Hirohama 2023)의
제목이 곧 "MMP7 을 신장질환 바이오마커로 식별한다" 입니다. 우리가 그 코호트에서
MMP7 을 유의하다고 재는 것은 독립 발견이 아니라 같은 자료에서 같은 결과를 다시
계산한 것입니다. 새로운 것은 **두 번째 코호트(질량분석)에서도 재현된다**는 쪽입니다.
이 구분을 표에 넣습니다.
"""
import os
import sys

import pandas as pd

OUT = 'results/proteome_meta'

# 값은 전부 아래 출처의 문장에서 옮긴 것입니다. 계산으로 얻은 값이 아닙니다.
#   Hirohama D et al. J Am Soc Nephrol 2023;34:1279-1291.  PMID 37022120
#   Schwab SK et al. Proteomics Clin Appl 2024;18:e202400018.  PMID 38923810
#   PRIDE PXD041884 sampleProcessingProtocol (저장소 메타데이터)
#   Li D et al. Diabetes 2024;73:1188-1195.  PMID 38394643  (자료를 재사용한 논문)
ROWS = [
    ('Repository / accession',
     'Mendeley Data 83k89shdx5',
     'PRIDE PXD041884'),
    ('Source publication',
     'Hirohama 2023 (JASN)',
     'Schwab 2024 (Proteomics Clin Appl)'),
    ('Measurement principle',
     'SomaScan aptamer panel, 1,305 proteins',
     'Label-free LC-MS/MS, Q-TOF, shotgun'),
    ('Case group as defined by the source',
     '23 individuals with diabetic kidney disease',
     '5 clinically validated diabetic cases, all with nodular sclerosis'),
    ('Control group as defined by the source',
     '10 healthy controls',
     '7 non-diabetic control cases'),
    ('Sample labels in the released matrix',
     'Con_01..Con_10, DKD_01..DKD_23',
     'C1, C2d, C3, C4, C6, C7, C9 / D1, D3, D4, D5, D6'),
    ('Tissue',
     'Kidney cortex',
     'Kidney FFPE tissue curls, over 60% cortex per case'),
    ('Procurement route stated for the cases',
     'Kidney tissue; route not stated in the accessible text',
     'Explant tissue, archival FFPE, one institution'),
    ('Procurement route stated for the controls',
     'Kidney tissue; route not stated in the accessible text',
     'Explant tissue, archival FFPE, same institution'),
    ('Biopsy-versus-nephrectomy asymmetry by design',
     'Cannot be established from the accessible text',
     'Absent: both arms are archival explant cases from one source'),
    ('Pre-analytical exposure reported',
     'No ischaemia, fixation or storage time reported',
     'No ischaemia, fixation or storage time reported'),
    ('Already-published protein finding in this cohort',
     'MMP7 is the source publication\'s headline biomarker',
     'None of our candidates is a named finding of the source'),
]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    os.makedirs(OUT, exist_ok=True)

    d = pd.DataFrame(ROWS, columns=['attribute', 'somascan_cohort', 'massspec_cohort'])
    d.to_csv(os.path.join(OUT, 'cohort_audit.tsv'), sep='\t', index=False)

    log('=' * 78)
    log('두 단백체 코호트 감사 — 문헌에 적힌 것만')
    log('=' * 78)
    for a, b, c in ROWS:
        log('  %s' % a)
        log('     SOMAscan  %s' % b)
        log('     질량분석  %s' % c)
    log('')
    log('여기서 나오는 결론')
    log('  1  대조군 정의가 다릅니다. "같은 case/control 대조" 라고 쓰면 안 됩니다.')
    log('     "DKD 대 비DKD" 가 두 코호트를 함께 덮는 정확한 표현입니다.')
    log('  2  질량분석 코호트는 설계상 생검 대 신절제 대비가 없습니다. 다만 전분석')
    log('     노출을 잰 것은 아니므로 "조달이 맞춰졌다" 가 아니라 "그 대비가 설계에')
    log('     없다" 까지만 말할 수 있습니다.')
    log('  3  SOMAscan 코호트에서의 MMP7 은 그 자료의 원논문이 이미 낸 결과입니다.')
    log('     독립 확인은 질량분석 코호트 쪽입니다.')
    log('')
    log('  %s 에 표를 썼습니다.' % os.path.join(OUT, 'cohort_audit.tsv'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
