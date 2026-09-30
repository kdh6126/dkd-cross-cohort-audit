#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""표에 이름이 오른 단백질을 정말로 그 유전자의 단백질로 잰 것인가.

단백체 심사자가 가장 먼저 묻는 것입니다. 두 기술 모두 유전자에 이름을 붙이는 방식이
간접적이기 때문입니다.

    SOMAscan   앱타머 하나가 표적 단백질에 붙는다고 **가정**합니다. 결합 특이성은
               제조사 주석에 달려 있고, 우리가 자료에서 확인할 수 있는 것은 그
               앱타머가 어떤 UniProt 항목에 배정되었는가까지입니다.
    질량분석   펩타이드를 단백질에 배정합니다. 어떤 펩타이드는 여러 단백질에 공통이라
               하나로 좁혀지지 않습니다. 그런 경우 accession 칸에 여러 항목이 파이프로
               묶여 나옵니다. 그 묶음을 한 유전자로 읽으면 안 됩니다.

그래서 유전자마다 다음을 적습니다.

    앱타머 식별자와 표적 설명, 배정된 UniProt
    질량분석의 accession 묶음 크기, 펩타이드 행 수, 고유 서열 수
    묶음이 둘 이상이면 그 사실과 함께 묶인 상대

이 감사에서 실제로 하나가 걸렸습니다. ACTN1 의 묶음에는 ACTN4 가 함께 들어 있습니다.
ACTN4 는 발세포 유전자로 신장에서 잘 알려져 있으므로, 그 신호를 ACTN1 것이라고 단정할
수 없습니다. 표에서 지우지 않고 단서를 붙여 남깁니다.
"""
import os
import sys

import pandas as pd

SOMA_X = 'data/raw/proteomics/mendeley_83k89shdx5/SOMAscan.xlsx'
MS_X = ('data/raw/proteomics/pride_PXD041884/'
        'Secocnd_Extraction_GPQ_Norm_TS_All_Peptides_FFPE_Protein_Rollup_Workbook.xlsx')
MAP = 'data/raw/proteomics/pride_PXD041884/uniprot_to_gene.tsv'
COV = 'results/proteome_meta/candidate_protein_coverage.tsv'
OUT = 'results/proteome_meta'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def extra_symbols(pep, known):
    """펩타이드 그룹에만 나타나는 accession 의 기호를 UniProt 에서 받아 온다.

    매핑 파일은 롤업된 단백질만 담고 있어, 그룹에만 나타나는 accession 은 비어 있습니다.
    그대로 두면 표에 "Q08043 (Q08043)" 처럼 accession 이 기호 자리에도 찍힙니다.
    조회에 실패하면 비워 두고 넘어갑니다. 이 감사가 네트워크 때문에 멈출 이유는 없습니다.
    """
    import io as _io
    import urllib.parse
    import urllib.request
    want = sorted({x for w in pep['Accession # Whole'].astype(str)
                   for x in w.split('|') if x and x not in known})
    out = {}
    for i in range(0, len(want), 40):
        q = ' OR '.join('accession:%s' % a for a in want[i:i + 40])
        u = ('https://rest.uniprot.org/uniprotkb/search?'
             + urllib.parse.urlencode({'query': q, 'fields': 'accession,gene_primary',
                                       'format': 'tsv', 'size': '500'}))
        try:
            t = pd.read_csv(_io.StringIO(
                urllib.request.urlopen(u, timeout=60).read().decode()), sep='\t')
        except Exception:
            continue
        for _, r in t.iterrows():
            g = str(r.get('Gene Names (primary)') or '').strip()
            if g and g.lower() != 'nan':
                out[str(r['Entry']).strip()] = g.split()[0]
    return out


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    cov = pd.read_csv(COV, sep='\t')

    soma = pd.read_excel(SOMA_X)
    soma['gene'] = soma['Entrez Gene Symbol'].astype(str).str.strip()
    ms_roll = pd.read_excel(MS_X, sheet_name='All Protien Rolled Up')
    ms_pep = pd.read_excel(MS_X, sheet_name='All')
    u2g = dict(pd.read_csv(MAP, sep='\t').values)
    ms_roll['gene'] = ms_roll['Accession #'].astype(str).str.strip().map(u2g)
    u2g.update(extra_symbols(ms_pep, u2g))

    rows = []
    for _, r in cov.iterrows():
        g = r['gene']
        d = dict(gene=g)

        # ---- 앱타머 쪽
        s = soma[soma['gene'] == g]
        d['on_somascan'] = bool(r['on_somascan'])
        d['somamer_id'] = '|'.join(s['SOMAmer seqID'].astype(str)) if len(s) else ''
        d['somascan_uniprot'] = '|'.join(s['UniProt ID'].astype(str)) if len(s) else ''
        d['somascan_target'] = '; '.join(s['SomaLogic Target Description'].astype(str)) \
            if len(s) else ''
        d['n_somamers'] = len(s)

        # ---- 질량분석 쪽
        m = ms_roll[ms_roll['gene'] == g]
        d['on_massspec'] = bool(r['on_massspec'])
        if len(m):
            acc = str(m['Accession #'].iloc[0]).strip()
            d['ms_accession'] = acc
            pep = ms_pep[ms_pep['Accession # Delimited'].astype(str).str.strip() == acc]
            # 'Accession # Whole' 는 그 펩타이드가 배정될 수 있는 단백질 전부를 파이프로
            # 잇는다. 값이 이 accession 하나뿐인 행만 그 단백질에 고유한 펩타이드다.
            whole = pep['Accession # Whole'].astype(str).str.strip()
            uniq = pep[whole == acc]
            partners = set()
            for w in whole[whole != acc]:
                partners.update(x for x in w.split('|') if x and x != acc)
            d['ms_accession'] = acc
            d['ms_peptides_total'] = int(pep['Sequence'].nunique())
            d['ms_peptides_unique'] = int(uniq['Sequence'].nunique())
            d['ms_peptides_shared'] = d['ms_peptides_total'] - d['ms_peptides_unique']
            # 파이프로 이으면 LaTeX 가 끊을 곳을 찾지 못해 표가 넘친다. 쉼표로 잇는다.
            d['ms_shared_with'] = ', '.join(
                '%s (%s)' % (a, u2g.get(a, a)) for a in sorted(partners))
            d['ms_rollup_peptide_count'] = float(m['Peptide Count'].iloc[0])
        else:
            for k in ('ms_accession', 'ms_shared_with'):
                d[k] = ''
            for k in ('ms_peptides_total', 'ms_peptides_unique', 'ms_peptides_shared',
                      'ms_rollup_peptide_count'):
                d[k] = 0

        flags = []
        if d['n_somamers'] > 1:
            flags.append('%d aptamers map to this gene' % d['n_somamers'])
        if d['ms_peptides_shared'] > d['ms_peptides_unique']:
            flags.append('most peptides are shared with %s'
                         % (d['ms_shared_with'].split(',')[0] or 'other proteins'))
        elif d['ms_peptides_shared']:
            flags.append('%d of %d peptides shared'
                         % (d['ms_peptides_shared'], d['ms_peptides_total']))
        if d['on_massspec'] and d['ms_peptides_unique'] < 3:
            flags.append('fewer than 3 unique peptides')
        d['caveat'] = '; '.join(flags)
        rows.append(d)

    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'protein_identity_audit.tsv'), sep='\t', index=False)

    log('=' * 78)
    log('단백질 식별성 감사 — 표에 오른 유전자 %d개' % len(t))
    log('')
    log('  %-9s %-15s %-9s %6s %6s  %s'
        % ('유전자', '앱타머', 'MS acc', '고유펩', '공유펩', '단서'))
    for _, r in t.iterrows():
        log('  %-9s %-15s %-9s %6d %6d  %s'
            % (r['gene'], (r['somamer_id'] or '-')[:15], (r['ms_accession'] or '-')[:9],
               r['ms_peptides_unique'], r['ms_peptides_shared'], r['caveat'][:40]))
    log('')
    bad = t[t['caveat'] != '']
    log('  단서가 붙은 유전자 %d개' % len(bad))
    for _, r in bad.iterrows():
        log('    %-9s %s' % (r['gene'], r['caveat']))
        if r['ms_shared_with']:
            log('              공유 상대: %s' % r['ms_shared_with'])
    log('')
    log('  펩타이드가 공유되면 그 신호를 한 유전자 것이라고 단정할 수 없습니다.')
    log('  지우지 않고 단서를 붙여 남깁니다 — 심사자가 물을 것을 먼저 적는 편이 낫습니다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
