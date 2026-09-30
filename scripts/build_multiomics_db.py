#!/usr/bin/env python
"""다중 오믹스 확장 스키마를 적용하고, 지금 채울 수 있는 것을 채운다. 그리고 엑셀로 낸다.

지금 채울 수 있는 것은 전사체와 대사체뿐입니다. 나머지 칸은 비워두되, "왜 비어 있는가"를
modality_inventory에 기록합니다. 그것이 이 대장의 핵심 용도입니다 — 다음 사람이 같은 검색을
반복하고 같은 막다른 길에 도달하지 않게 하는 것.

엑셀은 DDL에서 자동 생성합니다. 손으로 만든 엑셀과 실제 스키마는 반드시 어긋나므로,
한쪽에서만 정의하고 다른 쪽은 생성물로 둡니다.
"""
import os
import sqlite3
import sys
from datetime import date

import pandas as pd

DB = 'db/dkd.sqlite'
DDL = 'db/schema_multiomics.sql'
XLSX = 'db/multiomics_schema.xlsx'
TODAY = date.today().isoformat()


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ------------------------------------------------------------------ 채울 내용

def modality_rows():
    """지금까지 조사한 모달리티별 가용성. 못 쓴 것도 이유와 함께 남긴다."""
    return [
        # (modality, disease, resource, status, n_case, n_control, matched, blocker, url)
        ('transcriptomics', 'DKD', 'GSE30528/96804/104948/142025/294519', 'used',
         112, 81, 0, None, 'https://www.ncbi.nlm.nih.gov/geo/'),
        ('transcriptomics', 'DKD', 'GSE104948/104954 (ERCB other-CKD arms)', 'used',
         29, 310, 1, None, 'https://www.ncbi.nlm.nih.gov/geo/'),
        ('transcriptomics', 'DKD', 'GSE175759', 'used', 3, 84, 1, None,
         'https://www.ncbi.nlm.nih.gov/geo/'),
        ('single-nucleus', 'DKD', 'KPMP v1.5 snRNA', 'used', 27, 28, 1, None,
         'https://atlas.kpmp.org/'),
        ('proteomics', 'DKD', 'KPMP regional proteomics', 'unusable', 0, None, 0,
         '당뇨병 참여자 0명. 대조 자체가 성립하지 않음.', 'https://atlas.kpmp.org/'),
        ('metabolomics', 'DKD', 'ST003255 (Seoul National University)', 'used',
         64, 222, 1, None, 'https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST003255'),
        ('metabolomics', 'DKD', 'ST004442 / ST004483', 'unusable', None, None, 0,
         '원시 스펙트럼만 등록. 정량 매트릭스 없음.',
         'https://www.metabolomicsworkbench.org/'),
        ('metabolomics', 'DKD', 'ST000691', 'unusable', 17, 0, 0,
         'n=17이고 비교군 없음.', 'https://www.metabolomicsworkbench.org/'),
        ('metabolomics', 'CKD', 'Metabolomics Workbench "chronic kidney" 검색', 'available_unused',
         None, None, None, '16건 확인. 아직 검토 안 함.',
         'https://www.metabolomicsworkbench.org/'),
        ('metabolomics', 'T2D', 'ST003390', 'available_unused', 100, 200, 0,
         '신장 표현형 없음. DKD 질문에 직접 답하지 못함.',
         'https://www.metabolomicsworkbench.org/'),
        ('metabolomics', 'diabetic_neuropathy', 'ST001411 (Univ. of Michigan)',
         'available_unused', 48, 58, 1,
         '"당뇨는 같고 합병증만 다른" 대조가 내장됨. 설계가 가장 좋음.',
         'https://www.metabolomicsworkbench.org/'),
        ('genomics', 'DKD', 'GWAS Catalog (GCST90018612 외)', 'used', None, None, 0,
         None, 'https://www.ebi.ac.uk/gwas/'),
        ('proteomics', 'DKD', 'PRIDE 인간 DKD 소변/혈장', 'not_found', None, None, None,
         '검토하지 않음. 큐레이션 필요.', 'https://www.ebi.ac.uk/pride/'),
        ('multi-omics (paired)', 'DKD', '동일 환자 다층 공개 데이터', 'not_found',
         None, None, None,
         '공개 데이터에는 없음. UK Biobank류가 현실적 경로이며 별도 승인 필요.', None),
    ]


def contrast_rows():
    """이 프로젝트가 실제로 계산한 대조들. 모달리티가 달라도 같은 구조로 적는다."""
    return [
        ('ST003255', 'DKD vs healthy', 'DKD', 'healthy', 0, 64, 66, 'plasma',
         '기존 방식. 대조군이 질병 자체가 없음.'),
        ('ST003255', 'DKD vs other kidney disease', 'DKD', 'other_kidney_disease', 1,
         64, 156, 'plasma', '교란 매칭. 단, 세 질환을 묶은 것이 결과를 만들었음.'),
        ('ST003255', 'DKD vs membranous nephropathy', 'DKD', 'other_kidney_disease', 1,
         64, 66, 'plasma', '단일 질환 대조. 88개 유의 — 통합(29개)보다 많음.'),
        ('ST003255', 'DKD vs IgA nephropathy', 'DKD', 'other_kidney_disease', 1,
         64, 66, 'plasma', '단일 질환 대조. 16개 유의.'),
        ('GSE104948', 'DKD vs other CKD (glomerulus)', 'DKD', 'other_kidney_disease', 1,
         12, 158, 'kidney_glomerulus', '특이성 게이트의 근거.'),
        ('GSE104954', 'DKD vs other CKD (tubulointerstitium)', 'DKD',
         'other_kidney_disease', 1, 17, 152, 'kidney_tubulointerstitium', None),
    ]


def feature_link_rows():
    """지금 성립하는 모달리티 간 연결. 하나뿐이다."""
    genes = ['COL1A2', 'LUM', 'FMOD', 'THBS2', 'VCAN', 'MMP2', 'MMP7']
    return [(g, '4-Hydroxyproline', 'pathway', 'degrades' if g.startswith('MMP') else 'member_of',
             'manual_curation',
             '하이드록시프롤린은 콜라겐에만 있는 아미노산이라 혈장 농도가 콜라겐 회전율을 '
             '반영한다. 이 유전자군은 콜라겐 자체이거나 그 결합·분해 단백질이다.',
             'medium') for g in genes]


# ------------------------------------------------------------------ 적용

def main():
    if not os.path.exists(DB):
        log('  %s 가 없습니다. 먼저 run_all.py --stage database 를 실행하세요.' % DB)
        return 1
    con = sqlite3.connect(DB)
    con.executescript(open(DDL, encoding='utf-8').read())
    cur = con.cursor()

    # -------- modality_inventory
    cur.execute('DELETE FROM modality_inventory')
    cur.executemany(
        'INSERT INTO modality_inventory (modality, disease, resource, status, n_case, '
        'n_control, has_matched_control, blocker, url, checked_at) '
        'VALUES (?,?,?,?,?,?,?,?,?,?)',
        [r + (TODAY,) for r in modality_rows()])

    # -------- feature: 전사체 후보 + 대사체 패널
    cur.execute('DELETE FROM feature')
    genes = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    cur.executemany(
        'INSERT OR IGNORE INTO feature (modality, id_type, id_value, symbol, entrez_id) '
        'VALUES (?,?,?,?,?)',
        [('transcriptomics', 'entrez', r.entrez_id, r.symbol, int(r.entrez_id))
         for r in genes.itertuples() if str(r.entrez_id).isdigit()])
    met_file = 'data/raw/metabolomics/ST003255/AN005337_datatable.tsv'
    n_met = 0
    if os.path.exists(met_file):
        mets = open(met_file, encoding='utf-8').readline().rstrip('\n').split('\t')[2:]
        cur.executemany(
            'INSERT OR IGNORE INTO feature (modality, id_type, id_value, symbol) '
            'VALUES (?,?,?,?)',
            [('metabolomics', 'name', m, m) for m in mets])
        n_met = len(mets)

    # -------- study: 전사체 파이프라인이 등록하지 않은 비-GEO 자원을 먼저 넣는다.
    # 이것을 빠뜨려 ST003255의 대조 4건이 조용히 누락된 적이 있다.
    cur.execute("INSERT OR IGNORE INTO source_repository (name, base_url) VALUES (?,?)",
                ('MetabolomicsWorkbench', 'https://www.metabolomicsworkbench.org/'))
    rid = cur.execute("SELECT repository_id FROM source_repository WHERE name=?",
                      ('MetabolomicsWorkbench',)).fetchone()[0]
    cur.execute(
        "INSERT OR IGNORE INTO study (repository_id, accession, title, omics_type, assay, "
        "n_samples, local_path, downloaded_at) VALUES (?,?,?,?,?,?,?,?)",
        (rid, 'ST003255',
         'Metabolomics study to identify cause-specific biomarkers for Chronic Kidney Disease',
         'metabolomics', 'lc-ms', 286,
         'data/raw/metabolomics/ST003255/AN005337_datatable.tsv', TODAY))

    # -------- contrast
    cur.execute('DELETE FROM contrast')
    acc2id = dict(cur.execute('SELECT accession, study_id FROM study').fetchall())
    missing_acc = []
    for acc, name, case, ctrl, matched, nc, nk, tissue, note in contrast_rows():
        sid = acc2id.get(acc)
        if sid is None:
            missing_acc.append(acc)
            continue
        cur.execute(
            'INSERT OR IGNORE INTO contrast (study_id, name, case_group, control_group, '
            'control_is_matched, n_case, n_control, tissue, note) VALUES (?,?,?,?,?,?,?,?,?)',
            (sid, name, case, ctrl, matched, nc, nk, tissue, note))

    # -------- feature_link
    cur.execute('DELETE FROM feature_link')
    n_link = 0
    for gsym, msym, level, rel, src, detail, conf in feature_link_rows():
        a = cur.execute("SELECT feature_id FROM feature WHERE modality='transcriptomics' "
                        "AND symbol=?", (gsym,)).fetchone()
        b = cur.execute("SELECT feature_id FROM feature WHERE modality='metabolomics' "
                        "AND symbol=?", (msym,)).fetchone()
        if a and b:
            cur.execute(
                'INSERT OR IGNORE INTO feature_link (from_feature_id, to_feature_id, '
                'link_level, relation, evidence_source, evidence_detail, confidence, '
                'created_at) VALUES (?,?,?,?,?,?,?,?)',
                (a[0], b[0], level, rel, src, detail, conf, TODAY))
            n_link += 1

    con.commit()

    # ------------------------------------------------------------------ 보고
    log('=' * 74)
    log('다중 오믹스 스키마 적용')
    log('=' * 74)
    tables = ['feature', 'feature_link', 'pathway', 'pathway_member', 'contrast',
              'contrast_result', 'crossmodal_support', 'modality_inventory']
    for t in tables:
        n = cur.execute('SELECT COUNT(*) FROM %s' % t).fetchone()[0]
        log('  %-22s %6d 행%s' % (t, n, '   (비어 있음 — 나중에 채울 자리)' if n == 0 else ''))
    if missing_acc:
        log('')
        log('  study 테이블에 없어 건너뛴 accession: %s' % ', '.join(sorted(set(missing_acc))))

    log('')
    log('  모달리티 가용성 요약')
    for row in cur.execute(
            'SELECT modality, disease, n_used, n_available, n_unusable, n_not_found, '
            'n_with_matched_control FROM v_modality_gaps ORDER BY modality, disease'):
        log('    %-22s %-20s 사용 %d · 미검토 %d · 불가 %d · 미발견 %d · 매칭대조 %s'
            % (row[0], row[1], row[2], row[3], row[4], row[5], row[6] or 0))

    # ------------------------------------------------------------------ 엑셀
    with pd.ExcelWriter(XLSX, engine='openpyxl') as xw:
        readme = pd.DataFrame({
            '항목': ['목적', '작성 규칙', '채워진 시트', '빈 시트', '주의'],
            '내용': [
                '공개 데이터에서 여러 오믹스가 같은 환자에서 나오지 않으므로, sample 단위로 '
                '조인할 수 없습니다. 이 파일은 조인 가능한 세 층위(분자/경로/대조)를 각각 '
                '시트로 만든 것입니다.',
                '이 엑셀은 db/schema_multiomics.sql 에서 자동 생성됩니다. 여기서 고친 내용은 '
                '다시 생성하면 사라집니다. 구조를 바꾸려면 DDL을 고치세요. 데이터를 채우려면 '
                '시트에 적은 뒤 build_multiomics_db.py 에 로더를 추가하세요.',
                'modality_inventory, contrast, feature_link 는 현재까지 확인된 내용이 '
                '들어 있습니다.',
                'pathway, pathway_member, contrast_result, crossmodal_support 는 비어 '
                '있습니다. 구조만 잡아둔 것입니다.',
                'crossmodal_support 의 support_level 은 네 값만 허용됩니다. 등급을 임의로 '
                '만들면 DB가 거부합니다.',
            ]})
        readme.to_excel(xw, sheet_name='README', index=False)
        for t in tables:
            df = pd.read_sql_query('SELECT * FROM %s' % t, con)
            if df.empty:
                cols = [r[1] for r in cur.execute('PRAGMA table_info(%s)' % t)]
                df = pd.DataFrame(columns=cols)
            df.to_excel(xw, sheet_name=t[:31], index=False)
        for v in ('v_modality_gaps', 'v_matched_contrast', 'v_crossmodal_summary'):
            pd.read_sql_query('SELECT * FROM %s' % v, con).to_excel(
                xw, sheet_name=v[:31], index=False)

    con.close()
    log('')
    log('  wrote %s' % XLSX)
    log('  전사체 %d개 + 대사체 %d개 분자 등록, 모달리티 간 연결 %d건'
        % (len(genes), n_met, n_link))
    return 0


if __name__ == '__main__':
    sys.exit(main())
