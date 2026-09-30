#!/usr/bin/env python
"""What was acquired, what was actually used, and what is sitting unused.

Written because "6.4 GB downloaded" says nothing about how much of it entered an analysis.
Usage is determined from the pipeline, not from intent: a dataset counts as USED only if a
harmonised matrix or a result file derives from it.
"""
import os, sys, glob, json, argparse
import pandas as pd

H = 'data/processed/harmonized'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def du(path):
    if not os.path.exists(path):
        return 0
    if os.path.isfile(path):
        return os.path.getsize(path)
    t = 0
    for r, _, fs in os.walk(path):
        for f in fs:
            try:
                t += os.path.getsize(os.path.join(r, f))
            except OSError:
                pass
    return t


def gb(n):
    return n / 1e9


# asset -> (path, role, used_by)
ASSETS = [
    # --- transcriptomics actually in the pipeline
    ('GSE142025', 'data/raw/geo/GSE142025', 'RNA-seq, whole cortex, 27 DKD / 9 ctrl',
     'LODO member (all4), harmonised'),
    ('GSE30528', 'data/raw/geo/GSE30528', 'array, glomeruli, 9 DKD / 13 ctrl',
     'LODO member, harmonised'),
    ('GSE96804', 'data/raw/geo/GSE96804', 'array, glomeruli, 41 DKD / 20 ctrl',
     'LODO member, harmonised'),
    ('GSE104948', 'data/raw/geo/GSE104948', 'array, glomeruli, 12 DKD / 26 ctrl / 158 other CKD',
     'LODO member + specificity gate, harmonised'),
    ('GSE30529', 'data/raw/geo/GSE30529', 'array, tubulointerstitium, 10 DKD / 12 ctrl',
     'compartment validation'),
    ('GSE104954', 'data/raw/geo/GSE104954', 'array, tubulointerstitium, 17 DKD / 26 ctrl / 152 other CKD',
     'compartment validation + specificity gate'),
    ('GSE294519', 'data/raw/geo/GSE294519',
     'RNA-seq, tubulointerstitium, 23 DKD / 13 tumour-nephrectomy control',
     'fifth cohort, mapped into the frozen gene space at 98.5% coverage'),
    # --- excluded by the audit
    ('GSE30122', 'data/raw/geo/GSE30122', 'array, superseries', 'EXCLUDED - duplicates GSE30528+GSE30529'),
    ('GSE47183', 'data/raw/geo/GSE47183', 'array, ERCB glomeruli', 'EXCLUDED - DN subjects identical to GSE104948'),
    ('GSE99340', 'data/raw/geo/GSE99340', 'array, ERCB + cell lines', 'EXCLUDED - ERCB superset'),
    # --- downloaded, not used
    ('GSE131882', 'data/raw/geo/GSE131882', 'snRNA-seq, 3 DKD / 3 ctrl', 'UNUSED'),
    ('GSE175759', 'data/raw/geo/GSE175759', 'RNA-seq tubulointerstitium, 90 samples',
     'procurement replication (independent platform)'),
    ('GSE142153', 'data/raw/geo/GSE142153', 'array, PBMC blood, 23 DKD / 10 ctrl / 7 ESRD', 'UNUSED'),
    ('GSE1009', 'data/raw/geo/GSE1009', 'array, 3v3, HG-U95Av2', 'UNUSED - too small / obsolete platform'),
    ('GSE111154', 'data/raw/geo/GSE111154', 'array, early DN, 4v4', 'UNUSED - too small'),
    ('GSE20602', 'data/raw/geo/GSE20602', 'array, nephrosclerosis', 'UNUSED - not DKD'),
    ('GSE163603', 'data/raw/geo/GSE163603',
     'RNA-seq, 레이저 미세절단 구획별, DKD 6명 / Reference 9명',
     'preanalytical test - the only cohort besides GSE162830 that reports storage'),
    # --- confounder controls acquired later
    ('GSE162830', 'data/raw/geo/GSE162830', 'RNA-seq, histology-matched DN vs ING',
     'histology-matched contrast (the one dataset with no procurement effect)'),
    ('GSE5406', 'data/raw/geo/GSE5406', 'array, human heart failure vs donor',
     'external organ control - bounds the claim to kidney procurement'),
    # --- KPMP
    ('KPMP v1.5 snRNA', 'data/raw/kpmp/KPMP_v1.5_snRNA_human_kidney.h5ad',
     'single-nucleus atlas, 304,989 nuclei',
     'donor-level pseudobulk negative control (DKD vs non-diabetic CKD biopsy)'),
    ('KPMP v1.5 scRNA', 'data/raw/kpmp/KPMP_v1.5_scRNA_human_kidney.h5ad',
     'single-cell atlas, 225,177 cells', 'UNUSED (the GraphQL API was used instead)'),
    ('KPMP Atlas API', None, 'per-gene cell type + regional proteomics',
     'USED - 30 candidates corroborated'),
    # --- other omics
    ('GWAS Catalog', 'data/raw/gwas', 'GCST90018612/832, GCST90179152, GCST005886, GCST90134333-4',
     'germline layer, SNP-count-matched nulls'),
    # --- 단백체 (전사체 후보의 교차 오믹스 확인. 두 코호트를 층화 결합한다)
    ('Mendeley 83k89shdx5', 'data/raw/proteomics/mendeley_83k89shdx5',
     'SOMAscan 신장 피질, DKD 23 / 정상 10',
     'USED - 후보 7개를 쟀고 4개가 q<0.05'),
    ('PRIDE PXD041884', 'data/raw/proteomics/pride_PXD041884',
     'LC-MS/MS FFPE 신장 피질, DKD 5 / 비당뇨 7',
     'USED - 후보 7개를 쟀고 2개가 q<0.05. 양쪽 조달 경로가 같다'),
    # --- 대사체 (개별 등록: 디스크 대조에서 누락이 잡혀 추가)
    ('ST003255', 'data/raw/metabolomics/ST003255',
     'LC-MS 혈장, DKD 64 / IgAN 66 / MN 66 / HN 24 / 정상 66',
     'second-modality contrast; the confounder-matched design in metabolomics'),
    ('ST001411', 'data/raw/metabolomics/ST001411',
     'LC-MS 혈장, 당뇨+신경병증 48 / 당뇨만 49 / 정상 9',
     'generalisation test: does the matched-comparator recommendation hold elsewhere'),
    ('ST000691', 'data/raw/metabolomics/ST000691', 'DKD 대사체 n=17',
     'UNUSED - no comparison group'),
    ('ST004442', 'data/raw/metabolomics/ST004442', 'DKD 혈장 대사체 2025',
     'UNUSED - raw spectra only, no quantified matrix'),
    ('ST004483', 'data/raw/metabolomics/ST004483', 'DKD 소변 대사체 2025',
     'UNUSED - raw spectra only, no quantified matrix'),
    ('ST002145', 'data/raw/metabolomics/ST002145', 'ChREBP 관련 대사체',
     'UNUSED - not a human DKD case/control design'),
    # --- annotation
    ('GPL571 / GPL17586 / gene_info', 'data/raw/annotation', 'probe -> Entrez mapping',
     'USED - gene space construction'),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='results/data_inventory.tsv')
    args = ap.parse_args()

    rows = []
    for name, path, role, usage in ASSETS:
        size = du(path) if path else 0
        # Default to 'used' and name the exceptions. The reverse -- whitelisting phrases that
        # count as used -- silently demoted assets whenever a usage string was reworded.
        state = ('excluded' if usage.startswith('EXCLUDED') else
                 'not_downloaded' if usage.startswith('NOT DOWNLOADED') else
                 'unused' if usage.startswith('UNUSED') else 'used')
        rows.append(dict(asset=name, size_gb=round(gb(size), 3), role=role,
                         usage=usage, state=state, path=path or ''))
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)

    # 디스크에 있는데 이 표에 없는 자산을 찾는다. GSE163603과 대사체 6건이 이런 식으로
    # 조용히 빠져 있었고, 그중 하나는 논문의 한계 서술을 바꾸는 자료였다.
    import glob as _glob
    import re as _re
    listed = set()
    for a in df['asset'].astype(str):
        listed |= set(_re.findall('(GSE[0-9]+|ST[0-9]+)', a))
    disk = {os.path.basename(d) for d in _glob.glob('data/raw/geo/GSE*') if os.path.isdir(d)}
    disk |= {os.path.basename(d) for d in _glob.glob('data/raw/metabolomics/ST*')
             if os.path.isdir(d)}
    missing = sorted(disk - listed)
    if missing:
        log('')
        log('  경고 - 디스크에 있으나 인벤토리에 없는 자산 %d건: %s'
            % (len(missing), ', '.join(missing)))
        log('  ASSETS 목록에 추가하세요. 등록되지 않은 자산은 집계에 잡히지 않습니다.')

    log('%-24s %9s  %s' % ('asset', 'size(GB)', 'state / usage'))
    for st in ('used', 'excluded', 'unused', 'not_downloaded'):
        d = df[df['state'] == st]
        log('\n--- %s : %d assets, %.2f GB ---' % (st.upper(), len(d), d['size_gb'].sum()))
        for _, r in d.iterrows():
            log('  %-24s %9.2f  %s' % (r['asset'], r['size_gb'], r['usage']))

    tot = df['size_gb'].sum()
    used = df[df['state'] == 'used']['size_gb'].sum()
    unused = df[df['state'] == 'unused']['size_gb'].sum()
    excl = df[df['state'] == 'excluded']['size_gb'].sum()
    log('\n' + '=' * 78)
    log('downloaded total        %6.2f GB' % tot)
    log('  entered the analysis  %6.2f GB  (%.0f%%)' % (used, 100 * used / max(tot, 1e-9)))
    log('  downloaded but unused %6.2f GB  (%.0f%%)' % (unused, 100 * unused / max(tot, 1e-9)))
    log('  excluded by the audit %6.2f GB  (%.0f%%)' % (excl, 100 * excl / max(tot, 1e-9)))

    # sample-level accounting from the harmonised matrices
    log('\n=== samples actually modelled ===')
    tot_s = tot_c = 0
    for f in sorted(glob.glob(os.path.join(H, '*_pheno.tsv'))):
        p = pd.read_csv(f, sep='\t')
        n = len(p)
        lab = int((p['label'] == 1).sum())
        ctl = int((p['label'] == 0).sum())
        oth = int((p['label'] == -1).sum())
        tot_s += lab + ctl
        tot_c += lab
        log('  %-11s %3d samples: %2d DKD, %2d control, %3d other-CKD'
            % (os.path.basename(f).split('_')[0], n, lab, ctl, oth))
    log('  ---> %d labelled samples in the models, of which %d DKD cases' % (tot_s, tot_c))
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
