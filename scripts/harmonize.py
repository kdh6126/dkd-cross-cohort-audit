#!/usr/bin/env python
"""Collapse probes to genes and emit harmonised gene x sample matrices + label tables.

Collapse rule: MAX-MEAN. When several probes map to the same Entrez gene, keep the single
probe with the highest mean expression across that cohort's samples (WGCNA collapseRows
"MaxMean"). Fixed up front and applied identically in every cohort so that LODO folds all
operate on the same feature definition.

Outputs (data/processed/harmonized/):
  <cohort>_expr.tsv    genes (Entrez) x samples, restricted to the common gene space
  <cohort>_pheno.tsv   sample, group, label, platform, subject_id, control_type
"""
import gzip, os, sys, csv, re
import numpy as np
import pandas as pd

OUT = 'data/processed/harmonized'
os.makedirs(OUT, exist_ok=True)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def read_series_matrix(path):
    """Return (DataFrame probes x samples, dict of header lines)."""
    hdr = {}
    rows, cols = [], None
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!'):
                p = line.rstrip('\n').split('\t')
                hdr.setdefault(p[0][1:], []).append([x.strip('"') for x in p[1:]])
        cols = [x.strip('"') for x in f.readline().rstrip('\n').split('\t')][1:]
        for line in f:
            if line.startswith('!series_matrix_table_end'):
                break
            p = line.rstrip('\n').split('\t')
            rows.append((p[0].strip('"'), p[1:]))
    idx = [r[0] for r in rows]
    data = np.array([[np.nan if v in ('', 'null', 'NA') else float(v) for v in r[1]] for r in rows],
                    dtype=float)
    return pd.DataFrame(data, index=idx, columns=cols), hdr


def load_probe_map():
    m = {}
    with open('data/processed/probe2gene.tsv', encoding='utf-8') as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            m.setdefault(row['cohort'], {})[row['probe_id']] = row['entrez_id']
    return m


def collapse_max_mean(df, probe2gene, gene_space):
    """probes x samples -> genes x samples using the max-mean probe per gene."""
    mapped = [p for p in df.index if probe2gene.get(p) in gene_space]
    sub = df.loc[mapped]
    genes = np.array([probe2gene[p] for p in mapped])
    means = sub.mean(axis=1, skipna=True).values
    order = np.lexsort((-means, genes))          # per gene, highest mean first
    keep, seen = [], set()
    for i in order:
        g = genes[i]
        if g not in seen:
            seen.add(g)
            keep.append(i)
    out = sub.iloc[keep].copy()
    out.index = genes[keep]
    out.index.name = 'entrez_id'
    return out.sort_index(key=lambda ix: ix.astype(int))


def pheno_from_header(hdr, n):
    """sample titles + flattened characteristics per sample."""
    def flat(key):
        v = hdr.get(key, [[]])
        return v[0] if v else [''] * n
    titles = flat('Sample_title')
    gsms = flat('Sample_geo_accession')
    srcs = flat('Sample_source_name_ch1')
    chars = hdr.get('Sample_characteristics_ch1', [])
    out = []
    for i in range(n):
        c = ' | '.join(x[i] for x in chars if i < len(x) and x[i])
        out.append({'gsm': gsms[i] if i < len(gsms) else '',
                    'title': titles[i] if i < len(titles) else '',
                    'source': srcs[i] if i < len(srcs) else '',
                    'characteristics': c})
    return out


# ---------------------------------------------------------------- label rules
def label_GSE142025(p):
    g = re.search(r'group:\s*(\S+)', p['characteristics'])
    g = g.group(1) if g else ''
    return {'group': g, 'label': 0 if g == 'Control' else 1,
            'stage': g, 'control_type': 'tumor_nephrectomy' if g == 'Control' else '',
            'subject_id': p['title']}


def label_woroniecka(p):
    c = p['characteristics'].lower()
    dkd = 'diabetic kidney disease' in c or 'dkd' in p['source'].lower()
    sid = re.search(r'Kidney (\d+)', p['title'])
    return {'group': 'DKD' if dkd else 'control', 'label': 1 if dkd else 0, 'stage': '',
            'control_type': '' if dkd else 'tumor_nephrectomy',
            'subject_id': sid.group(1) if sid else p['title']}


def label_GSE96804(p):
    dn = 'diabetic' in p['source'].lower()
    return {'group': 'DN' if dn else 'control', 'label': 1 if dn else 0, 'stage': '',
            'control_type': '' if dn else 'tumor_nephrectomy', 'subject_id': p['title']}


def label_ercb(p):
    t = p['title']
    diag = re.search(r'diagnosis:\s*([^|]+)', p['characteristics'])
    diag = diag.group(1).strip() if diag else ''
    if re.search(r'-LD', t, re.I) or (not diag and '-TN' not in t.upper()):
        return {'group': 'living_donor', 'label': 0, 'stage': '',
                'control_type': 'living_donor', 'subject_id': t}
    if 'tumor nephrectomy' in diag.lower():
        return {'group': 'tumor_nephrectomy', 'label': 0, 'stage': '',
                'control_type': 'tumor_nephrectomy', 'subject_id': t}
    if 'diabetic' in diag.lower():
        return {'group': 'DN', 'label': 1, 'stage': '', 'control_type': '', 'subject_id': t}
    return {'group': diag or 'other', 'label': -1, 'stage': '', 'control_type': '', 'subject_id': t}


COHORTS = {
    'GSE142025': dict(kind='tsv', path='data/processed/GSE142025_expr_matrix.tsv',
                      labeller=label_GSE142025, compartment='whole_cortex'),
    'GSE30528':  dict(kind='sm', paths=[('GPL571', 'data/raw/geo/GSE30528/GSE30528_series_matrix.txt.gz')],
                      labeller=label_woroniecka, compartment='glomerulus'),
    'GSE96804':  dict(kind='sm', paths=[('GPL17586', 'data/raw/geo/GSE96804/GSE96804_series_matrix.txt.gz')],
                      labeller=label_GSE96804, compartment='glomerulus'),
    'GSE104948': dict(kind='sm', paths=[('GPL22945', 'data/raw/geo/GSE104948/GSE104948-GPL22945_series_matrix.txt.gz'),
                                        ('GPL24120', 'data/raw/geo/GSE104948/GSE104948-GPL24120_series_matrix.txt.gz')],
                      labeller=label_ercb, compartment='glomerulus'),
    'GSE30529':  dict(kind='sm', paths=[('GPL571', 'data/raw/geo/GSE30529/GSE30529_series_matrix.txt.gz')],
                      labeller=label_woroniecka, compartment='tubulointerstitium'),
    'GSE104954': dict(kind='sm', paths=[('GPL22945', 'data/raw/geo/GSE104954/GSE104954-GPL22945_series_matrix.txt.gz'),
                                        ('GPL24120', 'data/raw/geo/GSE104954/GSE104954-GPL24120_series_matrix.txt.gz')],
                      labeller=label_ercb, compartment='tubulointerstitium'),
}


def main():
    gene_space = set(pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t',
                                 dtype=str)['entrez_id'])
    log('gene space: %d genes' % len(gene_space))
    p2g = load_probe_map()
    summary = []

    for name, cfg in COHORTS.items():
        pm = p2g[name]
        if cfg['kind'] == 'tsv':
            df = pd.read_csv(cfg['path'], sep='\t', index_col=0)
            phenos = [{'gsm': '', 'title': c, 'source': '', 'characteristics': ''} for c in df.columns]
            inv = {r['title']: r for r in csv.DictReader(
                open('data/metadata/sample_inventory.tsv', encoding='utf-8'), delimiter='\t')
                if r['gse'] == name}
            for p in phenos:
                r = inv.get(p['title'])
                if r:
                    p.update({'gsm': r['gsm'], 'source': r['source'], 'characteristics': r['characteristics']})
            expr = collapse_max_mean(df, pm, gene_space)
            plats = ['Illumina_HiSeq4000'] * df.shape[1]
        else:
            parts, phenos, plats = [], [], []
            for plat, path in cfg['paths']:
                d, hdr = read_series_matrix(path)
                parts.append(collapse_max_mean(d, pm, gene_space))
                ph = pheno_from_header(hdr, d.shape[1])
                phenos += ph
                plats += [plat] * d.shape[1]
            expr = pd.concat(parts, axis=1) if len(parts) > 1 else parts[0]

        rows = []
        for p, pl in zip(phenos, plats):
            lab = cfg['labeller'](p)
            lab.update({'sample': p['gsm'] or p['title'], 'title': p['title'],
                        'platform': pl, 'compartment': cfg['compartment'], 'cohort': name})
            rows.append(lab)
        pheno = pd.DataFrame(rows)[['cohort', 'sample', 'title', 'subject_id', 'group', 'label',
                                    'stage', 'control_type', 'platform', 'compartment']]
        expr.columns = pheno['sample'].values

        expr.to_csv(os.path.join(OUT, '%s_expr.tsv' % name), sep='\t')
        pheno.to_csv(os.path.join(OUT, '%s_pheno.tsv' % name), sep='\t', index=False)

        vals = expr.values[np.isfinite(expr.values)]
        n_case = int((pheno['label'] == 1).sum())
        n_ctrl = int((pheno['label'] == 0).sum())
        n_excl = int((pheno['label'] == -1).sum())
        summary.append(dict(cohort=name, genes=expr.shape[0], samples=expr.shape[1],
                            case=n_case, control=n_ctrl, other=n_excl,
                            vmin=round(float(vals.min()), 2), vmed=round(float(np.median(vals)), 2),
                            vmax=round(float(vals.max()), 2), compartment=cfg['compartment']))
        log('%-11s %5d genes x %3d samples | case=%-3d ctrl=%-3d other=%-3d | range %.2f..%.2f'
            % (name, expr.shape[0], expr.shape[1], n_case, n_ctrl, n_excl, vals.min(), vals.max()))

    pd.DataFrame(summary).to_csv('data/processed/harmonized_summary.tsv', sep='\t', index=False)
    log('\nwrote %s/ and data/processed/harmonized_summary.tsv' % OUT)


if __name__ == '__main__':
    main()
