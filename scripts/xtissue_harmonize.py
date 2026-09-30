#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""조직 확장 코호트를 받아 Entrez 유전자 공간의 행렬과 표본 표로 만든다.

신장 쪽 조화(build_gene_space.py, harmonize.py)와 같은 규칙을 쓴다. 여러 유전자에 걸리는
프로브는 버리고, 한 유전자에 프로브가 여럿이면 평균 발현이 가장 높은 것 하나를 남긴다
(max-mean). 다른 점은 플랫폼이 많아 주석을 해석하는 경로를 자동으로 고른다는 것뿐이다.

    1. GEO 의 .annot.gz 가 있으면 'Gene ID' 열을 쓴다.
    2. 없으면 플랫폼 전체 표에서 Entrez 열, gene_assignment, 유전자 기호, RefSeq 순으로 찾는다.
       기호는 NCBI gene_info 로, RefSeq 는 refGene 의 기호를 거쳐 Entrez 로 바꾼다.

값의 척도가 로그가 아니면(99번째 백분위수가 100 초과) log2(x+1) 로 바꾼다. 분석은 어차피
코호트 안에서 유전자별 z 로 가므로, 여기서는 척도를 로그로 맞추는 데까지만 한다.

출력은 data/processed/xtissue/<domain>/<GSE>_expr.tsv, _pheno.tsv 와
results/xtissue/<domain>/cohorts.tsv 이다. 제외된 코호트도 이유와 함께 표에 남는다.
"""
import argparse
import gzip
import io
import os
import re
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xtissue_config import COHORTS, arm_match   # noqa: E402

HDR = {'User-Agent': 'dkd-xtissue/1.0 (mailto:etri.jhdh@gmail.com)'}
RAW = 'data/raw/geo'
ANN = 'data/raw/annotation'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url, dest=None):
    for k in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=HDR), timeout=300) as r:
                data = r.read()
            if dest:
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                with open(dest + '.part', 'wb') as f:
                    f.write(data)
                os.replace(dest + '.part', dest)
            return data
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(3 + 5 * k)
    raise err


# ---------------------------------------------------------------- series matrix
def matrix_path(gse, platform):
    d = os.path.join(RAW, gse)
    stub = gse[:-3] + 'nnn'
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/' % (stub, gse)
    names = sorted(set(re.findall(r'href="(%s[^"]*_series_matrix\.txt\.gz)"' % gse,
                                  get(url).decode('utf-8', 'replace'))))
    if len(names) > 1:
        names = [n for n in names if platform in n] or names
    name = names[0]
    dest = os.path.join(d, name)
    if not os.path.exists(dest):
        log('  download %s' % name)
        get(url + name, dest)
    return dest


def parse_matrix(path):
    hdr, rows, cols = {}, [], None
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        intable = False
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                intable = True
                cols = next(f).rstrip('\n').split('\t')
                cols = [c.strip('"') for c in cols]
                continue
            if line.startswith('!series_matrix_table_end'):
                break
            if intable:
                rows.append(line.rstrip('\n').split('\t'))
            elif line.startswith('!'):
                p = line.rstrip('\n').split('\t')
                hdr.setdefault(p[0][1:], []).append([x.strip('"') for x in p[1:]])
    ids = hdr['Sample_geo_accession'][0]
    titles = hdr.get('Sample_title', [[''] * len(ids)])[0]
    source = hdr.get('Sample_source_name_ch1', [[''] * len(ids)])[0]
    desc = hdr.get('Sample_description', [[''] * len(ids)])[0]
    chars = [{} for _ in ids]
    for row in hdr.get('Sample_characteristics_ch1', []):
        for i, v in enumerate(row):
            if ':' in v:
                k, val = v.split(':', 1)
                chars[i][k.strip()] = val.strip()
    samples = [dict(gsm=g, title=t, source=s, chars=c, description=d)
               for g, t, s, c, d in zip(ids, titles, source, chars, desc)]
    df = pd.DataFrame(rows, columns=cols)
    df['ID_REF'] = df['ID_REF'].str.strip('"')
    df = df.set_index('ID_REF')
    df = df.apply(pd.to_numeric, errors='coerce')
    return samples, df


# ---------------------------------------------------------------- annotation
_GI = None


def gene_info():
    global _GI
    if _GI is None:
        sym, syn = {}, {}
        with gzip.open(os.path.join(ANN, 'Homo_sapiens.gene_info.gz'), 'rt',
                       encoding='utf-8', errors='replace') as f:
            hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
            ix = {c: i for i, c in enumerate(hdr)}
            for line in f:
                p = line.rstrip('\n').split('\t')
                sym[p[ix['Symbol']]] = p[ix['GeneID']]
                for a in p[ix['Synonyms']].split('|'):
                    if a and a != '-' and a not in syn:
                        syn[a] = p[ix['GeneID']]
        _GI = (sym, syn)
    return _GI


def sym2entrez(s):
    sym, syn = gene_info()
    s = s.strip()
    return sym.get(s) or syn.get(s)


_ENS = None


def ensembl2entrez(ens):
    """gene_info 의 dbXrefs 에 적힌 Ensembl:ENSG... 로 맞춘다. 한 ENSG 가 여러 Entrez 에
    걸리면 버린다."""
    global _ENS
    if _ENS is None:
        m, bad = {}, set()
        with gzip.open(os.path.join(ANN, 'Homo_sapiens.gene_info.gz'), 'rt',
                       encoding='utf-8', errors='replace') as f:
            hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
            ix = {c: i for i, c in enumerate(hdr)}
            for line in f:
                p = line.rstrip('\n').split('\t')
                for x in p[ix['dbXrefs']].split('|'):
                    if x.startswith('Ensembl:'):
                        e = x.split(':', 1)[1]
                        if e in m and m[e] != p[ix['GeneID']]:
                            bad.add(e)
                        m[e] = p[ix['GeneID']]
        _ENS = {k: v for k, v in m.items() if k not in bad}
    return _ENS.get(ens.split('.')[0])


_RG = None


def refseq2entrez(acc):
    global _RG
    if _RG is None:
        _RG = {}
        with gzip.open(os.path.join(ANN, 'refGene_hg19.txt.gz'), 'rt') as f:
            for line in f:
                p = line.split('\t')
                _RG.setdefault(p[1], p[12])
    s = _RG.get(acc.split('.')[0])
    return sym2entrez(s) if s else None


def platform_table(gpl):
    """GEO .annot.gz 가 있으면 그것을, 없으면 SOFT 플랫폼 표를 받아 DataFrame 으로."""
    stub = gpl[:-3] + 'nnn' if len(gpl) > 6 else 'GPLnnn'
    annot = os.path.join(ANN, '%s.annot.gz' % gpl)
    if not os.path.exists(annot):
        try:
            get('https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s.annot.gz'
                % (stub, gpl, gpl), annot)
        except Exception:  # noqa: BLE001
            annot = None
    if annot and os.path.exists(annot):
        with gzip.open(annot, 'rt', encoding='utf-8', errors='replace') as f:
            txt = f.read()
        body = txt.split('!platform_table_begin\n', 1)[-1].split('!platform_table_end', 1)[0]
        return pd.read_csv(io.StringIO(body), sep='\t', dtype=str, low_memory=False), 'annot'
    soft = os.path.join(ANN, '%s.soft.txt' % gpl)
    if not os.path.exists(soft):
        get('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=data'
            % gpl, soft)
    txt = open(soft, encoding='utf-8', errors='replace').read()
    body = txt.split('!platform_table_begin\n', 1)[-1].split('!platform_table_end', 1)[0]
    return pd.read_csv(io.StringIO(body), sep='\t', dtype=str, low_memory=False), 'soft'


def probe_map(gpl):
    cache = os.path.join(ANN, '%s.entrez.tsv' % gpl)
    if os.path.exists(cache):
        m = pd.read_csv(cache, sep='\t', dtype=str)
        return dict(zip(m['probe'], m['entrez'])), open(cache + '.route').read()
    t, kind = platform_table(gpl)
    t.columns = [c.strip() for c in t.columns]
    low = {c.lower(): c for c in t.columns}
    idcol = low.get('id')

    def single(ids):
        ids = {i for i in ids if i}
        return ids.pop() if len(ids) == 1 else None

    route, out = None, {}
    for key in ('gene id', 'entrez_gene_id', 'entrez_gene', 'gene'):
        if key in low and t[low[key]].dropna().str.match(r'^\d+(\s*///\s*\d+)*$').mean() > 0.5:
            c = low[key]
            for p, v in zip(t[idcol], t[c]):
                if isinstance(v, str):
                    g = single(x.strip() for x in v.split('///') if x.strip().isdigit())
                    if g:
                        out[p] = g
            route = '%s:%s' % (kind, c)
            break
    if not out and 'gene_assignment' in low:
        c = low['gene_assignment']
        for p, v in zip(t[idcol], t[c]):
            if isinstance(v, str) and v != '---':
                ids = set()
                for rec in v.split('///'):
                    fld = [x.strip() for x in rec.split('//')]
                    if len(fld) >= 5 and fld[4].isdigit():
                        ids.add(fld[4])
                g = single(ids)
                if g:
                    out[p] = g
        route = '%s:gene_assignment' % kind
    if not out:
        for key in ('gene symbol', 'symbol', 'gene_symbol', 'ilmn_gene', 'genesymbol'):
            if key in low:
                c = low[key]
                for p, v in zip(t[idcol], t[c]):
                    if isinstance(v, str) and v.strip():
                        g = single(sym2entrez(x) for x in re.split(r'\s*///\s*', v))
                        if g:
                            out[p] = g
                route = '%s:%s->gene_info' % (kind, c)
                break
    if not out:
        for key in ('gb_acc', 'refseq', 'gb_list'):
            if key in low:
                c = low[key]
                for p, v in zip(t[idcol], t[c]):
                    if isinstance(v, str) and v.strip():
                        g = single(refseq2entrez(x.strip()) for x in re.split(r'[,;/ ]+', v)
                                   if x.strip())
                        if g:
                            out[p] = g
                route = '%s:%s->refGene->gene_info' % (kind, c)
                break
    if not out:
        raise RuntimeError('%s: no usable gene column in %s' % (gpl, list(t.columns)[:20]))
    pd.DataFrame(dict(probe=list(out), entrez=list(out.values()))).to_csv(cache, sep='\t',
                                                                          index=False)
    open(cache + '.route', 'w').write(route)
    return out, route


# ---------------------------------------------------------------- per cohort
def counts_matrix(gse, cfg, samples):
    """보충 파일의 원시 계수를 받아 표본 GSM 열로 바꾼 DataFrame 을 돌려준다."""
    name = cfg['counts']
    dest = os.path.join(RAW, gse, name)
    if not os.path.exists(dest):
        get('https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s'
            % (gse[:-3] + 'nnn', gse, name), dest)
    sep = ',' if name.endswith('.csv.gz') else '\t'
    c = pd.read_csv(dest, sep=sep, index_col=0)
    c.columns = [str(x).strip('"') for x in c.columns]
    colmap = {}
    for s in samples:
        if cfg.get('count_column'):
            key = cfg['count_column'](s)
        else:
            key = s['description'].strip().zfill(4)
        if key in c.columns:
            colmap[key] = s['gsm']
    c = c[list(colmap)].rename(columns=colmap)
    return c


def harmonise_counts(gse, cfg, samples):
    c = counts_matrix(gse, cfg, samples)
    lib = c.sum(axis=0)
    cpm = c / lib * 1e6
    keep = (cpm >= 1).mean(axis=1) >= 0.2          # 표본의 20% 이상에서 CPM 1 이상
    logcpm = np.log2(cpm.loc[keep] + 1)
    ids = [ensembl2entrez(i) for i in logcpm.index]
    logcpm = logcpm.loc[[i is not None for i in ids]]
    logcpm.insert(0, 'entrez', [i for i in ids if i is not None])
    route = 'counts:Ensembl->gene_info dbXrefs; filter CPM>=1 in >=20%; log2(CPM+1)'
    return logcpm, route


def harmonise(domain, gse, cfg, outdir):
    if cfg['label'] is None:
        return dict(gse=gse, status=cfg['status'], n_case=0, n_control=0)
    path = matrix_path(gse, cfg['platform'])
    samples, df = parse_matrix(path)
    if cfg.get('counts'):
        df, route = harmonise_counts(gse, cfg, samples)
        pmap = None
    else:
        pmap, route = probe_map(cfg['platform'])

    # 표본 규칙
    rows = []
    seen = set()
    for s in samples:
        lab = cfg['label'](s)
        if lab is None:
            continue
        subj = cfg['subject'](s) if cfg['subject'] else s['gsm']
        if (subj, lab) in seen:
            continue
        seen.add((subj, lab))
        rows.append(dict(sample=s['gsm'], label=lab, subject_id=subj, title=s['title'],
                         source=s['source']))
    ph = pd.DataFrame(rows)
    if pmap is None:
        ent = df['entrez']
        df = df[[c for c in ph['sample'] if c in df.columns]]
        ph = ph[ph['sample'].isin(df.columns)]
        scale = 'log2(CPM+1)'
        df.insert(0, 'entrez', ent)
    else:
        df = df[[c for c in ph['sample'] if c in df.columns]]
        ph = ph[ph['sample'].isin(df.columns)]

        # 척도
        vals = df.values[np.isfinite(df.values)]
        p99 = float(np.percentile(vals, 99)) if len(vals) else float('nan')
        scale = 'as deposited'
        if p99 > 100:
            df = np.log2(df.clip(lower=0) + 1)
            scale = 'log2(x+1) applied (p99=%.0f)' % p99

        # 프로브 -> 유전자
        df = df.loc[[i for i in df.index if i in pmap]]
        df.insert(0, 'entrez', [pmap[i] for i in df.index])
    # max-mean
    df['_mean'] = df.drop(columns='entrez').mean(axis=1)
    df = df.sort_values('_mean', ascending=False).drop_duplicates('entrez')
    expr = df.drop(columns='_mean').set_index('entrez')
    expr = expr.loc[expr.notna().mean(axis=1) >= 0.8]

    os.makedirs(outdir, exist_ok=True)
    expr.to_csv(os.path.join(outdir, '%s_expr.tsv' % gse), sep='\t')
    ph.to_csv(os.path.join(outdir, '%s_pheno.tsv' % gse), sep='\t', index=False)
    pc, pk, ev = cfg['procurement']
    return dict(gse=gse, platform=cfg['platform'], status=cfg['status'],
                n_case=int((ph['label'] == 1).sum()), n_control=int((ph['label'] == 0).sum()),
                n_genes=int(expr.shape[0]), annotation_route=route, scale=scale,
                procurement_case=pc, procurement_control=pk, arms=arm_match((pc, pk)),
                evidence=ev)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--domain', required=True, choices=sorted(COHORTS))
    ap.add_argument('--only', default='')
    args = ap.parse_args()
    outdir = 'data/processed/xtissue/%s' % args.domain
    resdir = 'results/xtissue/%s' % args.domain
    os.makedirs(resdir, exist_ok=True)
    only = set(x for x in args.only.split(',') if x)
    rows = []
    for gse, cfg in COHORTS[args.domain].items():
        if only and gse not in only:
            continue
        log('== %s %s' % (gse, cfg['platform']))
        try:
            r = harmonise(args.domain, gse, cfg, outdir)
        except Exception as e:  # noqa: BLE001
            r = dict(gse=gse, platform=cfg['platform'], status='failed: %s' % str(e)[:200])
        log('   %s' % {k: v for k, v in r.items() if k != 'evidence'})
        rows.append(r)
    tab = pd.DataFrame(rows)
    path = os.path.join(resdir, 'cohorts.tsv')
    if only and os.path.exists(path):
        old = pd.read_csv(path, sep='\t')
        tab = pd.concat([old[~old['gse'].isin(tab['gse'])], tab], ignore_index=True)
    tab.to_csv(path, sep='\t', index=False)
    log('wrote %s' % path)


if __name__ == '__main__':
    main()
