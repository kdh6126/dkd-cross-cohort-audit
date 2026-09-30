#!/usr/bin/env python
"""Novelty assessment for the candidate genes, via PubMed E-utilities.

Why not a web search: "no results" from a search engine is not evidence of absence, and it
already misled us once - FMOD looked novel until a targeted query surfaced a 2021 paper naming
FMOD and LUM together as DKD hub genes. PubMed counts are auditable and reproducible.

Three queries per gene, from broad to specific:
  ANY        gene name anywhere                       -> is the gene studied at all?
  KIDNEY     gene AND kidney/renal/nephropathy        -> is it studied in kidney?
  DKD        gene AND diabetic kidney disease terms   -> is it already a DKD gene?
  BIOMARKER  gene AND DKD AND biomarker               -> is it already proposed as a marker?

A gene is only a defensible novel candidate if DKD count is 0 (or the hits are incidental).
PMIDs are saved so every claim can be checked.
"""
import os, sys, time, json, argparse, urllib.request, urllib.parse
import pandas as pd

import datetime as _dt

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
TODAY = _dt.date.today().isoformat()

DKD_TERMS = ('("diabetic nephropathy"[Title/Abstract] OR "diabetic kidney disease"[Title/Abstract] '
             'OR "diabetic kidney"[Title/Abstract])')
KIDNEY_TERMS = ('(kidney[Title/Abstract] OR renal[Title/Abstract] OR nephropathy[Title/Abstract] '
                'OR glomerul*[Title/Abstract] OR podocyte[Title/Abstract])')
BIOMARKER_TERMS = '(biomarker*[Title/Abstract] OR marker*[Title/Abstract])'


# Symbols that collide with a common non-gene abbreviation, so PubMed counts for the bare
# symbol are meaningless. Verified by reading the returned titles: every "MSC + diabetic
# nephropathy" hit is about mesenchymal stem/stromal cells, not the gene musculin.
AMBIGUOUS = {
    'MSC': 'collides with "mesenchymal stem/stromal cell"',
    'CFD': 'collides with "computational fluid dynamics"; DKD hits are complement panels, indirect',
    'CA10': 'collides with carbonic-anhydrase numbering and "CA 10" measurements',
    'MS4A6A': 'bare-symbol count inconsistent with the AND query; counts unreliable',
}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def count(term, retmax=8):
    url = EUTILS + '?' + urllib.parse.urlencode(
        {'db': 'pubmed', 'term': term, 'retmax': retmax, 'retmode': 'json'})
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                d = json.loads(r.read().decode('utf-8', 'replace'))
            es = d.get('esearchresult', {})
            return int(es.get('count', 0)), es.get('idlist', [])
        except Exception:
            time.sleep(2)
    return -1, []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', default='results/final_candidates_gated.tsv')
    ap.add_argument('--topn', type=int, default=30)
    ap.add_argument('--out', default='results/literature_review.tsv')
    args = ap.parse_args()

    cand = pd.read_csv(args.candidates, sep='\t').head(args.topn)
    genes = [g for g in cand['symbol'].dropna().astype(str) if g]
    log('querying PubMed for %d genes\n' % len(genes))
    log('%-10s %8s %8s %8s %8s   verdict' % ('gene', 'any', 'kidney', 'DKD', 'DKD+bm'))

    rows = []
    for g in genes:
        gq = '"%s"[Title/Abstract]' % g
        n_any, _ = count(gq); time.sleep(0.35)
        n_kid, _ = count('%s AND %s' % (gq, KIDNEY_TERMS)); time.sleep(0.35)
        n_dkd, ids = count('%s AND %s' % (gq, DKD_TERMS)); time.sleep(0.35)
        n_bm, bm_ids = count('%s AND %s AND %s' % (gq, DKD_TERMS, BIOMARKER_TERMS)); time.sleep(0.35)

        if g in AMBIGUOUS:
            v = 'AMBIGUOUS SYMBOL - counts unreliable'
        elif n_bm > 0:
            v = 'ALREADY PROPOSED as DKD marker'
        elif n_dkd >= 3:
            v = 'established in DKD'
        elif n_dkd > 0:
            v = 'reported in DKD (few papers)'
        elif n_kid > 0:
            v = 'kidney-studied, NOT in DKD'
        else:
            v = 'NOVEL - no kidney literature'
        log('%-10s %8d %8d %8d %8d   %s' % (g, n_any, n_kid, n_dkd, n_bm, v))
        # 조회한 날짜를 결과에 남긴다. PubMed 는 시간이 지나면 늘어나므로, 날짜가 없는
        # 건수는 나중에 다시 돌렸을 때 왜 달라졌는지 설명할 수 없다.
        rows.append(dict(gene=g, queried_on=TODAY, n_any=n_any, n_kidney=n_kid, n_dkd=n_dkd,
                         n_dkd_biomarker=n_bm, verdict=v,
                         ambiguity_note=AMBIGUOUS.get(g, ''),
                         dkd_pmids=';'.join(ids), dkd_biomarker_pmids=';'.join(bm_ids)))

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)

    log('\n=== summary ===')
    for v, d in df.groupby('verdict'):
        log('  %-32s %2d  %s' % (v, len(d), ', '.join(d['gene'])))
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
