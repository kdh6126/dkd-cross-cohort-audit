#!/usr/bin/env python
"""Cross-omics corroboration of the final candidates against KPMP.

Scope is set by what KPMP actually has (see docs/worklog/KPMP_ACCESS.md): regional proteomics and
regional transcriptomics contain ZERO diabetic participants, and the API's `dmr` slice
returns cell counts with no statistics. So KPMP cannot validate a DKD-vs-control claim.

What it CAN do, and what is asked of it here:
  * cell-type attribution   which kidney cell type expresses the gene (sn/sc, pooled cohort)
  * compartment attribution glomerulus vs tubulointerstitium at RNA level, and - the valuable
                            part - at PROTEIN level, from an independent cohort and assay

The second one is real cross-omics evidence for the compartment analysis, which is why it is
worth doing even though the DKD-specific contrast is unavailable.
"""
import os, sys, time, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kpmp_client import (data_types_for_gene, single_cell, regional_proteomics,
                         regional_transcriptomics_by_structure)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', default='results/final_candidates_gated.tsv')
    ap.add_argument('--topn', type=int, default=30)
    ap.add_argument('--out', default='results/kpmp_corroboration.tsv')
    args = ap.parse_args()

    cand = pd.read_csv(args.candidates, sep='\t').head(args.topn)
    genes = [g for g in cand['symbol'].dropna().astype(str) if g]
    log('corroborating %d candidates against KPMP\n' % len(genes))

    log('fetching regional transcriptomics for the whole glomerulus (one call, ~26k genes) ...')
    rt = {}
    try:
        for r in regional_transcriptomics_by_structure('Glomerulus'):
            rt[r['geneSymbol']] = r
        log('  got %d genes\n' % len(rt))
    except Exception as e:
        log('  ! failed: %s\n' % str(e)[:120])

    rows = []
    log('%-10s %-14s %-34s %-9s %s'
        % ('gene', 'KPMP layers', 'top cell type (sn, all)', 'sn FC', 'protein Glom-vs-TI'))
    for g in genes:
        rec = dict(symbol=g)
        try:
            rec['layers'] = ','.join(data_types_for_gene(g) or [])
        except Exception:
            rec['layers'] = ''
        top_ct, top_fc = '', np.nan
        try:
            sc = [r for r in single_cell(g, 'sn', 'all')
                  if r.get('pValAdj') is not None and r['pValAdj'] < 0.05]
            if sc:
                b = max(sc, key=lambda r: r.get('foldChange') or -9e9)
                top_ct, top_fc = (b['clusterName'] or '')[:34], b['foldChange']
                rec['sn_n_sig_clusters'] = len(sc)
        except Exception:
            pass
        rec['sn_top_celltype'] = top_ct
        rec['sn_top_foldchange'] = top_fc

        prot = ''
        try:
            rp = regional_proteomics(g)
            glom = [r for r in rp if r.get('region') == 'Glom']
            if glom:
                r0 = glom[0]
                rec['protein_glom_vs_TI_FC'] = r0['foldChange']
                rec['protein_adjP'] = r0['adjPVal']
                rec['protein_peptides'] = r0['numPeptides']
                prot = 'FC=%+.2f adjP=%.1e' % (r0['foldChange'], r0['adjPVal'])
        except Exception:
            pass

        if g in rt:
            rec['rt_glom_FC'] = rt[g].get('foldChange')
            rec['rt_glom_adjP'] = rt[g].get('adjPVal')

        log('%-10s %-14s %-34s %-9s %s'
            % (g, rec['layers'], top_ct, ('%+.2f' % top_fc) if not np.isnan(top_fc) else '-', prot))
        rows.append(rec)
        time.sleep(0.2)

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    n_prot = int(df.get('protein_glom_vs_TI_FC', pd.Series(dtype=float)).notna().sum())
    n_cell = int((df['sn_top_celltype'].astype(str) != '').sum())
    log('\n%d/%d candidates have KPMP protein-level evidence; %d have a significant cell-type assignment'
        % (n_prot, len(df), n_cell))
    log('wrote %s' % args.out)


if __name__ == '__main__':
    main()
