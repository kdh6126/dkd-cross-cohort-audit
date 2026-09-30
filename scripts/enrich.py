#!/usr/bin/env python
"""GO / KEGG / Reactome enrichment for a candidate signature, via Enrichr.

Uses Enrichr's **custom background** endpoint (speedrichr). This matters: the analysis space
here is the 9,900-gene cross-cohort intersection, not the whole genome. Testing a signature
drawn from 9,900 genes against a ~20,000-gene background inflates every p-value, because the
background contains genes that could never have been selected. The standard Enrichr endpoint
does exactly that, so it is used only as a fallback.
"""
import os, sys, json, time, argparse
import urllib.request
import pandas as pd

SPEED = 'https://maayanlab.cloud/speedrichr/api'
CLASSIC = 'https://maayanlab.cloud/Enrichr'
LIBRARIES = ['GO_Biological_Process_2023', 'KEGG_2021_Human',
             'Reactome_2022', 'GO_Molecular_Function_2023']


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def _post_form(url, fields):
    """multipart/form-data POST - Enrichr rejects urlencoded bodies on addList."""
    boundary = '----dkdboundary1724'
    parts = []
    for k, v in fields.items():
        parts.append('--%s\r\nContent-Disposition: form-data; name="%s"\r\n\r\n%s\r\n' % (boundary, k, v))
    body = (''.join(parts) + '--%s--\r\n' % boundary).encode('utf-8')
    req = urllib.request.Request(url, data=body, headers={
        'Content-Type': 'multipart/form-data; boundary=%s' % boundary})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.loads(r.read().decode('utf-8', 'replace'))


def _post_urlencoded(url, fields):
    """speedrichr's /backgroundenrich takes form-encoded fields, NOT JSON - posting JSON
    there returns HTTP 500 with no message."""
    import urllib.parse
    body = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(url, data=body, headers={
        'Content-Type': 'application/x-www-form-urlencoded'})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.loads(r.read().decode('utf-8', 'replace'))


def enrich(signature, background, libraries=LIBRARIES):
    """Returns {library: DataFrame}. Falls back to the genome-wide background if the
    custom-background service is unavailable, and says so."""
    out = {}
    try:
        ul = _post_form(SPEED + '/addList',
                        {'list': '\n'.join(signature), 'description': 'DKD signature'})
        bg = _post_form(SPEED + '/addbackground', {'background': '\n'.join(background)})
        for lib in libraries:
            d = _post_urlencoded(SPEED + '/backgroundenrich',
                                 {'userListId': ul['userListId'],
                                  'backgroundid': bg['backgroundid'],
                                  'backgroundType': lib})
            out[lib] = _to_frame(d.get(lib, []), custom_bg=True)
            time.sleep(0.5)
        return out, True
    except Exception as e:
        log('  ! custom-background enrichment failed (%s); falling back to genome-wide' % str(e)[:120])

    ul = _post_form(CLASSIC + '/addList',
                    {'list': '\n'.join(signature), 'description': 'DKD signature'})
    for lib in libraries:
        u = '%s/enrich?userListId=%s&backgroundType=%s' % (CLASSIC, ul['userListId'], lib)
        with urllib.request.urlopen(u, timeout=180) as r:
            d = json.loads(r.read().decode('utf-8', 'replace'))
        out[lib] = _to_frame(d.get(lib, []), custom_bg=False)
        time.sleep(0.5)
    return out, False


def _to_frame(rows, custom_bg):
    recs = []
    for r in rows:
        # rank, term, pval, zscore, combined, overlapping genes, adj p, old p, old adj p
        recs.append(dict(term=r[1], p_value=r[2], z=r[3], combined_score=r[4],
                         overlap=len(r[5]), genes=';'.join(r[5]), adj_p=r[6]))
    df = pd.DataFrame(recs)
    if not df.empty:
        df = df.sort_values('adj_p')
    df.attrs['custom_background'] = custom_bg
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--signature', default='results/validation/consensus_signature.tsv')
    ap.add_argument('--topk', type=int, default=50)
    ap.add_argument('--background', default='data/processed/gene_space_all6.tsv')
    ap.add_argument('--out', default='results/enrichment')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    sig = pd.read_csv(args.signature, sep='\t', dtype={'entrez_id': str})
    sig = sig.head(args.topk)
    genes = [s for s in sig['symbol'].dropna().astype(str) if s and s != '?']
    bg = pd.read_csv(args.background, sep='\t', dtype=str)
    background = [s for s in bg['symbol'].dropna().astype(str) if s]
    log('signature: %d symbols | background: %d symbols' % (len(genes), len(background)))
    log('  %s' % ', '.join(genes[:25]))

    res, custom = enrich(genes, background)
    log('\nbackground used: %s' % ('custom 9,900-gene space' if custom else 'GENOME-WIDE (inflated)'))
    for lib, df in res.items():
        path = os.path.join(args.out, '%s.tsv' % lib)
        df.to_csv(path, sep='\t', index=False)
        sig_rows = df[df['adj_p'] < 0.05] if not df.empty else df
        log('\n=== %s: %d terms at adj p < 0.05 ===' % (lib, len(sig_rows)))
        for _, r in sig_rows.head(10).iterrows():
            log('   adj p=%.2e  n=%-3d %s' % (r['adj_p'], r['overlap'], r['term'][:78]))
            log('        %s' % r['genes'][:110])
    log('\nwrote %s/' % args.out)


if __name__ == '__main__':
    main()
