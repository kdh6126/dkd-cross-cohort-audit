#!/usr/bin/env python
"""Extract sample-level metadata from GEO series_matrix files -> TSV inventory."""
import gzip, glob, os, re, sys, csv

def parse(path):
    op = gzip.open if path.endswith('.gz') else open
    hdr = {}
    with op(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'): break
            if not line.startswith('!'): continue
            parts = line.rstrip('\n').split('\t')
            key = parts[0][1:]
            vals = [p.strip('"') for p in parts[1:]]
            hdr.setdefault(key, []).append(vals)
    return hdr

def flat(hdr, key):
    return hdr.get(key, [[]])[0]

def rows(path):
    hdr = parse(path)
    gse = (flat(hdr,'Series_geo_accession') or ['NA'])[0]
    plat = ';'.join(sorted(set(flat(hdr,'Sample_platform_id')))) or 'NA'
    gsms = flat(hdr,'Sample_geo_accession')
    titles = flat(hdr,'Sample_title')
    src = flat(hdr,'Sample_source_name_ch1')
    chars = hdr.get('Sample_characteristics_ch1', [])
    out = []
    for i, g in enumerate(gsms):
        ch = ' | '.join(c[i] for c in chars if i < len(c) and c[i])
        out.append({
            'gse': gse, 'platform': plat, 'gsm': g,
            'title': titles[i] if i < len(titles) else '',
            'source': src[i] if i < len(src) else '',
            'characteristics': ch,
            'file': os.path.basename(path),
        })
    return out

if __name__ == '__main__':
    files = sorted(glob.glob('data/raw/geo/*/*series_matrix.txt.gz'))
    allrows = []
    for p in files:
        try:
            r = rows(p); allrows += r
            print("  %-45s %3d samples" % (os.path.basename(p), len(r)), file=sys.stderr)
        except Exception as e:
            print("  ! %s: %s" % (p, e), file=sys.stderr)
    os.makedirs('data/metadata', exist_ok=True)
    with open('data/metadata/sample_inventory.tsv','w',newline='',encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=['gse','platform','gsm','title','source','characteristics','file'], delimiter='\t')
        w.writeheader(); w.writerows(allrows)
    print("\nTOTAL %d samples -> data/metadata/sample_inventory.tsv" % len(allrows), file=sys.stderr)
