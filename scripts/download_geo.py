#!/usr/bin/env python
"""DKD public dataset downloader - GEO series matrix + selected supplementary files."""
import urllib.request, os, sys, re, time, hashlib

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "raw", "geo")

def geo_dir(gse, kind):
    stub = gse[:-3] + "nnn"
    return "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/%s/" % (stub, gse, kind)

def listdir(url):
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            html = r.read().decode('utf-8', 'replace')
    except Exception as e:
        print("   ! list failed: %s" % e); return []
    return [f for f in re.findall(r'href="([^"?/][^"]*)"', html) if not f.startswith('http')]

def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print("   = skip (exists) %s [%.1f MB]" % (os.path.basename(dest), os.path.getsize(dest)/1e6)); return True
    tmp = dest + ".part"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=300) as r, open(tmp, 'wb') as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk: break
                    f.write(chunk)
            os.replace(tmp, dest)
            print("   + %s [%.1f MB]" % (os.path.basename(dest), os.path.getsize(dest)/1e6)); return True
        except Exception as e:
            print("   ! retry %d: %s" % (attempt+1, e)); time.sleep(3)
    return False

def download(gse, kinds=("matrix",), suppl_filter=None):
    print("== %s ==" % gse)
    out = os.path.join(BASE, gse); os.makedirs(out, exist_ok=True)
    for kind in kinds:
        url = geo_dir(gse, kind)
        files = listdir(url)
        if kind == "suppl" and suppl_filter:
            files = [f for f in files if suppl_filter(f)]
        if not files: print("   (no files in %s)" % kind)
        for f in files:
            fetch(url + f, os.path.join(out, f))

if __name__ == "__main__":
    for gse in sys.argv[1:]:
        download(gse)
