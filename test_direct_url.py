import json
import urllib.request
import zipfile
import io
import gzip
import sys

sys.stdout.reconfigure(encoding='utf-8')

with open('organisms.json', 'r', encoding='utf-8') as f:
    org_data = json.load(f)

org_list = []
for cat, grps in org_data.items():
    if isinstance(grps, dict):
        for grp, items in grps.items():
            if isinstance(items, list):
                for item in items:
                    if item.get('bestPrimers', 0) > 0:
                        clean_key = item['orgKey'].replace('___', '__')
                        url = f"https://qprimerdb.biodb.org/tmp/qprimerdb/downloads/{clean_key}.best.primers.gz.zip"
                        org_list.append({
                            'name': item['orgName'],
                            'key': item['orgKey'],
                            'url': url
                        })

print(f"[+] Total best primer URLs constructed: {len(org_list)}")
print("[*] Testing direct download of first 3 best primer zip files...")

for item in org_list[:3]:
    print(f"[*] Fetching: {item['name']} ({item['url']})...")
    req = urllib.request.Request(item['url'], headers={
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)',
        'Referer': 'https://qprimerdb.biodb.org/downloads'
    })
    try:
        data = urllib.request.urlopen(req, timeout=15).read()
        print(f"    [+] Downloaded {len(data)} bytes!")
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            name = z.namelist()[0]
            with gzip.GzipFile(fileobj=io.BytesIO(z.read(name))) as gz:
                lines = gz.read().decode('utf-8', errors='ignore').splitlines()
                print(f"    [+] Rows: {len(lines)}, Header: {lines[0][:80]}")
                print(f"    [+] Sample Pair: {lines[1][:100] if len(lines)>1 else 'N/A'}")
    except Exception as e:
        print(f"    [!] Error: {e}")
