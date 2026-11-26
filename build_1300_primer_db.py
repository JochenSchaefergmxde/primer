import json
import urllib.request
import zipfile
import io
import csv
import re
import sys
import concurrent.futures

sys.stdout.reconfigure(encoding='utf-8')

# Load all 1315 organisms
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
                        full_key = f"{cat}___{grp}___{item['orgKey']}"
                        url = f"https://qprimerdb.biodb.org/tmp/qprimerdb/downloads/{clean_key}.best.primers.xls.zip"
                        org_list.append({
                            'name': item['orgName'],
                            'key': item['orgKey'],
                            'full_key': full_key,
                            'url': url
                        })

print(f"[+] Total organisms to process: {len(org_list)}")

def fetch_top_primer(org):
    # Call DownloadFile API to trigger file generation if needed
    api_url = "https://qprimerdb.biodb.org/bioapi/"
    data = urllib.parse.urlencode({
        'database': 'qprimerdb2',
        'method': 'DownloadFile',
        'orgKey': org['full_key'],
        'fType': 'best'
    }).encode()
    
    req_api = urllib.request.Request(api_url, data=data, headers={
        'Content-Type': 'application/x-www-form-urlencoded',
        'User-Agent': 'Mozilla/5.0'
    })
    
    try:
        res_api = urllib.request.urlopen(req_api, timeout=10).read().decode('utf-8', errors='ignore')
        gen_url = res_api.strip().replace('"', '').replace('\\', '')
        if not gen_url.startswith('http'):
            gen_url = org['url']
    except Exception:
        gen_url = org['url']

    # Download zip file
    req_zip = urllib.request.Request(gen_url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        zip_bytes = urllib.request.urlopen(req_zip, timeout=15).read()
        with zipfile.ZipFile(io.BytesIO(zip_bytes)) as z:
            fname = z.namelist()[0]
            content = z.read(fname).decode('utf-8', errors='ignore')
            lines = content.splitlines()
            if len(lines) > 1:
                parts = lines[1].split('\t')
                if len(parts) >= 10:
                    gene_id = parts[0]
                    pair_name = parts[1]
                    fwd = parts[3].upper().strip()
                    rev = parts[8].upper().strip()
                    if fwd and rev and re.match(r'^[ACGTRYSWKMBDHVN]+$', fwd) and re.match(r'^[ACGTRYSWKMBDHVN]+$', rev):
                        return {
                            'Organism': org['name'],
                            'PairID': f"{org['name']}_{pair_name}",
                            'Forward': fwd,
                            'Reverse': rev
                        }
    except Exception as e:
        pass
    return None

def build_db():
    results = []
    print("[*] Fetching top best primer pair for all 1,315 organisms using parallel threads...")
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
        future_to_org = {executor.submit(fetch_top_primer, org): org for org in org_list}
        completed = 0
        for future in concurrent.futures.as_completed(future_to_org):
            res = future.result()
            if res:
                results.append(res)
            completed += 1
            if completed % 100 == 0 or completed == len(org_list):
                print(f"    - Processed {completed}/{len(org_list)} organisms ({len(results)} valid primer pairs found)...")

    print(f"\n[+] Successfully collected {len(results)} best primer pairs across organisms!")

    # Write out best_primers_db.csv
    csv_path = r"C:\Users\internet\Downloads\best_primers_db.csv"
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
        writer.writeheader()
        writer.writerows(results)

    print(f"[+] Saved complete ~1,300 Organisms Primer Database to: {csv_path}")

if __name__ == "__main__":
    build_db()
