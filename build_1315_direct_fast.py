import json
import urllib.request
import zipfile
import io
import csv
import re
import sys
import concurrent.futures

sys.stdout.reconfigure(encoding='utf-8')

# Load all 1315 organisms metadata
with open('organisms.json', 'r', encoding='utf-8') as f:
    org_data = json.load(f)

org_list = []
for cat, grps in org_data.items():
    if isinstance(grps, dict):
        for grp, items in grps.items():
            if isinstance(items, list):
                for item in items:
                    if item.get('bestPrimers', 0) > 0:
                        # Construct direct URL using orgKey
                        clean_key = item['orgKey'].replace('___', '__')
                        url = f"https://qprimerdb.biodb.org/tmp/qprimerdb/downloads/{clean_key}.best.primers.xls.zip"
                        org_list.append({
                            'name': item['orgName'],
                            'key': item['orgKey'],
                            'url': url
                        })

print(f"[+] Total organisms with best primers: {len(org_list)}")

def fetch_organism_best_primer(org):
    req = urllib.request.Request(org['url'], headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
    try:
        data = urllib.request.urlopen(req, timeout=12).read()
        if len(data) > 1000:
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                fname = z.namelist()[0]
                content = z.read(fname).decode('utf-8', errors='ignore')
                lines = content.splitlines()
                if len(lines) > 1:
                    parts = lines[1].split('\t')
                    if len(parts) >= 9:
                        fwd = parts[3].upper().strip()
                        rev = parts[8].upper().strip()
                        if fwd and rev and re.match(r'^[ACGTRYSWKMBDHVN]+$', fwd) and re.match(r'^[ACGTRYSWKMBDHVN]+$', rev):
                            return {
                                'Organism': org['name'],
                                'PairID': f"{org['name']}_{parts[1]}",
                                'Forward': fwd,
                                'Reverse': rev
                            }
    except Exception:
        pass
    return None

def main():
    print("[*] Fetching top best primer pair for all 1,315 organisms via 30 parallel threads...")
    results = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=30) as executor:
        future_to_org = {executor.submit(fetch_organism_best_primer, org): org for org in org_list}
        completed = 0
        for future in concurrent.futures.as_completed(future_to_org):
            res = future.result()
            if res:
                results.append(res)
            completed += 1
            if completed % 100 == 0 or completed == len(org_list):
                print(f"    - Processed {completed}/{len(org_list)} organisms ({len(results)} valid primer pairs collected)...")

    print(f"\n[+] Total Valid Best Primer Pairs Collected: {len(results)}")
    
    csv_path = r"C:\Users\internet\Downloads\best_primers_db.csv"
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
        writer.writeheader()
        writer.writerows(results)

    print(f"[+] Saved complete ~1,300 Organisms Primer Database to: {csv_path}")

if __name__ == "__main__":
    main()
