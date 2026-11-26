import json
import os
import csv
import urllib.request
import zipfile
import io
import gzip
import sys
import re
from playwright.sync_api import sync_playwright

sys.stdout.reconfigure(encoding='utf-8')

# Load organism metadata
with open('organisms.json', 'r', encoding='utf-8') as f:
    org_data = json.load(f)

org_list = []
for cat, grps in org_data.items():
    if isinstance(grps, dict):
        for grp, items in grps.items():
            if isinstance(items, list):
                for item in items:
                    if item.get('bestPrimers', 0) > 0:
                        full_key = f"{cat}___{grp}___{item['orgKey']}"
                        org_list.append({
                            'cat': cat,
                            'grp': grp,
                            'orgName': item['orgName'],
                            'orgVersion': item['orgVersion'],
                            'orgKey': item['orgKey'],
                            'fullKey': full_key
                        })

print(f"[+] Total organisms with best primers: {len(org_list)}")

def run():
    all_rows = []
    
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        print("[*] Connecting to qprimerdb downloads page...")
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        batch_size = 50
        all_keys = [o['fullKey'] for o in org_list]
        
        script = """
        async (keys) => {
            let resList = [];
            for (let k of keys) {
                try {
                    let r = await fetch('https://qprimerdb.biodb.org/bioapi/', {
                        method: 'POST',
                        headers: {'Content-Type': 'application/x-www-form-urlencoded'},
                        body: 'database=qprimerdb2&method=DownloadFile&orgKey=' + encodeURIComponent(k) + '&fType=best'
                    });
                    let rawText = await r.text();
                    let url = rawText.trim().replace(/^"|"$/g, '').split('\\\\').join('');
                    resList.push({fullKey: k, url: url});
                } catch(e) {
                    resList.push({fullKey: k, url: null});
                }
            }
            return resList;
        }
        """

        url_results = {}
        for i in range(0, len(all_keys), batch_size):
            chunk = all_keys[i:i+batch_size]
            print(f"[*] Fetching download links {i+1}..{i+len(chunk)} / {len(all_keys)}...")
            res = page.evaluate(script, chunk)
            for item in res:
                if item['url'] and item['url'].startswith('http'):
                    url_results[item['fullKey']] = item['url']

        browser.close()

    print(f"[+] Successfully fetched {len(url_results)} download URLs!")

    print("[*] Downloading primer zip files and parsing sequences...")
    parsed_count = 0
    
    for idx, org in enumerate(org_list):
        fk = org['fullKey']
        url = url_results.get(fk)
        if not url:
            continue

        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            data = urllib.request.urlopen(req, timeout=10).read()
            
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                for name in z.namelist():
                    with gzip.GzipFile(fileobj=io.BytesIO(z.read(name))) as gz:
                        tsv_content = gz.read().decode('utf-8', errors='ignore')
                        lines = tsv_content.splitlines()
                        for line in lines:
                            if not line or line.startswith("#") or line.startswith("Gene_ID"):
                                continue
                            parts = line.split('\t')
                            if len(parts) >= 11:
                                gene_id = parts[0]
                                pair_id = parts[1]
                                fwd = parts[9].upper().strip()
                                rev = parts[10].upper().strip()
                                if fwd and rev and re.match(r'^[ACGTRYSWKMBDHVN]+$', fwd) and re.match(r'^[ACGTRYSWKMBDHVN]+$', rev):
                                    all_rows.append({
                                        'Organism': org['orgName'],
                                        'PairID': f"{org['orgName']}_{gene_id}_{pair_id}",
                                        'Forward': fwd,
                                        'Reverse': rev
                                    })
            parsed_count += 1
            if (idx + 1) % 100 == 0:
                print(f"    - Parsed {idx + 1}/{len(org_list)} organism datasets ({len(all_rows)} primer pairs)...")
        except Exception as e:
            pass

    print(f"[+] Total primer pairs parsed: {len(all_rows)}")

    # Write out comprehensive best_primers_db.csv
    csv_file = r"C:\Users\internet\Downloads\best_primers_db.csv"
    with open(csv_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"[+] Written full database to {csv_file}")

if __name__ == "__main__":
    run()
