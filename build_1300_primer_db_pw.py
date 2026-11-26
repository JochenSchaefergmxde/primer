import json
import zipfile
import io
import csv
import re
import sys
from playwright.sync_api import sync_playwright

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
                        full_key = f"{cat}___{grp}___{item['orgKey']}"
                        org_list.append({
                            'name': item['orgName'],
                            'key': item['orgKey'],
                            'full_key': full_key
                        })

print(f"[+] Total organisms to process: {len(org_list)}")

def build_db():
    results = []
    
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        print("[*] Opening qprimerdb downloads page...")
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        # Batch request generated URLs in JS
        batch_size = 50
        print("[*] Generating download links for all organisms...")
        
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
                    let t = await r.text();
                    let url = t.trim().replace(/^"|"$/g, '').split('\\\\').join('');
                    resList.push({key: k, url: url});
                } catch(e) {
                    resList.push({key: k, url: null});
                }
            }
            return resList;
        }
        """

        url_map = {}
        all_keys = [o['full_key'] for o in org_list]
        for i in range(0, len(all_keys), batch_size):
            chunk = all_keys[i:i+batch_size]
            print(f"    - Fetching links {i+1}..{i+len(chunk)} / {len(all_keys)}...")
            res = page.evaluate(script, chunk)
            for item in res:
                if item['url'] and item['url'].startswith('http'):
                    url_map[item['key']] = item['url']

        print(f"[+] Generated {len(url_map)} file URLs!")

        print("[*] Downloading zip files and parsing top primer pair for each organism...")
        for idx, org in enumerate(org_list):
            fk = org['full_key']
            url = url_map.get(fk)
            if not url:
                continue

            try:
                r = page.request.get(url, timeout=30000)
                if r.status == 200:
                    data = r.body()
                    with zipfile.ZipFile(io.BytesIO(data)) as z:
                        fname = z.namelist()[0]
                        txt = z.read(fname).decode('utf-8', errors='ignore')
                        lines = txt.splitlines()
                        if len(lines) > 1:
                            parts = lines[1].split('\t')
                            if len(parts) >= 9:
                                fwd = parts[3].upper().strip()
                                rev = parts[8].upper().strip()
                                if fwd and rev and re.match(r'^[ACGTRYSWKMBDHVN]+$', fwd) and re.match(r'^[ACGTRYSWKMBDHVN]+$', rev):
                                    results.append({
                                        'Organism': org['name'],
                                        'PairID': f"{org['name']}_{parts[1]}",
                                        'Forward': fwd,
                                        'Reverse': rev
                                    })
            except Exception:
                pass

            if (idx + 1) % 100 == 0 or (idx + 1) == len(org_list):
                print(f"    - Downloaded & parsed {idx + 1}/{len(org_list)} organisms ({len(results)} valid primer pairs)...")

        browser.close()

    print(f"\n[+] Total valid best primer pairs collected: {len(results)}")

    csv_path = r"C:\Users\internet\Downloads\best_primers_db.csv"
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
        writer.writeheader()
        writer.writerows(results)

    print(f"[+] Written complete database to {csv_path}")

if __name__ == "__main__":
    build_db()
