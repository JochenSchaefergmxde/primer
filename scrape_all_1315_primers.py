import json
import zipfile
import io
import csv
import re
import sys
import time
import concurrent.futures
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

print(f"[+] Total organisms with best primers: {len(org_list)}")

def scrape_worker(chunk):
    results = []
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        for org in chunk:
            fk = org['full_key']
            org_name = org['name']
            try:
                # Trigger Best Primers in UI state via page evaluate
                url = page.evaluate("""
                async (key) => {
                    window.g_app._store.dispatch({
                        type: 'fetchData/downloadFile',
                        payload: {
                            params: {
                                database: 'qprimerdb2',
                                method: 'DownloadFile',
                                orgKey: key,
                                fType: 'best'
                            }
                        }
                    });
                    for (let i = 0; i < 20; i++) {
                        await new Promise(r => setTimeout(r, 250));
                        let s = window.g_app._store.getState();
                        if (s.fetchData && s.fetchData.dlLink) {
                            return s.fetchData.dlLink;
                        }
                    }
                    return null;
                }
                """, fk)

                if url:
                    # Fix extension if returned .gz.zip -> .xls.zip
                    xls_url = url.replace(".gz.zip", ".xls.zip")
                    res = page.request.get(xls_url, timeout=15000)
                    if res.status == 200:
                        with zipfile.ZipFile(io.BytesIO(res.body())) as z:
                            fname = z.namelist()[0]
                            content = z.read(fname).decode('utf-8', errors='ignore')
                            lines = content.splitlines()
                            if len(lines) > 1:
                                parts = lines[1].split('\t')
                                if len(parts) >= 9:
                                    pair_id = parts[1]
                                    fwd = parts[3].upper().strip()
                                    rev = parts[8].upper().strip()
                                    if fwd and rev and re.match(r'^[ACGTRYSWKMBDHVN]+$', fwd) and re.match(r'^[ACGTRYSWKMBDHVN]+$', rev):
                                        results.append({
                                            'Organism': org_name,
                                            'PairID': f"{org_name}_{pair_id}",
                                            'Forward': fwd,
                                            'Reverse': rev
                                        })
            except Exception:
                pass

        browser.close()
    return results

def main():
    num_threads = 8
    chunk_size = (len(org_list) + num_threads - 1) // num_threads
    chunks = [org_list[i:i+chunk_size] for i in range(0, len(org_list), chunk_size)]
    
    print(f"[*] Starting {len(chunks)} parallel browser workers to scrape all 1,315 organisms...")
    
    all_results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(scrape_worker, chunk) for chunk in chunks]
        for f in concurrent.futures.as_completed(futures):
            res = f.result()
            all_results.extend(res)
            print(f"    [+] Worker finished batch ({len(res)} primer pairs found, total: {len(all_results)})...")

    print(f"\n[+] Total Best Primer Pairs Scraped: {len(all_results)}")
    
    csv_path = r"C:\Users\internet\Downloads\best_primers_db.csv"
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
        writer.writeheader()
        writer.writerows(all_results)

    print(f"[+] Saved complete ~1,300 Organisms Primer Database to: {csv_path}")

if __name__ == "__main__":
    main()
