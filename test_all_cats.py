from playwright.sync_api import sync_playwright
import json
import zipfile
import io
import re
import sys

sys.stdout.reconfigure(encoding='utf-8')

with open('organisms.json', 'r', encoding='utf-8') as f:
    org_data = json.load(f)

test_orgs = []
for cat, grps in org_data.items():
    if isinstance(grps, dict):
        for grp, items in grps.items():
            if isinstance(items, list):
                for item in items[:2]: # Pick 2 from each group
                    if item.get('bestPrimers', 0) > 0:
                        full_key = f"{cat}___{grp}___{item['orgKey']}"
                        test_orgs.append({
                            'cat': cat,
                            'grp': grp,
                            'name': item['orgName'],
                            'full_key': full_key
                        })

print(f"[+] Total test organisms across categories: {len(test_orgs)}")

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        results = []
        for org in test_orgs:
            fk = org['full_key']
            print(f"[*] Calling DownloadFile for: {fk}...")
            
            raw_url = page.evaluate("""
            async (key) => {
                let r = await fetch('https://qprimerdb.biodb.org/bioapi/', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/x-www-form-urlencoded'},
                    body: 'database=qprimerdb2&method=DownloadFile&orgKey=' + encodeURIComponent(key) + '&fType=best'
                });
                let t = await r.text();
                return t.trim().replace(/^"|"$/g, '').split('\\\\').join('');
            }
            """, fk)

            print(f"    [+] API returned: '{raw_url}'")
            xls_url = raw_url.replace(".gz.zip", ".xls.zip")

            try:
                res = page.request.get(xls_url, timeout=30000)
                print(f"    [+] GET Status: {res.status}, Bytes: {len(res.body())}")
                if res.status == 200 and len(res.body()) > 1000:
                    with zipfile.ZipFile(io.BytesIO(res.body())) as z:
                        fname = z.namelist()[0]
                        txt = z.read(fname).decode('utf-8', errors='ignore')
                        lines = txt.splitlines()
                        if len(lines) > 1:
                            parts = lines[1].split('\t')
                            fwd = parts[3].upper().strip()
                            rev = parts[8].upper().strip()
                            print(f"    [+] SUCCESS! FWD: {fwd[:10]}... REV: {rev[:10]}... (Total rows: {len(lines)})")
                            results.append(org['name'])
            except Exception as e:
                print(f"    [!] Error: {e}")

        browser.close()

    print(f"\n[+] Successfully parsed {len(results)}/{len(test_orgs)} test organisms!")

if __name__ == "__main__":
    run()
