from playwright.sync_api import sync_playwright
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
                        full_key = f"{cat}___{grp}___{item['orgKey']}"
                        org_list.append({
                            'name': item['orgName'],
                            'full_key': full_key
                        })

print(f"[+] Total organisms: {len(org_list)}")

with sync_playwright() as p:
    browser = p.chromium.launch(headless=True)
    page = browser.new_page()
    page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
    page.wait_for_selector("table", timeout=15000)

    for item in org_list[:3]:
        fk = item['full_key']
        print(f"[*] Calling DownloadFile for: {fk}...")
        url = page.evaluate("""
        async (k) => {
            let r = await fetch('https://qprimerdb.biodb.org/bioapi/', {
                method: 'POST',
                headers: {'Content-Type': 'application/x-www-form-urlencoded'},
                body: 'database=qprimerdb2&method=DownloadFile&orgKey=' + encodeURIComponent(k) + '&fType=best'
            });
            let t = await r.text();
            return t.trim().replace(/^"|"$/g, '').split('\\\\').join('');
        }
        """, fk)

        print(f"    [+] API returned URL: '{url}'")

        if url.startswith('http'):
            try:
                res = page.request.get(url, timeout=30000)
                body = res.body()
                print(f"    [+] Downloaded {len(body)} bytes!")
                with zipfile.ZipFile(io.BytesIO(body)) as z:
                    name = z.namelist()[0]
                    with gzip.GzipFile(fileobj=io.BytesIO(z.read(name))) as gz:
                        lines = gz.read().decode('utf-8', errors='ignore').splitlines()
                        print(f"    [+] Rows: {len(lines)}, Header: {lines[0][:60]}")
            except Exception as e:
                print(f"    [!] Error: {e}")

    browser.close()
