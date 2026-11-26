from playwright.sync_api import sync_playwright
import zipfile
import io
import sys

sys.stdout.reconfigure(encoding='utf-8')

with sync_playwright() as p:
    browser = p.chromium.launch(headless=True)
    page = browser.new_page()
    page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
    page.wait_for_selector("table", timeout=15000)

    key = "animal___invertebrate___Acyrthosiphon_pisum___GCF_005508785.1_pea_aphid_22Mar2018_4r6ur"
    
    # 1. Trigger DownloadFile API
    raw_url = page.evaluate("""
    async (k) => {
        let r = await fetch('https://qprimerdb.biodb.org/bioapi/', {
            method: 'POST',
            headers: {'Content-Type': 'application/x-www-form-urlencoded'},
            body: 'database=qprimerdb2&method=DownloadFile&orgKey=' + encodeURIComponent(k) + '&fType=best'
        });
        let t = await r.text();
        return t.trim().replace(/^"|"$/g, '').split('\\\\').join('');
    }
    """, key)

    print("API returned:", raw_url)

    # 2. Fix extension to .xls.zip if needed
    xls_url = raw_url.replace(".gz.zip", ".xls.zip")
    print("Testing XLS Zip URL:", xls_url)

    res = page.request.get(xls_url)
    print("Status:", res.status)
    print("Bytes:", len(res.body()))

    with zipfile.ZipFile(io.BytesIO(res.body())) as z:
        print("Zip entries:", z.namelist())
        txt = z.read(z.namelist()[0]).decode('utf-8', errors='ignore')
        lines = txt.splitlines()
        print("Total rows:", len(lines))
        print("Row 1:", lines[1].split('\t'))

    browser.close()
