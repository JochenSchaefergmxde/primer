from playwright.sync_api import sync_playwright
import zipfile
import io
import gzip
import urllib.request
import sys

sys.stdout.reconfigure(encoding='utf-8')

def test_fetch():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        url = page.evaluate("""
        async () => {
            let key = "animal___invertebrate___Acyrthosiphon_pisum___GCF_005508785.1_pea_aphid_22Mar2018_4r6ur";
            let r = await fetch('https://qprimerdb.biodb.org/bioapi/', {
                method: 'POST',
                headers: {'Content-Type': 'application/x-www-form-urlencoded'},
                body: 'database=qprimerdb2&method=DownloadFile&orgKey=' + encodeURIComponent(key) + '&fType=best'
            });
            let rawText = await r.text();
            return rawText.trim().replace(/^"|"$/g, '').split('\\\\').join('');
        }
        """)

        print(f"[+] Returned URL from browser API: '{url}'")

        cookies = page.context.cookies()
        cookie_str = "; ".join([f"{c['name']}={c['value']}" for c in cookies])

        req = urllib.request.Request(url, headers={
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)',
            'Referer': 'https://qprimerdb.biodb.org/downloads',
            'Cookie': cookie_str
        })
        
        data = urllib.request.urlopen(req).read()
        print(f"[+] Downloaded bytes: {len(data)}")

        with zipfile.ZipFile(io.BytesIO(data)) as z:
            print("Zip files inside:", z.namelist())
            for fname in z.namelist():
                with gzip.GzipFile(fileobj=io.BytesIO(z.read(fname))) as gz:
                    txt = gz.read().decode('utf-8', errors='ignore')
                    lines = txt.splitlines()
                    print("Header:", lines[0])
                    print("Row 1:", lines[1] if len(lines) > 1 else 'EMPTY')

        browser.close()

if __name__ == "__main__":
    test_fetch()
