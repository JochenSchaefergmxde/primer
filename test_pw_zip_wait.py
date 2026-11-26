from playwright.sync_api import sync_playwright
import zipfile
import io
import gzip
import urllib.request
import time
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

        print("[*] Waiting 2 seconds for server file creation...")
        time.sleep(2)

        # Download via page evaluate fetch arraybuffer
        b64 = page.evaluate("""
        async (fileUrl) => {
            let res = await fetch(fileUrl);
            let buf = await res.arrayBuffer();
            let bytes = new Uint8Array(buf);
            let binary = '';
            for (let i = 0; i < bytes.byteLength; i++) {
                binary += String.fromCharCode(bytes[i]);
            }
            return btoa(binary);
        }
        """, url)

        import base64
        data = base64.b64decode(b64)
        print(f"[+] Downloaded bytes via browser fetch: {len(data)}")

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
