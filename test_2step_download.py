from playwright.sync_api import sync_playwright
import zipfile
import io
import gzip
import sys

sys.stdout.reconfigure(encoding='utf-8')

def run_test():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        test_keys = [
            "animal___invertebrate___Acyrthosiphon_pisum___GCF_005508785.1_pea_aphid_22Mar2018_4r6ur",
            "microorganism___bacteria___Escherichia_coli___GCF_000005845.2_ASM584v2",
            "microorganism___bacteria___Staphylococcus_aureus___GCF_000013425.1_ASM1342v1"
        ]

        for key in test_keys:
            print(f"[*] Calling DownloadFile API for key: {key}...")
            raw_url = page.evaluate("""
            async (k) => {
                let r = await fetch('https://qprimerdb.biodb.org/bioapi/', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/x-www-form-urlencoded'},
                    body: 'database=qprimerdb2&method=DownloadFile&orgKey=' + encodeURIComponent(k) + '&fType=best'
                });
                let text = await r.text();
                return text.trim().replace(/^"|"$/g, '').split('\\\\').join('');
            }
            """, key)

            print(f"    [+] Got generated file URL: '{raw_url}'")
            print(f"    [*] Fetching file content...")
            
            res = page.request.get(raw_url, timeout=60000)
            data = res.body()
            print(f"    [+] Status: {res.status}, Downloaded {len(data)} bytes!")

            with zipfile.ZipFile(io.BytesIO(data)) as z:
                print("    [+] Zip contents:", z.namelist())
                for name in z.namelist():
                    with gzip.GzipFile(fileobj=io.BytesIO(z.read(name))) as gz:
                        lines = gz.read().decode('utf-8', errors='ignore').splitlines()
                        print(f"    [+] Total rows: {len(lines)}")
                        print(f"    [+] Sample row: {lines[1][:100] if len(lines)>1 else 'EMPTY'}")

        browser.close()

if __name__ == "__main__":
    run_test()
