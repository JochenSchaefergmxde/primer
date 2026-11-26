from playwright.sync_api import sync_playwright
import urllib.request
import zipfile
import io
import gzip
import csv
import sys
import re

sys.stdout.reconfigure(encoding='utf-8')

def run_download():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        print("[*] Opening qprimerdb downloads...")
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("tr.ant-table-row", timeout=15000)

        # Execute DownloadFile API call for best primers
        test_keys = [
            "animal___invertebrate___Acyrthosiphon_pisum___GCF_005508785.1_pea_aphid_22Mar2018_4r6ur",
            "microorganism___bacteria___Escherichia_coli___GCF_000005845.2_ASM584v2",
            "microorganism___bacteria___Staphylococcus_aureus___GCF_000013425.1_ASM1342v1",
            "microorganism___fungi___Candida_albicans___GCF_000182965.3_ASM18296v3"
        ]

        results = []
        for key in test_keys:
            print(f"[*] Triggering Best Primers generation for: {key}...")
            url = page.evaluate("""
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
            
            print(f"    [+] Generated URL: '{url}'")
            results.append((key, url))

        browser.close()

    print("\n[*] Downloading files using urllib...")
    all_rows = []
    
    for key, url in results:
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            data = urllib.request.urlopen(req, timeout=30).read()
            print(f"[+] Downloaded {len(data)} bytes for {key}!")
            
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                for fname in z.namelist():
                    with gzip.GzipFile(fileobj=io.BytesIO(z.read(fname))) as gz:
                        txt = gz.read().decode('utf-8', errors='ignore')
                        lines = txt.splitlines()
                        header = lines[0].split('\t')
                        print(f"    [+] TSV rows: {len(lines)}")
                        
                        org_name = key.split('___')[-1]
                        for line in lines[1:]:
                            parts = line.split('\t')
                            if len(parts) >= 11:
                                fwd = parts[9].upper().strip()
                                rev = parts[10].upper().strip()
                                if fwd and rev and re.match(r'^[ACGTRYSWKMBDHVN]+$', fwd) and re.match(r'^[ACGTRYSWKMBDHVN]+$', rev):
                                    all_rows.append({
                                        'Organism': org_name,
                                        'PairID': f"{org_name}_{parts[0]}_{parts[1]}",
                                        'Forward': fwd,
                                        'Reverse': rev
                                    })
        except Exception as e:
            print(f"[!] Error downloading/parsing {key}: {e}")

    print(f"\n[+] Total parsed Best Primers: {len(all_rows)}")

if __name__ == "__main__":
    run_download()
