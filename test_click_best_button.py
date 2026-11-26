from playwright.sync_api import sync_playwright
import urllib.request
import zipfile
import io
import gzip
import sys

sys.stdout.reconfigure(encoding='utf-8')

def test_click_best():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        print("[*] Opening qprimerdb downloads...")
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("tr.ant-table-row", timeout=15000)

        # Inspect all buttons in the first row
        info = page.evaluate("""
        async () => {
            let row = document.querySelector("tr.ant-table-row");
            let btns = row.querySelectorAll("button");
            let btnTexts = Array.from(btns).map(b => b.innerText || b.innerHTML);
            
            // Click second button (best primers)
            btns[0].click();
            await new Promise(r => setTimeout(r, 2000));
            let pop1 = document.querySelector(".ant-popover a") ? document.querySelector(".ant-popover a").href : null;

            if (btns.length > 1) {
                btns[1].click();
                await new Promise(r => setTimeout(r, 2000));
            }
            let pop2 = document.querySelector(".ant-popover a") ? document.querySelector(".ant-popover a").href : null;

            return {
                btnTexts: btnTexts,
                pop1: pop1,
                pop2: pop2
            };
        }
        """)

        print("[+] Button Info:", info)

        target_url = info['pop1'] or info['pop2']
        if target_url:
            print(f"[*] Downloading: {target_url} using urllib...")
            req = urllib.request.Request(target_url, headers={
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)',
                'Referer': 'https://qprimerdb.biodb.org/downloads'
            })
            data = urllib.request.urlopen(req, timeout=30).read()
            print(f"[+] Downloaded {len(data)} bytes!")
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                name = z.namelist()[0]
                with gzip.GzipFile(fileobj=io.BytesIO(z.read(name))) as gz:
                    lines = gz.read().decode('utf-8', errors='ignore').splitlines()
                    print(f"[+] Rows in TSV: {len(lines)}")
                    print(f"[+] Sample row: {lines[1][:100] if len(lines)>1 else 'EMPTY'}")

        browser.close()

if __name__ == "__main__":
    test_click_best()
