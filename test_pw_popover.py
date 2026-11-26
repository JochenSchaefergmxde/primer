from playwright.sync_api import sync_playwright
import zipfile
import io
import gzip
import sys

sys.stdout.reconfigure(encoding='utf-8')

def test_popover():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        print("[*] Navigating to qprimerdb downloads...")
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        print("[*] Clicking Best Primers button...")
        page.evaluate("""
        () => {
            let btns = document.querySelectorAll("tr.ant-table-row button");
            if (btns.length >= 2) btns[1].click();
        }
        """)

        page.wait_for_timeout(3000)

        link_href = page.evaluate("""
        () => {
            let a = document.querySelector(".ant-popover a");
            return a ? a.href : null;
        }
        """)

        print(f"[+] Popover Link Href: '{link_href}'")

        if link_href:
            # Download file using page.request.get
            print("[*] Fetching file via page.request.get...")
            res = page.request.get(link_href)
            print(f"[+] Status: {res.status}, Content Length: {len(res.body())}")
            
            with zipfile.ZipFile(io.BytesIO(res.body())) as z:
                print("Zip files inside:", z.namelist())
                for fname in z.namelist():
                    with gzip.GzipFile(fileobj=io.BytesIO(z.read(fname))) as gz:
                        txt = gz.read().decode('utf-8', errors='ignore')
                        lines = txt.splitlines()
                        print("Header:", lines[0])
                        print("Row 1:", lines[1] if len(lines) > 1 else 'EMPTY')

        browser.close()

if __name__ == "__main__":
    test_popover()
