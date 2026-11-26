from playwright.sync_api import sync_playwright
import zipfile
import io
import gzip
import sys

sys.stdout.reconfigure(encoding='utf-8')

def test_download():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        print("[*] Opening qprimerdb downloads...")
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        # Trigger download via expectation
        print("[*] Triggering download via button click...")
        with page.expect_download(timeout=15000) as download_info:
            page.evaluate("""
            () => {
                let btns = document.querySelectorAll("tr.ant-table-row button");
                if (btns.length >= 2) btns[1].click(); // Click best primers
            }
            """)
        
        download = download_info.value
        path = download.save_as("test_downloaded.zip")
        print(f"[+] Successfully downloaded file: test_downloaded.zip (suggested filename: {download.suggested_filename})")

        with zipfile.ZipFile("test_downloaded.zip") as z:
            print("Zip files inside:", z.namelist())
            for name in z.namelist():
                with gzip.GzipFile(fileobj=io.BytesIO(z.read(name))) as gz:
                    txt = gz.read().decode('utf-8', errors='ignore')
                    lines = txt.splitlines()
                    print("Header:", lines[0])
                    print("Row 1:", lines[1] if len(lines) > 1 else 'EMPTY')

        browser.close()

if __name__ == "__main__":
    test_download()
