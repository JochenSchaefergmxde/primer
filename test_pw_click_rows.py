from playwright.sync_api import sync_playwright
import zipfile
import io
import gzip
import sys

sys.stdout.reconfigure(encoding='utf-8')

def test_click_rows():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        print("[*] Navigating to qprimerdb downloads...")
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("tr.ant-table-row", timeout=15000)

        # Get first 3 rows
        print("[*] Interacting with UI table rows...")
        
        script = """
        async () => {
            let rows = document.querySelectorAll("tr.ant-table-row");
            let results = [];
            for (let i = 0; i < Math.min(3, rows.length); i++) {
                let row = rows[i];
                let orgName = row.children[0] ? row.children[0].innerText : "Unknown";
                let btns = row.querySelectorAll("button");
                if (btns.length >= 2) {
                    btns[1].click(); // Click best primers
                    await new Promise(r => setTimeout(r, 2000));
                    let popoverA = document.querySelector(".ant-popover a");
                    let url = popoverA ? popoverA.href : null;
                    results.push({org: orgName, url: url});
                }
            }
            return results;
        }
        """

        res = page.evaluate(script)
        print("[+] Click results:", res)

        for item in res:
            if item['url']:
                print(f"[*] Fetching: {item['org']} ({item['url']})...")
                r = page.request.get(item['url'], timeout=60000)
                body = r.body()
                print(f"    [+] Status: {r.status}, Bytes: {len(body)}")
                if len(body) > 1000:
                    with zipfile.ZipFile(io.BytesIO(body)) as z:
                        print("    [+] Zip contents:", z.namelist())

        browser.close()

if __name__ == "__main__":
    test_click_rows()
