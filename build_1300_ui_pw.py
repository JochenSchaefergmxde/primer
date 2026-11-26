from playwright.sync_api import sync_playwright
import zipfile
import io
import csv
import re
import sys

sys.stdout.reconfigure(encoding='utf-8')

def build_1315_database():
    print("[*] Launching browser to interact with qPrimerDB UI table...")
    results = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("tr.ant-table-row", timeout=15000)

        # 1. Expand pagination to maximum size (100 per page)
        print("[*] Setting page size to 100 items per page...")
        page.evaluate("""
        () => {
            let select = document.querySelector(".ant-pagination-options-size-changer");
            if (select) select.click();
        }
        """)
        page.wait_for_timeout(1000)
        page.evaluate("""
        () => {
            let options = document.querySelectorAll(".ant-select-item-option");
            for (let opt of options) {
                if (opt.innerText.includes("100")) {
                    opt.click();
                    break;
                }
            }
        }
        """)
        page.wait_for_timeout(3000)

        page_num = 1
        while True:
            print(f"[*] Processing UI Page {page_num}...")
            rows_data = page.evaluate("""
            async () => {
                let rows = document.querySelectorAll("tr.ant-table-row");
                let items = [];
                for (let row of rows) {
                    let nameCell = row.children[0];
                    let orgName = nameCell ? nameCell.innerText.trim() : "";
                    let btns = row.querySelectorAll("button");
                    if (btns.length >= 2) {
                        btns[0].click(); // Best Primers button
                        await new Promise(r => setTimeout(r, 600));
                        let popA = document.querySelector(".ant-popover a");
                        let url = popA ? popA.href : null;
                        if (url) {
                            items.append ? items.append({org: orgName, url: url}) : items.push({org: orgName, url: url});
                        }
                    }
                }
                return items;
            }
            """)

            print(f"    [+] Collected {len(rows_data)} download links from page {page_num}.")

            # Fetch files for this page
            for item in rows_data:
                org_name = item['org']
                url = item['url']
                try:
                    res = page.request.get(url, timeout=30000)
                    if res.status == 200:
                        with zipfile.ZipFile(io.BytesIO(res.body())) as z:
                            fname = z.namelist()[0]
                            content = z.read(fname).decode('utf-8', errors='ignore')
                            lines = content.splitlines()
                            if len(lines) > 1:
                                parts = lines[1].split('\t')
                                if len(parts) >= 10:
                                    pair_id = parts[1]
                                    fwd = parts[3].upper().strip()
                                    rev = parts[8].upper().strip()
                                    if fwd and rev and re.match(r'^[ACGTRYSWKMBDHVN]+$', fwd) and re.match(r'^[ACGTRYSWKMBDHVN]+$', rev):
                                        results.append({
                                            'Organism': org_name,
                                            'PairID': f"{org_name}_{pair_id}",
                                            'Forward': fwd,
                                            'Reverse': rev
                                        })
                except Exception:
                    pass

            # Next page
            next_btn = page.query_selector(".ant-pagination-next:not(.ant-pagination-disabled)")
            if next_btn:
                print("[*] Clicking Next Page...")
                next_btn.click()
                page.wait_for_timeout(3000)
                page_num += 1
            else:
                break

        browser.close()

    print(f"\n[+] Total parsed Best Primers: {len(results)}")
    csv_path = r"C:\Users\internet\Downloads\best_primers_db.csv"
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
        writer.writeheader()
        writer.writerows(results)

    print(f"[+] Saved complete ~1,300 Organisms Primer Database to: {csv_path}")

if __name__ == "__main__":
    build_1315_database()
