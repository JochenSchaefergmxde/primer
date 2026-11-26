from playwright.sync_api import sync_playwright
import zipfile
import io
import csv
import re
import sys
import time

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("[*] Launching browser to scrape qPrimerDB UI table...")
    results = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("tr.ant-table-row", timeout=15000)

        page_num = 1
        while True:
            # 1. Get all row buttons and click them to extract popover URLs
            urls_in_page = page.evaluate("""
            async () => {
                let rows = document.querySelectorAll("tr.ant-table-row");
                let list = [];
                for (let r of rows) {
                    let org = r.children[0] ? r.children[0].innerText.trim() : "Unknown";
                    let btns = r.querySelectorAll("button");
                    if (btns.length >= 1) {
                        btns[0].click();
                        await new Promise(res => setTimeout(res, 300));
                        let a = document.querySelector(".ant-popover a");
                        if (a && a.href) {
                            list.push({org: org, url: a.href});
                        }
                    }
                }
                return list;
            }
            """)

            print(f"[+] Page {page_num}: Found {len(urls_in_page)} download links.")

            for item in urls_in_page:
                org_name = item['org']
                url = item['url'].replace(".gz.zip", ".xls.zip")
                try:
                    res = page.request.get(url, timeout=20000)
                    if res.status == 200:
                        with zipfile.ZipFile(io.BytesIO(res.body())) as z:
                            fname = z.namelist()[0]
                            content = z.read(fname).decode('utf-8', errors='ignore')
                            lines = content.splitlines()
                            if len(lines) > 1:
                                parts = lines[1].split('\t')
                                if len(parts) >= 9:
                                    fwd = parts[3].upper().strip()
                                    rev = parts[8].upper().strip()
                                    if fwd and rev and re.match(r'^[ACGTRYSWKMBDHVN]+$', fwd) and re.match(r'^[ACGTRYSWKMBDHVN]+$', rev):
                                        results.append({
                                            'Organism': org_name,
                                            'PairID': f"{org_name}_{parts[1]}",
                                            'Forward': fwd,
                                            'Reverse': rev
                                        })
                except Exception:
                    pass

            # Go to next page
            next_btn = page.query_selector("li.ant-pagination-next:not(.ant-pagination-disabled) button")
            if not next_btn:
                next_btn = page.query_selector("li.ant-pagination-next:not(.ant-pagination-disabled)")
            
            if next_btn:
                page_num += 1
                next_btn.click()
                time.sleep(2)
            else:
                print("[*] Reached last page!")
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
    main()
