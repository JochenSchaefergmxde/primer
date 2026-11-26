from playwright.sync_api import sync_playwright
import json

with open('organisms.json', 'r', encoding='utf-8') as f:
    org_data = json.load(f)

org_list = []
for cat, grps in org_data.items():
    if isinstance(grps, dict):
        for grp, items in grps.items():
            if isinstance(items, list):
                for item in items:
                    if item.get('bestPrimers', 0) > 0:
                        full_key = f"{cat}___{grp}___{item['orgKey']}"
                        org_list.append({
                            'cat': cat,
                            'grp': grp,
                            'orgName': item['orgName'],
                            'orgVersion': item['orgVersion'],
                            'orgKey': item['orgKey'],
                            'fullKey': full_key
                        })

print(f"[+] Found {len(org_list)} organisms with Best Primers.")

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        print("[*] Opening qprimerdb...")
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        # Batch request 20 URLs via fetch inside page
        test_keys = [o['fullKey'] for o in org_list[:20]]
        
        script = """
        async (keys) => {
            let resList = [];
            for (let k of keys) {
                let r = await fetch('https://qprimerdb.biodb.org/bioapi/', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/x-www-form-urlencoded'},
                    body: 'database=qprimerdb2&method=DownloadFile&orgKey=' + encodeURIComponent(k) + '&fType=best'
                });
                let rawText = await r.text();
                let url = rawText.trim().replace(/^"|"$/g, '').split('\\\\').join('');
                resList.push({fullKey: k, url: url});
            }
            return resList;
        }
        """

        urls = page.evaluate(script, test_keys)
        print("[+] Returned URLs:")
        for u in urls[:5]:
            print(u)

        browser.close()

if __name__ == "__main__":
    run()
