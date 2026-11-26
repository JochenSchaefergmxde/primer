from playwright.sync_api import sync_playwright
import json
import os
import csv
import urllib.request
import zipfile
import io
import gzip
import time

def fetch_primers():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        print("[*] Navigating to qprimerdb downloads page to establish session...")
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        # 1. Get complete organism list with orgKey and categories
        print("[*] Fetching complete org list...")
        org_list_json = page.evaluate("""
            async () => {
                let res = await fetch('https://qprimerdb.biodb.org/bioapi/', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/x-www-form-urlencoded'},
                    body: 'database=qprimerdb2&method=GetOrgList'
                });
                return await res.json();
            }
        """)

        organisms = []
        for cat, groups in org_list_json.items():
            if isinstance(groups, dict):
                for grp, items in groups.items():
                    if isinstance(items, list):
                        for item in items:
                            if item.get('bestPrimers', 0) > 0 or item.get('allPrimers', 0) > 0:
                                full_key = f"{cat}___{grp}___{item['orgKey']}"
                                organisms.append({
                                    'name': item['orgName'],
                                    'version': item['orgVersion'],
                                    'orgKey': item['orgKey'],
                                    'fullKey': full_key
                                })

        print(f"[+] Found {len(organisms)} organisms with primers.")

        # Batch fetch download URLs via page.evaluate
        all_primers_data = []
        print(f"[*] Fetching primer download links for {len(organisms)} organisms...")
        
        # Test first 10 for rapid verification
        test_orgs = organisms[:10]
        
        script = """
        async (orgKeys) => {
            let results = [];
            for (let key of orgKeys) {
                try {
                    let res = await fetch('https://qprimerdb.biodb.org/bioapi/', {
                        method: 'POST',
                        headers: {'Content-Type': 'application/x-www-form-urlencoded'},
                        body: 'database=qprimerdb2&method=DownloadFile&orgKey=' + encodeURIComponent(key) + '&fType=best'
                    });
                    let url = await res.json();
                    results.push({key: key, url: url});
                } catch(e) {
                    results.push({key: key, url: null});
                }
            }
            return results;
        }
        """

        urls = page.evaluate(script, [o['fullKey'] for o in test_orgs])
        print("[+] Download URLs sample:", urls[:3])

        browser.close()

if __name__ == "__main__":
    fetch_primers()
