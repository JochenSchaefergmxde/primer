from playwright.sync_api import sync_playwright
import time
import json

def run(playwright):
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page()
    print("Navigating to downloads page...")
    page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
    
    # wait a bit for table to render if needed
    page.wait_for_selector("table", timeout=10000)
    
    data = page.evaluate("""
    () => {
        let results = [];
        let rows = document.querySelectorAll("tr");
        rows.forEach(r => {
            let cols = r.querySelectorAll("td");
            if(cols.length > 0) {
                let textData = Array.from(cols).map(c => c.innerText);
                let links = Array.from(cols).map(c => {
                    let a = c.querySelector("a");
                    return a ? a.href : null;
                });
                results.push({text: textData, links: links});
            }
        });
        return results;
    }
    """)
    print("Row data:")
    for d in data[:5]:
        print(d)
        
    print(f"Total rows found: {len(data)}")
    browser.close()

with sync_playwright() as p:
    run(p)
