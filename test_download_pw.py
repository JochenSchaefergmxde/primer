from playwright.sync_api import sync_playwright
import json

def run(playwright):
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page()
    
    def log_response(response):
        if "bioapi" in response.url and response.request.method == "POST":
            print(">>> POST URL:", response.url)
            print(">>> Post Data:", response.request.post_data)
            print(">>> Response:", response.text()[:500])

    page.on("response", log_response)

    print("Navigating to downloads...")
    page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
    page.wait_for_selector("table", timeout=10000)
    
    print("Clicking first 'Best Primers' button...")
    page.evaluate("""
    () => {
        let btns = document.querySelectorAll("tr.ant-table-row button");
        if(btns.length >= 2) {
            btns[1].click(); // Best primers button
        }
    }
    """)
    page.wait_for_timeout(3000)
    browser.close()

with sync_playwright() as p:
    run(p)
