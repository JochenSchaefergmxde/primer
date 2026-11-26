from playwright.sync_api import sync_playwright

def run(playwright):
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page()
    print("Navigating...")
    page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
    page.wait_for_selector("table", timeout=10000)
    
    html = page.evaluate('document.querySelector("tbody").innerHTML')
    print("HTML:")
    print(html[:2000])
    browser.close()

with sync_playwright() as p:
    run(p)
