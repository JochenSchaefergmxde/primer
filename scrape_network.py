from playwright.sync_api import sync_playwright

def run(playwright):
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page()
    
    # Listen to all network requests
    page.on("request", lambda request: print(f"> {request.method} {request.url}"))
    page.on("response", lambda response: print(f"< {response.status} {response.url}"))

    print("Navigating to downloads page...")
    page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
    page.wait_for_selector("table", timeout=10000)
    browser.close()

with sync_playwright() as p:
    run(p)
