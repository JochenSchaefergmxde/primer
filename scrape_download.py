from playwright.sync_api import sync_playwright

def run(playwright):
    browser = playwright.chromium.launch(headless=True)
    context = browser.new_context(accept_downloads=True)
    page = context.new_page()
    page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
    
    page.wait_for_selector("table", timeout=10000)
    
    # We click the first "Best Primer" button (which is the third column, i.e., second button)
    print("Clicking download...")
    with page.expect_download() as download_info:
        # 3rd column in the first data row
        page.locator('tr.ant-table-row:nth-child(2) td:nth-child(3) button').click()
    
    download = download_info.value
    print(f"Download URL: {download.url}")
    print(f"Suggested filename: {download.suggested_filename}")
    browser.close()

with sync_playwright() as p:
    run(p)
