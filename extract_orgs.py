from playwright.sync_api import sync_playwright
import json

def run(playwright):
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page()
    
    def log_response(response):
        if "bioapi" in response.url and response.request.method == "POST":
            try:
                data = response.json()
                if "animal" in data or "plant" in data or "fungi" in data:
                    with open(r"C:\Users\internet\Downloads\repository\primer\organisms.json", "w", encoding="utf-8") as f:
                        json.dump(data, f)
                    print("Successfully saved organisms.json")
            except Exception as e:
                pass

    page.on("response", log_response)

    print("Navigating...")
    page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
    page.wait_for_selector("table", timeout=10000)
    browser.close()

with sync_playwright() as p:
    run(p)
