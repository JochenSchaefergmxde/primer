from playwright.sync_api import sync_playwright

def run(playwright):
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page()
    
    def log_request(request):
        if "bioapi" in request.url and request.method == "POST":
            print(f"> POST {request.url}")
            print(f"Payload: {request.post_data}")
            
    def log_response(response):
        if "bioapi" in response.url and response.request.method == "POST":
            print(f"< POST {response.url} Status: {response.status}")
            try:
                print(f"Response: {response.text()[:500]}")
            except:
                pass

    page.on("request", log_request)
    page.on("response", log_response)

    print("Navigating...")
    page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
    page.wait_for_selector("table", timeout=10000)
    browser.close()

with sync_playwright() as p:
    run(p)
