from playwright.sync_api import sync_playwright
import zipfile
import io
import gzip
import sys

sys.stdout.reconfigure(encoding='utf-8')

def test_redux():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        print("[*] Opening qprimerdb downloads...")
        page.goto("https://qprimerdb.biodb.org/downloads", wait_until="networkidle")
        page.wait_for_selector("table", timeout=15000)

        # Dispatch Redux action for best primers of E. coli / Acyrthosiphon pisum
        test_key = "animal___invertebrate___Acyrthosiphon_pisum___GCF_005508785.1_pea_aphid_22Mar2018_4r6ur"
        print(f"[*] Dispatching Redux action for key: {test_key}...")

        url = page.evaluate("""
        async (key) => {
            window.g_app._store.dispatch({
                type: 'fetchData/downloadFile',
                payload: {
                    params: {
                        database: 'qprimerdb2',
                        method: 'DownloadFile',
                        orgKey: key,
                        fType: 'best'
                    }
                }
            });

            // Wait for dlLink in state
            for (let i = 0; i < 30; i++) {
                await new Promise(r => setTimeout(r, 500));
                let state = window.g_app._store.getState();
                if (state.fetchData && state.fetchData.dlLink) {
                    return state.fetchData.dlLink;
                }
            }
            return null;
        }
        """, test_key)

        print(f"[+] Redux returned dlLink: '{url}'")

        if url:
            print("[*] Downloading file content via page.request...")
            res = page.request.get(url, timeout=60000)
            body = res.body()
            print(f"[+] Downloaded {len(body)} bytes!")
            with zipfile.ZipFile(io.BytesIO(body)) as z:
                print("Zip files inside:", z.namelist())
                for fname in z.namelist():
                    with gzip.GzipFile(fileobj=io.BytesIO(z.read(fname))) as gz:
                        lines = gz.read().decode('utf-8', errors='ignore').splitlines()
                        print("Header:", lines[0])
                        print("Sample Row:", lines[1] if len(lines)>1 else 'EMPTY')

        browser.close()

if __name__ == "__main__":
    test_redux()
