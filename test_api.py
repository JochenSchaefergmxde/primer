import urllib.request
import json
import traceback

endpoints = ['organisms', 'download', 'downloads', 'species', 'primer', 'primers', 'data', 'v1/organisms']

for ep in endpoints:
    url = f"https://qprimerdb.biodb.org/bioapi/{ep}"
    print(f"Trying {url} ...")
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            print("SUCCESS!")
            if isinstance(data, list):
                print(f"List of {len(data)} items. First item keys: {list(data[0].keys())}")
                print(data[0])
            elif isinstance(data, dict):
                print(f"Dict keys: {list(data.keys())}")
                if 'data' in data and isinstance(data['data'], list) and len(data['data']) > 0:
                    print(f"data[0] keys: {list(data['data'][0].keys())}")
            break
    except Exception as e:
        print("Failed:", str(e))
