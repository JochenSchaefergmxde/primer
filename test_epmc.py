import urllib.request
import urllib.parse
import json

query = 'Babesia AND PCR AND primer AND OPEN_ACCESS:y'
url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={urllib.parse.quote(query)}&format=json&pageSize=5&resultType=core"
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
with urllib.request.urlopen(req) as resp:
    data = json.loads(resp.read().decode('utf-8'))

hit_count = data.get("hitCount", 0)
results = data.get("resultList", {}).get("result", [])
print(f"Total hits: {hit_count}")
for r in results:
    print(f"ID: {r.get('id')}, PMCID: {r.get('pmcid')}, Title: {r.get('title')[:80]}...")
