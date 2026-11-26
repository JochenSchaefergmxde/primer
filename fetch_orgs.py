import urllib.request
import urllib.parse
import json

data = urllib.parse.urlencode({'database': 'qprimerdb2', 'method': 'GetOrgList'}).encode('utf-8')
headers = {
    'Content-Type': 'application/x-www-form-urlencoded',
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36'
}
req = urllib.request.Request("https://qprimerdb.biodb.org/bioapi/", data=data, headers=headers)
with urllib.request.urlopen(req) as f:
    resp = json.loads(f.read().decode('utf-8'))

orgs = []
for category, groups in resp.items():
    if isinstance(groups, dict):
        for group_name, items in groups.items():
            if isinstance(items, list):
                orgs.extend(items)

print(f"Total organisms found: {len(orgs)}")
if len(orgs) > 0:
    print(f"First organism: {orgs[0]}")
