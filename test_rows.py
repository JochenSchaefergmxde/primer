import sys
import urllib.request
import re

sys.stdout.reconfigure(encoding='utf-8')

url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13517960/fullTextXML'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req) as resp:
    text = resp.read().decode('utf-8', errors='ignore')

tables = re.findall(r'<table-wrap.*?</table-wrap>', text, re.DOTALL)
if tables:
    rows = re.findall(r'<tr.*?>.*?</tr>', tables[0], re.DOTALL)
    print(f"Total rows in table 1: {len(rows)}")
    for r in rows:
        cells = re.findall(r'<t[dh].*?>(.*?)</t[dh]>', r, re.DOTALL)
        clean_cells = [re.sub(r'<[^>]+>', '', c).strip() for c in cells]
        print("ROW:", clean_cells)
