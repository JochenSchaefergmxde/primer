import sys
import re
import urllib.request
import urllib.parse
import json

sys.stdout.reconfigure(encoding='utf-8')

IUPAC_REGEX = re.compile(r'^[ACGTNRYMKSwbdhvnACGTURYMKSWBDHVN]+$')

def clean_seq(s):
    if not s:
        return ""
    # remove 5', 3', hyphens, spaces, digits
    s = re.sub(r"5['′’]?", "", s)
    s = re.sub(r"3['′’]?", "", s)
    s = re.sub(r"[-–—\s\d\(\)]", "", s)
    s = s.upper().replace('U', 'T')
    if 15 <= len(s) <= 45 and IUPAC_REGEX.match(s):
        return s
    return ""

def classify_primer(name, text_context=""):
    combined = (name + " " + text_context).lower()
    if any(p in combined for p in ['probe', 'sonde', 'taqman', '_pr_', '-p', '_p']):
        return 'probe'
    if any(f in combined for f in ['forward', 'sense', '_f', '-f', ' f1', ' f2', ' f3', 'forw', '(f)', '_for']):
        return 'forward'
    if any(r in combined for r in ['reverse', 'antisense', '_r', '-r', ' r1', ' r2', ' r3', 'rev', '(r)', '_rev']):
        return 'reverse'
    return 'unknown'

def search_epmc(query, max_results=10):
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={urllib.parse.quote(query)}&format=json&pageSize={max_results}&resultType=core"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            return data.get("resultList", {}).get("result", [])
    except Exception as e:
        print(f"Error searching {query}: {e}")
        return []

def get_xml(pmcid):
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            return resp.read().decode('utf-8', errors='ignore')
    except Exception as e:
        return None

# Test on Borrelia and Mycobacterium as well
queries = [
    'Babesia PCR primer OPEN_ACCESS:y',
    'Borrelia burgdorferi PCR primer OPEN_ACCESS:y',
    'Mycobacterium tuberculosis PCR primer OPEN_ACCESS:y'
]

for q in queries:
    results = search_epmc(q, max_results=3)
    print(f"\n================ Query: {q} ================")
    for r in results:
        pmcid = r.get("pmcid")
        if not pmcid:
            continue
        title = r.get("title", "")
        print(f"\nFetching {pmcid}: {title[:70]}...")
        xml = get_xml(pmcid)
        if not xml:
            continue
        
        tables = re.findall(r'<table-wrap.*?</table-wrap>', xml, re.DOTALL)
        found_primers = 0
        for t in tables:
            rows = re.findall(r'<tr.*?>.*?</tr>', t, re.DOTALL)
            for row in rows:
                cells = re.findall(r'<t[dh].*?>(.*?)</t[dh]>', row, re.DOTALL)
                clean_cells = [re.sub(r'<[^>]+>', '', c).strip() for c in cells]
                for c in clean_cells:
                    seq = clean_seq(c)
                    if seq:
                        # Validate that seq literally exists in the article xml!
                        assert seq in xml.upper()
                        found_primers += 1
                        print(f"  [Found seq in table] {seq} (len={len(seq)}) in row: {clean_cells[:3]}")
        print(f"Total verified primers found in {pmcid}: {found_primers}")
