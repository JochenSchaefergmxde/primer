import urllib.request
import re

pmcids = ['PMC13517960', 'PMC13541529', 'PMC13603127', 'PMC13610672']

for pmcid in pmcids:
    url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req) as resp:
            text = resp.read().decode('utf-8', errors='ignore')
        
        has_table = "<table" in text.lower()
        # Find DNA sequences: 16 to 35 base pairs
        # Often written like: 5'-GACGATCAGATACCGTCGTAGTC-3' or in table cells <td>GACGATCAGATACCGTCGTAGTC</td>
        raw_seqs = re.findall(r'>\s*([ACGTNRYMKSwbdhvn]{16,35})\s*<', text)
        p53 = re.findall(r"5['′’]?\s*[-–—]?\s*([ACGTNRYMKSwbdhvn]{16,35})\s*[-–—]?\s*3['′’]?", text, re.I)
        print(f"{pmcid}: length={len(text)}, tables={has_table}, in_tags={len(raw_seqs)}, 5to3={len(p53)}")
        if raw_seqs:
            print(f"  Sample in_tags: {raw_seqs[:5]}")
        if p53:
            print(f"  Sample 5to3: {p53[:5]}")
    except Exception as e:
        print(f"{pmcid} error: {e}")
