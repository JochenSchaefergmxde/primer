import urllib.request
import re
import xml.etree.ElementTree as ET

url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13479736/fullTextXML'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req) as resp:
    xml_data = resp.read()

root = ET.fromstring(xml_data)

# Look for table-wrap elements
tables = root.findall('.//table-wrap')
print(f"Tables found: {len(tables)}")
for i, t in enumerate(tables):
    title = "".join(t.itertext())
    print(f"\n--- Table {i+1} ---")
    print(title[:300].strip())

# Search for primer sequence patterns in XML text
# Usually 15-35 bp DNA sequences with standard bases A, C, G, T and degenerate R, Y, M, K, S, W, H, B, V, D, N
primer_pattern = re.compile(r'\b[ACGTURYMKSWBDHVNacgturymkswbdhvn]{15,35}\b')
matches = primer_pattern.findall(xml_data.decode('utf-8'))
print(f"\nPotential DNA sequences found: {len(matches)}")
print("Sample matches:", matches[:10])
