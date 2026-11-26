import json
import urllib.request
import zipfile
import gzip
import io
import csv
import re
import sys
import concurrent.futures

sys.stdout.reconfigure(encoding='utf-8')

# Load all 1,315 organisms metadata
with open('organisms.json', 'r', encoding='utf-8') as f:
    org_data = json.load(f)

org_list = []
keywords = ['candida', 'plasmodium', 'toxoplasma', 'leishmania', 'trypanosoma', 'babesia', 'cryptosporidium', 'giardia', 'entamoeba', 'trichomonas', 'aspergillus', 'cryptococcus', 'pneumocystis', 'histoplasma', 'coccidioides', 'blastomyces', 'paracoccidioides', 'sporothrix', 'microsporidia', 'acanthamoeba', 'naegleria', 'balantidium']

for cat, grps in org_data.items():
    if isinstance(grps, dict):
        for grp, items in grps.items():
            if isinstance(items, list):
                for item in items:
                    org_name = item['orgName']
                    if any(k in org_name.lower() for k in keywords):
                        clean_key = item['orgKey'].replace('___', '__')
                        url_all = f"https://qprimerdb.biodb.org/tmp/qprimerdb/downloads/{clean_key}.all.primers.gz.zip"
                        url_best = f"https://qprimerdb.biodb.org/tmp/qprimerdb/downloads/{clean_key}.best.primers.xls.zip"
                        org_list.append({
                            'name': org_name,
                            'key': item['orgKey'],
                            'url_all': url_all,
                            'url_best': url_best
                        })

print(f"[+] Total organisms to process for ALL primers database (Filtered for Pathogens): {len(org_list)}")

def fetch_organism_all_primers(org):
    results = []
    # 1. Try downloading .all.primers.gz.zip
    req = urllib.request.Request(org['url_all'], headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
    try:
        data = urllib.request.urlopen(req, timeout=60).read()
        if len(data) > 1000:
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                fname = z.namelist()[0]
                gz_bytes = z.read(fname)
                content = gzip.GzipFile(fileobj=io.BytesIO(gz_bytes)).read().decode('utf-8', errors='ignore')
                lines = content.splitlines()
                for line in lines[1:]: # Skip header
                    parts = line.split('\t')
                    if len(parts) < 11:
                        parts = line.split(',')
                    if len(parts) >= 11:
                        fwd = parts[9].upper().strip()
                        rev = parts[10].upper().strip()
                        pair_id = parts[1].strip() if len(parts) > 1 else "pair"
                        if fwd and rev and re.match(r'^[ACGTRYSWKMBDHVN]+$', fwd) and re.match(r'^[ACGTRYSWKMBDHVN]+$', rev):
                            results.append({
                                'Organism': org['name'],
                                'PairID': f"{org['name']}_{pair_id}",
                                'Forward': fwd,
                                'Reverse': rev
                            })
                if results:
                    return results
    except Exception:
        pass

    # 2. Fallback to .best.primers.xls.zip if .all.primers not generated
    req_best = urllib.request.Request(org['url_best'], headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
    try:
        data_best = urllib.request.urlopen(req_best, timeout=30).read()
        if len(data_best) > 1000:
            with zipfile.ZipFile(io.BytesIO(data_best)) as z:
                fname = z.namelist()[0]
                content = z.read(fname).decode('utf-8', errors='ignore')
                lines = content.splitlines()
                for line in lines[1:]:
                    parts = line.split('\t')
                    if len(parts) >= 9:
                        fwd = parts[3].upper().strip()
                        rev = parts[8].upper().strip()
                        pair_id = parts[1].strip()
                        if fwd and rev and re.match(r'^[ACGTRYSWKMBDHVN]+$', fwd) and re.match(r'^[ACGTRYSWKMBDHVN]+$', rev):
                            results.append({
                                'Organism': org['name'],
                                'PairID': f"{org['name']}_{pair_id}",
                                'Forward': fwd,
                                'Reverse': rev
                            })
    except Exception:
        pass

    return results

import os
import threading

def main():
    print("[*] Fetching ALL candidate primers across all 1,315 organisms using 10 parallel threads...")
    
    out_csv_path = "all_primers_db.csv"
    best_csv_path = "best_primers_db.csv"
    
    # Initialize the CSV with headers
    with open(out_csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
        writer.writeheader()
        
    total_primers = 0
    write_lock = threading.Lock()
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        future_to_org = {executor.submit(fetch_organism_all_primers, org): org for org in org_list}
        completed = 0
        for future in concurrent.futures.as_completed(future_to_org):
            res = future.result()
            if res:
                # Write to disk immediately to save memory
                with write_lock:
                    with open(out_csv_path, 'a', newline='', encoding='utf-8') as f:
                        writer = csv.DictWriter(f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
                        writer.writerows(res)
                
                # Let's save memory by not storing all PairIDs in RAM.
                total_primers += len(res)
                
                # Delete `res` so Python GC can reclaim memory
                del res
                
            completed += 1
            if completed % 10 == 0 or completed == len(org_list):
                print(f"    - Processed {completed}/{len(org_list)} organisms (Total primer pairs collected: {total_primers:,})...", flush=True)

    print(f"\n[+] Total Primer Pairs Collected for FULL Database: {total_primers:,}")
    
    # Also load best_primers_db.csv to include diagnostic primers (viruses, Babesia, Borrelia, etc.)
    if os.path.exists(best_csv_path):
        with open(best_csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            added_extra = 0
            with open(out_csv_path, 'a', newline='', encoding='utf-8') as out_f:
                writer = csv.DictWriter(out_f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
                for row in reader:
                    # We just append them. We don't check for duplicates to save RAM.
                    # CD-HIT or downstream tools can handle redundant identical primers easily.
                    writer.writerow(row)
                    added_extra += 1
            print(f"[+] Appended {added_extra} diagnostic primer pairs (viruses/pathogens) from best_primers_db.csv")

    print(f"[+] Successfully saved FULL Primer Database to: {out_csv_path}")

if __name__ == "__main__":
    main()
