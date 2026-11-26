#!/usr/bin/env python3
"""
Europe PMC Primer Harvester
Automated extraction, strict in-text validation, and assembly of pathogen PCR primers from Open Access literature.
"""

import sys
import os
import re
import shutil
import argparse
import urllib.request
import urllib.parse
import json
import csv
from datetime import datetime

sys.stdout.reconfigure(encoding='utf-8')

IUPAC_PATTERN = re.compile(r'^[ACGTNRYMKSwbdhvnACGTURYMKSWBDHVN]+$')
DOWNLOADS_DIR = r"C:\Users\internet\Downloads"

def clean_seq(s):
    if not s or not isinstance(s, str):
        return ""
    # Strip orientation labels (5', 3', etc.)
    s = re.sub(r"5['′’]?", "", s)
    s = re.sub(r"3['′’]?", "", s)
    # Strip spaces, dashes, parentheses, non-letters
    s = re.sub(r"[-–—\s\d\(\)\[\]:/,\.]", "", s)
    s = s.upper().replace('U', 'T')
    if 15 <= len(s) <= 45 and IUPAC_PATTERN.match(s):
        return s
    return ""

def classify_role(name, context=""):
    combined = (name + " " + context).lower()
    if any(p in combined for p in ['probe', 'sonde', 'taqman', '_pr_', '-p', '_p', 'probe:']):
        return 'P'
    if any(f in combined for f in ['forward', 'sense', '_f', '-f', ' f1', ' f2', ' f3', 'forw', '(f)', '_for', 'fw', '-fw']):
        return 'F'
    if any(r in combined for r in ['reverse', 'antisense', '_r', '-r', ' r1', ' r2', ' r3', 'rev', '(r)', '_rev', 'rv', '-rv']):
        return 'R'
    return 'UNKNOWN'

def fetch_search_results(query, page_size=25, max_pages=1):
    all_results = []
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) PrimerHarvester/1.0'}
    
    cursor = "*"
    for page in range(max_pages):
        params = {
            'query': f"{query} AND OPEN_ACCESS:y",
            'format': 'json',
            'pageSize': page_size,
            'resultType': 'core',
            'cursorMark': cursor
        }
        url = f"{base_url}?{urllib.parse.urlencode(params)}"
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                results = data.get("resultList", {}).get("result", [])
                if not results:
                    break
                all_results.extend(results)
                cursor = data.get("nextCursorMark")
                if not cursor:
                    break
        except Exception as e:
            print(f"[!] Search query error: {e}")
            break
    return all_results

def fetch_fulltext_xml(pmcid):
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    headers = {'User-Agent': 'Mozilla/5.0 PrimerHarvester/1.0'}
    try:
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=25) as resp:
            return resp.read().decode('utf-8', errors='ignore')
    except Exception as e:
        return None

def parse_primer_table(table_html, xml_full, article_meta):
    extracted_pairs = []
    rows = re.findall(r'<tr.*?>.*?</tr>', table_html, re.DOTALL)
    if not rows:
        return extracted_pairs

    current_organism = article_meta.get("target_organism", "")
    current_gene = ""
    current_size = ""
    
    # Track pending forward primer to pair with upcoming reverse
    pending_forward = None

    for r in rows:
        cells = re.findall(r'<t[dh].*?>(.*?)</t[dh]>', r, re.DOTALL)
        clean_cells = [re.sub(r'<[^>]+>', ' ', c).strip() for c in cells]
        clean_cells = [" ".join(c.split()) for c in clean_cells if c]
        if not clean_cells:
            continue

        row_text = " | ".join(clean_cells)
        # Check for potential target gene or organism in row
        gene_match = re.search(r'\b(16S|18S|23S|ITS\d?|cytb|cox\d?|osp[A-C]|gltA|rpoB|hsp65|IS6110|groEL|flaS|p\d+|dsb|purB|gap)\b', row_text, re.I)
        if gene_match:
            current_gene = gene_match.group(0)

        # Check for amplicon size (e.g. "232 bp" or "232.00")
        size_match = re.search(r'\b(\d{2,4})\s*(?:bp)?\b', row_text)
        if size_match and int(size_match.group(1)) > 50 and int(size_match.group(1)) < 3000:
            current_size = size_match.group(1)

        # Find cells that have valid primer sequences
        seq_cells = []
        for idx, c in enumerate(clean_cells):
            seq = clean_seq(c)
            if seq:
                # Anti-hallucination check: Sequence must exist verbatim in the original full-text XML
                if seq in xml_full.upper():
                    seq_cells.append((idx, c, seq))

        if not seq_cells:
            continue

        # Scenario 1: Row has both Forward and Reverse in distinct columns
        if len(seq_cells) >= 2:
            f_seq = None
            r_seq = None
            f_name = "Forward"
            r_name = "Reverse"
            for idx, raw_val, s in seq_cells:
                role = classify_role(raw_val, row_text)
                if role == 'F' and not f_seq:
                    f_seq = s
                    f_name = clean_cells[idx-1] if idx > 0 else "F"
                elif role == 'R' and not r_seq:
                    r_seq = s
                    r_name = clean_cells[idx-1] if idx > 0 else "R"

            if not f_seq or not r_seq:
                # Default to first = F, second = R
                f_seq = seq_cells[0][2]
                r_seq = seq_cells[1][2]
                f_name = clean_cells[seq_cells[0][0]-1] if seq_cells[0][0] > 0 else "F"
                r_name = clean_cells[seq_cells[1][0]-1] if seq_cells[1][0] > 0 else "R"

            extracted_pairs.append({
                "Organism": current_organism or article_meta.get("target_organism", ""),
                "Target_Gene": current_gene or "PCR",
                "PairID": f"{f_name}_{r_name}".replace(" ", "_"),
                "Forward": f_seq,
                "Reverse": r_seq,
                "Amplicon_Size": current_size,
                "PMCID": article_meta.get("pmcid", ""),
                "DOI": article_meta.get("doi", ""),
                "Title": article_meta.get("title", "")
            })
            pending_forward = None

        # Scenario 2: Row has a single sequence (e.g. Forward on row 1, Reverse on row 2)
        elif len(seq_cells) == 1:
            idx, raw_val, s = seq_cells[0]
            name = clean_cells[idx-1] if idx > 0 else clean_cells[0]
            role = classify_role(name, row_text)

            if role == 'F':
                pending_forward = {
                    "name": name,
                    "seq": s,
                    "gene": current_gene,
                    "organism": current_organism or article_meta.get("target_organism", ""),
                    "size": current_size
                }
            elif role == 'R' and pending_forward:
                # Pair with pending forward!
                extracted_pairs.append({
                    "Organism": pending_forward["organism"],
                    "Target_Gene": pending_forward["gene"] or current_gene or "PCR",
                    "PairID": f"{pending_forward['name']}_{name}".replace(" ", "_"),
                    "Forward": pending_forward["seq"],
                    "Reverse": s,
                    "Amplicon_Size": current_size or pending_forward["size"],
                    "PMCID": article_meta.get("pmcid", ""),
                    "DOI": article_meta.get("doi", ""),
                    "Title": article_meta.get("title", "")
                })
                pending_forward = None

    return extracted_pairs

def harvest_pathogens(targets, pages_per_target=2, output_csv="harvested_pathogen_primers.csv"):
    all_pairs = []
    seen_pairs = set()

    for target in targets:
        print(f"\n=======================================================")
        print(f"[*] Querying Europe PMC for: {target}")
        print(f"=======================================================")
        query = f'"{target}" PCR primer'
        articles = fetch_search_results(query, page_size=15, max_pages=pages_per_target)
        print(f"[+] Found {len(articles)} open-access articles.")

        for art in articles:
            pmcid = art.get("pmcid")
            if not pmcid:
                continue
            title = art.get("title", "")
            doi = art.get("doi", "")
            print(f" -> Scanning {pmcid}: {title[:75]}...")

            xml = fetch_fulltext_xml(pmcid)
            if not xml:
                continue

            art_meta = {
                "target_organism": target,
                "pmcid": pmcid,
                "doi": doi,
                "title": title
            }

            tables = re.findall(r'<table-wrap.*?</table-wrap>', xml, re.DOTALL)
            for t in tables:
                pairs = parse_primer_table(t, xml, art_meta)
                for p in pairs:
                    pair_key = (p["Forward"], p["Reverse"])
                    if pair_key not in seen_pairs:
                        seen_pairs.add(pair_key)
                        all_pairs.append(p)
                        print(f"    [+] Discovered Pair: {p['Target_Gene']} ({p['PairID']}) -> F: {p['Forward']} | R: {p['Reverse']}")

    # Write output CSV
    fieldnames = ["Organism", "Target_Gene", "PairID", "Forward", "Reverse", "Amplicon_Size", "PMCID", "DOI", "Title"]
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for p in all_pairs:
            writer.writerow(p)

    print(f"\n[DONE] Successfully harvested and verified {len(all_pairs)} unique primer pairs!")
    print(f"[+] Saved locally to: {os.path.abspath(output_csv)}")

    # Copy to user's Downloads directory per global instructions
    dl_path = os.path.join(DOWNLOADS_DIR, os.path.basename(output_csv))
    shutil.copy2(output_csv, dl_path)
    print(f"[+] Copied to Downloads directory: {dl_path}")
    return dl_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Europe PMC Primer Harvester")
    parser.add_argument("--targets", nargs="+", default=["Babesia", "Borrelia", "Bartonella", "Mycobacterium"], help="Pathogens to search for")
    parser.add_argument("--pages", type=int, default=2, help="Pages per target (15 articles/page)")
    parser.add_argument("--out", default="harvested_pathogen_primers.csv", help="Output CSV path")
    args = parser.parse_args()

    harvest_pathogens(args.targets, pages_per_target=args.pages, output_csv=args.out)
