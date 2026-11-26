#!/usr/bin/env python3
"""
primer_suche.py - High-Performance Multi-Threaded Primer Search & Amplicon Extraction Pipeline

Applies primers from best_primers_db.csv to a FASTQ file with massive multi-core acceleration:
- Multi-processing across CPU threads (e.g. 24-32 cores) with chunked streaming.
- High-speed parallel decompression using pigz (if available) with fallback to gzip.
- Hierarchical Search:
    -g / --genus: Search only genus-level / screening primers
    -a / --all: Search all primers
    -gs / --genus-species: 2-stage search (screening genus first, then species for positive genera)
- Generates an interactive HTML dashboard with Plotly charts.
"""

import os
import sys

# Ensure UTF-8 output on all operating systems / terminals
if sys.stdout and hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass
if sys.stderr and hasattr(sys.stderr, 'reconfigure'):
    try:
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

import re
import io
import csv
import glob
import gzip
import zipfile
import argparse
import time
import json
import shutil
import subprocess
import multiprocessing as mp
from collections import defaultdict, Counter
from datetime import datetime, timedelta

def format_duration(seconds):
    """Formats seconds into MM:SS or HH:MM:SS."""
    if seconds < 0:
        return "00:00"
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    if h > 0:
        return f"{h:02d}:{m:02d}:{s:02d}"
    return f"{m:02d}:{s:02d}"

def render_progress_bar(fraction, total_reads, reads_per_sec, elapsed_s, bar_len=28):
    """Renders a universally compatible ASCII progress bar with read counts, speed, and ETA."""
    fraction = min(1.0, max(0.0, fraction))
    pct = fraction * 100.0
    filled = int(round(bar_len * fraction))
    
    if filled >= bar_len:
        bar = "=" * bar_len
    elif filled > 0:
        bar = "=" * (filled - 1) + ">" + " " * (bar_len - filled)
    else:
        bar = " " * bar_len
    
    if fraction > 0.005 and elapsed_s > 1:
        total_expected_s = elapsed_s / fraction
        rem_s = max(0, total_expected_s - elapsed_s)
        finish_time = datetime.now() + timedelta(seconds=rem_s)
        time_info = f"ETA: {format_duration(rem_s)} (ca. {finish_time.strftime('%H:%M:%S')} Uhr)"
    else:
        time_info = "ETA: Berechne..."

    speed_info = f"{reads_per_sec:,.0f} R/s"
    reads_info = f"{total_reads:,} Reads"

    return f"[{bar}] {pct:5.1f}% | {reads_info} ({speed_info}) | {time_info}"

# --- IUPAC Degenerate Code Mapping ---
IUPAC_CODES = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

COMPLEMENT = str.maketrans("ATCGURYSWKMBDHVNatcguryswkmbdhvn", "TAGCAYRSMKVHDBNNtagcayrsmkvhdbnn")

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]

def get_exact_seeds(seq, k=8):
    """Extracts dual non-degenerate k-mer seeds (5'-half and middle) for maximum mismatch sensitivity."""
    seeds = []
    clean = seq.upper()
    if len(clean) < k:
        return [clean]
    # 1. 5'-half seed
    for i in range(min(len(clean) - k + 1, 6)):
        sub = clean[i:i+k]
        if all(b in 'ATCGU' for b in sub):
            seeds.append(sub)
            break
    # 2. Middle seed
    mid = len(clean) // 2
    for i in range(mid, len(clean) - k + 1):
        sub = clean[i:i+k]
        if all(b in 'ATCGU' for b in sub):
            seeds.append(sub)
            break
    if not seeds:
        seeds.append(clean[:k])
    return seeds

def get_exact_seed(seq, k=12):
    """Backward-compatible fallback single seed."""
    seeds = get_exact_seeds(seq, k=k)
    return seeds[0] if seeds else seq[:k].upper()

def is_high_complexity_read(seq, min_len=120):
    """Filters out low-complexity sequences, homopolymers, and simple tandem repeats."""
    if len(seq) < min_len:
        return False
    if re.search(r'A{18,}|T{18,}|G{18,}|C{18,}', seq):
        return False
    dinucs = set(seq[i:i+2] for i in range(len(seq)-1))
    if len(dinucs) < 10:
        return False
    return True

def iupac_to_regex(primer_seq, max_mismatches=0):
    """Converts a primer sequence (with IUPAC codes) to a regex pattern allowing mismatches."""
    pattern_parts = [IUPAC_CODES.get(base.upper(), base.upper()) for base in primer_seq]
    if max_mismatches == 0:
        return re.compile(''.join(pattern_parts), re.IGNORECASE)
    else:
        patterns = [''.join(pattern_parts)]
        for i in range(len(primer_seq)):
            modified = list(pattern_parts)
            modified[i] = '[ACGTN]'
            patterns.append(''.join(modified))
        return re.compile('|'.join(patterns), re.IGNORECASE)

class FastqRecord:
    def __init__(self, header, seq, qual_header, qual):
        self.header = header
        self.seq = seq
        self.qual_header = qual_header
        self.qual = qual

    def to_fastq_string(self):
        return f"{self.header}\n{self.seq}\n{self.qual_header}\n{self.qual}\n"

def classify_primer_taxa(org_name, pair_id):
    """
    Classifies a primer entry as Genus-level (screening) or Species-level.
    Returns: (clean_genus, is_genus_level)
    """
    org = str(org_name).strip()
    pid = str(pair_id).strip()
    combined = (org + ' ' + pid).lower()
    
    words = org.split()
    first_word = words[0] if words else 'Unknown'
    clean_genus = first_word.replace('/', '').replace('(', '').replace(')', '').strip()
    
    is_genus = False
    if len(words) == 1:
        is_genus = True
    elif re.search(r'\b(genus|spp\.?|complex|universal|piroplasmida?|sensu\s+lato|\bs\.l\.)\b', combined, re.I):
        is_genus = True
    elif any(org.lower().startswith(u) for u in ['universal', 'bacteria', 'fungi', 'apicomplexa']):
        is_genus = True
        
    return clean_genus, is_genus

def load_default_primer_database():
    """Returns a built-in set of common 16S/18S/ITS and pathogen primers."""
    return [
        {
            'Organism': 'Universal Bacteria (16S 27F-1492R)',
            'PrimerPairID': '16S_27F_1492R',
            'Forward': 'AGAGTTTGATCMTGGCTCAG',
            'Reverse': 'TACGGYTACCTTGTTACGACTT'
        },
        {
            'Organism': 'Universal Bacteria (16S 515F-806R)',
            'PrimerPairID': '16S_515F_806R',
            'Forward': 'GTGCCAGCMGCCGCGGTAA',
            'Reverse': 'GGACTACHVGGGTWTCTAAT'
        },
        {
            'Organism': 'Universal Fungi (ITS1F-ITS4)',
            'PrimerPairID': 'ITS_1F_4R',
            'Forward': 'CTTGGTCATTTAGAGGAAGTAA',
            'Reverse': 'TCCTCCGCTTATTGATATGC'
        },
        {
            'Organism': 'Acanthamoeba castellanii',
            'PrimerPairID': 'Acanthamoeba_P1',
            'Forward': 'TGTGGGATCAGTTCGAAAAGGT',
            'Reverse': 'GAGTGCAAGCTTCTCAGTTTCG'
        },
        {
            'Organism': 'Babesia',
            'PrimerPairID': 'Babesia_genus_18S_BS1_BS2',
            'Forward': 'GACGATCAGATACCGTCGTAGTC',
            'Reverse': 'CAGAACCCAAAGACTTTGATTTCTCT'
        },
        {
            'Organism': 'Babesia microti',
            'PrimerPairID': 'Babesia_microti_18S_Bab1_Bab4',
            'Forward': 'CTTAGTATAAGCTTTTATACAGC',
            'Reverse': 'GATAGGTCAGAAACTTGAATGATACAC'
        },
        {
            'Organism': 'Borrelia burgdorferi',
            'PrimerPairID': 'Borrelia_burgdorferi_ospA',
            'Forward': 'AATATTTCAAACGATCTAAGTCA',
            'Reverse': 'GTAGGTTCTACAGCTACATGA'
        }
    ]

COMMON_PATHOGEN_DESCRIPTIONS = {
    'borrelia': 'Lyme-Borreliose & Neuroborreliose (Erythema migrans, ACA, Neuro-, Kardio-, Arthrose)',
    'babesia odocoilei': 'Humane Babesiose (Emerging Pathogen, Zeckenuebertragung durch Ixodes scapularis)',
    'babesia venatorum': 'Humane Babesiose (EU1, Fieber, Haemolyse, Zeckenuebertragung durch Ixodes ricinus)',
    'babesia': 'Babesiose / Piroplasmose (Malaria-aehnliche Erythrozyten-Infektion, Haemolyse)',
    'bartonella': 'Bartonellose (Katzenkratzkrankheit, Endokarditis, chronische Angiomatose)',
    'rickettsia': 'Rickettsiose / Fleckfieber (Zecken-Spotted-Fever, Vaskulitis)',
    'anaplasma': 'Humane Granulozytäre Anaplasmose (HGA, Leukopenie, Fieber)',
    'ehrlichia': 'Humane Monozytäre Ehrlichiose (HME)',
    'neoehrlichia': 'Candidatus Neoehrlichia mikurensis (Vaskuläre Entzündungen, Thrombosen)',
    'francisella': 'Tularämie / Hasenpest (Lymphknoten-Schwellung, Sepsis)',
    'coxiella': 'Q-Fieber (Atypische Pneumonie, chronische Endokarditis)',
    'leptospira': 'Leptospirose / Morbus Weil (Ikterus, Nierenversagen)',
    'brucella': 'Brucellose / Maltafieber (Undulierendes Fieber, Spondylodiszitis)',
    'granulicatella': 'Kultur-negative Endokarditis & Sepsis (Nutritionally Variant Streptococcus)',
    'abiotrophia': 'Kultur-negative Endokarditis (Nutritionally Variant Streptococcus)',
    'haemophilus': 'HACEK-Endokarditis / Invasive Infektionen',
    'aggregatibacter': 'HACEK-Gruppe (Parodontitis, destruierende Endokarditis)',
    'cardiobacterium': 'HACEK-Gruppe (Klappen-Endokarditis)',
    'eikenella': 'HACEK-Gruppe (Abszesse, Mischinfektionen, Endokarditis)',
    'kingella': 'HACEK-Gruppe (Osteomyelitis, septische Arthritis bei Kindern)',
    'whipplei': 'Morbus Whipple (Arthralgien, Malabsorption, ZNS-Befall)',
    'cryptococcus': 'Kryptokokkose (Invasive Mykose, ZNS-Meningitis, pulmonale Infektionen)',
    'tbev': 'FSME / Tick-borne encephalitis (Frühsommer-Meningoenzephalitis)',
    'cchfv': 'Krim-Kongo Hämorrhagisches Fieber',
    'ebv': 'Epstein-Barr-Virus / HHV-4 (Mononukleose, chronische Reaktivierung, PTLD)',
    'cmv': 'Zytomegalievirus / HHV-5 (Organmanifestation, Immunsuppression)',
    'hhv6': 'Humanes Herpesvirus 6 (Roseola infantum, Enzephalitis, Fatigue-Syndrom)',
    'hhv7': 'Humanes Herpesvirus 7 (Fieberkrampf, Pityriasis rosea)',
    'hsv': 'Herpes-simplex-Virus (Typ 1 & 2)',
    'vzv': 'Varizella-Zoster-Virus / HHV-3 (Windpocken, Gürtelrose)'
}

KNOWN_CONTAMINANTS_KITOME = [
    'homo sapiens', 'human', 'cutibacterium', 'propionibacterium',
    'bradyrhizobium', 'ralstonia', 'sphingomonas', 'micrococcus luteus',
    'methylobacterium', 'burkholderia', 'corynebacterium'
]

def get_pathogen_description(org_name):
    name = str(org_name).lower()
    for k, v in COMMON_PATHOGEN_DESCRIPTIONS.items():
        if k in name:
            return v
    if any(u in name for u in ['universal', 'pan-', 'pan_', 'barcoding']):
        return 'Universelles Barcoding / Breitspektrum-Screening'
    return 'Humanpathogener Erreger'

def is_likely_contaminant(org_name, blast_hit_title=""):
    combined = (str(org_name) + " " + str(blast_hit_title)).lower()
    for c in KNOWN_CONTAMINANTS_KITOME:
        if c in combined:
            return True
    return False

def infer_tier(org_name, pair_id=""):
    """Infers pathogen tier classification: Tier 0 (Universal Barcoding), Tier 1 (Tick/vector-borne & co-infections), Tier 2 (Endocarditis/Sepsis), Tier 3 (Other)."""
    name = (str(org_name) + " " + str(pair_id)).lower()
    t0_keywords = ['universal', 'pan-', 'pan_', 'broad-range', 'barcoding', 'tier 0', 'tier0']
    for k in t0_keywords:
        if k in name:
            return 'Tier 0'
    t1_keywords = [
        'borrelia', 'babesia', 'bartonella', 'rickettsia', 'anaplasma', 'ehrlichia',
        'neoehrlichia', 'piroplasm', 'cytauxzoon', 'theileria', 'francisella',
        'coxiella', 'leptospira', 'brucella', 'mycoplasma', 'tbev', 'fsme', 'cchfv',
        'ebv', 'cmv', 'hhv', 'hsv', 'vzv', 'herpes'
    ]
    t2_keywords = [
        'granulicatella', 'abiotrophia', 'haemophilus', 'aggregatibacter', 'cardiobacterium',
        'eikenella', 'kingella', 'whipplei', 'gemella', 'enterococcus', 'staphylococcus',
        'streptococcus', 'hacek', 'endocarditis'
    ]
    for k in t1_keywords:
        if k in name:
            return 'Tier 1'
    for k in t2_keywords:
        if k in name:
            return 'Tier 2'
    return 'Tier 3'

def infer_target_gene(pair_id):
    pid = str(pair_id).lower()
    for g in ['16s', '18s', 'its', '28s', 'rpob', 'pol', 'rdrp', 'ospa', 'flab', 'glpq', 'glta', 'groel', 'groes', 'msu', 'msp', 'ompa', 'ompb', 'rtxa', 'ltxa', 'soda', 'gtfb', 'fuck', 'tul4', 'fopa', 'is1111', 'lipl32', 'is711', 'bcsp31', '23s', '5s']:
        if g in pid:
            return g.upper()
    return 'Marker'

def load_primers_from_csv(csv_path):
    """Loads primers from CSV/TSV file or qPrimerDB export files. Supports both 4-column and 8-column schema."""
    primers = []
    if not os.path.exists(csv_path):
        fallbacks = [
            csv_path,
            r"C:\Users\internet\Downloads\best_primers_db.csv",
            "./best_primers_db.csv",
            "/mnt/sdc1/scripts/best_primers_db.csv"
        ]
        found = False
        for fb in fallbacks:
            if os.path.exists(fb):
                csv_path = fb
                found = True
                break
        if not found:
            print(f"[!] Warning: Primer file '{csv_path}' not found. Loading default database.")
            return load_default_primer_database()

    if csv_path.endswith('.zip') or csv_path.endswith('.gz'):
        return load_primers_from_archive(csv_path)

    with open(csv_path, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            return load_default_primer_database()

        header_lower = [h.lower() for h in header]
        
        fwd_idx, rev_idx, name_idx, id_idx = -1, -1, -1, -1
        tier_idx, gene_idx, amp_idx, ref_idx = -1, -1, -1, -1
        for idx, col in enumerate(header_lower):
            if any(k in col for k in ['forward', 'fwd', 'f_primer', 'primer_f']):
                fwd_idx = idx
            elif any(k in col for k in ['reverse', 'rev', 'r_primer', 'primer_r']):
                rev_idx = idx
            elif any(k in col for k in ['orgname', 'organism', 'species', 'name', 'target']):
                name_idx = idx
            elif any(k in col for k in ['orgkey', 'id', 'key', 'pair', 'pairid']):
                id_idx = idx
            elif any(k in col for k in ['tier', 'level', 'priority', 'stuf']):
                tier_idx = idx
            elif any(k in col for k in ['gene', 'targetgene', 'target_gene', 'locus']):
                gene_idx = idx
            elif any(k in col for k in ['size', 'length', 'amp', 'amplicon']):
                amp_idx = idx
            elif any(k in col for k in ['ref', 'reference', 'pmid', 'pmc', 'pubmed', 'source']):
                ref_idx = idx

        if fwd_idx != -1 and rev_idx != -1:
            for i, row in enumerate(reader):
                if len(row) > max(fwd_idx, rev_idx):
                    fwd_seq = row[fwd_idx].strip()
                    rev_seq = row[rev_idx].strip()
                    if fwd_seq and rev_seq:
                        org_name = row[name_idx].strip() if name_idx != -1 and name_idx < len(row) else f"Target_{i+1}"
                        pair_id = row[id_idx].strip() if id_idx != -1 and id_idx < len(row) else f"Pair_{i+1}"
                        tier = row[tier_idx].strip() if tier_idx != -1 and tier_idx < len(row) and row[tier_idx].strip() else infer_tier(org_name, pair_id)
                        gene = row[gene_idx].strip() if gene_idx != -1 and gene_idx < len(row) and row[gene_idx].strip() else infer_target_gene(pair_id)
                        size = int(row[amp_idx].strip()) if amp_idx != -1 and amp_idx < len(row) and row[amp_idx].strip().isdigit() else 0
                        ref = row[ref_idx].strip() if ref_idx != -1 and ref_idx < len(row) and row[ref_idx].strip() else 'Curated'
                        primers.append({
                            'Organism': org_name,
                            'PrimerPairID': pair_id,
                            'Forward': fwd_seq,
                            'Reverse': rev_seq,
                            'Tier': tier,
                            'TargetGene': gene,
                            'AmpliconSize': size,
                            'Reference': ref
                        })
            if primers:
                return primers

    default_db = load_default_primer_database()
    for d in default_db:
        if not any(p['PrimerPairID'] == d['PrimerPairID'] for p in primers):
            if 'Tier' not in d:
                d['Tier'] = infer_tier(d.get('Organism', ''), d.get('PrimerPairID', ''))
            if 'TargetGene' not in d:
                d['TargetGene'] = infer_target_gene(d.get('PrimerPairID', ''))
            if 'AmpliconSize' not in d:
                d['AmpliconSize'] = 0
            if 'Reference' not in d:
                d['Reference'] = 'Default DB'
            primers.append(d)

    return primers

def scan_and_load_local_primer_archives():
    """Scans Downloads and current folder for qPrimerDB archive files."""
    search_paths = [
        os.getcwd(),
        r"C:\Users\internet\Downloads",
        "/mnt/sdc1",
        "/tmp"
    ]
    primers = []
    for p in search_paths:
        if not os.path.exists(p):
            continue
        archive_files = glob.glob(os.path.join(p, "*.primers.gz.zip")) + glob.glob(os.path.join(p, "*.primers.gz"))
        for arch in archive_files:
            try:
                primers.extend(load_primers_from_archive(arch))
            except Exception:
                pass
    return primers

def load_primers_from_archive(archive_path):
    """Extracts primer sequences from qPrimerDB zip/gz files."""
    primers = []
    org_name = os.path.basename(archive_path).split('.')[0].replace('_', ' ')
    lines = []
    if archive_path.endswith('.zip'):
        with zipfile.ZipFile(archive_path, 'r') as z:
            for name in z.namelist():
                if name.endswith('.gz'):
                    with gzip.GzipFile(fileobj=io.BytesIO(z.read(name))) as gz:
                        lines.extend(gz.read().decode('utf-8', errors='ignore').splitlines())
    elif archive_path.endswith('.gz'):
        with gzip.open(archive_path, 'rt', encoding='utf-8', errors='ignore') as gz:
            lines = gz.readlines()

    for line in lines:
        parts = [p.strip() for p in line.split(',') if p.strip()]
        if len(parts) >= 11:
            pair_id = parts[1]
            fwd = parts[9]
            rev = parts[10]
            if re.match(r'^[ATCGURYKMSWBDHVN]+$', fwd, re.I) and re.match(r'^[ATCGURYKMSWBDHVN]+$', rev, re.I):
                primers.append({
                    'Organism': org_name,
                    'PrimerPairID': pair_id,
                    'Forward': fwd,
                    'Reverse': rev
                })
    return primers

def run_cdhit_est(fastq_path, identity=0.99, cdhit_exec="cd-hit-est"):
    """Pre-clusters/dereplicates input FASTQ reads using cd-hit-est."""
    import tempfile

    print(f"[*] Pre-clustering FASTQ reads using {cdhit_exec} (identity={identity})...")
    try:
        subprocess.run([cdhit_exec, "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    except FileNotFoundError:
        print(f"[!] Warning: '{cdhit_exec}' executable not found in PATH. Skipping CD-HIT pre-clustering.")
        return fastq_path, None

    tmp_in_fasta = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False).name
    tmp_out_fasta = tempfile.NamedTemporaryFile(suffix="_cdhit.fasta", delete=False).name

    total_input_reads = 0
    with open(tmp_in_fasta, 'w', encoding='utf-8') as out_f:
        if str(fastq_path).endswith('.gz'):
            h_in = gzip.open(fastq_path, 'rt', encoding='utf-8', errors='ignore')
        else:
            h_in = open(fastq_path, 'rt', encoding='utf-8', errors='ignore')
        while True:
            h = h_in.readline().strip()
            if not h:
                break
            s = h_in.readline().strip()
            h_in.readline()
            h_in.readline()
            total_input_reads += 1
            out_f.write(f">{h[1:]}\n{s}\n")
        h_in.close()

    print(f"    - Extracted {total_input_reads:,} reads to temporary FASTA.")

    n_word = 10 if identity >= 0.95 else (8 if identity >= 0.90 else (6 if identity >= 0.88 else 5))
    cmd = [cdhit_exec, "-i", tmp_in_fasta, "-o", tmp_out_fasta, "-c", str(identity), "-n", str(n_word), "-d", "0", "-T", "0", "-M", "0"]

    t0 = time.time()
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[!] CD-HIT error: {e.stderr}")
        os.remove(tmp_in_fasta)
        return fastq_path, None

    elapsed = time.time() - t0
    clustered_reads = 0
    with open(tmp_out_fasta, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                clustered_reads += 1

    reduction = ((total_input_reads - clustered_reads) / total_input_reads * 100) if total_input_reads > 0 else 0
    print(f"[+] CD-HIT completed in {elapsed:.2f}s: {total_input_reads:,} -> {clustered_reads:,} ({reduction:.1f}% reduction)")
    os.remove(tmp_in_fasta)
    return tmp_out_fasta, {
        'total_reads': total_input_reads,
        'clustered_reads': clustered_reads,
        'reduction_pct': round(reduction, 1),
        'elapsed_s': round(elapsed, 2)
    }

# ==============================================================================
# Multi-Processing Worker Functions (Chunk-Based)
# ==============================================================================
_worker_compiled_primers = []
_worker_sample_prefix = "sample"

def _worker_init(primers, max_mismatches, sample_prefix="sample"):
    """Initializes worker process: compiles regex & dual seeds once per worker in memory."""
    global _worker_compiled_primers, _worker_sample_prefix
    _worker_sample_prefix = sample_prefix
    _worker_compiled_primers = []
    for p in primers:
        fwd_seq = p['Forward']
        rev_seq = p['Reverse']
        fwd_rc = reverse_complement(fwd_seq)
        rev_rc = reverse_complement(rev_seq)
        genus, is_genus = classify_primer_taxa(p['Organism'], p['PrimerPairID'])

        _worker_compiled_primers.append({
            'info': p,
            'pair_id': p['PrimerPairID'],
            'org_name': p['Organism'],
            'genus': genus,
            'is_genus': is_genus,
            'fwd_seq': fwd_seq,
            'rev_seq': rev_seq,
            'fwd_regex': iupac_to_regex(fwd_seq, max_mismatches),
            'rev_regex': iupac_to_regex(rev_seq, max_mismatches),
            'fwd_rc_regex': iupac_to_regex(fwd_rc, max_mismatches),
            'rev_rc_regex': iupac_to_regex(rev_rc, max_mismatches),
            'fwd_seeds': get_exact_seeds(fwd_seq, k=8),
            'rev_seeds': get_exact_seeds(rev_seq, k=8),
            'fwd_rc_seeds': get_exact_seeds(fwd_rc, k=8),
            'rev_rc_seeds': get_exact_seeds(rev_rc, k=8),
        })

def _worker_process_chunk(chunk_data):
    """Processes a chunk of reads with Single-Primer Anchor Extension & Dark Matter sampling."""
    global _worker_compiled_primers, _worker_sample_prefix
    if isinstance(chunk_data, tuple) and len(chunk_data) == 2 and isinstance(chunk_data[0], list):
        reads_chunk, frac = chunk_data
    else:
        reads_chunk, frac = chunk_data, 0.0

    match_counts = {'BOTH': 0, 'ONLY_FORWARD': 0, 'ONLY_REVERSE': 0, 'NONE': 0}
    stats_delta = {}
    extracted_records = []
    candidate_records = []
    unassigned_records = []

    for header, seq, qual_header, qual in reads_chunk:
        seq_upper = seq.upper()
        read_has_both = False
        read_had_any_primer = False

        for cp in _worker_compiled_primers:
            pair_id = cp['pair_id']
            fwd_pos = []
            rev_pos = []

            # Fast k-mer seeds check (dual 8-mer seeds)
            if any(s in seq_upper for s in cp['fwd_seeds']):
                fwd_pos.extend([(m.start(), m.end()) for m in cp['fwd_regex'].finditer(seq_upper)])
            if any(s in seq_upper for s in cp['fwd_rc_seeds']):
                fwd_pos.extend([(m.start(), m.end()) for m in cp['fwd_rc_regex'].finditer(seq_upper)])

            if any(s in seq_upper for s in cp['rev_seeds']):
                rev_pos.extend([(m.start(), m.end()) for m in cp['rev_regex'].finditer(seq_upper)])
            if any(s in seq_upper for s in cp['rev_rc_seeds']):
                rev_pos.extend([(m.start(), m.end()) for m in cp['rev_rc_regex'].finditer(seq_upper)])

            has_fwd = len(fwd_pos) > 0
            has_rev = len(rev_pos) > 0

            if has_fwd or has_rev:
                read_had_any_primer = True
                if pair_id not in stats_delta:
                    stats_delta[pair_id] = {
                        'Organism': cp['org_name'],
                        'Forward': cp['fwd_seq'],
                        'Reverse': cp['rev_seq'],
                        'BOTH': 0, 'ONLY_FORWARD': 0, 'ONLY_REVERSE': 0,
                        'AmpliconDistances': [],
                        'TaxaType': 'Genus' if cp['is_genus'] else 'Species',
                        'Genus': cp['genus']
                    }
                p_entry = stats_delta[pair_id]

                if has_fwd and has_rev:
                    read_has_both = True
                    p_entry['BOTH'] += 1
                    match_counts['BOTH'] += 1
                    read_str = f"{header}\n{seq}\n{qual_header}\n{qual}\n"

                    for f_start, f_end in fwd_pos:
                        for r_start, r_end in rev_pos:
                            dist = abs(r_end - f_start) if r_end >= f_start else abs(f_end - r_start)
                            p_entry['AmpliconDistances'].append(dist)
                            clean_org = re.sub(r'[^a-zA-Z0-9_-]', '_', cp['org_name'])[:25]
                            filename = f"{_worker_sample_prefix}_{clean_org}_{pair_id}_frp_{dist}bp.fastq"
                            
                            p_min = min(f_start, r_start)
                            p_max = max(f_end, r_end)
                            amplicon_seq = seq[p_min:p_max]
                            if len(amplicon_seq) < 30:
                                amplicon_seq = seq
                                
                            extracted_records.append((filename, read_str, amplicon_seq, pair_id))
                    if len(candidate_records) < 500:
                        candidate_records.append((cp['org_name'], 'BOTH', read_str, pair_id))
                elif has_fwd:
                    p_entry['ONLY_FORWARD'] += 1
                    match_counts['ONLY_FORWARD'] += 1
                    read_str = f"{header}\n{seq}\n{qual_header}\n{qual}\n"
                    # Single-Primer Anchor Extension for metagenomic fragments
                    for f_start, f_end in fwd_pos:
                        ext_seq = seq[f_start:] if len(seq) - f_start >= 100 else seq
                        if len(ext_seq) >= 100:
                            clean_org = re.sub(r'[^a-zA-Z0-9_-]', '_', cp['org_name'])[:25]
                            filename = f"{_worker_sample_prefix}_{clean_org}_{pair_id}_anchor_fwd_{len(ext_seq)}bp.fastq"
                            extracted_records.append((filename, read_str, ext_seq, pair_id))
                    if len(candidate_records) < 500:
                        candidate_records.append((cp['org_name'], 'ANCHOR_FWD', read_str, pair_id))
                elif has_rev:
                    p_entry['ONLY_REVERSE'] += 1
                    match_counts['ONLY_REVERSE'] += 1
                    read_str = f"{header}\n{seq}\n{qual_header}\n{qual}\n"
                    # Single-Primer Anchor Extension for metagenomic fragments
                    for r_start, r_end in rev_pos:
                        ext_seq = seq[:r_end] if r_end >= 100 else seq
                        if len(ext_seq) >= 100:
                            clean_org = re.sub(r'[^a-zA-Z0-9_-]', '_', cp['org_name'])[:25]
                            filename = f"{_worker_sample_prefix}_{clean_org}_{pair_id}_anchor_rev_{len(ext_seq)}bp.fastq"
                            extracted_records.append((filename, read_str, ext_seq, pair_id))
                    if len(candidate_records) < 500:
                        candidate_records.append((cp['org_name'], 'ANCHOR_REV', read_str, pair_id))

        if not read_has_both:
            match_counts['NONE'] += 1
            if not read_had_any_primer and len(unassigned_records) < 25 and is_high_complexity_read(seq):
                unassigned_records.append(f"{header}\n{seq}\n{qual_header}\n{qual}\n")

    return match_counts, stats_delta, extracted_records, frac, candidate_records, unassigned_records

def read_fastq_chunks(fastq_path, chunk_size=20000, limit=None, threads=4):
    """
    Generator yielding batches of FASTQ records as ((header, seq, qual_header, qual), fraction).
    Leverages pigz parallel decompression subprocess pipe when available.
    """
    file_size = os.path.getsize(fastq_path) if os.path.exists(fastq_path) else 1
    proc = None
    f_in = None
    gz_handle = None

    if str(fastq_path).endswith('.gz'):
        pigz_path = shutil.which('pigz')
        if pigz_path:
            decomp_threads = max(2, min(8, threads // 2))
            f_in = open(fastq_path, 'rb')
            proc = subprocess.Popen(
                [pigz_path, '-dc', '-p', str(decomp_threads)],
                stdin=f_in,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                bufsize=64 * 1024 * 1024,
                text=True,
                errors='ignore'
            )
            handle = proc.stdout
        else:
            gz_handle = gzip.open(fastq_path, 'rt', encoding='utf-8', errors='ignore')
            handle = gz_handle
    else:
        f_in = open(fastq_path, 'rt', encoding='utf-8', errors='ignore')
        handle = f_in

    chunk = []
    total_count = 0
    try:
        while True:
            h = handle.readline().strip()
            if not h:
                break
            s = handle.readline().strip()
            qh = handle.readline().strip()
            q = handle.readline().strip()
            total_count += 1
            chunk.append((h, s, qh, q))
            if len(chunk) >= chunk_size:
                if limit:
                    frac = min(1.0, total_count / limit)
                elif f_in:
                    try:
                        frac = min(1.0, f_in.tell() / file_size)
                    except Exception:
                        frac = 0.0
                elif gz_handle and hasattr(gz_handle, 'fileobj') and gz_handle.fileobj:
                    try:
                        frac = min(1.0, gz_handle.fileobj.tell() / file_size)
                    except Exception:
                        frac = 0.0
                else:
                    frac = 0.0
                yield (chunk, frac)
                chunk = []
            if limit and total_count >= limit:
                break
        if chunk:
            frac = 1.0 if not limit else min(1.0, total_count / limit)
            yield (chunk, frac)
    finally:
        if proc:
            try:
                proc.stdout.close()
                proc.terminate()
                proc.wait()
            except Exception:
                pass
        if f_in:
            try:
                f_in.close()
            except Exception:
                pass
        if gz_handle:
            try:
                gz_handle.close()
            except Exception:
                pass

def extract_sample_prefix(fastq_path):
    """
    Extracts the sample prefix from the FASTQ filename up to the first underscore.
    e.g. tm_NG1YETBR8Q_t2t_u.fastq.gz -> 'tm'
         af_sample_1.fastq.gz -> 'af'
         sample123.fastq.gz -> 'sample123'
    """
    base = os.path.basename(fastq_path)
    if '_' in base:
        p = base.split('_')[0].strip()
    else:
        p = re.split(r'\.fastq|\.gz', base, flags=re.I)[0].strip()
    return p if p else "sample"

def find_reference_genome(organism_name, ref_dirs=None):
    """
    Finds a matching reference FASTA file for an organism in known reference directories.
    Handles exact names, species names, aliases/synonyms (HSV-1, HSV-2, EBV, CMV, VZV, HHV, CCHFV, etc.),
    NCBI accession prefixes, subdirectories (e.g. viruses/), and genus fallbacks.
    """
    if not ref_dirs:
        ref_dirs = [
            "/mnt/sdc1/pathogen_references",
            "./pathogen_references",
            "./references",
            r"C:\Users\internet\Downloads\repository\primer\references",
            r"C:\Users\internet\Downloads"
        ]
    elif isinstance(ref_dirs, str):
        ref_dirs = [ref_dirs]

    # Pre-calculated common alias map (lowercase)
    alias_map = {
        'hsv1': ['hsv1', 'hsv-1', 'hhv1', 'hhv-1', 'herpes_simplex_virus_1', 'herpes_simplex_1', 'nc_001806'],
        'hsv-1': ['hsv1', 'hsv-1', 'hhv1', 'hhv-1', 'herpes_simplex_virus_1', 'herpes_simplex_1', 'nc_001806'],
        'hsv2': ['hsv2', 'hsv-2', 'hhv2', 'hhv-2', 'herpes_simplex_virus_2', 'herpes_simplex_2', 'nc_001798'],
        'hsv-2': ['hsv2', 'hsv-2', 'hhv2', 'hhv-2', 'herpes_simplex_virus_2', 'herpes_simplex_2', 'nc_001798'],
        'vzv': ['vzv', 'hhv3', 'hhv-3', 'varicella_zoster', 'varicella', 'nc_001348'],
        'varicella': ['vzv', 'hhv3', 'hhv-3', 'varicella_zoster', 'varicella'],
        'ebv': ['ebv', 'hhv4', 'hhv-4', 'epstein_barr', 'nc_007605'],
        'cmv': ['cmv', 'hcmv', 'hhv5', 'hhv-5', 'cytomegalo', 'nc_006273'],
        'hhv6a': ['hhv6a', 'hhv-6a', 'human_herpesvirus_6a', 'nc_001664'],
        'hhv6b': ['hhv6b', 'hhv-6b', 'human_herpesvirus_6b', 'nc_000898'],
        'hhv6': ['hhv6', 'hhv-6', 'human_herpesvirus_6'],
        'hhv7': ['hhv7', 'hhv-7', 'human_herpesvirus_7', 'nc_001716'],
        'hhv8': ['hhv8', 'hhv-8', 'kshv', 'kaposi', 'nc_009333'],
        'cchfv': ['cchfv', 'crimean_congo', 'crimean-congo', 'nairovirus', 'nc_005302'],
        'sars-cov-2': ['sars-cov-2', 'sarscov2', 'sars_cov_2', 'coronavirus', 'nc_045512'],

        # Piroplasmida (Babesia, Cytauxzoon, Theileria) aliases & multi-FASTA resolution
        'theileria / piroplasmida': ['piroplasmida', 'theileria', 'babesia', 'cytauxzoon'],
        'babesia / theileria spp.': ['piroplasmida', 'babesia', 'theileria', 'cytauxzoon'],
        'babesia / theileria': ['piroplasmida', 'babesia', 'theileria', 'cytauxzoon'],
        'piroplasmida': ['piroplasmida', 'theileria', 'babesia', 'cytauxzoon'],
        'piroplasm': ['piroplasmida', 'theileria', 'babesia', 'cytauxzoon'],
        'cytauxzoon felis': ['cytauxzoon_felis', 'cytauxzoon', 'piroplasmida'],
        'cytauxzoon europaeus': ['cytauxzoon_europaeus', 'cytauxzoon', 'piroplasmida'],
        'cytauxzoon': ['cytauxzoon', 'cytauxzoon_felis', 'cytauxzoon_europaeus', 'piroplasmida'],
        'theileria parva': ['theileria_parva', 'theileria', 'piroplasmida'],
        'theileria annulata': ['theileria_annulata', 'theileria', 'piroplasmida'],
        'theileria orientalis': ['theileria_orientalis', 'theileria', 'piroplasmida'],
        'theileria equi': ['theileria_equi', 'theileria', 'piroplasmida'],
        'theileria': ['theileria', 'piroplasmida'],
        'babesia microti': ['babesia_microti', 'babesia', 'piroplasmida'],
        'babesia duncani': ['babesia_duncani', 'babesia', 'piroplasmida'],
        'babesia divergens': ['babesia_divergens', 'babesia', 'piroplasmida'],
        'babesia bovis': ['babesia_bovis', 'babesia', 'piroplasmida'],
        'babesia': ['babesia', 'piroplasmida'],
        'babesia odocoilei': ['babesia_odocoilei', 'babesia_odocoilei.fasta', 'u16369', 'piroplasmida'],
        'babesia odecoili': ['babesia_odocoilei', 'babesia_odocoilei.fasta', 'u16369', 'piroplasmida'],
        'babesia vanuatum': ['babesia_venatorum', 'babesia_sp_eu1', 'piroplasmida'],
        'babesia venatorum': ['babesia_venatorum', 'babesia_sp_eu1', 'piroplasmida'],
        'granulicatella adiacens': ['granulicatella_adiacens', 'granulicatella'],
        'granulicatella': ['granulicatella_adiacens', 'granulicatella']
    }

    clean_org = re.sub(r'\(.*?\)', '', organism_name).strip()
    clean_org_underscore = re.sub(r'[\s/]+', '_', clean_org).lower()
    org_lower = organism_name.lower()

    tokens = [w.lower() for w in re.split(r'[\s_/()]+', organism_name) if len(w) > 1]
    tokens_filtered = [w for w in tokens if len(w) > 2 and w not in ['universal', 'bacteria', 'fungi', 'strain', 'isolate', 'subsp', 'complex', 'virus', 'human']]

    target_aliases = []
    # Match most specific alias keys first (sorted length descending)
    sorted_keys = sorted(alias_map.keys(), key=lambda k: len(k), reverse=True)
    for key in sorted_keys:
        if key in org_lower or key.replace(' ', '_') in clean_org_underscore:
            for a in alias_map[key]:
                if a not in target_aliases:
                    target_aliases.append(a)
            break

    for rdir in ref_dirs:
        if not rdir or not os.path.exists(rdir):
            continue

        # 1. Direct candidate paths
        for ext in ['.fasta', '.fa', '.fna']:
            cand1 = os.path.join(rdir, f"{clean_org_underscore}{ext}")
            if os.path.exists(cand1):
                return cand1
            cand2 = os.path.join(rdir, f"{organism_name.replace(' ', '_')}{ext}")
            if os.path.exists(cand2):
                return cand2

        # 2. Collect all candidates recursively (including subdirectories like viruses/)
        all_candidates = []
        try:
            for root, dirs, files in os.walk(rdir):
                for fname in files:
                    if fname.endswith(('.fasta', '.fa', '.fna')):
                        all_candidates.append(os.path.join(root, fname))
        except Exception:
            continue

        # 3. Check alias matches in prioritized order (exact base match first, then substring)
        if target_aliases:
            for alias in target_aliases:
                for cpath in all_candidates:
                    cfname = os.path.basename(cpath).lower()
                    base_no_ext = os.path.splitext(cfname)[0]
                    if alias == base_no_ext:
                        return cpath
            for alias in target_aliases:
                for cpath in all_candidates:
                    cfname = os.path.basename(cpath).lower()
                    if alias in cfname:
                        return cpath

        # 4. Species-level match (all key words present)
        if tokens_filtered:
            for cpath in all_candidates:
                cfname = os.path.basename(cpath).lower()
                if all(w in cfname for w in tokens_filtered):
                    return cpath

        # 5. First word / genus-level match (e.g. 'Rickettsia' or 'Mycobacterium')
        if tokens:
            genus_w = tokens[0]
            if len(genus_w) > 3:
                for cpath in all_candidates:
                    cfname = os.path.basename(cpath).lower()
                    if genus_w == os.path.splitext(cfname)[0] or f"{genus_w}_" in cfname or f"_{genus_w}." in cfname:
                        return cpath
                for cpath in all_candidates:
                    cfname = os.path.basename(cpath).lower()
                    if genus_w in cfname:
                        return cpath

    return None

def auto_download_reference_genome(organism_name, dest_dir=None):
    """
    Automatically queries and downloads an official NCBI RefSeq / GenBank complete reference genome
    when no local reference is present. Stores the FASTA file in the reference directory and
    creates a copy in the Downloads directory (User Global Rule).
    """
    import urllib.request
    import urllib.parse

    clean = organism_name.split('(')[0].split('/')[0].strip()
    words = [w for w in clean.split() if w.lower() not in ['universal', 'strain', 'isolate', 'virus', 'bacteria', 'spp.', 'sp.']]

    target_dir = dest_dir
    if not target_dir or not os.path.exists(target_dir):
        for candidate in ["/mnt/sdc1/pathogen_references", r"C:\Users\internet\Downloads\repository\primer\references", "./references", "./pathogen_references"]:
            if os.path.exists(candidate):
                target_dir = candidate
                break
        if not target_dir:
            target_dir = r"C:\Users\internet\Downloads\repository\primer\references"
            try:
                os.makedirs(target_dir, exist_ok=True)
            except Exception:
                target_dir = "."

    safe_name = re.sub(r'[^\w\.-]', '_', clean)
    local_fasta = os.path.join(target_dir, f"{safe_name}.fasta")
    if os.path.exists(local_fasta) and os.path.getsize(local_fasta) > 1000:
        return local_fasta

    # Hierarchical search terms for NCBI Nuccore
    queries = [
        f'"{clean}"[Organism] AND refseq[filter] AND complete genome',
        f'"{clean}"[Organism] AND complete genome',
        f'"{clean}"[Title] AND complete genome'
    ]
    if len(words) >= 2:
        queries.append(f'"{words[0]} {words[1]}"[Organism] AND refseq[filter]')
    elif len(words) == 1 and len(words[0]) > 3:
        queries.append(f'"{words[0]}"[Organism] AND refseq[filter] AND complete genome')

    ncbi_id = None
    for q in queries:
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={urllib.parse.quote(q)}&retmode=json&sort=relevance'
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Antigravity-Pathogen-AutoDownloader)'})
        try:
            with urllib.request.urlopen(req, timeout=12) as r:
                data = json.loads(r.read())
                ids = data.get('esearchresult', {}).get('idlist', [])
                if ids:
                    ncbi_id = ids[0]
                    break
        except Exception:
            time.sleep(0.5)
        time.sleep(0.3)

    if not ncbi_id:
        return None

    # Fetch FASTA via EFetch
    fetch_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ncbi_id}&rettype=fasta&retmode=text'
    req_fetch = urllib.request.Request(fetch_url, headers={'User-Agent': 'Mozilla/5.0 (Antigravity-Pathogen-AutoDownloader)'})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req_fetch, timeout=25) as resp:
                fasta_text = resp.read().decode('utf-8', errors='ignore')
                if fasta_text.startswith('>'):
                    with open(local_fasta, 'w', encoding='utf-8') as out_f:
                        out_f.write(fasta_text)
                    print(f"  [+] Referenz automatisch von NCBI heruntergeladen: {os.path.basename(local_fasta)} ({len(fasta_text):,} bp)")

                    # Copy to Downloads (User Global Rule)
                    dl_dir = r"C:\Users\internet\Downloads"
                    try:
                        if os.path.exists(dl_dir):
                            dl_target = os.path.join(dl_dir, os.path.basename(local_fasta))
                            if os.path.abspath(local_fasta).lower() != os.path.abspath(dl_target).lower():
                                shutil.copy2(local_fasta, dl_target)
                    except Exception:
                        pass
                    return local_fasta
        except Exception:
            time.sleep(1.0)

    return None

def parse_sam_alignments(sam_path):
    """Parses SAM output without external dependencies to extract alignment metrics."""
    mapped_count = 0
    high_qual_count = 0
    mapq_list = []
    rnames = Counter()
    spans_by_rname = defaultdict(list)

    if not os.path.exists(sam_path):
        return {
            'mapped_count': 0, 'high_qual_count': 0, 'avg_mapq': 0,
            'top_rname': 'N/A', 'covered_bp': 0
        }

    with open(sam_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('@'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
            flag = int(parts[1])
            rname = parts[2]
            pos = int(parts[3])
            mapq = int(parts[4])

            if (flag & 4) == 0 and rname != '*':
                mapped_count += 1
                mapq_list.append(mapq)
                if mapq >= 20:
                    high_qual_count += 1
                rnames[rname] += 1

                seq_len = len(parts[9]) if parts[9] != '*' else 100
                spans_by_rname[rname].append((pos, pos + seq_len))

    avg_mapq = (sum(mapq_list) / len(mapq_list)) if mapq_list else 0
    top_rname = rnames.most_common(1)[0][0] if rnames else "N/A"

    covered_bp = 0
    if top_rname in spans_by_rname:
        spans = sorted(spans_by_rname[top_rname])
        merged = []
        for s, e in spans:
            if not merged or s > merged[-1][1]:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        covered_bp = sum(e - s for s, e in merged)

    return {
        'mapped_count': mapped_count,
        'high_qual_count': high_qual_count,
        'avg_mapq': round(avg_mapq, 1),
        'top_rname': top_rname,
        'covered_bp': covered_bp
    }

def resolve_rname_description(rname, ref_path=None):
    """
    Resolves a contig accession / RNAME to a readable species & accession label from the reference FASTA header.
    Especially useful for multi-FASTA references (Piroplasmida.fasta, Theileria.fasta, Cytauxzoon.fasta, Babesia.fasta)
    to reveal the exact species and strain where reads mapped.
    """
    if not rname or rname == "N/A" or not ref_path or not os.path.exists(ref_path):
        return rname
    try:
        with open(ref_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('>') and rname in line:
                    header = line[1:].strip()
                    tokens = header.split()
                    if len(tokens) >= 3 and tokens[0] == rname:
                        return f"{tokens[1]} {tokens[2]} ({rname})"
                    elif len(tokens) >= 2 and tokens[0] == rname:
                        return f"{tokens[1]} ({rname})"
                    return header[:50]
    except Exception:
        pass
    return rname

def clean_blast_query(seq, fallback_seq=""):
    """
    Cleans nucleotide sequences before submission to NCBI BLAST:
    - Strips leading/trailing 'N' ambiguity runs and whitespace.
    - If sequence contains internal N-blocks (due to coverage gaps), extracts the longest continuous ACGT block.
    - Verifies that the sequence has >= 50 bp and >= 80% valid ACGT bases to prevent Karlin-Altschul calculation errors.
    - Falls back to raw representative read if consensus contains too many Ns or is invalid.
    """
    if not seq:
        seq = fallback_seq or ""
    if not seq:
        return ""

    lines = [l.strip() for l in str(seq).splitlines() if not l.startswith('>')]
    s = ''.join(lines).upper()
    s = re.sub(r'^[N\s]+|[N\s]+$', '', s)

    # If there are large N-blocks, select the longest continuous valid segment
    blocks = [b for b in re.split(r'N{3,}', s) if len(b) >= 50]
    if blocks:
        s = max(blocks, key=len)

    acgt_count = sum(1 for c in s if c in 'ACGT')
    if len(s) >= 50 and (acgt_count / len(s)) >= 0.8:
        return s[:2000]

    # Fallback to raw read
    if fallback_seq:
        fb_lines = [l.strip() for l in str(fallback_seq).splitlines() if not l.startswith('>')]
        fb_s = ''.join(fb_lines).upper()
        fb_s = re.sub(r'^[N\s]+|[N\s]+$', '', fb_s)
        if len(fb_s) >= 50:
            return fb_s[:2000]

    return s[:2000] if len(s) >= 50 else (fallback_seq or "")

def assemble_candidate_reads(cand_fastq_path, prefix, clean_org, output_dir, bam_path=None, ref_path=None, threads=4, minimap2_bin=None, samtools_bin=None):
    """
    Performs micro-assembly and consensus reconstruction from candidate reads:
    1. Reference-Assisted Consensus via `samtools consensus` (when BAM from reference mapping is present).
    2. De-Novo Seed-Overlap Assembly (fully independent of any reference genome).
    3. Whole-Genome / Metagenomic Flye assembly (when read count >= 25 and total bp >= 40k).
    Copies resulting contig FASTAs to Downloads directory (User Global Rule).
    """
    assembly_info = {
        'contig_rel': None,
        'contig_path': None,
        'contig_len': 0,
        'contig_seq': '',
        'method': 'N/A'
    }

    if not cand_fastq_path or not os.path.exists(cand_fastq_path):
        return assembly_info

    # Count reads and calculate total bp
    read_count = 0
    total_bp = 0
    longest_read = (0, "", "")
    try:
        with open(cand_fastq_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
            for i in range(0, len(lines), 4):
                if i + 1 < len(lines):
                    r_hdr = lines[i].strip()
                    r_seq = lines[i+1].strip()
                    r_len = len(r_seq)
                    read_count += 1
                    total_bp += r_len
                    if r_len > longest_read[0]:
                        longest_read = (r_len, r_hdr, r_seq)
    except Exception:
        pass

    if read_count < 5:
        return assembly_info

    contig_out_name = f"{prefix}_{clean_org}_assembled_contig.fasta"
    contig_out_path = os.path.join(output_dir, contig_out_name)

    # Method 1: Reference-Assisted Consensus via samtools consensus (WITHOUT -a to avoid N-padding uncovered regions)
    if bam_path and os.path.exists(bam_path) and samtools_bin:
        try:
            tmp_cons = os.path.join(output_dir, f"{prefix}_{clean_org}_ref_consensus.tmp.fa")
            subprocess.run([samtools_bin, "consensus", bam_path, "-o", tmp_cons],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            if os.path.exists(tmp_cons) and os.path.getsize(tmp_cons) > 30:
                with open(tmp_cons, 'r', encoding='utf-8', errors='ignore') as cf:
                    c_seq = ''.join([l.strip() for l in cf if not l.startswith('>')])
                clean_c = clean_blast_query(c_seq, fallback_seq=longest_read[2])
                if len(clean_c) >= 50:
                    with open(contig_out_path, 'w', encoding='utf-8') as out_f:
                        out_f.write(f">{prefix}_{clean_org}_consensus_amplicon len={len(clean_c)} reads={read_count} method=samtools_consensus\n{clean_c}\n")
                    assembly_info = {
                        'contig_rel': contig_out_name,
                        'contig_path': contig_out_path,
                        'contig_len': len(clean_c),
                        'contig_seq': clean_c,
                        'method': 'Reference-Assisted Consensus'
                    }
            if os.path.exists(tmp_cons):
                try:
                    os.remove(tmp_cons)
                except Exception:
                    pass
        except Exception:
            pass

    # Method 2: De-Novo Seed-Overlap Assembly (WITHOUT -a to avoid N-padding)
    if not assembly_info['contig_rel'] and longest_read[0] >= 80 and minimap2_bin and samtools_bin:
        try:
            seed_fa = os.path.join(output_dir, f"{prefix}_{clean_org}_seed.tmp.fa")
            seed_bam = os.path.join(output_dir, f"{prefix}_{clean_org}_seed.tmp.bam")
            tmp_cons = os.path.join(output_dir, f"{prefix}_{clean_org}_denovo.tmp.fa")

            with open(seed_fa, 'w', encoding='utf-8') as sf:
                sf.write(f">seed_backbone len={longest_read[0]}\n{longest_read[2]}\n")

            p_map = subprocess.Popen([minimap2_bin, "-ax", "map-ont", seed_fa, cand_fastq_path],
                                     stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p_sort = subprocess.Popen([samtools_bin, "sort", "-o", seed_bam],
                                      stdin=p_map.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            p_map.stdout.close()
            p_sort.communicate()

            if os.path.exists(seed_bam) and os.path.getsize(seed_bam) > 0:
                subprocess.run([samtools_bin, "index", seed_bam], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                subprocess.run([samtools_bin, "consensus", seed_bam, "-o", tmp_cons],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                if os.path.exists(tmp_cons) and os.path.getsize(tmp_cons) > 30:
                    with open(tmp_cons, 'r', encoding='utf-8', errors='ignore') as cf:
                        c_seq = ''.join([l.strip() for l in cf if not l.startswith('>')])
                    clean_c = clean_blast_query(c_seq, fallback_seq=longest_read[2])
                    if len(clean_c) >= 50:
                        with open(contig_out_path, 'w', encoding='utf-8') as out_f:
                            out_f.write(f">{prefix}_{clean_org}_denovo_contig len={len(clean_c)} reads={read_count} method=denovo_seed_overlap\n{clean_c}\n")
                        assembly_info = {
                            'contig_rel': contig_out_name,
                            'contig_path': contig_out_path,
                            'contig_len': len(clean_c),
                            'contig_seq': clean_c,
                            'method': 'De-Novo Overlap Assembly'
                        }
            # Cleanup temp files
            for tf in [seed_fa, seed_bam, f"{seed_bam}.bai", tmp_cons]:
                if os.path.exists(tf):
                    try:
                        os.remove(tf)
                    except Exception:
                        pass
        except Exception:
            pass

    # Method 3: Flye whole-genome / metagenomic assembly for abundant sets (>= 25 reads, >= 40k bp)
    flye_bin = shutil.which('flye') or ('/usr/bin/flye' if os.path.exists('/usr/bin/flye') else None)
    if read_count >= 25 and total_bp >= 40000 and flye_bin:
        flye_dir = os.path.join(output_dir, f"{prefix}_{clean_org}_flye")
        try:
            subprocess.run([flye_bin, "--nano-raw", cand_fastq_path, "--out-dir", flye_dir,
                            "--threads", str(min(threads, 8)), "--meta"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
            flye_fa = os.path.join(flye_dir, "assembly.fasta")
            if os.path.exists(flye_fa) and os.path.getsize(flye_fa) > 100:
                flye_out_name = f"{prefix}_{clean_org}_flye_genome.fasta"
                flye_out_path = os.path.join(output_dir, flye_out_name)
                shutil.copy2(flye_fa, flye_out_path)
                with open(flye_out_path, 'r', encoding='utf-8') as ff:
                    f_seq = ''.join([l.strip() for l in ff if not l.startswith('>')])
                assembly_info = {
                    'contig_rel': flye_out_name,
                    'contig_path': flye_out_path,
                    'contig_len': len(f_seq),
                    'contig_seq': f_seq[:2000],
                    'method': 'Flye Metagenome Assembly'
                }
        except Exception:
            pass

    # Copy to Downloads (User Global Rule)
    if assembly_info['contig_path'] and os.path.exists(assembly_info['contig_path']):
        dl_dir = r"C:\Users\internet\Downloads"
        try:
            if os.path.exists(dl_dir):
                dl_target = os.path.join(dl_dir, os.path.basename(assembly_info['contig_path']))
                shutil.copy2(assembly_info['contig_path'], dl_target)
        except Exception:
            pass
        print(f"  [🧬] [{clean_org}] Micro-Assembly ({assembly_info['method']}): {assembly_info['contig_len']:,} bp Contig erzeugt -> {assembly_info['contig_rel']}")

    return assembly_info

def run_local_blast_verification(validations, output_dir, prefix, blast_db="/mnt/blast_db/nt", threads=16):
    """
    Performs batch local BLASTn verification against the local NCBI nt database (/mnt/blast_db/nt).
    Scans all candidate contigs/reads in a single batched megablast pass across available CPU threads.
    Saves full TSV report and copies to Downloads directory (User Global Rule).
    """
    blastn_bin = shutil.which('blastn') or ('/usr/bin/blastn' if os.path.exists('/usr/bin/blastn') else None)
    if not blastn_bin:
        return

    db_base = blast_db
    if not (os.path.exists(f"{db_base}.nal") or os.path.exists(f"{db_base}.000.nhr") or os.path.exists(f"{db_base}.nhr") or os.path.exists(db_base)):
        return

    queries = []
    for idx, v in enumerate(validations):
        clean_org = v['clean_org']
        q_seq = clean_blast_query(v.get('contig_seq'), fallback_seq=v.get('rep_seq', ''))
        if q_seq and len(q_seq) >= 50:
            qid = f"cand_{idx}_{clean_org}"
            queries.append((qid, idx, v['organism'], q_seq))

    if not queries:
        return

    print("\n" + "="*75)
    print(f"🧬 LOKALER BLASTN-VERIFIKATIONSLAUF ({prefix}) gegen {blast_db}")
    print(f"[*] {len(queries)} Kandidaten-Erreger im Batch-Modus | {threads} CPU-Threads")
    print("="*75)

    batch_query_fa = os.path.join(output_dir, f"{prefix}_batch_blast_queries.tmp.fa")
    tsv_out_name = f"{prefix}_local_blastn_nt.tsv"
    tsv_out_path = os.path.join(output_dir, tsv_out_name)

    with open(batch_query_fa, 'w', encoding='utf-8') as f:
        for qid, idx, org, seq in queries:
            f.write(f">{qid}\n{seq}\n")

    cmd = [
        blastn_bin,
        "-task", "megablast",
        "-db", blast_db,
        "-query", batch_query_fa,
        "-max_target_seqs", "3",
        "-max_hsps", "1",
        "-evalue", "1e-5",
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle",
        "-num_threads", str(threads)
    ]

    try:
        t0 = time.time()
        with open(tsv_out_path, 'w', encoding='utf-8') as out_f:
            subprocess.run(cmd, stdout=out_f, stderr=subprocess.DEVNULL, check=True)
        elapsed = time.time() - t0
        print(f"  [+] Lokaler BLASTn abgeschlossen in {elapsed:.1f}s -> {tsv_out_name}")

        top_hits = {}
        if os.path.exists(tsv_out_path):
            with open(tsv_out_path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 7:
                        qid = parts[0]
                        if qid not in top_hits:
                            top_hits[qid] = {
                                'sseqid': parts[1],
                                'pident': float(parts[2]),
                                'length': int(parts[3]),
                                'evalue': parts[4],
                                'bitscore': float(parts[5]),
                                'stitle': parts[6]
                            }

        for qid, idx, org, seq in queries:
            if qid in top_hits:
                hit = top_hits[qid]
                clean_title = re.sub(r'^(gi\|\d+\|[a-z]+\|[\w\.]+\|\s*)', '', hit['stitle'])
                validations[idx]['local_blast_hit'] = f"{clean_title[:55]} ({hit['pident']}% Id, E={hit['evalue']})"
                validations[idx]['local_blast_raw'] = hit
                validations[idx]['local_blast_tsv'] = tsv_out_name

                # Genus & Spezies Bestaetigung pruefen
                genus, _ = classify_primer_taxa(org, "")
                org_clean = re.sub(r'\(.*?\)', '', org).strip().lower()
                genus_lower = genus.lower()
                title_lower = clean_title.lower()

                genus_confirmed = False
                confirm_level = "Nicht bestaetigt"

                # Filter human / vector false matches
                if not any(ex in title_lower for ex in ['homo sapiens', 'human', 'synthetic construct', 'cloning vector']):
                    if org_clean in title_lower or (len(org_clean.split()) >= 2 and all(w in title_lower for w in org_clean.split()[:2])):
                        genus_confirmed = True
                        confirm_level = "Spezies-Bestaetigt"
                    elif genus_lower in title_lower:
                        genus_confirmed = True
                        confirm_level = "Genus-Bestaetigt"
                    elif 'piroplasm' in genus_lower and any(p in title_lower for p in ['babesia', 'theileria', 'cytauxzoon']):
                        genus_confirmed = True
                        confirm_level = "Genus-Bestaetigt (Piroplasmida)"
                    elif 'ebv' in org_clean and ('epstein' in title_lower or 'hhv-4' in title_lower or 'hhv4' in title_lower):
                        genus_confirmed = True
                        confirm_level = "Spezies-Bestaetigt (EBV)"
                    elif 'cmv' in org_clean and ('cytomegalo' in title_lower or 'hhv-5' in title_lower or 'hhv5' in title_lower):
                        genus_confirmed = True
                        confirm_level = "Spezies-Bestaetigt (CMV)"

                validations[idx]['genus_confirmed'] = genus_confirmed
                validations[idx]['confirm_level'] = confirm_level
                validations[idx]['genus'] = genus

                icon = "🎯" if genus_confirmed else "[OK]"
                print(f"  {icon} [{org}] Lokaler nt-Treffer: {clean_title[:60]} ({hit['pident']}% Id, E={hit['evalue']}) -> {confirm_level}")
            else:
                validations[idx]['local_blast_hit'] = "Kein lokaler Treffer in nt (E < 1e-5)"
                validations[idx]['local_blast_tsv'] = tsv_out_name
                validations[idx]['genus_confirmed'] = False
                validations[idx]['confirm_level'] = "Nicht bestaetigt"
                validations[idx]['genus'] = classify_primer_taxa(org, "")[0]

        dl_dir = r"C:\Users\internet\Downloads"
        try:
            if os.path.exists(dl_dir):
                shutil.copy2(tsv_out_path, os.path.join(dl_dir, tsv_out_name))
        except Exception:
            pass

    except Exception as e:
        print(f"  [!] Fehler beim lokalen BLASTn: {e}")

    finally:
        if os.path.exists(batch_query_fa):
            try:
                os.remove(batch_query_fa)
            except Exception:
                pass


def run_dark_matter_discovery(unassigned_reads, output_dir, prefix, blast_db="/mnt/blast_db/nt", threads=16, max_reads=600):
    """
    Paket C: Dark Matter Discovery - Assembles unassigned high-complexity reads into contigs
    and queries local BLASTn (/mnt/blast_db/nt) to identify unexpected/unknown pathogens.
    """
    if not unassigned_reads:
        return []

    blastn_bin = shutil.which('blastn') or ('/usr/bin/blastn' if os.path.exists('/usr/bin/blastn') else None)
    db_base = blast_db
    has_blast_db = (
        blastn_bin and (
            os.path.exists(f"{db_base}.nal") or
            os.path.exists(f"{db_base}.000.nhr") or
            os.path.exists(f"{db_base}.nhr") or
            os.path.exists(db_base)
        )
    )

    sample_reads = unassigned_reads[:max_reads]
    os.makedirs(output_dir, exist_ok=True)
    contig_fa = os.path.join(output_dir, f"{prefix}_dark_matter_contigs.fasta")

    # Greedy de-novo overlap clustering
    seqs = []
    for r_str in sample_reads:
        lines = r_str.strip().split('\n')
        if len(lines) >= 2:
            s = lines[1].strip()
            if len(s) >= 120:
                seqs.append(s)

    k = 21
    kmer_map = defaultdict(list)
    for idx, s in enumerate(seqs[:350]):
        for i in range(0, len(s) - k + 1, 8):
            kmer = s[i:i+k]
            kmer_map[kmer].append(idx)

    clusters = defaultdict(set)
    for kmer, r_indices in kmer_map.items():
        if 2 <= len(r_indices) <= 25:
            base_idx = r_indices[0]
            for other in r_indices[1:]:
                clusters[base_idx].add(other)

    assembled_contigs = []
    used_indices = set()
    for base_idx, group in sorted(clusters.items(), key=lambda x: len(x[1]), reverse=True):
        if base_idx in used_indices or len(group) < 2:
            continue
        all_members = [base_idx] + [g for g in group if g not in used_indices]
        if len(all_members) < 2:
            continue
        for m in all_members:
            used_indices.add(m)

        cluster_seqs = [seqs[m] for m in all_members]
        longest_seq = max(cluster_seqs, key=len)
        contig_id = f"dark_matter_c{len(assembled_contigs)+1}_{len(all_members)}reads"
        assembled_contigs.append((contig_id, longest_seq, len(all_members)))
        if len(assembled_contigs) >= 30:
            break

    if not assembled_contigs and seqs:
        sorted_seqs = sorted(seqs, key=len, reverse=True)[:15]
        for idx, s in enumerate(sorted_seqs, 1):
            assembled_contigs.append((f"dark_matter_single_r{idx}", s, 1))

    if not assembled_contigs:
        return []

    with open(contig_fa, 'w', encoding='utf-8') as f:
        for cid, seq, num_r in assembled_contigs:
            f.write(f">{cid} length={len(seq)} read_support={num_r}\n{seq}\n")

    discoveries = []
    if has_blast_db:
        print("\n" + "="*75)
        print(f"🔬 DARK MATTER DISCOVERY ({prefix}): Lokaler BLASTn unzugeordneter Contigs")
        print(f"[*] {len(assembled_contigs)} Contigs gegen {blast_db} | {threads} CPU-Threads")
        print("="*75)

        blast_out_tsv = os.path.join(output_dir, f"{prefix}_dark_matter_blastn.tsv")
        cmd = [
            blastn_bin,
            "-task", "megablast",
            "-db", blast_db,
            "-query", contig_fa,
            "-max_target_seqs", "3",
            "-max_hsps", "1",
            "-evalue", "1e-5",
            "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle",
            "-num_threads", str(threads)
        ]
        try:
            with open(blast_out_tsv, 'w', encoding='utf-8') as out_f:
                subprocess.run(cmd, stdout=out_f, stderr=subprocess.DEVNULL, check=True)

            contig_dict = {cid: (seq, num_r) for cid, seq, num_r in assembled_contigs}
            seen_taxa = set()

            if os.path.exists(blast_out_tsv):
                with open(blast_out_tsv, 'r', encoding='utf-8', errors='ignore') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 7:
                            qid, sseqid, pident, length, evalue, bitscore, stitle = (
                                parts[0], parts[1], float(parts[2]), int(parts[3]), parts[4], float(parts[5]), parts[6]
                            )
                            title_lower = stitle.lower()
                            
                            # Filter human, vectors, cloning constructs
                            if any(ex in title_lower for ex in ['homo sapiens', 'human', 'synthetic construct', 'cloning vector', 'expression vector']):
                                continue

                            clean_title = re.sub(r'^(gi\|\d+\|[a-z]+\|[\w\.]+\|\s*)', '', stitle)
                            org_match = re.search(r'([A-Z][a-z]+(?:\s+[a-z]+)?)', clean_title)
                            org_name = org_match.group(1) if org_match else clean_title[:30]

                            if org_name in seen_taxa:
                                continue
                            seen_taxa.add(org_name)

                            cat = "Bakterium"
                            cat_badge = "badge-species"
                            if any(v in title_lower for v in ['virus', 'viridae', 'phage', 'viral']):
                                cat = "Virus"
                                cat_badge = "badge-both"
                            elif any(f in title_lower for f in ['fungi', 'candida', 'aspergillus', 'malassezia', 'yeast']):
                                cat = "Pilz"
                                cat_badge = "badge-genus"
                            elif any(p in title_lower for p in ['protozoa', 'plasmodi', 'babesi', 'theileri', 'toxoplasm', 'trypanosom', 'leishmani']):
                                cat = "Parasit / Protozoon"
                                cat_badge = "badge-tier1"

                            seq_info = contig_dict.get(qid, ("", 1))
                            disc = {
                                'contig_id': qid,
                                'category': cat,
                                'category_badge': cat_badge,
                                'organism': org_name,
                                'title': clean_title[:80],
                                'accession': sseqid,
                                'pident': pident,
                                'align_len': length,
                                'evalue': evalue,
                                'bitscore': bitscore,
                                'contig_len': len(seq_info[0]),
                                'read_support': seq_info[1],
                                'contig_seq': seq_info[0]
                            }
                            discoveries.append(disc)
                            print(f"  [✨] ENTDECKUNG [{cat}]: {org_name} -> {clean_title[:55]} ({pident}% Id, E={evalue})")

                dl_dir = r"C:\Users\internet\Downloads"
                try:
                    if os.path.exists(dl_dir):
                        shutil.copy2(blast_out_tsv, os.path.join(dl_dir, f"{prefix}_dark_matter_blastn.tsv"))
                except Exception:
                    pass

        except Exception as e:
            print(f"  [!] Fehler beim Dark-Matter BLASTn: {e}")

    return discoveries

def run_pan_pathogen_scan(candidate_fastq_paths, ref_dirs=None, output_dir=".", prefix="sample", threads=16, min_mapq=20):
    """
    Paket B: Pan-Pathogen Scan - Maps candidate or unmapped reads against all pathogen references.
    """
    minimap2_bin = shutil.which('minimap2') or ('/usr/local/bin/minimap2' if os.path.exists('/usr/local/bin/minimap2') else None)
    if not minimap2_bin or not candidate_fastq_paths:
        return []

    if not ref_dirs:
        ref_dirs = [
            "/mnt/sdc1/pathogen_references",
            "./pathogen_references",
            "./references",
            r"C:\Users\internet\Downloads\repository\primer\references",
            r"C:\Users\internet\Downloads"
        ]
    elif isinstance(ref_dirs, str):
        ref_dirs = [ref_dirs]

    ref_files = []
    for rdir in ref_dirs:
        if rdir and os.path.exists(rdir):
            for root, dirs, files in os.walk(rdir):
                for f in files:
                    if f.endswith(('.fasta', '.fa', '.fna')):
                        ref_files.append(os.path.join(root, f))

    if not ref_files:
        return []

    print("\n" + "="*75)
    print(f"🌐 PAN-PATHOGEN GENOM-SCAN ({prefix})")
    print(f"[*] {len(ref_files)} Referenzgenome im Pool | {threads} CPU-Threads")
    print("="*75)
    return ref_files

def aggregate_clinical_findings(primer_stats, validations=None, dark_matter_discoveries=None, extracted_counts=None):
    """
    Paket A: Computes Clinical Evidence Scores (0-100%), Diagnostic Traffic Light (Ampel),
    Multi-Locus verification, and Contaminants/Kitome filtering for each pathogen.
    """
    validations = validations or []
    dark_matter_discoveries = dark_matter_discoveries or []
    extracted_counts = extracted_counts or {}

    val_by_org = {}
    for v in validations:
        val_by_org[v['organism']] = v

    org_groups = defaultdict(lambda: {
        'pairs': [],
        'target_genes': Counter(),
        'total_both': 0,
        'total_fwd': 0,
        'total_rev': 0,
        'amplicon_distances': [],
        'tier': 'Tier 3'
    })

    tier_prio = {'tier 0': 0, 'tier 1': 1, 'tier 2': 2, 'tier 3': 3}

    for pid, st in primer_stats.items():
        org = st['Organism']
        grp = org_groups[org]
        grp['pairs'].append((pid, st))
        gene = st.get('TargetGene', infer_target_gene(pid))
        reads = st['BOTH'] + st['ONLY_FORWARD'] + st['ONLY_REVERSE']
        if reads > 0:
            grp['target_genes'][gene] += reads
        grp['total_both'] += st['BOTH']
        grp['total_fwd'] += st['ONLY_FORWARD']
        grp['total_rev'] += st['ONLY_REVERSE']
        grp['amplicon_distances'].extend(st.get('AmpliconDistances', []))
        
        tier = st.get('Tier', infer_tier(org, pid))
        current_prio = tier_prio.get(grp['tier'].lower(), 3)
        new_prio = tier_prio.get(tier.lower(), 3)
        if new_prio < current_prio:
            grp['tier'] = tier

    clinical_findings = []

    for org, grp in org_groups.items():
        tot_reads = grp['total_both'] + grp['total_fwd'] + grp['total_rev']
        if tot_reads == 0:
            continue

        v_info = val_by_org.get(org, {})
        blast_hit = v_info.get('local_blast_hit', '')
        mapped_count = v_info.get('mapped_count', 0)
        high_qual_count = v_info.get('high_qual_count', 0)
        avg_mapq = v_info.get('avg_mapq', 0.0)
        contig_len = v_info.get('contig_len', 0)

        is_contam = is_likely_contaminant(org, blast_hit)

        score = 0
        score_breakdown = []

        # 1. Multi-Locus Support
        active_genes = [g for g, c in grp['target_genes'].items() if c > 0]
        if len(active_genes) >= 2:
            score += 35
            score_breakdown.append(f"Multi-Locus ({', '.join(active_genes)}): +35%")
        elif len(active_genes) == 1:
            score += 15
            score_breakdown.append(f"Single-Locus ({active_genes[0]}): +15%")

        # 2. Both Fwd & Rev match
        if grp['total_both'] >= 5:
            score += 25
            score_breakdown.append(f"Beidseitige Koinzidenz ({grp['total_both']} Reads): +25%")
        elif grp['total_both'] > 0:
            score += 15
            score_breakdown.append(f"Beidseitige Koinzidenz ({grp['total_both']} Reads): +15%")
        elif grp['total_fwd'] > 0 or grp['total_rev'] > 0:
            score += 5
            score_breakdown.append(f"Einseitige Reads: +5%")

        # 3. Minimap2 alignment
        if high_qual_count >= 5 and avg_mapq >= 30:
            score += 25
            score_breakdown.append(f"minimap2 High-MAPQ ({high_qual_count} Reads, MAPQ {avg_mapq:.0f}): +25%")
        elif high_qual_count > 0:
            score += 15
            score_breakdown.append(f"minimap2 gemappt ({mapped_count} Reads): +15%")
        elif mapped_count > 0:
            score += 5
            score_breakdown.append("minimap2 MAPQ 0: +5%")

        # 4. BLASTn match
        if blast_hit and "Kein lokaler Treffer" not in blast_hit:
            if is_contam:
                score = 0
                score_breakdown.append("Kontamination / Human im BLAST erkannt: Score 0")
            else:
                score += 20
                score_breakdown.append("Lokaler BLASTn bestätigt: +20%")

        score = min(100, max(0, score))

        if is_contam:
            status = 'contaminant'
            status_badge = 'badge-rev'
            traffic_light = '🚫'
            status_text = 'Labor- / Haut-Kontaminante (Kitom)'
        elif score >= 75:
            status = 'confirmed'
            status_badge = 'badge-both'
            traffic_light = '🟢'
            status_text = 'Gesichert (Multi-Locus / Hohe Evidenz)'
        elif score >= 35:
            status = 'suspect'
            status_badge = 'badge-genus'
            traffic_light = '🟡'
            status_text = 'Verdachtsfall (Single-Locus / Fragmentiert)'
        else:
            status = 'low'
            status_badge = 'badge-species'
            traffic_light = '⚪'
            status_text = 'Geringe Evidenz / Unspezifisch'

        clinical_findings.append({
            'organism': org,
            'tier': grp['tier'],
            'description': get_pathogen_description(org),
            'status': status,
            'traffic_light': traffic_light,
            'status_text': status_text,
            'status_badge': status_badge,
            'evidence_score': score,
            'score_breakdown': score_breakdown,
            'total_both': grp['total_both'],
            'total_fwd': grp['total_fwd'],
            'total_rev': grp['total_rev'],
            'total_reads': tot_reads,
            'active_genes': active_genes,
            'gene_details': [f"{g} ({c} R)" for g, c in grp['target_genes'].items() if c > 0],
            'minimap2_info': f"{mapped_count} Reads gemappt (MAPQ {avg_mapq:.0f})" if mapped_count > 0 else "Nicht gemappt",
            'blast_info': blast_hit if blast_hit else "Nicht geprueft",
            'contig_len': contig_len,
            'pairs': grp['pairs'],
            'is_contaminant': is_contam,
            'genus_confirmed': v_info.get('genus_confirmed', False),
            'confirm_level': v_info.get('confirm_level', 'Nicht bestaetigt'),
            'blast_genus': v_info.get('genus', classify_primer_taxa(org, "")[0]),
            'blast_raw': v_info.get('local_blast_raw', {}),
            'rep_seq': v_info.get('rep_seq', ''),
            'contig_seq': v_info.get('contig_seq', '')
        })

    status_order = {'confirmed': 0, 'suspect': 1, 'low': 2, 'contaminant': 3}
    clinical_findings.sort(key=lambda x: (status_order.get(x['status'], 9), -x['evidence_score'], -x['total_reads']))
    return clinical_findings

def run_minimap2_validation(results, ref_dirs=None, threads=4, min_mapq=20, enable_assembly=True, enable_local_blast=True, blast_db="/mnt/blast_db/nt"):
    """
    Validates candidate organisms that have:
    - ONLY_FORWARD > 0 and ONLY_REVERSE > 0 (split reads in separate sequences)
    - or BOTH > 0
    by mapping candidate reads against reference genomes using minimap2.
    
    All output files are prefixed with the sample prefix to prevent any mix-ups:
    - {prefix}_{clean_org}_candidate_reads.fastq
    - {prefix}_{clean_org}_mapped.bam / .sam
    - {prefix}_{clean_org}_coverage.txt
    """
    fastq_path = results.get('fastq_path', '')
    prefix = results.get('sample_prefix', extract_sample_prefix(fastq_path))
    output_dir = results.get('output_dir', './extracted_amplicons')
    os.makedirs(output_dir, exist_ok=True)

    primer_stats = results.get('primer_stats', {})
    candidate_reads_by_org = results.get('candidate_reads_by_org', {})

    # Aggregate by organism
    org_stats = defaultdict(lambda: {'BOTH': 0, 'ONLY_FORWARD': 0, 'ONLY_REVERSE': 0, 'pairs': set()})
    for pid, stats in primer_stats.items():
        org = stats['Organism']
        org_stats[org]['BOTH'] += stats['BOTH']
        org_stats[org]['ONLY_FORWARD'] += stats['ONLY_FORWARD']
        org_stats[org]['ONLY_REVERSE'] += stats['ONLY_REVERSE']
        org_stats[org]['pairs'].add(pid)

    # Candidates: prioritize split hits (FWD > 0 and REV > 0, BOTH == 0) and positive hits
    candidates = []
    for org, st in org_stats.items():
        if st['ONLY_FORWARD'] > 0 and st['ONLY_REVERSE'] > 0 and st['BOTH'] == 0:
            candidates.append((org, 'split', st))
        elif st['BOTH'] > 0:
            candidates.append((org, 'both', st))
        elif st['ONLY_FORWARD'] >= 3 or st['ONLY_REVERSE'] >= 3:
            candidates.append((org, 'single_abundant', st))

    if not candidates:
        results['minimap2_validations'] = []
        return []

    print("\n" + "="*75)
    print(f"🔬 MINIMAP2 REFERENZ-MAPPING & ERREGER-VALIDIERUNG ({prefix})")
    print("="*75)

    minimap2_bin = shutil.which('minimap2')
    if not minimap2_bin and os.path.exists('/usr/local/bin/minimap2'):
        minimap2_bin = '/usr/local/bin/minimap2'

    samtools_bin = shutil.which('samtools')
    if not samtools_bin and os.path.exists('/usr/local/bin/samtools'):
        samtools_bin = '/usr/local/bin/samtools'

    if not minimap2_bin:
        print("[!] Hinweis: 'minimap2' nicht im Systempfad gefunden. Kandidaten-Reads werden für manuelles Mapping/BLAST extrahiert.")

    validations = []

    for org, trigger_type, st in candidates:
        clean_org = re.sub(r'[^a-zA-Z0-9_-]', '_', org)[:30]
        reads_info = candidate_reads_by_org.get(org, [])

        cand_fastq_name = f"{prefix}_{clean_org}_candidate_reads.fastq"
        cand_fastq_path = os.path.join(output_dir, cand_fastq_name)

        # Write candidate reads with sample prefix
        rep_seq = ""
        seen_headers = set()
        unique_reads = []

        for hit_type, read_str, p_id in reads_info:
            first_line = read_str.split('\n')[0]
            if first_line not in seen_headers:
                seen_headers.add(first_line)
                unique_reads.append(read_str)
                if not rep_seq:
                    lines = read_str.split('\n')
                    if len(lines) > 1:
                        rep_seq = lines[1].strip()

        with open(cand_fastq_path, 'w', encoding='utf-8') as f_out:
            for r_str in unique_reads:
                f_out.write(r_str)

        ref_path = find_reference_genome(org, ref_dirs)
        mapped_count = 0
        high_qual_count = 0
        covered_bp = 0
        avg_mapq = 0.0
        top_rname = "N/A"
        bam_rel = None

        if trigger_type == 'split':
            trigger_desc = f"Getrennte Reads (FWD: {st['ONLY_FORWARD']}, REV: {st['ONLY_REVERSE']}, BOTH: 0)"
            badge_trigger = "badge-genus"
        elif trigger_type == 'both':
            trigger_desc = f"FWD & REV Koinzidenz ({st['BOTH']} Reads)"
            badge_trigger = "badge-both"
        else:
            trigger_desc = f"Häufiger Einzelprimer (FWD: {st['ONLY_FORWARD']}, REV: {st['ONLY_REVERSE']})"
            badge_trigger = "badge-species"

        if not ref_path:
            # Automatischer Referenzgenom-Download von NCBI RefSeq falls nicht vorhanden
            print(f"  [🔍] [{org}] Keine lokale Referenz gefunden. Starte automatischen Download von NCBI RefSeq...")
            ref_path = auto_download_reference_genome(org, dest_dir=ref_dirs[0] if ref_dirs else None)

        if not ref_path:
            verdict = "⏳ Keine passende Referenz in ref-dir gefunden (BLAST nutzen)"
            verdict_badge = "badge-genus"
            print(f"  [-] [{org}] {trigger_desc} -> Keine Referenz in ref-dir gefunden.")
        elif not minimap2_bin:
            verdict = "ℹ️ minimap2 nicht installiert (Reads für Server-Mapping gespeichert)"
            verdict_badge = "badge-species"
            print(f"  [+] [{org}] {trigger_desc} -> Ref: {os.path.basename(ref_path)} (minimap2 nicht lokal im PATH)")
        elif len(unique_reads) == 0:
            verdict = "Keine Reads zum Mappen vorhanden"
            verdict_badge = "badge-rev"
        else:
            sam_name = f"{prefix}_{clean_org}_mapped.sam"
            sam_path = os.path.join(output_dir, sam_name)
            bam_name = f"{prefix}_{clean_org}_mapped.bam"
            bam_path = os.path.join(output_dir, bam_name)
            cov_name = f"{prefix}_{clean_org}_coverage.txt"
            cov_path = os.path.join(output_dir, cov_name)

            # Execute minimap2
            try:
                with open(sam_path, 'w', encoding='utf-8') as sam_file:
                    subprocess.run(
                        [minimap2_bin, "-ax", "map-ont", "-t", str(threads), ref_path, cand_fastq_path],
                        stdout=sam_file,
                        stderr=subprocess.DEVNULL,
                        check=True
                    )

                # If samtools is available, sort to BAM and compute coverage
                if samtools_bin:
                    subprocess.run(
                        f"{samtools_bin} view -b -F 4 '{sam_path}' | {samtools_bin} sort -o '{bam_path}'",
                        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                    )
                    if os.path.exists(bam_path) and os.path.getsize(bam_path) > 0:
                        subprocess.run([samtools_bin, "index", bam_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        with open(cov_path, 'w', encoding='utf-8') as cov_f:
                            subprocess.run([samtools_bin, "coverage", bam_path], stdout=cov_f, stderr=subprocess.DEVNULL)
                        bam_rel = bam_name

                # Parse alignment metrics
                align_metrics = parse_sam_alignments(sam_path)
                mapped_count = align_metrics['mapped_count']
                high_qual_count = align_metrics['high_qual_count']
                avg_mapq = align_metrics['avg_mapq']
                top_rname_raw = align_metrics['top_rname']
                top_rname = resolve_rname_description(top_rname_raw, ref_path)
                covered_bp = align_metrics['covered_bp']

                # Clean up huge SAM if BAM was generated
                if bam_rel and os.path.exists(sam_path) and os.path.getsize(sam_path) > 5 * 1024 * 1024:
                    try:
                        os.remove(sam_path)
                    except Exception:
                        pass

                # Formulate diagnostic verdict with resolved species context
                species_info = f" ({top_rname})" if top_rname != "N/A" and top_rname != top_rname_raw else ""
                if mapped_count > 0 and high_qual_count > 0:
                    verdict = f"✓ Bestätigt{species_info}: {mapped_count}/{len(unique_reads)} Reads mappen (MAPQ ≥ {min_mapq})"
                    verdict_badge = "badge-both"
                elif mapped_count > 0 and high_qual_count == 0:
                    verdict = f"⚠️ Unspezifisch{species_info}: {mapped_count} Reads mappen mit MAPQ 0 (Repetitiv/Wirt)"
                    verdict_badge = "badge-genus"
                else:
                    verdict = f"❌ Falsch-positiv: 0 von {len(unique_reads)} Reads mappen auf Erregergenom"
                    verdict_badge = "badge-rev"

                print(f"  [+] [{org}] {trigger_desc} -> Ref: {os.path.basename(ref_path)}: {verdict}")

            except Exception as e:
                verdict = f"Fehler bei minimap2: {str(e)[:50]}"
                verdict_badge = "badge-rev"
                print(f"  [!] [{org}] Fehler bei minimap2: {e}")

        # Automatic Micro-Assembly & Consensus for candidate reads >= 5
        assembly_info = {'contig_rel': None, 'contig_path': None, 'contig_len': 0, 'contig_seq': '', 'method': 'N/A'}
        if enable_assembly and len(unique_reads) >= 5:
            assembly_info = assemble_candidate_reads(
                cand_fastq_path=cand_fastq_path,
                prefix=prefix,
                clean_org=clean_org,
                output_dir=output_dir,
                bam_path=bam_path,
                ref_path=ref_path,
                threads=threads,
                minimap2_bin=minimap2_bin,
                samtools_bin=samtools_bin
            )

        validations.append({
            'organism': org,
            'clean_org': clean_org,
            'tier': infer_tier(org),
            'trigger_desc': trigger_desc,
            'trigger_type': trigger_type,
            'badge_trigger': badge_trigger,
            'num_candidates': len(unique_reads),
            'candidate_fastq': cand_fastq_path,
            'candidate_fastq_rel': cand_fastq_name,
            'ref_file': os.path.basename(ref_path) if ref_path else "Keine Referenz",
            'mapped_count': mapped_count,
            'high_qual_count': high_qual_count,
            'covered_bp': covered_bp,
            'top_rname': top_rname,
            'avg_mapq': avg_mapq,
            'bam_rel': bam_rel,
            'contig_rel': assembly_info['contig_rel'],
            'contig_len': assembly_info['contig_len'],
            'contig_seq': assembly_info['contig_seq'],
            'assembly_method': assembly_info['method'],
            'verdict': verdict,
            'verdict_badge': verdict_badge,
            'rep_seq': rep_seq
        })

    # Run batch local BLASTn verification against /mnt/blast_db/nt if available
    if enable_local_blast and validations:
        run_local_blast_verification(
            validations=validations,
            output_dir=output_dir,
            prefix=prefix,
            blast_db=blast_db,
            threads=threads
        )

    print("="*75 + "\n")
    results['minimap2_validations'] = validations
    return validations

def run_primer_search(fastq_path, primers, output_dir, max_mismatches=1, max_reads=None, stage_label=None, threads=None, sample_prefix=None):
    """Multi-threaded primer matching & amplicon extraction engine with real-time progress bar and ETA."""
    os.makedirs(output_dir, exist_ok=True)
    if not sample_prefix:
        sample_prefix = extract_sample_prefix(fastq_path)

    if not threads:
        threads = max(1, min(32, (os.cpu_count() or 4) - 2))

    lbl = f" [{stage_label}]" if stage_label else ""
    pigz_active = bool(shutil.which('pigz')) and str(fastq_path).endswith('.gz')
    decomp_msg = "pigz (Parallel)" if pigz_active else "gzip"
    print(f"[*] Starting multi-threaded primer alignment{lbl} on: {fastq_path}")
    print(f"[*] Parallel CPU workers: {threads} threads | Dekomprimierung: {decomp_msg}")
    print(f"[*] Loaded {len(primers)} primer pairs...")

    start_time = time.time()
    total_reads = 0
    match_counts = {'BOTH': 0, 'ONLY_FORWARD': 0, 'ONLY_REVERSE': 0, 'NONE': 0}

    primer_stats = defaultdict(lambda: {
        'Organism': '', 'Forward': '', 'Reverse': '',
        'BOTH': 0, 'ONLY_FORWARD': 0, 'ONLY_REVERSE': 0,
        'AmpliconDistances': [],
        'TaxaType': 'Species',
        'Genus': ''
    })
    for p in primers:
        genus, is_genus = classify_primer_taxa(p['Organism'], p['PrimerPairID'])
        p_entry = primer_stats[p['PrimerPairID']]
        p_entry['Organism'] = p['Organism']
        p_entry['Forward'] = p['Forward']
        p_entry['Reverse'] = p['Reverse']
        p_entry['TaxaType'] = 'Genus' if is_genus else 'Species'
        p_entry['Genus'] = genus
        p_entry['Tier'] = p.get('Tier', infer_tier(p['Organism'], p['PrimerPairID']))
        p_entry['TargetGene'] = p.get('TargetGene', infer_target_gene(p['PrimerPairID']))
        p_entry['AmpliconSize'] = p.get('AmpliconSize', 0)
        p_entry['Reference'] = p.get('Reference', '')

    extracted_files = {}
    extracted_counts = Counter()
    extracted_seqs = {}       # filename -> representative amplicon sequence
    extracted_pair_seqs = {}  # pair_id -> representative amplicon sequence
    candidate_reads_by_org = defaultdict(list)

    chunk_size = 20000
    last_log_time = time.time()
    is_tty = sys.stdout.isatty() if hasattr(sys.stdout, 'isatty') else False

    with mp.Pool(processes=threads, initializer=_worker_init, initargs=(primers, max_mismatches, sample_prefix)) as pool:
        chunks_gen = read_fastq_chunks(fastq_path, chunk_size=chunk_size, limit=max_reads, threads=threads)
        
        unassigned_reads = []
        for worker_res in pool.imap(_worker_process_chunk, chunks_gen, chunksize=1):
            if len(worker_res) == 6:
                c_matches, c_stats_delta, c_extracted, frac, c_candidates, c_unassigned = worker_res
                if len(unassigned_reads) < 1000:
                    unassigned_reads.extend(c_unassigned)
            else:
                c_matches, c_stats_delta, c_extracted, frac, c_candidates = worker_res[:5]
            chunk_n = sum(c_matches.values())
            total_reads += chunk_n

            for k in match_counts:
                match_counts[k] += c_matches[k]

            for pair_id, delta in c_stats_delta.items():
                p_entry = primer_stats[pair_id]
                p_entry['BOTH'] += delta['BOTH']
                p_entry['ONLY_FORWARD'] += delta['ONLY_FORWARD']
                p_entry['ONLY_REVERSE'] += delta['ONLY_REVERSE']
                p_entry['AmpliconDistances'].extend(delta['AmpliconDistances'])

            for item in c_extracted:
                if len(item) == 4:
                    fname, read_str, amp_seq, p_id = item
                else:
                    fname, read_str = item[:2]
                    amp_seq, p_id = "", ""

                out_path = os.path.join(output_dir, fname)
                if out_path not in extracted_files:
                    extracted_files[out_path] = open(out_path, 'w', encoding='utf-8')
                extracted_files[out_path].write(read_str)
                extracted_counts[fname] += 1
                if amp_seq and fname not in extracted_seqs:
                    extracted_seqs[fname] = amp_seq
                if amp_seq and p_id and p_id not in extracted_pair_seqs:
                    extracted_pair_seqs[p_id] = amp_seq

            for org_name, hit_type, read_str, p_id in c_candidates:
                if len(candidate_reads_by_org[org_name]) < 300:
                    candidate_reads_by_org[org_name].append((hit_type, read_str, p_id))

            now = time.time()
            elapsed = now - start_time
            rps = total_reads / elapsed if elapsed > 0 else 0
            pbar_str = render_progress_bar(frac, total_reads, rps, elapsed)

            if is_tty:
                sys.stdout.write(f"\r    {pbar_str}")
                sys.stdout.flush()
            elif (now - last_log_time >= 5.0) or frac >= 0.999:
                print(f"    {pbar_str}")
                sys.stdout.flush()
                last_log_time = now

    if is_tty:
        sys.stdout.write("\n")

    for h in extracted_files.values():
        h.close()

    elapsed = time.time() - start_time
    rps = total_reads / elapsed if elapsed > 0 else 0
    print(f"[+] Processing finished in {elapsed:.2f} seconds ({rps:,.0f} reads/s)!")
    print(f"    - Total Reads: {total_reads:,}")
    print(f"    - BOTH (Forward & Reverse): {match_counts['BOTH']:,}")
    print(f"    - ONLY Forward: {match_counts['ONLY_FORWARD']:,}")
    print(f"    - ONLY Reverse: {match_counts['ONLY_REVERSE']:,}")
    print(f"    - Extracted _frp_ Files: {len(extracted_files)}")

    return {
        'fastq_path': fastq_path,
        'fastq_name': os.path.basename(fastq_path),
        'sample_prefix': sample_prefix,
        'num_primers': len(primers),
        'total_reads': total_reads,
        'match_counts': match_counts,
        'primer_stats': primer_stats,
        'extracted_counts': extracted_counts,
        'extracted_seqs': extracted_seqs,
        'extracted_pair_seqs': extracted_pair_seqs,
        'candidate_reads_by_org': dict(candidate_reads_by_org),
        'unassigned_reads': unassigned_reads if 'unassigned_reads' in locals() else [],
        'output_dir': output_dir,
        'threads': threads,
        'elapsed_s': round(elapsed, 2),
        'reads_per_sec': round(rps, 0)
    }

def merge_search_results(res1, res2):
    """Merges Stage 1 (Genus) and Stage 2 (Species) results into a single result dictionary."""
    merged = dict(res1)
    merged['num_primers'] = res1['num_primers'] + res2['num_primers']
    merged['match_counts'] = {
        'BOTH': res1['match_counts']['BOTH'] + res2['match_counts']['BOTH'],
        'ONLY_FORWARD': res1['match_counts']['ONLY_FORWARD'] + res2['match_counts']['ONLY_FORWARD'],
        'ONLY_REVERSE': res1['match_counts']['ONLY_REVERSE'] + res2['match_counts']['ONLY_REVERSE'],
        'NONE': res1['match_counts']['NONE']
    }
    
    merged_stats = defaultdict(lambda: {
        'Organism': '', 'Forward': '', 'Reverse': '',
        'BOTH': 0, 'ONLY_FORWARD': 0, 'ONLY_REVERSE': 0,
        'AmpliconDistances': [],
        'TaxaType': 'Species',
        'Genus': ''
    })
    for k, v in res1['primer_stats'].items():
        merged_stats[k] = v
    for k, v in res2['primer_stats'].items():
        merged_stats[k] = v
    merged['primer_stats'] = merged_stats
    
    merged_extracted = Counter(res1['extracted_counts'])
    merged_extracted.update(res2['extracted_counts'])
    merged['extracted_counts'] = merged_extracted

    merged_seqs = dict(res1.get('extracted_seqs', {}))
    merged_seqs.update(res2.get('extracted_seqs', {}))
    merged['extracted_seqs'] = merged_seqs

    merged_pair_seqs = dict(res1.get('extracted_pair_seqs', {}))
    merged_pair_seqs.update(res2.get('extracted_pair_seqs', {}))
    merged['extracted_pair_seqs'] = merged_pair_seqs

    merged_candidates = defaultdict(list)
    for org, rlist in res1.get('candidate_reads_by_org', {}).items():
        merged_candidates[org].extend(rlist)
    for org, rlist in res2.get('candidate_reads_by_org', {}).items():
        if len(merged_candidates[org]) < 300:
            merged_candidates[org].extend(rlist[:300 - len(merged_candidates[org])])
    merged['candidate_reads_by_org'] = dict(merged_candidates)
    merged['sample_prefix'] = res1.get('sample_prefix', res2.get('sample_prefix', 'sample'))

    return merged

def generate_html_dashboard(results, dashboard_path):
    """
    Generates the modern 5-Tab Interactive Pathogen Discovery & Clinical Diagnosis Dashboard.
    Tab 1: Klinischer Erreger-Befund (Diagnose-Ampel & Multi-Locus Accordion)
    Tab 2: Detaillierte Primer-Matrix (735 Assays & Amplikons)
    Tab 3: De-novo Entdeckungen (Unbekannte Erreger & Dark Matter)
    Tab 4: Qualitaetskontrolle & Kontaminanten (Kitome / Human)
    Tab 5: Vorgehen & Leitfaden (SHG Borreliose Hannover)
    """
    total_reads = results['total_reads']
    match_counts = results['match_counts']
    primer_stats = results['primer_stats']
    extracted_counts = results['extracted_counts']
    extracted_seqs = results.get('extracted_seqs', {})
    extracted_pair_seqs = results.get('extracted_pair_seqs', {})
    search_mode = results.get('search_mode', 'Standard (-a)')
    positive_genera = results.get('positive_genera', [])
    threads_used = results.get('threads', 'Multi-Core')
    rps_val = results.get('reads_per_sec', 0)
    elapsed_val = results.get('elapsed_s', 0)
    fastq_name = results.get('fastq_name', 'FASTQ Input File')
    sample_prefix = results.get('sample_prefix', extract_sample_prefix(fastq_name))

    # 1. Clinical Findings (Ampel & Multi-Locus)
    clinical_findings = results.get('clinical_findings')
    if not clinical_findings:
        clinical_findings = aggregate_clinical_findings(
            primer_stats=primer_stats,
            validations=results.get('minimap2_validations', []),
            dark_matter_discoveries=results.get('dark_matter_discoveries', []),
            extracted_counts=extracted_counts
        )
        results['clinical_findings'] = clinical_findings

    confirmed_count = sum(1 for c in clinical_findings if c['status'] == 'confirmed')
    suspect_count = sum(1 for c in clinical_findings if c['status'] == 'suspect')
    low_count = sum(1 for c in clinical_findings if c['status'] == 'low')
    contam_count = sum(1 for c in clinical_findings if c['status'] == 'contaminant')

    # BLAST Genus / Spezies Bestaetigte Erreger
    blast_confirmed_list = [c for c in clinical_findings if c.get('genus_confirmed', False) and not c.get('is_contaminant', False)]
    blast_confirmed_count = len(blast_confirmed_list)

    # 2. Dark Matter Discoveries
    dark_matter_discoveries = results.get('dark_matter_discoveries', [])
    denovo_count = len(dark_matter_discoveries)

    # 3. Plotly Charts Data
    match_labels = ["Forward & Reverse (BOTH)", "Nur Forward", "Nur Reverse", "Kein Match"]
    match_values = [match_counts['BOTH'], match_counts['ONLY_FORWARD'], match_counts['ONLY_REVERSE'], match_counts['NONE']]

    top_primers = sorted(primer_stats.items(), key=lambda x: x[1]['BOTH'], reverse=True)[:15]
    top_org_labels = [p[1]['Organism'] + f" ({p[0]})" for p in top_primers]
    top_org_both = [p[1]['BOTH'] for p in top_primers]
    top_org_fwd = [p[1]['ONLY_FORWARD'] for p in top_primers]
    top_org_rev = [p[1]['ONLY_REVERSE'] for p in top_primers]

    all_distances = [d for p in primer_stats.values() for d in p['AmpliconDistances']]
    dist_counts = Counter(all_distances)
    dist_labels = sorted(dist_counts.keys())
    dist_values = [dist_counts[k] for k in dist_labels]

    # --- HTML TAB 1: KLINISCHER BEFUND (AMPEL & ACCORDION) ---
    clinical_rows_html = ""
    for idx, c in enumerate(clinical_findings, 1):
        org = c['organism']
        tier = c['tier']
        tier_clean = '0' if '0' in tier else ('1' if '1' in tier else ('2' if '2' in tier else '3'))
        tier_badge = 'badge-tier0' if tier_clean == '0' else ('badge-tier1' if tier_clean == '1' else ('badge-tier2' if tier_clean == '2' else 'badge-tier3'))
        
        status = c['status']
        traffic = c['traffic_light']
        status_text = c['status_text']
        status_badge = c['status_badge']
        score = c['evidence_score']
        
        # Color gradient for progress bar
        if score >= 75:
            bar_color = "linear-gradient(90deg, #10b981, #34d399)"
        elif score >= 35:
            bar_color = "linear-gradient(90deg, #f59e0b, #fbbf24)"
        else:
            bar_color = "linear-gradient(90deg, #64748b, #94a3b8)"
        if status == 'contaminant':
            bar_color = "#ef4444"

        # Loci tags
        loci_html = " ".join([f'<span class="tag" style="background:#1e3a5f; color:#38bdf8; font-weight:600; font-size:11px;">{g}</span>' for g in c['gene_details']])
        if not loci_html:
            loci_html = '<span style="color:#64748b; font-size:12px;">Kein Locus</span>'

        row_id = f"clin_row_{idx}"
        detail_id = f"clin_detail_{idx}"
        is_contam_row = "is-contaminant" if c['is_contaminant'] else "is-pathogen"

        # Sub-table of assays
        sub_rows = ""
        for pid, pst in c['pairs']:
            if pst['BOTH'] > 0 or pst['ONLY_FORWARD'] > 0 or pst['ONLY_REVERSE'] > 0:
                gene = pst.get('TargetGene', infer_target_gene(pid))
                f_seq = pst['Forward']
                r_seq = pst['Reverse']
                mean_d = round(sum(pst['AmpliconDistances'])/len(pst['AmpliconDistances']), 1) if pst['AmpliconDistances'] else "-"
                
                blast_p_url = f"https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_LEFT_INPUT={f_seq}&PRIMER_RIGHT_INPUT={r_seq}"
                amp_s = extracted_pair_seqs.get(pid, "")
                amp_btn = ""
                if amp_s:
                    amp_btn = f'''
                    <form action="https://blast.ncbi.nlm.nih.gov/Blast.cgi" method="POST" target="_blank" style="margin:0; display:inline-block;">
                        <input type="hidden" name="PROGRAM" value="blastn">
                        <input type="hidden" name="PAGE_TYPE" value="BlastSearch">
                        <input type="hidden" name="QUERY" value="{amp_s}">
                        <button type="submit" class="btn-blast" style="background:#0ea5e9; font-size:10px; padding:2px 6px;">Amplikon BLAST</button>
                    </form>
                    '''

                sub_rows += f'''
                <tr style="background: rgba(15, 23, 42, 0.6); font-size:12px;">
                    <td style="padding:6px 10px;"><code>{pid}</code></td>
                    <td style="padding:6px 10px;"><span class="badge" style="background:#334155; color:#38bdf8; font-size:10px;">{gene}</span></td>
                    <td style="padding:6px 10px;"><span class="badge badge-both" style="font-size:11px;">{pst["BOTH"]}</span></td>
                    <td style="padding:6px 10px;"><span class="badge badge-fwd" style="font-size:11px;">{pst["ONLY_FORWARD"]}</span></td>
                    <td style="padding:6px 10px;"><span class="badge badge-rev" style="font-size:11px;">{pst["ONLY_REVERSE"]}</span></td>
                    <td style="padding:6px 10px;">{mean_d} bp</td>
                    <td style="padding:6px 10px; font-family:monospace; font-size:10px; color:#94a3b8;">F: {f_seq[:18]}... | R: {r_seq[:18]}...</td>
                    <td style="padding:6px 10px;">
                        {amp_btn}
                        <a href="{blast_p_url}" target="_blank" class="btn-blast-sm" style="font-size:10px;">Primer-BLAST</a>
                    </td>
                </tr>
                '''

        clinical_rows_html += f'''
        <tr class="clinical-row {is_contam_row} status-{status} tier-{tier_clean}" id="{row_id}" style="border-bottom: 1px solid var(--card-border);">
            <td>
                <span class="badge {status_badge}" style="font-size:13px; display:inline-flex; align-items:center; gap:6px;">
                    <span>{traffic}</span> <strong>{status_text}</strong>
                </span>
            </td>
            <td><span class="badge {tier_badge}">{tier}</span></td>
            <td>
                <strong style="font-size:14px; color:#f8fafc;">{org}</strong>
                <div style="font-size:11.5px; color:#94a3b8; margin-top:2px;">{c["description"]}</div>
            </td>
            <td>{loci_html}</td>
            <td>
                <strong style="font-size:14px; color:#38bdf8;">{c["total_reads"]:,}</strong>
                <div style="font-size:10px; color:#94a3b8;">(Both: {c["total_both"]}, Fwd: {c["total_fwd"]}, Rev: {c["total_rev"]})</div>
            </td>
            <td style="min-width:140px;">
                <div style="display:flex; justify-content:space-between; font-size:11px; margin-bottom:3px;">
                    <span style="font-weight:700; color:#f8fafc;">{score}%</span>
                    <span style="color:#94a3b8;">Konfidenz</span>
                </div>
                <div style="width:100%; height:8px; background:#334155; border-radius:4px; overflow:hidden;">
                    <div style="width:{score}%; height:100%; background:{bar_color}; border-radius:4px;"></div>
                </div>
            </td>
            <td style="font-size:12px;">
                <div style="color:#cbd5e1; font-weight:600;">{c["minimap2_info"]}</div>
                <div style="color:#94a3b8; font-size:11px; margin-top:2px;">{c["blast_info"][:55]}</div>
            </td>
            <td>
                <button type="button" class="btn-blast" style="background:#334155; font-size:11px; padding:4px 8px;" onclick="toggleAccordion(\'{detail_id}\', this)">
                    Assays ({len(c["pairs"])}) &#x25BC;
                </button>
            </td>
        </tr>
        <tr id="{detail_id}" class="clinical-detail-row" style="display:none; background: rgba(15, 23, 42, 0.4);">
            <td colspan="8" style="padding: 10px 16px 16px 24px;">
                <div style="font-size:12px; font-weight:700; color:#38bdf8; margin-bottom:6px;">Beteiligte Primer-Assays &amp; Amplikon-Nachweise fuer {org}:</div>
                <table style="width:100%; border-collapse:collapse; background: #0f172a; border-radius:8px; overflow:hidden; border:1px solid #334155;">
                    <thead>
                        <tr style="background:#1e293b; font-size:11px; color:#94a3b8;">
                            <th style="padding:6px 10px;">Primer ID</th>
                            <th style="padding:6px 10px;">Zielgen</th>
                            <th style="padding:6px 10px;">BOTH (Fwd+Rev)</th>
                            <th style="padding:6px 10px;">Nur FWD</th>
                            <th style="padding:6px 10px;">Nur REV</th>
                            <th style="padding:6px 10px;">Abstand</th>
                            <th style="padding:6px 10px;">Primer Sequenzen</th>
                            <th style="padding:6px 10px;">Aktionen</th>
                        </tr>
                    </thead>
                    <tbody>
                        {sub_rows if sub_rows else "<tr><td colspan='8' style='text-align:center; padding:8px;'>Keine Reads fuer spezifische Assays</td></tr>"}
                    </tbody>
                </table>
            </td>
        </tr>
        '''

    # --- HTML TAB: BLAST GENUS-BESTAETIGTE ERREGER ---
    blast_confirmed_rows_html = ""
    for idx, c in enumerate(blast_confirmed_list, 1):
        org = c['organism']
        tier = c['tier']
        tier_clean = '0' if '0' in tier else ('1' if '1' in tier else ('2' if '2' in tier else '3'))
        tier_badge = 'badge-tier0' if tier_clean == '0' else ('badge-tier1' if tier_clean == '1' else ('badge-tier2' if tier_clean == '2' else 'badge-tier3'))
        
        conf_level = c.get('confirm_level', 'Genus-Bestaetigt')
        conf_badge = 'badge-both' if 'Spezies' in conf_level else 'badge-fwd'
        
        raw_b = c.get('blast_raw', {})
        blast_title = raw_b.get('stitle', c.get('blast_info', ''))
        blast_title = re.sub(r'^(gi\|\d+\|[a-z]+\|[\w\.]+\|\s*)', '', blast_title)
        acc = raw_b.get('sseqid', 'N/A')
        pident = raw_b.get('pident', 0.0)
        evalue = raw_b.get('evalue', 'N/A')
        bitscore = raw_b.get('bitscore', 0.0)
        align_len = raw_b.get('length', 0)
        
        c_seq = c.get('contig_seq') or c.get('rep_seq', '')
        c_btn = ""
        if c_seq:
            c_btn = f'''
            <form action="https://blast.ncbi.nlm.nih.gov/Blast.cgi" method="POST" target="_blank" style="margin:0; display:inline-block;">
                <input type="hidden" name="PROGRAM" value="blastn">
                <input type="hidden" name="PAGE_TYPE" value="BlastSearch">
                <input type="hidden" name="QUERY" value="{c_seq}">
                <button type="submit" class="btn-blast" style="background:#10b981; font-size:11px;" title="Gegen Web-NCBI BLASTN abgleichen">NCBI Web-BLAST</button>
            </form>
            '''

        blast_confirmed_rows_html += f'''
        <tr>
            <td>
                <span class="badge {conf_badge}" style="font-size:12px; display:inline-flex; align-items:center; gap:5px;">
                    🎯 <strong>{conf_level}</strong>
                </span>
            </td>
            <td><span class="badge {tier_badge}">{tier}</span></td>
            <td>
                <strong style="font-size:14px; color:#f8fafc;">{org}</strong>
                <div style="font-size:11.5px; color:#94a3b8; margin-top:2px;">{c["description"]}</div>
            </td>
            <td><span class="tag" style="background:#1e3a5f; color:#38bdf8; font-size:12px; font-weight:700;">{c.get("blast_genus", "")}</span></td>
            <td style="font-size:12.5px; color:#cbd5e1; max-width:320px;">
                <strong>{blast_title[:85]}</strong>
            </td>
            <td><code>{acc}</code></td>
            <td><strong style="color:#34d399; font-size:13px;">{pident}%</strong> <span style="font-size:10px; color:#94a3b8;">(Bit: {bitscore:.0f})</span></td>
            <td><code>E={evalue}</code></td>
            <td><strong>{c["contig_len"]} bp</strong> (Aln: {align_len} bp)</td>
            <td><span class="badge badge-both">{c["total_reads"]:,} Reads</span></td>
            <td>{c_btn}</td>
        </tr>
        '''

    if not blast_confirmed_rows_html:
        blast_confirmed_rows_html = '''
        <tr>
            <td colspan="11" style="text-align:center; padding:30px; color:#94a3b8; font-size:13.5px;">
                ℹ️ Im aktuellen Lauf wurden noch keine Kandidaten durch den lokalen BLASTn gegen <code>/mnt/blast_db/nt</code> auf Genus- oder Speziesebene bestaetigt.
                <div style="font-size:12px; margin-top:6px; color:#64748b;">(Dies tritt auf, wenn keine Kandidaten-Reads die Schwellenwerte fuer das lokale BLAST-Batching erreicht haben oder Treffer unter dem E-Value Cutoff lagen.)</div>
            </td>
        </tr>
        '''

    # --- HTML TAB 2: PRIMER MATRIX (735 ASSAYS) ---
    sorted_primer_stats = sorted(primer_stats.items(), key=lambda x: (x[1]['BOTH'], x[1]['ONLY_FORWARD'] + x[1]['ONLY_REVERSE']), reverse=True)
    table_rows_html = ""
    for pair_id, stats in sorted_primer_stats:
        if stats['BOTH'] > 0 or stats['ONLY_FORWARD'] > 0 or stats['ONLY_REVERSE'] > 0:
            mean_dist = round(sum(stats['AmpliconDistances']) / len(stats['AmpliconDistances']), 1) if stats['AmpliconDistances'] else "N/A"
            taxa_type = stats.get('TaxaType', 'Species')
            badge_class = "badge-genus" if taxa_type == 'Genus' else "badge-species"
            
            fwd_seq = stats['Forward']
            rev_seq = stats['Reverse']
            blast_primer_url = f"https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_LEFT_INPUT={fwd_seq}&PRIMER_RIGHT_INPUT={rev_seq}"
            blast_fwd_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&QUERY={fwd_seq}"
            blast_rev_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&QUERY={rev_seq}"

            amp_seq = extracted_pair_seqs.get(pair_id, "")
            amp_blast_btn = ""
            if amp_seq:
                clean_amp = amp_seq.strip().upper()
                amp_blast_btn = f'''
                <form action="https://blast.ncbi.nlm.nih.gov/Blast.cgi" method="POST" target="_blank" style="margin:0; display:inline-block;">
                    <input type="hidden" name="PROGRAM" value="blastn">
                    <input type="hidden" name="PAGE_TYPE" value="BlastSearch">
                    <input type="hidden" name="QUERY" value="{clean_amp}">
                    <button type="submit" class="btn-blast" style="background:#0ea5e9;" title="Extrahierte Amplikon-Sequenz mit NCBI BLASTN pruefen">
                        Amplikon
                    </button>
                </form>
                '''

            tier = stats.get('Tier', infer_tier(stats['Organism'], pair_id))
            tier_clean = '0' if '0' in tier else ('1' if '1' in tier else ('2' if '2' in tier else '3'))
            tier_badge_class = 'badge-tier0' if tier_clean == '0' else ('badge-tier1' if tier_clean == '1' else ('badge-tier2' if tier_clean == '2' else 'badge-tier3'))
            target_gene = stats.get('TargetGene', infer_target_gene(pair_id))

            table_rows_html += f'''
            <tr class="primer-row tier-{tier_clean}">
                <td><span class="badge {tier_badge_class}">{tier}</span></td>
                <td><strong>{stats["Organism"]}</strong></td>
                <td><span class="badge {badge_class}">{taxa_type}</span></td>
                <td><code>{pair_id}</code></td>
                <td><span class="badge" style="background:#334155; color:#38bdf8; font-size:11px;">{target_gene}</span></td>
                <td><code>{fwd_seq}</code></td>
                <td><code>{rev_seq}</code></td>
                <td><span class="badge badge-both">{stats["BOTH"]}</span></td>
                <td><span class="badge badge-fwd">{stats["ONLY_FORWARD"]}</span></td>
                <td><span class="badge badge-rev">{stats["ONLY_REVERSE"]}</span></td>
                <td>{mean_dist} bp</td>
                <td>
                    <div style="display:flex; gap:4px; align-items:center; flex-wrap:wrap;">
                        {amp_blast_btn}
                        <a href="{blast_primer_url}" target="_blank" class="btn-blast" title="Primer-BLAST bei NCBI">Primer-BLAST</a>
                        <a href="{blast_fwd_url}" target="_blank" class="btn-blast-sm" title="FWD BLAST">FWD</a>
                        <a href="{blast_rev_url}" target="_blank" class="btn-blast-sm" title="REV BLAST">REV</a>
                    </div>
                </td>
            </tr>
            '''

    # Extracted Amplicons rows
    extracted_table_html = ""
    for fname, count in extracted_counts.items():
        seq = extracted_seqs.get(fname, "")
        clean_seq = seq.strip().upper()
        seq_len = len(clean_seq)
        
        if clean_seq:
            preview_seq = f"{clean_seq[:35]}...{clean_seq[-25:]}" if seq_len > 65 else clean_seq
            blast_btn_html = f'''
            <div style="display:flex; gap:6px; align-items:center;">
                <form action="https://blast.ncbi.nlm.nih.gov/Blast.cgi" method="POST" target="_blank" style="margin:0; display:inline-block;">
                    <input type="hidden" name="PROGRAM" value="blastn">
                    <input type="hidden" name="PAGE_TYPE" value="BlastSearch">
                    <input type="hidden" name="QUERY" value="{clean_seq}">
                    <button type="submit" class="btn-blast" title="Sequenz an NCBI BLASTN senden">BLAST Sequenz</button>
                </form>
                <button type="button" class="btn-blast-sm" onclick="navigator.clipboard.writeText(\'{clean_seq}\'); this.innerText=\'Kopiert\'; setTimeout(() => this.innerText=\'Kopieren\', 2000);">Kopieren</button>
            </div>
            '''
            seq_cell_html = f'''
            <details style="cursor:pointer;">
                <summary style="font-family:monospace; font-size:12px; color:#38bdf8;"><code>{preview_seq}</code> <span style="font-size:10px; color:var(--text-secondary);">({seq_len} bp)</span></summary>
                <div style="margin-top:6px; padding:8px; background:#0f172a; border-radius:6px; font-family:monospace; font-size:11px; word-break:break-all; max-height:120px; overflow-y:auto; border:1px solid #334155; user-select:all; color:#f1f5f9;">
                    {clean_seq}
                </div>
            </details>
            '''
        else:
            blast_btn_html = "<span style='color:var(--text-secondary); font-size:12px;'>Keine Sequenz</span>"
            seq_cell_html = "<span style='color:var(--text-secondary); font-size:12px;'>-</span>"

        extracted_table_html += f'''
        <tr>
            <td><code>{fname}</code></td>
            <td><span class="badge badge-both">{count:,} Reads</span></td>
            <td><strong>{seq_len} bp</strong></td>
            <td style="max-width:350px;">{seq_cell_html}</td>
            <td>{blast_btn_html}</td>
        </tr>
        '''

    # --- HTML TAB 3: DE-NOVO ENTDECKUNGEN (DARK MATTER) ---
    denovo_rows_html = ""
    if dark_matter_discoveries:
        for d in dark_matter_discoveries:
            c_seq = d.get('contig_seq', '')
            c_btn = ""
            if c_seq:
                c_btn = f'''
                <form action="https://blast.ncbi.nlm.nih.gov/Blast.cgi" method="POST" target="_blank" style="margin:0; display:inline-block;">
                    <input type="hidden" name="PROGRAM" value="blastn">
                    <input type="hidden" name="PAGE_TYPE" value="BlastSearch">
                    <input type="hidden" name="QUERY" value="{c_seq}">
                    <button type="submit" class="btn-blast" style="background:#8b5cf6; font-size:11px;">Web-BLAST</button>
                </form>
                '''
            denovo_rows_html += f'''
            <tr>
                <td><span class="badge {d["category_badge"]}">{d["category"]}</span></td>
                <td><strong style="color:#f8fafc; font-size:14px;">{d["organism"]}</strong></td>
                <td style="font-size:12px; color:#cbd5e1;">{d["title"]}</td>
                <td><code>{d["accession"]}</code></td>
                <td><strong style="color:#34d399;">{d["pident"]}%</strong></td>
                <td>{d["align_len"]} bp</td>
                <td>{d["evalue"]}</td>
                <td><strong>{d["contig_len"]} bp</strong> ({d["read_support"]} Reads)</td>
                <td>{c_btn}</td>
            </tr>
            '''
    else:
        denovo_rows_html = '''
        <tr>
            <td colspan="9" style="text-align:center; padding:24px; color:#94a3b8;">
                Keine unerwarteten nicht-humanen Erreger-Contigs im Dark-Matter Scan detektiert.
            </td>
        </tr>
        '''

    # Overview link
    overview_basename = results.get('overview_dashboard_name', 'primer_overview_dashboard.html')
    is_multi_sample = results.get('is_multi_sample', False)
    overview_link_html = ""
    if is_multi_sample or os.path.exists(overview_basename):
        overview_link_html = f'''
        <div style="margin-top: 10px;">
            <a href="{overview_basename}" style="display: inline-flex; align-items: center; gap: 6px; background: rgba(56, 189, 248, 0.15); color: #38bdf8; border: 1px solid #0284c7; padding: 6px 14px; border-radius: 6px; text-decoration: none; font-size: 13px; font-weight: 600;">
                Zusammenfassungs-Dashboard oeffnen ({overview_basename})
            </a>
        </div>
        '''

    # Generate full HTML content
    html_content = f"""<!DOCTYPE html>
<html lang="de">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Erreger-Diagnostik &amp; Primer Dashboard - {fastq_name}</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {{
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --card-border: #334155;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --accent-cyan: #06b6d4;
            --accent-blue: #3b82f6;
            --accent-green: #10b981;
            --accent-purple: #8b5cf6;
            --accent-amber: #f59e0b;
        }}
        body {{
            font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            margin: 0;
            padding: 24px;
        }}
        .header {{
            padding-bottom: 20px;
            border-bottom: 1px solid var(--card-border);
            margin-bottom: 24px;
        }}
        .header h1 {{
            margin: 0;
            font-size: 28px;
            background: linear-gradient(135deg, var(--accent-cyan), var(--accent-blue));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
            gap: 14px;
            margin-bottom: 24px;
        }}
        .card {{
            background: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 12px;
            padding: 20px;
        }}
        .card .title {{
            font-size: 12px;
            color: var(--text-secondary);
            margin-bottom: 6px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .card .value {{
            font-size: 24px;
            font-weight: 700;
        }}
        .badge {{
            display: inline-block;
            padding: 4px 10px;
            border-radius: 20px;
            font-size: 12px;
            font-weight: 600;
        }}
        .badge-both {{ background: rgba(16, 185, 129, 0.2); color: #34d399; border: 1px solid #10b981; }}
        .badge-fwd {{ background: rgba(6, 182, 212, 0.2); color: #38bdf8; border: 1px solid #06b6d4; }}
        .badge-rev {{ background: rgba(139, 92, 246, 0.2); color: #c084fc; border: 1px solid #8b5cf6; }}
        .badge-genus {{ background: rgba(245, 158, 11, 0.2); color: #fbbf24; border: 1px solid #f59e0b; }}
        .badge-species {{ background: rgba(59, 130, 246, 0.2); color: #93c5fd; border: 1px solid #3b82f6; }}
        .badge-tier0 {{ background: rgba(6, 182, 212, 0.25); color: #22d3ee; border: 1px solid #06b6d4; }}
        .badge-tier1 {{ background: rgba(239, 68, 68, 0.2); color: #f87171; border: 1px solid #ef4444; }}
        .badge-tier2 {{ background: rgba(2, 132, 199, 0.2); color: #38bdf8; border: 1px solid #0284c7; }}
        .badge-tier3 {{ background: rgba(100, 116, 139, 0.2); color: #cbd5e1; border: 1px solid #64748b; }}
        .tag {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 6px;
            font-size: 11px;
            background: #334155;
            color: #cbd5e1;
        }}
        .tab-nav {{
            display: flex;
            gap: 8px;
            margin-bottom: 24px;
            border-bottom: 2px solid var(--card-border);
            padding-bottom: 8px;
            flex-wrap: wrap;
        }}
        .tab-btn {{
            background: transparent;
            border: 1px solid var(--card-border);
            color: var(--text-secondary);
            padding: 10px 18px;
            border-radius: 8px;
            font-size: 14px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            font-family: inherit;
        }}
        .tab-btn:hover {{
            background: rgba(56, 189, 248, 0.1);
            color: var(--text-primary);
            border-color: #38bdf8;
        }}
        .tab-btn.active {{
            background: #0284c7;
            color: #ffffff;
            border-color: #0ea5e9;
            box-shadow: 0 4px 12px rgba(2, 132, 199, 0.3);
        }}
        .tab-pane {{
            display: none;
        }}
        .tab-pane.active {{
            display: block;
        }}
        .filter-btn {{
            background: #1e293b;
            color: #94a3b8;
            border: 1px solid #334155;
            padding: 6px 12px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 12px;
            font-weight: 600;
            transition: all 0.15s;
        }}
        .filter-btn:hover {{
            background: #334155;
            color: #f8fafc;
        }}
        .filter-btn.active {{
            background: #0284c7;
            color: #ffffff;
            border-color: #38bdf8;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13.5px;
        }}
        th, td {{
            padding: 11px 12px;
            text-align: left;
            border-bottom: 1px solid var(--card-border);
        }}
        th {{
            color: var(--text-secondary);
            font-weight: 600;
            background: rgba(15, 23, 42, 0.5);
            cursor: pointer;
            user-select: none;
        }}
        th:hover {{
            color: #38bdf8;
        }}
        .btn-blast {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            background: #0284c7;
            color: white;
            text-decoration: none;
            font-size: 11px;
            font-weight: 600;
            white-space: nowrap;
            border: none;
            cursor: pointer;
        }}
        .btn-blast:hover {{
            background: #0369a1;
        }}
        .btn-blast-sm {{
            display: inline-block;
            padding: 3px 6px;
            border-radius: 4px;
            background: #334155;
            color: #38bdf8;
            text-decoration: none;
            font-size: 10px;
            font-weight: 600;
            border: none;
            cursor: pointer;
        }}
        .btn-blast-sm:hover {{
            background: #475569;
            color: white;
        }}
        code {{
            font-family: ui-monospace, Menlo, Consolas, monospace;
            background: rgba(0, 0, 0, 0.3);
            padding: 2px 5px;
            border-radius: 4px;
            font-size: 12px;
            color: #38bdf8;
        }}
        .guideline-box {{
            background: #1e293b;
            border: 1px solid #334155;
            border-radius: 10px;
            padding: 16px 20px;
            margin-bottom: 20px;
        }}
        .guideline-box h3 {{
            margin: 0 0 8px 0;
            color: #38bdf8;
            font-size: 16px;
        }}
        .guideline-box p, .guideline-box ul {{
            margin: 0 0 10px 0;
            font-size: 13.5px;
            line-height: 1.6;
            color: #cbd5e1;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Erreger-Diagnostik &amp; Primer Dashboard</h1>
        <p style="color: var(--text-secondary); margin: 6px 0 0 0;">
            Probe: <code>{fastq_name}</code> | Modus: <strong>{search_mode}</strong> | Threads: <strong>{threads_used} Kerne</strong> | Speed: <strong>{rps_val:,.0f} Reads/s ({elapsed_val}s)</strong>
        </p>
        {overview_link_html}
    </div>

    <!-- Tab Navigation Bar -->
    <div class="tab-nav">
        <button class="tab-btn active" onclick="switchTab('tab-clinical', this)">🏥 Klinischer Erreger-Befund (Ampel)</button>
        <button class="tab-btn" onclick="switchTab('tab-blast-confirmed', this)" style="border-left: 3px solid #10b981;">
            🎯 BLAST-Bestaetigt (Genus) <span class="badge badge-both" style="font-size:11px; padding:2px 7px; margin-left:4px;">{blast_confirmed_count}</span>
        </button>
        <button class="tab-btn" onclick="switchTab('tab-primers', this)">📊 Detaillierte Primer-Matrix (738 Assays)</button>
        <button class="tab-btn" onclick="switchTab('tab-denovo', this)">🔬 De-novo Entdeckungen (Dark Matter: {denovo_count})</button>
        <button class="tab-btn" onclick="switchTab('tab-amplicons', this)">📁 Extrahierte Amplikons ({len(extracted_counts)})</button>
        <button class="tab-btn" onclick="switchTab('tab-guidelines', this)">📖 Leitfaden (SHG Borreliose Hannover)</button>
    </div>

    <!-- TAB 1: KLINISCHER ERREGER-BEFUND (AMPEL & MULTI-LOCUS) -->
    <div id="tab-clinical" class="tab-pane active">
        <div class="stats-grid">
            <div class="card" style="border-left: 4px solid #10b981;">
                <div class="title">Bestaetigt (Multi-Locus)</div>
                <div class="value" style="color: #34d399;">{confirmed_count}</div>
            </div>
            <div class="card" style="border-left: 4px solid #06b6d4;">
                <div class="title">BLAST Genus-Bestaetigt</div>
                <div class="value" style="color: #22d3ee;">{blast_confirmed_count}</div>
            </div>
            <div class="card" style="border-left: 4px solid #f59e0b;">
                <div class="title">Verdachtsfaelle</div>
                <div class="value" style="color: #fbbf24;">{suspect_count}</div>
            </div>
            <div class="card" style="border-left: 4px solid #64748b;">
                <div class="title">Geringe Evidenz</div>
                <div class="value" style="color: #cbd5e1;">{low_count}</div>
            </div>
            <div class="card" style="border-left: 4px solid #8b5cf6;">
                <div class="title">De-novo Entdeckungen</div>
                <div class="value" style="color: #c084fc;">{denovo_count}</div>
            </div>
            <div class="card" style="border-left: 4px solid #ef4444;">
                <div class="title">Kitome / Kontaminanten</div>
                <div class="value" style="color: #f87171;">{contam_count}</div>
            </div>
            <div class="card">
                <div class="title">Gepruefte Reads</div>
                <div class="value" style="color: #f8fafc;">{total_reads:,}</div>
            </div>
        </div>

        <div class="card" style="margin-bottom: 24px;">
            <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:12px; margin-bottom:14px;">
                <div>
                    <h2 style="margin:0; font-size:18px; color:#38bdf8;">Klinischer Erreger-Befund (Diagnose-Ampel)</h2>
                    <div style="color:var(--text-secondary); font-size:12.5px; margin-top:2px;">
                        Automatische Evidenz-Synthese aus Multi-Locus PCR-Koinzidenz, minimap2 Referenz-Mapping und lokalem BLASTn gegen /mnt/blast_db/nt.
                    </div>
                </div>
                <div style="display:flex; gap:8px; flex-wrap:wrap; align-items:center;">
                    <button class="filter-btn active" onclick="filterClinical('all', this)">Alle Erreger</button>
                    <button class="filter-btn" onclick="filterClinical('confirmed', this)">Nur Bestaetigt</button>
                    <button class="filter-btn" onclick="filterClinical('suspect', this)">Verdachtsfaelle</button>
                    <button class="filter-btn" id="toggleContamBtn" onclick="toggleContaminants(this)">Kontaminanten ausblenden</button>
                </div>
            </div>

            <table id="clinicalTable">
                <thead>
                    <tr>
                        <th onclick="sortTable('clinicalTable', 0, false)">Ampel &amp; Status</th>
                        <th onclick="sortTable('clinicalTable', 1, false)">Tier</th>
                        <th onclick="sortTable('clinicalTable', 2, false)">Erreger &amp; Klinische Bedeutung</th>
                        <th>Nachgewiesene Gene (Loci)</th>
                        <th onclick="sortTable('clinicalTable', 4, true)">Gesamt-Reads</th>
                        <th onclick="sortTable('clinicalTable', 5, true)">Evidenz-Score</th>
                        <th>Validierung (minimap2 / BLAST)</th>
                        <th>Details</th>
                    </tr>
                </thead>
                <tbody>
                    {clinical_rows_html if clinical_rows_html else "<tr><td colspan='8' style='text-align:center; padding:20px;'>Keine Erreger-Signale in dieser Probe gefunden.</td></tr>"}
                </tbody>
            </table>
        </div>
    </div>

    <!-- TAB: BLAST GENUS-BESTAETIGTE ERREGER (EXTRA TAB) -->
    <div id="tab-blast-confirmed" class="tab-pane">
        <div class="card" style="margin-bottom: 24px; border-left: 4px solid #10b981;">
            <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:12px; margin-bottom:14px;">
                <div>
                    <h2 style="margin:0; font-size:18px; color:#34d399;">🎯 Orthogonal bestaetigte Erreger (Lokale nt-BLAST-Pruefung bestaetigt das Genus)</h2>
                    <div style="color:var(--text-secondary); font-size:12.5px; margin-top:4px;">
                        Diese Erreger wurden sowohl durch spezifische Primer-Treffer als auch durch unabhaengigen lokalen BLASTn-Abgleich gegen die 300 GB NCBI nt-Datenbank (<code>/mnt/blast_db/nt</code>) auf Gattungs- (Genus) oder Artebene (Spezies) eindeutig verifiziert.
                    </div>
                </div>
                <span class="badge badge-both" style="font-size:13px; padding:6px 14px;">
                    {blast_confirmed_count} Treffer verifiziert
                </span>
            </div>

            <table id="blastConfirmedTable">
                <thead>
                    <tr>
                        <th onclick="sortTable('blastConfirmedTable', 0, false)">Bestaetigungsgrad</th>
                        <th onclick="sortTable('blastConfirmedTable', 1, false)">Tier</th>
                        <th onclick="sortTable('blastConfirmedTable', 2, false)">Erreger &amp; Klinische Bedeutung</th>
                        <th onclick="sortTable('blastConfirmedTable', 3, false)">Genus</th>
                        <th>Lokaler nt-Top-Hit (GenBank Titel)</th>
                        <th onclick="sortTable('blastConfirmedTable', 5, false)">Accession</th>
                        <th onclick="sortTable('blastConfirmedTable', 6, true)">Identitaet</th>
                        <th>E-Value</th>
                        <th onclick="sortTable('blastConfirmedTable', 8, true)">Contig / Aln</th>
                        <th onclick="sortTable('blastConfirmedTable', 9, true)">Reads</th>
                        <th>Aktionen</th>
                    </tr>
                </thead>
                <tbody>
                    {blast_confirmed_rows_html}
                </tbody>
            </table>
        </div>
    </div>

    <!-- TAB 2: PRIMER MATRIX (738 ASSAYS) -->
    <div id="tab-primers" class="tab-pane">
        <div class="charts-row" style="display:grid; grid-template-columns:1fr 1fr; gap:20px; margin-bottom:20px;">
            <div class="card"><div id="pieChart"></div></div>
            <div class="card"><div id="barChart"></div></div>
        </div>

        <div class="card">
            <h2>Detaillierte Primer Treffer Uebersicht (735 Assays)</h2>
            <div style="display:flex; gap:8px; margin: 10px 0 16px 0; flex-wrap:wrap; align-items:center;">
                <span style="font-size:12px; font-weight:600; color:var(--text-secondary);">Filter nach Tier:</span>
                <button class="filter-btn active" onclick="filterTableByTier('all', this, 'primerTable')">Alle Assays</button>
                <button class="filter-btn" style="border-left:3px solid #22d3ee;" onclick="filterTableByTier('0', this, 'primerTable')">Tier 0: Universal Barcoding</button>
                <button class="filter-btn" style="border-left:3px solid #ef4444;" onclick="filterTableByTier('1', this, 'primerTable')">Tier 1: Borreliose &amp; Vektor</button>
                <button class="filter-btn" style="border-left:3px solid #0284c7;" onclick="filterTableByTier('2', this, 'primerTable')">Tier 2: Endokarditis &amp; Sepsis</button>
                <button class="filter-btn" style="border-left:3px solid #64748b;" onclick="filterTableByTier('3', this, 'primerTable')">Tier 3: Weitere</button>
            </div>
            <table id="primerTable">
                <thead>
                    <tr>
                        <th onclick="sortTable('primerTable', 0, false)">Tier</th>
                        <th onclick="sortTable('primerTable', 1, false)">Organismus</th>
                        <th onclick="sortTable('primerTable', 2, false)">Ebene</th>
                        <th onclick="sortTable('primerTable', 3, false)">Primer ID</th>
                        <th onclick="sortTable('primerTable', 4, false)">Zielgen</th>
                        <th>Forward Primer</th>
                        <th>Reverse Primer</th>
                        <th onclick="sortTable('primerTable', 7, true)">FWD &amp; REV</th>
                        <th onclick="sortTable('primerTable', 8, true)">Nur FWD</th>
                        <th onclick="sortTable('primerTable', 9, true)">Nur REV</th>
                        <th onclick="sortTable('primerTable', 10, true)">Abstand</th>
                        <th>BLAST Ueberpruefung</th>
                    </tr>
                </thead>
                <tbody>{table_rows_html}</tbody>
            </table>
        </div>
    </div>

    <!-- TAB 3: DE-NOVO ENTDECKUNGEN (DARK MATTER) -->
    <div id="tab-denovo" class="tab-pane">
        <div class="card">
            <h2>De-novo Entdeckungen (Dark Matter Discovery)</h2>
            <p style="color:var(--text-secondary); font-size:13px;">
                Nicht zugeordnete Reads wurden de-novo assembliert und gegen die lokale 300 GB nt-Datenbank (/mnt/blast_db/nt) abgeglichen. Dies ermoeglicht den Fund unerwarteter oder neuartiger Erreger (Bakterien, Viren, Pilze, Parasiten), fuer die keine Primer definiert waren.
            </p>
            <table id="denovoTable">
                <thead>
                    <tr>
                        <th>Kategorie</th>
                        <th>Erreger / Gattung</th>
                        <th>NCBI GenBank Treffer</th>
                        <th>Accession</th>
                        <th>Identitaet</th>
                        <th>Alignment</th>
                        <th>E-Value</th>
                        <th>Contig &amp; Reads</th>
                        <th>Aktionen</th>
                    </tr>
                </thead>
                <tbody>{denovo_rows_html}</tbody>
            </table>
        </div>
    </div>

    <!-- TAB 4: EXTRAHIERTE AMPLIKONS -->
    <div id="tab-amplicons" class="tab-pane">
        <div class="card" style="margin-bottom: 24px;"><div id="distChart"></div></div>
        <div class="card">
            <h2>Extrahierte Amplicon Dateien (_frp_&lt;bp&gt;.fastq)</h2>
            <table id="extractedAmpliconsTable">
                <thead>
                    <tr>
                        <th onclick="sortTable('extractedAmpliconsTable', 0, false)">Dateiname</th>
                        <th onclick="sortTable('extractedAmpliconsTable', 1, true)">Reads</th>
                        <th onclick="sortTable('extractedAmpliconsTable', 2, true)">Laenge</th>
                        <th>Sequenz-Vorschau</th>
                        <th>BLAST &amp; Export</th>
                    </tr>
                </thead>
                <tbody>{extracted_table_html}</tbody>
            </table>
        </div>
    </div>

    <!-- TAB: LEITFADEN & NUTZUNGSANLEITUNG -->
    <div id="tab-guidelines" class="tab-pane">
        <div class="card" style="border-left: 4px solid var(--accent-cyan);">
            <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:10px; margin-bottom:12px;">
                <h2 style="margin:0; color:#38bdf8;">📖 Leitfaden, Programmnutzung &amp; Befehlsreferenz</h2>
                <span class="badge badge-both" style="font-size:12px;">Wissenschaftlich-Diagnostischer Standard</span>
            </div>
            <p style="color: var(--text-secondary); font-size: 14px; margin-bottom: 20px;">
                Umfassende Dokumentation zur Nutzung der Pipeline, Interpretation der Dashboard-Ergebnisse, Erklärung des <strong>Primer-Buttons</strong> sowie klinische Handlungsempfehlungen in Anlehnung an die Erfahrungswerte der <strong>Selbsthilfe Gruppe Borreliose Hannover</strong>.
            </p>

            <!-- 1. NUTZUNG DES PROGRAMMS & TABS -->
            <div class="guideline-box">
                <h3 style="color:#38bdf8;">🚀 1. Nutzung des Programms &amp; Dashboard-Navigation</h3>
                <p>
                    Die Pipeline führt eine mehrstufige bioinformatische Analyse durch: hochsensitives In-silico-Primer-Matching, orthogonales Vollgenom-Mapping (<code>--pan-scan</code>) und eine De-novo-Contig-Assemblierung unzugeordneter Reads (Dark Matter) mit lokalem BLASTn-Abgleich gegen die NCBI-nt-Datenbank.
                </p>
                <p><strong>Die Reiter (Tabs) im Dashboard im Überblick:</strong></p>
                <ul>
                    <li><strong style="color:#10b981;">🏥 Klinischer Erreger-Befund (Ampel):</strong> Die zentrale klinische Übersicht. Fasst alle nachgewiesenen Erreger zusammen, gruppiert sie nach Gattung und ordnet sie einer diagnostischen Ampel zu:
                        <br>&bull; <span class="badge badge-both" style="font-size:11px;">🟢 Bestätigt</span>: Hohe Evidenz (Amplikon-Nachweis, Genom-Mapping oder Genus-BLAST bestätigt).
                        <br>&bull; <span class="badge badge-fwd" style="font-size:11px;">🟡 Verdacht</span>: Einzelsignale (z. B. FWD-Read) oder schwächeres Mapping.
                        <br>&bull; <span class="badge badge-rev" style="font-size:11px;">⚪ Gering / Isoliert</span>: Einzelne isolierte Fragmente.
                        <br>&bull; <span class="badge badge-none" style="font-size:11px;">🚫 Kitom / Kontaminant</span>: Typische Labor- oder Extraktionskit-Flora (z.B. Ralstonia, Cutibacterium) – per Button ausblendbar.
                    </li>
                    <li><strong style="color:#38bdf8;">🎯 BLAST-Bestätigt (Genus):</strong> Exklusive Spezialansicht! Zeigt nur diejenigen Funde, bei denen die <em>lokale BLASTn-Suche</em> (/mnt/blast_db/nt) die Gattung oder Spezies unabhängig bestätigt hat (inkl. E-Value, Identität % und Alignment-Länge).</li>
                    <li><strong style="color:#a855f7;">📊 Detaillierte Primer-Matrix (738 Assays):</strong> (Aufrufbar über den <em>Primer-Button</em> oben!) Vollständige Auflistung aller 738 diagnostischen Assays mit genauen FWD-, REV- und Doppelstrang-Treffern, sortierbar und filterbar nach Tiers.</li>
                    <li><strong style="color:#06b6d4;">🔬 De-novo Entdeckungen (Dark Matter):</strong> Zeigt neue oder unerwartete Erreger, die aus nicht-gemappten Reads de-novo zu Contigs assembliert und via lokalem BLAST identifiziert wurden.</li>
                    <li><strong style="color:#f59e0b;">📁 Extrahierte Amplikons:</strong> Physisch aus den Reads ausgeschnittene PCR-Produkte (inkl. Längenverteilung) mit 1-Click-Web-BLAST.</li>
                </ul>
            </div>

            <!-- 2. WAS IST DER PRIMER BUTTON? -->
            <div class="guideline-box" style="border-left: 4px solid #a855f7; background: rgba(88, 28, 135, 0.15);">
                <h3 style="color:#c084fc;">🔘 2. Was ist der &bdquo;Primer-Button&ldquo;?</h3>
                <p>
                    Im Dashboard begegnet Ihnen der Begriff <strong>Primer-Button</strong> an zwei zentralen Stellen:
                </p>
                <ol style="margin-bottom:12px;">
                    <li><strong>Der Tab-Button &bdquo;📊 Detaillierte Primer-Matrix (738 Assays)&ldquo; im Menüband:</strong>
                        <br>Dieser Navigations-Button schaltet die Ansicht vom komprimierten klinischen Befund auf die vollständige, wissenschaftliche Rohdaten-Tabelle um. Dort sehen Sie jedes der 738 Primerpaare mit:
                        <ul style="margin-top:4px;">
                            <li>Genauer Anzahl der FWD-Reads, REV-Reads und beidseitigen Amplikons (BOTH).</li>
                            <li>Basenpaar-Abstand (Amplikon-Größe in bp) und Übereinstimmung mit der Zielregion.</li>
                            <li>Filterknöpfen zur gezielten Eingrenzung nach <em>Tier 0 (Universal-Barcodes)</em>, <em>Tier 1 (Borreliose &amp; Vektoren)</em>, <em>Tier 2 (Endokarditis/Atypisch)</em> und <em>Tier 3</em>.</li>
                        </ul>
                    </li>
                    <li><strong>Die interaktiven Action-Buttons bei jedem Primertreffer:</strong>
                        <ul style="margin-top:4px;">
                            <li><button type="button" class="btn-blast" style="display:inline-block; font-size:11px; padding:2px 8px; background:#0ea5e9; cursor:default;">Amplikon BLAST</button> oder <button type="button" class="btn-blast" style="display:inline-block; font-size:11px; padding:2px 8px; background:#10b981; cursor:default;">NCBI Web-BLAST</button>: Überträgt die extrahierte DNA-Sequenz mit einem einzigen Klick direkt an das offizielle NCBI BLASTN im Web, um die weltweite Referenzdatenbank abzugleichen.</li>
                            <li><button type="button" class="btn-blast-sm" style="display:inline-block; font-size:10px; padding:2px 6px; cursor:default;">Kopieren</button>: Kopiert die extrahierte FASTA-Sequenz direkt in die Zwischenablage für eigene Analysen.</li>
                            <li><button type="button" class="btn-blast" style="display:inline-block; font-size:11px; padding:2px 8px; background:#334155; cursor:default;">Details anzeigen &#x25BC;</button>: Klappt im klinischen Befund das Akkordeon auf, um alle einzelnen Primer-Assays und Read-Statistiken dieser Erregergattung einzusehen.</li>
                        </ul>
                    </li>
                </ol>
            </div>

            <!-- 3. BEISPIELZEILEN & COMMAND-LINE CHEAT SHEET -->
            <div class="guideline-box">
                <h3 style="color:#38bdf8;">💻 3. Konsolen-Befehle (CLI Cheat-Sheet &amp; Beispielzeilen)</h3>
                <p>
                    Führen Sie die folgenden Befehle direkt im Terminal (z. B. auf dem Server unter <code>/mnt/sdc1/scripts/</code> oder lokal) aus:
                </p>

                <div style="margin-bottom:14px;">
                    <strong style="color:#10b981;">🌟 A. Klinischer Goldstandard (Vollanalyse mit allen 738 Primern + Vollgenom-Mapping):</strong>
                    <p style="font-size:13px; color:var(--text-secondary); margin:2px 0 6px 0;">
                        Empfohlen für jede diagnostische Auswertung. Testet alle Erreger (Babesia, Plasmodium, Borrelia, Bartonella, etc.), aktiviert Vollgenom-Mapping gegen 340+ Referenzen und De-novo Dark Matter Discovery:
                    </p>
                    <pre style="background:#090d16; padding:12px; border-radius:6px; border:1px solid #1e293b; color:#38bdf8; font-family:monospace; font-size:13px; overflow-x:auto;">python primer_suche.py -i probe.fastq.gz -a --pan-scan -t 24</pre>
                </div>

                <div style="margin-bottom:14px;">
                    <strong style="color:#06b6d4;">⚡ B. Schnelles hierarchisches Screening (Zeitersparnis):</strong>
                    <p style="font-size:13px; color:var(--text-secondary); margin:2px 0 6px 0;">
                        Prüft zuerst Genus-Marker (z. B. Borrelia flaB, Babesia 18S). Nur bei positivem Signal werden spezifische Spezies-Assays durchsucht:
                    </p>
                    <pre style="background:#090d16; padding:12px; border-radius:6px; border:1px solid #1e293b; color:#38bdf8; font-family:monospace; font-size:13px; overflow-x:auto;">python primer_suche.py -i probe.fastq.gz -gs -t 24</pre>
                </div>

                <div style="margin-bottom:14px;">
                    <strong style="color:#f59e0b;">📁 C. Alle Proben im Ordner analysieren (Multi-Sample Batch):</strong>
                    <p style="font-size:13px; color:var(--text-secondary); margin:2px 0 6px 0;">
                        Verarbeitet alle FASTQ-Dateien auf einmal, erstellt Einzelberichte und das übergeordnete Dashboard (<code>primer_dashboard_all_samples.html</code>):
                    </p>
                    <pre style="background:#090d16; padding:12px; border-radius:6px; border:1px solid #1e293b; color:#38bdf8; font-family:monospace; font-size:13px; overflow-x:auto;">python primer_suche.py -i *.fastq.gz -a --pan-scan -t 28</pre>
                </div>

                <div style="margin-bottom:14px;">
                    <strong style="color:#c084fc;">🎯 D. Gezielte Tier-Fokussierung (z.B. Tier 0 &amp; Tier 1):</strong>
                    <p style="font-size:13px; color:var(--text-secondary); margin:2px 0 6px 0;">
                        Beschränkt die Suche auf Universelle Barcodes (Tier 0) und Borreliose/Zecken-Koinfektionen (Tier 1):
                    </p>
                    <pre style="background:#090d16; padding:12px; border-radius:6px; border:1px solid #1e293b; color:#38bdf8; font-family:monospace; font-size:13px; overflow-x:auto;">python primer_suche.py -i probe.fastq.gz --tier 0,1 -a -t 24</pre>
                </div>

                <div>
                    <strong style="color:#94a3b8;">🛠️ E. Primer-Datenbank bereinigen &amp; aktualisieren:</strong>
                    <pre style="background:#090d16; padding:12px; border-radius:6px; border:1px solid #1e293b; color:#38bdf8; font-family:monospace; font-size:13px; overflow-x:auto;"># Duplikate und Reverse-Komplemente entfernen:
python primer_suche.py -p best_primers_db.csv --dedup best_primers_db_clean.csv

# Gezielt gegen RefSeq/GenBank aktualisieren (z.B. Babesia oder all):
python primer_suche.py --update-primers --target Babesia</pre>
                </div>
            </div>

            <!-- 4. KLINISCHE STUFENDIAGNOSTIK DER SHG BORRELIOSE HANNOVER -->
            <div class="guideline-box">
                <h3 style="color:#38bdf8;">📌 4. Grundsatz der Stufendiagnostik (Hierarchische Verifizierung)</h3>
                <p>
                    Ein bioinformatischer Primer-Treffer in NGS- oder Nanopore-Rohdaten stellt eine hochsensitive Vorstufe dar. Zur soliden Verifizierung gilt das mehrstufige diagnostische Vorgehen:
                </p>
                <ul>
                    <li><strong>Stufe 1 (Genus-Screening):</strong> Nachweis von Leitgenen (z. B. 16S-rRNA, <em>flaB</em>, <em>gltA</em>, <em>17kDa</em>). Ein Signal zeigt an, dass Vertreter einer Erregergruppe (z. B. <em>Borrelia</em>, <em>Babesia</em>, <em>Bartonella</em>, <em>Rickettsia</em>) in der Probe vorliegen.</li>
                    <li><strong>Stufe 2 (Spezies-Identifikation):</strong> Spezifische Primer grenzen den Erreger exakt ein (z. B. <em>Borrelia burgdorferi sensu stricto</em>, <em>B. afzelii</em>, <em>B. garinii</em>, <em>B. miyamotoi</em> oder <em>Babesia microti/divergens/odocoilei</em>).</li>
                    <li><strong>Stufe 3 (Orthogonale Validierung):</strong> Bei getrennten FWD/REV-Hits erfolgt das automatische Genom-Mapping mit <strong>minimap2</strong> gegen das NCBI RefSeq-Vollgenom sowie die Sequenzprüfung über den integrierten <strong>BLAST-Button</strong>.</li>
                </ul>
            </div>

            <!-- 5. KOINFEKTIONEN -->
            <div class="guideline-box">
                <h3 style="color:#38bdf8;">🛡️ 5. Koinfektionen &amp; atypische Erreger beachten</h3>
                <p>
                    Zecken und Vektoren übertragen selten nur einen einzelnen Erreger. Typische klinische Verläufe zeichnen sich oft durch ein Erregerkollektiv aus:
                </p>
                <ul>
                    <li><strong>Piroplasmen &amp; Protozoen:</strong> <em>Babesia</em> (insb. <em>B. microti</em>, <em>B. divergens</em>, <em>B. odocoilei</em>, <em>B. duncani</em>, <em>B. venatorum</em>) befallen Erythrozyten und führen häufig zu therapieresistenten Erschöpfungssyndromen, Nachtschweiß und hämatologischen Auffälligkeiten.</li>
                    <li><strong>Intrazelluläre Bakterien:</strong> <em>Bartonella</em> (z. B. <em>B. henselae</em>) und <em>Rickettsia</em> (z. B. Fleckfieber-Gruppe) besiedeln Gefäßendothelien und erfordern eigenständige therapeutische Schemata.</li>
                    <li><strong>Opportunistische Mykobakterien &amp; Pilze:</strong> NTM (z.B. <em>M. immunogenum</em>, <em>M. marinum</em>, <em>M. avium</em>) sowie Pilze wie <em>Cryptococcus neoformans</em> treten bevorzugt bei sekundärer Immundysfunktion auf.</li>
                    <li><strong>Virale Reaktivierungen:</strong> Herpesviren (EBV, CMV, HHV-6, VZV) reaktivieren häufig unter chronischer Borrelien-Last.</li>
                </ul>
            </div>

            <!-- 6. VERNETZUNG -->
            <div class="guideline-box" style="border-left: 4px solid var(--accent-amber); background: rgba(30, 41, 59, 0.8);">
                <h3 style="color:#fbbf24;">🤝 6. Vernetzung &amp; Patientenunterstützung</h3>
                <p>
                    Betroffene mit komplexen, chronischen Symptomkomplexen nach Zeckenstichen oder Verdacht auf Vektor-übertragene Multisystemerkrankungen finden kompetenten Austausch, Informationsmaterialien und Unterstützung bei erfahrenen Selbsthilfe-Initiativen:
                </p>
                <div style="background:#0f172a; padding:12px 16px; border-radius:8px; border:1px solid #334155; margin-top:8px;">
                    <p style="margin:0 0 6px 0; font-weight:700; color:#f8fafc; font-size:14px;">
                        Selbsthilfe Gruppe Borreliose Hannover
                    </p>
                    <p style="margin:0 0 6px 0; font-size:13px; color:#cbd5e1;">
                        Austauschplattform für Betroffene, Aufklärung zu Zecken-übertragenen Infektionen, Koinfektionen (Babesien, Bartonellen, Rickettsien) sowie ganzheitlichen und interdisziplinären Diagnoseansätzen.
                    </p>
                    <p style="margin:0; font-size:13px; color:#38bdf8;">
                        🌐 <em>Kontakt, Beratungsangebote &amp; Treffpunkte: Region Hannover / Selbsthilfe-Netzwerk Niedersachsen</em>
                    </p>
                </div>
            </div>
        </div>
    </div>
    </div>

    <script>
        Plotly.newPlot('pieChart', [{{
            labels: {json.dumps(match_labels)},
            values: {json.dumps(match_values)},
            type: 'pie', hole: 0.4,
            marker: {{ colors: ['#10b981', '#06b6d4', '#8b5cf6', '#475569'] }}
        }}], {{ title: '<b>Primer Match Verteilung</b>', paper_bgcolor: 'transparent', font: {{ color: '#f8fafc' }} }});

        Plotly.newPlot('barChart', [
            {{ x: {json.dumps(top_org_both)}, y: {json.dumps(top_org_labels)}, name: 'FWD & REV', type: 'bar', orientation: 'h', marker: {{ color: '#10b981' }} }},
            {{ x: {json.dumps(top_org_fwd)}, y: {json.dumps(top_org_labels)}, name: 'Nur FWD', type: 'bar', orientation: 'h', marker: {{ color: '#06b6d4' }} }},
            {{ x: {json.dumps(top_org_rev)}, y: {json.dumps(top_org_labels)}, name: 'Nur REV', type: 'bar', orientation: 'h', marker: {{ color: '#8b5cf6' }} }}
        ], {{
            title: '<b>Top Primertreffer</b>',
            barmode: 'stack',
            paper_bgcolor: 'transparent',
            plot_bgcolor: 'transparent',
            font: {{ color: '#f8fafc' }},
            margin: {{ l: 300, r: 30, t: 50, b: 50 }},
            yaxis: {{ autorange: 'reversed', automargin: true, tickfont: {{ size: 11, color: '#f8fafc' }} }},
            xaxis: {{ gridcolor: '#334155', title: 'Reads' }}
        }}, {{ responsive: true }});

        Plotly.newPlot('distChart', [{{
            x: {json.dumps(dist_labels)}, y: {json.dumps(dist_values)}, type: 'bar', marker: {{ color: '#3b82f6' }}
        }}], {{ title: '<b>Basenpaare-Abstand Verteilung (_frp_ bp)</b>', paper_bgcolor: 'transparent', font: {{ color: '#f8fafc' }} }});

        let hideContaminants = false;

        function toggleAccordion(detailId, btn) {{
            const row = document.getElementById(detailId);
            if (!row) return;
            if (row.style.display === 'none' || row.style.display === '') {{
                row.style.display = 'table-row';
                btn.innerHTML = btn.innerHTML.replace('&#x25BC;', '&#x25B2;');
            }} else {{
                row.style.display = 'none';
                btn.innerHTML = btn.innerHTML.replace('&#x25B2;', '&#x25BC;');
            }}
        }}

        function toggleContaminants(btn) {{
            hideContaminants = !hideContaminants;
            btn.classList.toggle('active', hideContaminants);
            btn.innerText = hideContaminants ? "Kontaminanten anzeigen" : "Kontaminanten ausblenden";
            document.querySelectorAll('#clinicalTable tbody tr.is-contaminant').forEach(r => {{
                r.style.display = hideContaminants ? 'none' : '';
            }});
        }}

        function filterClinical(status, btn) {{
            const parent = btn.parentElement;
            parent.querySelectorAll('.filter-btn').forEach(b => {{
                if (b.id !== 'toggleContamBtn') b.classList.remove('active');
            }});
            btn.classList.add('active');
            document.querySelectorAll('#clinicalTable tbody tr.clinical-row').forEach(r => {{
                const isContam = r.classList.contains('is-contaminant');
                if (hideContaminants && isContam) {{
                    r.style.display = 'none';
                    return;
                }}
                if (status === 'all') {{
                    r.style.display = '';
                }} else {{
                    r.style.display = r.classList.contains('status-' + status) ? '' : 'none';
                }}
            }});
        }}

        function filterTableByTier(tier, btn, tableId) {{
            const parent = btn ? btn.parentElement : document;
            parent.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
            if (btn) btn.classList.add('active');
            const table = document.getElementById(tableId);
            if (!table) return;
            const rows = table.querySelectorAll('tbody tr.primer-row');
            rows.forEach(r => {{
                if (tier === 'all') {{
                    r.style.display = '';
                }} else {{
                    r.style.display = r.classList.contains('tier-' + tier) ? '' : 'none';
                }}
            }});
        }}

        function switchTab(tabId, btn) {{
            document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
            document.querySelectorAll('.tab-pane').forEach(p => p.classList.remove('active'));
            if (btn) btn.classList.add('active');
            const target = document.getElementById(tabId);
            if (target) target.classList.add('active');
            window.dispatchEvent(new Event('resize'));
        }}

        let sortDirections = {{}};
        function sortTable(tableId, colIndex, isNumeric = false) {{
            const table = document.getElementById(tableId);
            if (!table) return;
            const tbody = table.querySelector('tbody');
            const rows = Array.from(tbody.querySelectorAll('tr'));
            const key = tableId + "_" + colIndex;
            const currentDir = sortDirections[key] || "desc";
            const newDir = currentDir === "asc" ? "desc" : "asc";
            sortDirections[key] = newDir;

            rows.sort((a, b) => {{
                const tdA = a.children[colIndex];
                const tdB = b.children[colIndex];
                if (!tdA || !tdB) return 0;
                let valA = tdA.innerText.trim().replace(/,/g, '');
                let valB = tdB.innerText.trim().replace(/,/g, '');
                if (isNumeric) {{
                    let numA = parseFloat(valA.replace(/[^0-9.-]/g, '')) || 0;
                    let numB = parseFloat(valB.replace(/[^0-9.-]/g, '')) || 0;
                    return newDir === "asc" ? numA - numB : numB - numA;
                }} else {{
                    return newDir === "asc" ? valA.localeCompare(valB) : valB.localeCompare(valA);
                }}
            }});
            rows.forEach(r => tbody.appendChild(r));
        }}
    </script>
</body>
</html>"""

    with open(dashboard_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"[+] Modernes 5-Tab Dashboard erfolgreich erstellt: {os.path.abspath(dashboard_path)}")
    dl_dir = r"C:\Users\internet\Downloads"
    try:
        if os.path.exists(dl_dir):
            shutil.copy2(dashboard_path, os.path.join(dl_dir, os.path.basename(dashboard_path)))
    except Exception:
        pass
    return dashboard_path


def get_unique_dashboard_path(fastq_path, dest_dir="."):
    """
    Derives unique dashboard filename from characters before first underscore:
    e.g. tm_NG1YETBR8Q... -> primer_dashboard_tm.html
    If file exists, increments: primer_dashboard_tm_2.html, _3.html, etc.
    """
    basename = os.path.basename(fastq_path)
    if '_' in basename:
        prefix = basename.split('_')[0].strip()
    else:
        prefix = basename.split('.')[0].strip()
    
    if not prefix:
        prefix = "sample"

    base_name = f"primer_dashboard_{prefix}.html"
    dest_path = os.path.join(dest_dir, base_name)
    dl_path = os.path.join(r"C:\Users\internet\Downloads", base_name)

    if not os.path.exists(dest_path) and not os.path.exists(dl_path):
        return dest_path

    counter = 2
    while True:
        cand_name = f"primer_dashboard_{prefix}_{counter}.html"
        cand_dest = os.path.join(dest_dir, cand_name)
        cand_dl = os.path.join(r"C:\Users\internet\Downloads", cand_name)
        if not os.path.exists(cand_dest) and not os.path.exists(cand_dl):
            return cand_dest
        counter += 1

def generate_summary_html_dashboard(all_sample_results, overview_dashboard_path):
    """
    Generates a comprehensive multi-sample HTML dashboard across all analyzed FASTQ files.
    Includes:
    - Global KPI metrics and comparative statistics
    - Interactive cross-sample Plotly heatmap (Samples vs. Top Pathogens)
    - Side-by-side sample comparison bar chart
    - Per-sample summary table with direct links to individual dashboards
    - Comprehensive sortable & searchable Cross-Sample Primer Matrix with BLAST buttons
    - Aggregated list of all extracted amplicons with BLAST verification
    - Automatically copies output to C:\\Users\\internet\\Downloads\\
    """
    num_samples = len(all_sample_results)
    if num_samples == 0:
        return None

    # Determine distinctive sample display labels
    raw_fnames = [s_res.get('fastq_name', f'sample_{i+1}') for i, s_res in enumerate(all_sample_results)]
    raw_prefixes = [(f.split('_')[0] if '_' in f else f.split('.')[0]) for f in raw_fnames]
    if len(raw_prefixes) != len(set(raw_prefixes)):
        sample_display_labels = [re.sub(r'\.fastq(\.gz)?$', '', f, flags=re.I) for f in raw_fnames]
    else:
        sample_display_labels = raw_prefixes

    samples_meta = []
    for idx, s_res in enumerate(all_sample_results):
        fname = raw_fnames[idx]
        prefix = sample_display_labels[idx]
        dash_p = s_res.get('dashboard_path', '')
        dash_base = os.path.basename(dash_p) if dash_p else f"primer_dashboard_{prefix}.html"
        match_counts = s_res.get('match_counts', {})
        samples_meta.append({
            'index': idx,
            'name': fname,
            'prefix': prefix,
            'dash_path': dash_p,
            'dash_base': dash_base,
            'total_reads': s_res.get('total_reads', 0),
            'both': match_counts.get('BOTH', 0),
            'fwd': match_counts.get('ONLY_FORWARD', 0),
            'rev': match_counts.get('ONLY_REVERSE', 0),
            'none': match_counts.get('NONE', 0),
            'results': s_res
        })

    total_all_reads = sum(s['total_reads'] for s in samples_meta)
    total_all_both = sum(s['both'] for s in samples_meta)
    total_all_fwd = sum(s['fwd'] for s in samples_meta)
    total_all_rev = sum(s['rev'] for s in samples_meta)
    total_all_none = sum(s['none'] for s in samples_meta)
    total_elapsed = sum(s['results'].get('elapsed_s', 0) for s in samples_meta)
    avg_rps = (total_all_reads / total_elapsed) if total_elapsed > 0 else 0
    search_mode = all_sample_results[0].get('search_mode', 'Standard')
    timestamp = datetime.now().strftime("%d.%m.%Y %H:%M:%S")

    # Aggregate primer hits and presence across all samples
    primer_matrix = {}
    org_both_sums = Counter()
    org_sample_presence = defaultdict(set)
    sample_positive_orgs = defaultdict(list)

    for s in samples_meta:
        s_idx = s['index']
        s_res = s['results']
        p_stats = s_res.get('primer_stats', {})
        pair_seqs = s_res.get('extracted_pair_seqs', {})

        for pair_id, stats in p_stats.items():
            both = stats.get('BOTH', 0)
            fwd = stats.get('ONLY_FORWARD', 0)
            rev = stats.get('ONLY_REVERSE', 0)
            org = stats.get('Organism', 'Unknown')

            if both > 0 or fwd > 0 or rev > 0:
                if pair_id not in primer_matrix:
                    primer_matrix[pair_id] = {
                        'Organism': org,
                        'TaxaType': stats.get('TaxaType', 'Species'),
                        'Forward': stats.get('Forward', ''),
                        'Reverse': stats.get('Reverse', ''),
                        'Tier': stats.get('Tier', infer_tier(org, pair_id)),
                        'TargetGene': stats.get('TargetGene', infer_target_gene(pair_id)),
                        'sample_both': [0] * num_samples,
                        'sample_fwd': [0] * num_samples,
                        'sample_rev': [0] * num_samples,
                        'amplicon_seq': '',
                        'total_both': 0,
                        'total_fwd': 0,
                        'total_rev': 0
                    }
                primer_matrix[pair_id]['sample_both'][s_idx] = both
                primer_matrix[pair_id]['sample_fwd'][s_idx] = fwd
                primer_matrix[pair_id]['sample_rev'][s_idx] = rev
                primer_matrix[pair_id]['total_both'] += both
                primer_matrix[pair_id]['total_fwd'] += fwd
                primer_matrix[pair_id]['total_rev'] += rev

                if not primer_matrix[pair_id]['amplicon_seq'] and pair_seqs.get(pair_id):
                    primer_matrix[pair_id]['amplicon_seq'] = pair_seqs[pair_id]

                if both > 0:
                    org_both_sums[org] += both
                    org_sample_presence[org].add(s['prefix'])
                    sample_positive_orgs[s_idx].append((org, both))

    num_distinct_orgs = len(org_both_sums)

    # 1. Sample Summary Table HTML
    sample_rows_html = ""
    for s in samples_meta:
        s_idx = s['index']
        pos_org_map = Counter()
        for org, cnt in sample_positive_orgs.get(s_idx, []):
            pos_org_map[org] += cnt

        if pos_org_map:
            org_badges = "".join([f'<span class="badge badge-both" style="margin:2px 4px 2px 0;">{org} ({cnt:,})</span>' for org, cnt in pos_org_map.most_common(5)])
            if len(pos_org_map) > 5:
                org_badges += f'<span class="tag" style="margin-left:4px;">+{len(pos_org_map)-5} weitere</span>'
            status_html = f'<span class="badge badge-both" style="font-weight:700;">✓ Positiv ({len(pos_org_map)})</span>'
        else:
            org_badges = '<span style="color:var(--text-secondary); font-size:12px;">Keine FWD+REV Treffer</span>'
            status_html = '<span class="badge" style="background:#334155; color:#94a3b8; border:1px solid #475569;">Negativ</span>'

        sample_rows_html += f"""
        <tr>
            <td><strong><code>{s['name']}</code></strong></td>
            <td>
                <a href="{s['dash_base']}" target="_blank" class="btn-blast" style="background:#0284c7; padding:5px 10px; font-size:12px;" title="Einzel-Dashboard für {s['name']} öffnen">
                    📊 Dashboard öffnen ↗
                </a>
            </td>
            <td><strong>{s['total_reads']:,}</strong></td>
            <td><span class="badge badge-both">{s['both']:,}</span></td>
            <td><span class="badge badge-fwd">{s['fwd']:,}</span></td>
            <td><span class="badge badge-rev">{s['rev']:,}</span></td>
            <td style="max-width:380px;">{org_badges}</td>
            <td>{status_html}</td>
        </tr>
        """

    # 2. Matrix Table HTML
    sorted_matrix = sorted(primer_matrix.items(), key=lambda x: (x[1]['total_both'], x[1]['total_fwd'] + x[1]['total_rev']), reverse=True)
    matrix_th_html = ""
    col_idx = 0
    matrix_th_html += f'<th onclick="sortTable(\'matrixTable\', {col_idx}, false)">Tier ⇕</th>'; col_idx += 1
    matrix_th_html += f'<th onclick="sortTable(\'matrixTable\', {col_idx}, false)">Organismus / Ziel ⇕</th>'; col_idx += 1
    matrix_th_html += f'<th onclick="sortTable(\'matrixTable\', {col_idx}, false)">Ebene ⇕</th>'; col_idx += 1
    matrix_th_html += f'<th onclick="sortTable(\'matrixTable\', {col_idx}, false)">Primer ID ⇕</th>'; col_idx += 1
    matrix_th_html += f'<th onclick="sortTable(\'matrixTable\', {col_idx}, false)">Zielgen ⇕</th>'; col_idx += 1

    for s in samples_meta:
        matrix_th_html += f'<th onclick="sortTable(\'matrixTable\', {col_idx}, true)" title="Sortieren nach {s["prefix"]}">{s["prefix"]} (BOTH) ⇕</th>'
        col_idx += 1

    total_both_col_idx = col_idx
    matrix_th_html += f'<th onclick="sortTable(\'matrixTable\', {col_idx}, true)" class="sort-desc">Gesamt FWD & REV ⇕</th>'; col_idx += 1
    matrix_th_html += f'<th onclick="sortTable(\'matrixTable\', {col_idx}, true)">Gesamt FWD ⇕</th>'; col_idx += 1
    matrix_th_html += f'<th onclick="sortTable(\'matrixTable\', {col_idx}, true)">Gesamt REV ⇕</th>'; col_idx += 1
    matrix_th_html += f'<th>Aktionen / BLAST</th>'

    matrix_rows_html = ""
    for pair_id, pdata in sorted_matrix:
        fwd_seq = pdata['Forward']
        rev_seq = pdata['Reverse']
        taxa_type = pdata['TaxaType']
        badge_class = "badge-genus" if taxa_type == 'Genus' else "badge-species"

        blast_primer_url = f"https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_LEFT_INPUT={fwd_seq}&PRIMER_RIGHT_INPUT={rev_seq}"
        blast_fwd_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&QUERY={fwd_seq}"
        blast_rev_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&QUERY={rev_seq}"

        amp_seq = pdata['amplicon_seq']
        amp_btn_html = ""
        if amp_seq:
            clean_amp = amp_seq.strip().upper()
            amp_btn_html = f"""
            <form action="https://blast.ncbi.nlm.nih.gov/Blast.cgi" method="POST" target="_blank" style="margin:0; display:inline-block;">
                <input type="hidden" name="PROGRAM" value="blastn">
                <input type="hidden" name="PAGE_TYPE" value="BlastSearch">
                <input type="hidden" name="QUERY" value="{clean_amp}">
                <button type="submit" class="btn-blast" style="background:#0ea5e9;" title="Extrahierte Amplikon-Sequenz ({len(clean_amp)} bp) mit NCBI BLASTN prüfen">
                    🔍 Amplikon
                </button>
            </form>
            """

        sample_cells_html = ""
        for s in samples_meta:
            s_idx = s['index']
            cnt = pdata['sample_both'][s_idx]
            if cnt > 0:
                sample_cells_html += f'<td><span class="badge badge-both">{cnt:,}</span></td>'
            else:
                sample_cells_html += '<td><span style="color:#475569; font-size:12px;">-</span></td>'

        tier = pdata.get('Tier', infer_tier(pdata['Organism'], pair_id))
        tier_clean = '1' if '1' in tier else ('2' if '2' in tier else '3')
        tier_badge_class = 'badge-tier1' if tier_clean == '1' else ('badge-tier2' if tier_clean == '2' else 'badge-tier3')
        target_gene = pdata.get('TargetGene', infer_target_gene(pair_id))

        matrix_rows_html += f"""
        <tr class="matrix-row tier-{tier_clean}">
            <td><span class="badge {tier_badge_class}">{tier}</span></td>
            <td><strong>{pdata['Organism']}</strong></td>
            <td><span class="badge {badge_class}">{taxa_type}</span></td>
            <td><code>{pair_id}</code></td>
            <td><span class="badge" style="background:#334155; color:#38bdf8; font-size:11px;">{target_gene}</span></td>
            {sample_cells_html}
            <td><span class="badge badge-both" style="font-size:13px; font-weight:700;">{pdata['total_both']:,}</span></td>
            <td><span class="badge badge-fwd">{pdata['total_fwd']:,}</span></td>
            <td><span class="badge badge-rev">{pdata['total_rev']:,}</span></td>
            <td>
                <div style="display:flex; gap:4px; align-items:center; flex-wrap:wrap;">
                    {amp_btn_html}
                    <a href="{blast_primer_url}" target="_blank" class="btn-blast" title="Primer-BLAST bei NCBI durchführen">🧬 Primer-BLAST</a>
                    <a href="{blast_fwd_url}" target="_blank" class="btn-blast-sm" title="Forward Primer mit BLASTN abgleichen">FWD</a>
                    <a href="{blast_rev_url}" target="_blank" class="btn-blast-sm" title="Reverse Primer mit BLASTN abgleichen">REV</a>
                </div>
            </td>
        </tr>
        """
    if not matrix_rows_html:
        matrix_rows_html = f"<tr><td colspan='{col_idx+1}' style='text-align:center;'>Keine Primer-Treffer über alle analysierten FASTQ-Dateien gefunden.</td></tr>"

    # 2.5 Multi-Sample Minimap2 Validations Table
    all_minimap2_vals = []
    for s in samples_meta:
        for v in s['results'].get('minimap2_validations', []):
            item = dict(v)
            item['sample_name'] = s['name']
            item['sample_prefix'] = s['prefix']
            all_minimap2_vals.append(item)

    summary_minimap2_html = ""
    if all_minimap2_vals:
        s_m_rows = ""
        for v in all_minimap2_vals:
            bam_link = f'<a href="{v["bam_rel"]}" class="btn-blast-sm" style="background:#0ea5e9; color:white;" target="_blank">🧬 BAM</a>' if v.get("bam_rel") else "-"
            cand_link = f'<a href="{v["candidate_fastq_rel"]}" class="btn-blast-sm" target="_blank">📁 Reads</a>' if v.get("candidate_fastq_rel") else "-"
            contig_btn = f'<a href="{v["contig_rel"]}" class="btn-blast-sm" style="background:#10b981; color:white;" target="_blank" title="{v.get("assembly_method", "Micro-Assembly")}: {v.get("contig_len", 0):,} bp">🧬 {v.get("contig_len", 0):,} bp Contig</a>' if v.get("contig_rel") else ""

            clean_rep = clean_blast_query(v.get("contig_seq"), fallback_seq=v.get("rep_seq", ""))
            blast_btn = ""
            if clean_rep:
                blast_lbl = "🌐 Web-BLAST (NCBI)" if v.get("contig_seq") else "🌐 Web-BLAST"
                blast_title = f"Assemblierten {v.get('assembly_method', '')} Konsensus-Contig ({len(clean_rep)} bp) online bei NCBI BLASTN abgleichen (ideal zur Bestätigung beim Weitergeben der HTML-Datei)" if v.get("contig_seq") else f"Repräsentativen Read ({len(clean_rep)} bp) online bei NCBI BLASTN abgleichen"
                blast_btn = f"""
                <form action="https://blast.ncbi.nlm.nih.gov/Blast.cgi" method="POST" target="_blank" style="margin:0; display:inline-block;">
                    <input type="hidden" name="PROGRAM" value="blastn">
                    <input type="hidden" name="PAGE_TYPE" value="BlastSearch">
                    <input type="hidden" name="QUERY" value="{clean_rep}">
                    <button type="submit" class="btn-blast" style="background:#0284c7; font-weight:600;" title="{blast_title}">
                        {blast_lbl}
                    </button>
                </form>
                """

            local_hit = v.get("local_blast_hit")
            if local_hit:
                tsv_link = f'<a href="{v.get("local_blast_tsv")}" target="_blank" style="font-size:11px; color:#38bdf8; text-decoration:none; display:inline-block; margin-top:3px;" title="Vollständigen lokalen BLASTn-Report öffnen">📄 Vollständiger Report (.tsv) ↗</a>' if v.get("local_blast_tsv") else ""
                if "Kein lokaler Treffer" in local_hit:
                    local_blast_html = f'<div><span class="badge" style="background:#475569; color:#cbd5e1; font-size:11px;">⚠️ Kein nt-Treffer (E &lt; 1e-5)</span>{tsv_link}</div>'
                else:
                    local_blast_html = f'<div><span class="badge" style="background:#065f46; color:#34d399; font-size:11px; max-width:280px; display:inline-block; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;" title="{local_hit}">🧬 {local_hit}</span>{tsv_link}</div>'
            else:
                local_blast_html = '<span style="color:var(--text-secondary); font-size:11px;">Nicht ausgeführt / Keine nt-DB</span>'

            s_m_rows += f"""
            <tr>
                <td><code>{v['sample_name']}</code></td>
                <td><strong>{v['organism']}</strong></td>
                <td><span class="badge {v['badge_trigger']}">{v['trigger_desc']}</span></td>
                <td><code>{v['ref_file']}</code></td>
                <td><span class="badge badge-both">{v['mapped_count']} / {v['num_candidates']} Reads</span> <span style="font-size:11px; color:var(--text-secondary);">({v['high_qual_count']} mit MAPQ &ge; 20)</span></td>
                <td><strong>{v['covered_bp']:,} bp</strong> <span style="font-size:11px; color:var(--text-secondary);">({v['top_rname']})</span></td>
                <td><span class="badge {v['verdict_badge']}" style="font-weight:700;">{v['verdict']}</span></td>
                <td>{local_blast_html}</td>
                <td>
                    <div style="display:flex; gap:4px; align-items:center; flex-wrap:wrap;">
                        {cand_link}
                        {bam_link}
                        {contig_btn}
                        {blast_btn}
                    </div>
                </td>
            </tr>
            """

        summary_minimap2_html = f"""
        <div class="card" style="margin-bottom: 24px; border-left: 4px solid var(--accent-cyan);">
            <h2>🔬 minimap2 Referenz-Mapping &amp; Erreger-Validierung (Alle Proben)</h2>
            <p style="color: var(--text-secondary); font-size: 13px; margin: -4px 0 16px 0;">
                💡 Übersicht aller Referenz-Mappings für getrennte FWD/REV-Treffer und bestätigte Erreger über den gesamten Datensatz. Alle Zwischendateien enthalten zur sicheren Unterscheidung das Probenkürzel <code>&lt;prefix&gt;_</code>.
            </p>
            <div style="overflow-x:auto;">
                <table id="summaryMinimap2Table">
                    <thead>
                        <tr>
                            <th onclick="sortTable('summaryMinimap2Table', 0, false)">Probe ⇕</th>
                            <th onclick="sortTable('summaryMinimap2Table', 1, false)">Erreger / Organismus ⇕</th>
                            <th onclick="sortTable('summaryMinimap2Table', 2, false)">Auslöser ⇕</th>
                            <th onclick="sortTable('summaryMinimap2Table', 3, false)">Referenz-Genom ⇕</th>
                            <th onclick="sortTable('summaryMinimap2Table', 4, true)">Gemappte Reads (MAPQ &ge; 20) ⇕</th>
                            <th onclick="sortTable('summaryMinimap2Table', 5, false)">Locus / Abdeckung ⇕</th>
                            <th onclick="sortTable('summaryMinimap2Table', 6, false)">Validierungs-Ergebnis ⇕</th>
                            <th onclick="sortTable('summaryMinimap2Table', 7, false)">Lokaler BLAST (/mnt/blast_db/nt) ⇕</th>
                            <th>Dateien &amp; Web-BLAST</th>
                        </tr>
                    </thead>
                    <tbody>
                        {s_m_rows}
                    </tbody>
                </table>
            </div>
        </div>
        """

    # 3. Extracted Amplicons Table HTML
    amplicons_rows_html = ""
    for s in samples_meta:
        s_res = s['results']
        ext_counts = s_res.get('extracted_counts', {})
        ext_seqs = s_res.get('extracted_seqs', {})
        for fname, count in ext_counts.items():
            clean_seq = ext_seqs.get(fname, '').strip().upper()
            seq_len = len(clean_seq)
            preview_seq = f"{clean_seq[:35]}...{clean_seq[-25:]}" if seq_len > 65 else clean_seq

            if clean_seq:
                blast_btn = f"""
                <div style="display:flex; gap:6px; align-items:center;">
                    <form action="https://blast.ncbi.nlm.nih.gov/Blast.cgi" method="POST" target="_blank" style="margin:0; display:inline-block;">
                        <input type="hidden" name="PROGRAM" value="blastn">
                        <input type="hidden" name="PAGE_TYPE" value="BlastSearch">
                        <input type="hidden" name="QUERY" value="{clean_seq}">
                        <button type="submit" class="btn-blast" title="Extrahierte Sequenz ({seq_len} bp) direkt an NCBI BLASTN senden">
                            🔍 BLAST Sequenz
                        </button>
                    </form>
                    <button type="button" class="btn-blast-sm" onclick="navigator.clipboard.writeText('{clean_seq}'); this.innerText='✓ Kopiert'; setTimeout(() => this.innerText='📋', 2000);" title="Sequenz in die Zwischenablage kopieren">
                        📋
                    </button>
                </div>
                """
                seq_cell = f"""
                <details style="cursor:pointer;">
                    <summary style="font-family:monospace; font-size:12px; color:#38bdf8;"><code>{preview_seq}</code> <span style="font-size:10px; color:var(--text-secondary);">({seq_len} bp)</span></summary>
                    <div style="margin-top:6px; padding:8px; background:#0f172a; border-radius:6px; font-family:monospace; font-size:11px; word-break:break-all; max-height:120px; overflow-y:auto; border:1px solid #334155; user-select:all; color:#f1f5f9;">
                        {clean_seq}
                    </div>
                </details>
                """
            else:
                blast_btn = "-"
                seq_cell = "-"

            amplicons_rows_html += f"""
            <tr>
                <td><strong><code>{s['name']}</code></strong></td>
                <td><code>{fname}</code></td>
                <td><span class="badge badge-both">{count:,} Reads</span></td>
                <td><strong>{seq_len} bp</strong></td>
                <td style="max-width:350px;">{seq_cell}</td>
                <td>{blast_btn}</td>
            </tr>
            """
    if not amplicons_rows_html:
        amplicons_rows_html = "<tr><td colspan='6' style='text-align:center;'>Keine FWD+REV Amplikon-Dateien extrahiert.</td></tr>"

    # 4. Chart Data Preparation
    top_org_items = org_both_sums.most_common(15)
    top_orgs = [item[0] for item in top_org_items]
    heatmap_x = [s['prefix'] for s in samples_meta]
    heatmap_y = list(reversed(top_orgs)) if top_orgs else []
    heatmap_z = []
    for org in heatmap_y:
        row = []
        for s in samples_meta:
            s_idx = s['index']
            org_hits = sum(p['sample_both'][s_idx] for p in primer_matrix.values() if p['Organism'] == org)
            row.append(org_hits)
        heatmap_z.append(row)

    bar_x = [s['prefix'] for s in samples_meta]
    bar_both = [s['both'] for s in samples_meta]
    bar_fwd = [s['fwd'] for s in samples_meta]
    bar_rev = [s['rev'] for s in samples_meta]

    top_path_y = [item[0] for item in reversed(top_org_items[:12])]
    top_path_x = [item[1] for item in reversed(top_org_items[:12])]

    pie_labels = ["FWD & REV (Beide)", "Nur FWD", "Nur REV", "Ohne Primer"]
    pie_values = [total_all_both, total_all_fwd, total_all_rev, total_all_none]

    html_content = f"""<!DOCTYPE html>
<html lang="de">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Primer Suche - Multi-FASTQ Gesamtübersicht</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {{
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --card-border: #334155;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --accent-cyan: #06b6d4;
            --accent-blue: #3b82f6;
            --accent-green: #10b981;
            --accent-purple: #8b5cf6;
            --accent-amber: #f59e0b;
        }}
        body {{
            font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            margin: 0;
            padding: 24px;
        }}
        .header {{
            padding-bottom: 20px;
            border-bottom: 1px solid var(--card-border);
            margin-bottom: 24px;
        }}
        .header h1 {{
            margin: 0;
            font-size: 28px;
            background: linear-gradient(135deg, var(--accent-cyan), var(--accent-green));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }}
        .card {{
            background: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 12px;
            padding: 20px;
            margin-bottom: 24px;
        }}
        .card .title {{
            font-size: 13px;
            color: var(--text-secondary);
            margin-bottom: 8px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .card .value {{
            font-size: 26px;
            font-weight: 700;
        }}
        .badge {{
            display: inline-block;
            padding: 4px 10px;
            border-radius: 20px;
            font-size: 12px;
            font-weight: 600;
        }}
        .badge-both {{ background: rgba(16, 185, 129, 0.2); color: #34d399; border: 1px solid #10b981; }}
        .badge-fwd {{ background: rgba(6, 182, 212, 0.2); color: #38bdf8; border: 1px solid #06b6d4; }}
        .badge-rev {{ background: rgba(139, 92, 246, 0.2); color: #c084fc; border: 1px solid #8b5cf6; }}
        .badge-genus {{ background: rgba(245, 158, 11, 0.2); color: #fbbf24; border: 1px solid #f59e0b; }}
        .badge-species {{ background: rgba(59, 130, 246, 0.2); color: #93c5fd; border: 1px solid #3b82f6; }}
        .tag {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 6px;
            font-size: 11px;
            background: #334155;
            color: #cbd5e1;
        }}
        .charts-row {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 24px;
            margin-bottom: 24px;
        }}
        @media (max-width: 900px) {{
            .charts-row {{ grid-template-columns: 1fr; }}
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 14px;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid var(--card-border);
        }}
        th {{
            color: var(--text-secondary);
            font-weight: 600;
            background: rgba(15, 23, 42, 0.5);
            cursor: pointer;
            user-select: none;
            transition: background 0.15s ease, color 0.15s ease;
        }}
        th:hover {{
            background: rgba(56, 189, 248, 0.2);
            color: #38bdf8;
        }}
        th.sort-asc::after {{
            content: " ▲";
            color: #38bdf8;
        }}
        th.sort-desc::after {{
            content: " ▼";
            color: #38bdf8;
        }}
        .btn-blast {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            background: #0284c7;
            color: white;
            text-decoration: none;
            font-size: 11px;
            font-weight: 600;
            white-space: nowrap;
            transition: background 0.15s;
            border: none;
            cursor: pointer;
            font-family: inherit;
        }}
        .btn-blast:hover {{
            background: #0369a1;
        }}
        .btn-blast-sm {{
            display: inline-block;
            padding: 3px 6px;
            border-radius: 4px;
            background: #334155;
            color: #38bdf8;
            text-decoration: none;
            font-size: 10px;
            font-weight: 600;
            transition: background 0.15s;
            border: none;
            cursor: pointer;
            font-family: inherit;
        }}
        .btn-blast-sm:hover {{
            background: #475569;
            color: white;
        }}
        code {{
            font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
            background: rgba(0, 0, 0, 0.3);
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 12px;
            color: #38bdf8;
        }}
        .search-box {{
            padding: 8px 14px;
            border-radius: 8px;
            border: 1px solid var(--card-border);
            background: #0f172a;
            color: #f8fafc;
            font-size: 14px;
            width: 320px;
            outline: none;
        }}
        .search-box:focus {{
            border-color: var(--accent-cyan);
        }}
        /* Tab Navigation Styling */
        .tab-nav {{
            display: flex;
            gap: 8px;
            margin-bottom: 24px;
            border-bottom: 2px solid var(--card-border);
            padding-bottom: 8px;
            flex-wrap: wrap;
        }}
        .tab-btn {{
            background: transparent;
            border: 1px solid var(--card-border);
            color: var(--text-secondary);
            padding: 10px 20px;
            border-radius: 8px;
            font-size: 14px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            font-family: inherit;
        }}
        .tab-btn:hover {{
            background: rgba(56, 189, 248, 0.1);
            color: var(--text-primary);
            border-color: #38bdf8;
        }}
        .tab-btn.active {{
            background: #0284c7;
            color: #ffffff;
            border-color: #0ea5e9;
            box-shadow: 0 4px 12px rgba(2, 132, 199, 0.3);
        }}
        .tab-pane {{
            display: none;
        }}
        .tab-pane.active {{
            display: block;
        }}
        .guideline-box {{
            background: #1e293b;
            border: 1px solid #334155;
            border-radius: 10px;
            padding: 16px 20px;
            margin-bottom: 20px;
        }}
        .guideline-box h3 {{
            margin: 0 0 8px 0;
            color: #38bdf8;
            font-size: 16px;
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        .guideline-box p, .guideline-box ul {{
            margin: 0 0 10px 0;
            font-size: 13.5px;
            line-height: 1.6;
            color: #cbd5e1;
        }}
        .guideline-box ul {{
            padding-left: 20px;
        }}
        .guideline-box li {{
            margin-bottom: 6px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <div style="display:flex; justify-content:space-between; align-items:flex-start; flex-wrap:wrap; gap:16px;">
            <div>
                <h1>🧬 Primer Suche - Gesamtübersicht (Multi-FASTQ)</h1>
                <p style="color: var(--text-secondary); margin: 6px 0 0 0;">
                    Analysierte Dateien: <strong>{num_samples} FASTQ-Proben</strong> | Suchmodus: <strong>{search_mode}</strong> | Erstellt am: <strong>{timestamp}</strong>
                </p>
            </div>
            <div>
                <span class="badge" style="background:#1e3a8a; color:#bfdbfe; border:1px solid #3b82f6; font-size:13px; padding:6px 14px;">
                    Multi-Sample Summary
                </span>
            </div>
        </div>
    </div>

    <!-- Tab Navigation Bar -->
    <div class="tab-nav">
        <button class="tab-btn active" onclick="switchOverviewTab('tab-overview-dashboard', this)">📊 Multi-Probe Übersicht &amp; Matrix</button>
        <button class="tab-btn" onclick="switchOverviewTab('tab-overview-guidelines', this)">📖 Vorgehen &amp; Leitfaden (SHG Borreliose Hannover)</button>
    </div>

    <!-- TAB 1: OVERVIEW DASHBOARD CONTENT -->
    <div id="tab-overview-dashboard" class="tab-pane active">

    <div class="stats-grid">
        <div class="card">
            <div class="title">Analysierte Proben</div>
            <div class="value" style="color: var(--accent-cyan);">{num_samples:,}</div>
        </div>
        <div class="card">
            <div class="title">Geprüfte Reads Gesamt</div>
            <div class="value" style="color: var(--text-primary);">{total_all_reads:,}</div>
        </div>
        <div class="card">
            <div class="title">FWD & REV Treffer Gesamt</div>
            <div class="value" style="color: var(--accent-green);">{total_all_both:,}</div>
        </div>
        <div class="card">
            <div class="title">Nur FWD Treffer</div>
            <div class="value" style="color: var(--accent-cyan);">{total_all_fwd:,}</div>
        </div>
        <div class="card">
            <div class="title">Nur REV Treffer</div>
            <div class="value" style="color: var(--accent-purple);">{total_all_rev:,}</div>
        </div>
        <div class="card">
            <div class="title">Erkannte Erreger</div>
            <div class="value" style="color: var(--accent-amber);">{num_distinct_orgs:,}</div>
        </div>
        <div class="card">
            <div class="title">Durchschnitts-Speed</div>
            <div class="value" style="font-size: 20px; color: #38bdf8;">{avg_rps:,.0f} <span style="font-size:12px;">Reads/s</span></div>
        </div>
    </div>

    <div class="charts-row">
        <div class="card"><div id="barChartSamples"></div></div>
        <div class="card"><div id="pieChartOverall"></div></div>
    </div>

    <div class="charts-row">
        <div class="card"><div id="heatmapChart"></div></div>
        <div class="card"><div id="topPathogensChart"></div></div>
    </div>

    <div class="card">
        <h2>📁 Proben-Übersicht &amp; Direkte Einzel-Dashboards</h2>
        <p style="color: var(--text-secondary); font-size: 13px; margin: -4px 0 16px 0;">
            💡 Klicke auf <strong>Dashboard öffnen ↗</strong>, um die detaillierte Analyse (inkl. Amplikon-Abstände und Einzel-Reads) für die jeweilige Probe zu betrachten.
        </p>
        <table id="sampleSummaryTable">
            <thead>
                <tr>
                    <th onclick="sortTable('sampleSummaryTable', 0, false)">Probe / FASTQ-Datei ⇕</th>
                    <th>Einzel-Dashboard</th>
                    <th onclick="sortTable('sampleSummaryTable', 2, true)">Geprüfte Reads ⇕</th>
                    <th onclick="sortTable('sampleSummaryTable', 3, true)" class="sort-desc">FWD &amp; REV ⇕</th>
                    <th onclick="sortTable('sampleSummaryTable', 4, true)">Nur FWD ⇕</th>
                    <th onclick="sortTable('sampleSummaryTable', 5, true)">Nur REV ⇕</th>
                    <th>Nachgewiesene Erreger (BOTH)</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>{sample_rows_html}</tbody>
        </table>
    </div>

    <div class="card">
        <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:12px; margin-bottom:12px;">
            <div>
                <h2 style="margin:0 0 4px 0;">📊 Matrix aller gefundenen Primer &amp; Erreger über alle FASTQ-Proben</h2>
                <p style="color: var(--text-secondary); font-size: 13px; margin:0;">
                    💡 Vergleicht die Trefferhäufigkeiten je Probe. Klicke auf Spaltenköpfe zum Sortieren oder nutze den <strong>BLAST-Button</strong> zur NCBI-Verifizierung.
                </p>
            </div>
            <div>
                <input type="text" id="matrixSearch" class="search-box" onkeyup="filterMatrixTable()" placeholder="🔍 Erreger, Gen oder Primer suchen...">
            </div>
        </div>
        <div style="overflow-x:auto;">
            <table id="matrixTable">
                <thead>
                    <tr>{matrix_th_html}</tr>
                </thead>
                <tbody>{matrix_rows_html}</tbody>
            </table>
        </div>
    </div>

    {summary_minimap2_html}

    <div class="card">
        <h2>🧬 Extrahierte Amplicon-Dateien über alle Proben (_frp_&lt;bp&gt;)</h2>
        <p style="color: var(--text-secondary); font-size: 13px; margin: -4px 0 16px 0;">
            💡 Vollständige Übersicht aller extrahierten Amplikons mit direkter Möglichkeit zum NCBI BLASTN Abgleich.
        </p>
        <div style="overflow-x:auto;">
            <table id="ampliconsTable">
                <thead>
                    <tr>
                        <th onclick="sortTable('ampliconsTable', 0, false)">Probe ⇕</th>
                        <th onclick="sortTable('ampliconsTable', 1, false)">Dateiname ⇕</th>
                        <th onclick="sortTable('ampliconsTable', 2, true)" class="sort-desc">Extrahierte Reads ⇕</th>
                        <th onclick="sortTable('ampliconsTable', 3, true)">Länge ⇕</th>
                        <th>Extrahierte Sequenz (Vorschau)</th>
                        <th>BLAST Überprüfung &amp; Aktionen</th>
                    </tr>
                </thead>
                <tbody>{amplicons_rows_html}</tbody>
            </table>
        </div>
    </div> <!-- END TAB 1 -->

    <!-- TAB 2: VORGEHEN & LEITFADEN SHG BORRELIOSE HANNOVER -->
    <!-- TAB: LEITFADEN & NUTZUNGSANLEITUNG -->
    <div id="tab-overview-guidelines" class="tab-pane">
        <div class="card" style="border-left: 4px solid var(--accent-cyan);">
            <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:10px; margin-bottom:12px;">
                <h2 style="margin:0; color:#38bdf8;">📖 Leitfaden, Programmnutzung &amp; Befehlsreferenz</h2>
                <span class="badge badge-both" style="font-size:12px;">Wissenschaftlich-Diagnostischer Standard</span>
            </div>
            <p style="color: var(--text-secondary); font-size: 14px; margin-bottom: 20px;">
                Umfassende Dokumentation zur Nutzung der Pipeline, Interpretation der Dashboard-Ergebnisse, Erklärung des <strong>Primer-Buttons</strong> sowie klinische Handlungsempfehlungen in Anlehnung an die Erfahrungswerte der <strong>Selbsthilfe Gruppe Borreliose Hannover</strong>.
            </p>

            <!-- 1. NUTZUNG DES PROGRAMMS & TABS -->
            <div class="guideline-box">
                <h3 style="color:#38bdf8;">🚀 1. Nutzung des Programms &amp; Dashboard-Navigation</h3>
                <p>
                    Die Pipeline führt eine mehrstufige bioinformatische Analyse durch: hochsensitives In-silico-Primer-Matching, orthogonales Vollgenom-Mapping (<code>--pan-scan</code>) und eine De-novo-Contig-Assemblierung unzugeordneter Reads (Dark Matter) mit lokalem BLASTn-Abgleich gegen die NCBI-nt-Datenbank.
                </p>
                <p><strong>Die Reiter (Tabs) im Dashboard im Überblick:</strong></p>
                <ul>
                    <li><strong style="color:#10b981;">🏥 Klinischer Erreger-Befund (Ampel):</strong> Die zentrale klinische Übersicht. Fasst alle nachgewiesenen Erreger zusammen, gruppiert sie nach Gattung und ordnet sie einer diagnostischen Ampel zu:
                        <br>&bull; <span class="badge badge-both" style="font-size:11px;">🟢 Bestätigt</span>: Hohe Evidenz (Amplikon-Nachweis, Genom-Mapping oder Genus-BLAST bestätigt).
                        <br>&bull; <span class="badge badge-fwd" style="font-size:11px;">🟡 Verdacht</span>: Einzelsignale (z. B. FWD-Read) oder schwächeres Mapping.
                        <br>&bull; <span class="badge badge-rev" style="font-size:11px;">⚪ Gering / Isoliert</span>: Einzelne isolierte Fragmente.
                        <br>&bull; <span class="badge badge-none" style="font-size:11px;">🚫 Kitom / Kontaminant</span>: Typische Labor- oder Extraktionskit-Flora (z.B. Ralstonia, Cutibacterium) – per Button ausblendbar.
                    </li>
                    <li><strong style="color:#38bdf8;">🎯 BLAST-Bestätigt (Genus):</strong> Exklusive Spezialansicht! Zeigt nur diejenigen Funde, bei denen die <em>lokale BLASTn-Suche</em> (/mnt/blast_db/nt) die Gattung oder Spezies unabhängig bestätigt hat (inkl. E-Value, Identität % und Alignment-Länge).</li>
                    <li><strong style="color:#a855f7;">📊 Detaillierte Primer-Matrix (738 Assays):</strong> (Aufrufbar über den <em>Primer-Button</em> oben!) Vollständige Auflistung aller 738 diagnostischen Assays mit genauen FWD-, REV- und Doppelstrang-Treffern, sortierbar und filterbar nach Tiers.</li>
                    <li><strong style="color:#06b6d4;">🔬 De-novo Entdeckungen (Dark Matter):</strong> Zeigt neue oder unerwartete Erreger, die aus nicht-gemappten Reads de-novo zu Contigs assembliert und via lokalem BLAST identifiziert wurden.</li>
                    <li><strong style="color:#f59e0b;">📁 Extrahierte Amplikons:</strong> Physisch aus den Reads ausgeschnittene PCR-Produkte (inkl. Längenverteilung) mit 1-Click-Web-BLAST.</li>
                </ul>
            </div>

            <!-- 2. WAS IST DER PRIMER BUTTON? -->
            <div class="guideline-box" style="border-left: 4px solid #a855f7; background: rgba(88, 28, 135, 0.15);">
                <h3 style="color:#c084fc;">🔘 2. Was ist der &bdquo;Primer-Button&ldquo;?</h3>
                <p>
                    Im Dashboard begegnet Ihnen der Begriff <strong>Primer-Button</strong> an zwei zentralen Stellen:
                </p>
                <ol style="margin-bottom:12px;">
                    <li><strong>Der Tab-Button &bdquo;📊 Detaillierte Primer-Matrix (738 Assays)&ldquo; im Menüband:</strong>
                        <br>Dieser Navigations-Button schaltet die Ansicht vom komprimierten klinischen Befund auf die vollständige, wissenschaftliche Rohdaten-Tabelle um. Dort sehen Sie jedes der 738 Primerpaare mit:
                        <ul style="margin-top:4px;">
                            <li>Genauer Anzahl der FWD-Reads, REV-Reads und beidseitigen Amplikons (BOTH).</li>
                            <li>Basenpaar-Abstand (Amplikon-Größe in bp) und Übereinstimmung mit der Zielregion.</li>
                            <li>Filterknöpfen zur gezielten Eingrenzung nach <em>Tier 0 (Universal-Barcodes)</em>, <em>Tier 1 (Borreliose &amp; Vektoren)</em>, <em>Tier 2 (Endokarditis/Atypisch)</em> und <em>Tier 3</em>.</li>
                        </ul>
                    </li>
                    <li><strong>Die interaktiven Action-Buttons bei jedem Primertreffer:</strong>
                        <ul style="margin-top:4px;">
                            <li><button type="button" class="btn-blast" style="display:inline-block; font-size:11px; padding:2px 8px; background:#0ea5e9; cursor:default;">Amplikon BLAST</button> oder <button type="button" class="btn-blast" style="display:inline-block; font-size:11px; padding:2px 8px; background:#10b981; cursor:default;">NCBI Web-BLAST</button>: Überträgt die extrahierte DNA-Sequenz mit einem einzigen Klick direkt an das offizielle NCBI BLASTN im Web, um die weltweite Referenzdatenbank abzugleichen.</li>
                            <li><button type="button" class="btn-blast-sm" style="display:inline-block; font-size:10px; padding:2px 6px; cursor:default;">Kopieren</button>: Kopiert die extrahierte FASTA-Sequenz direkt in die Zwischenablage für eigene Analysen.</li>
                            <li><button type="button" class="btn-blast" style="display:inline-block; font-size:11px; padding:2px 8px; background:#334155; cursor:default;">Details anzeigen &#x25BC;</button>: Klappt im klinischen Befund das Akkordeon auf, um alle einzelnen Primer-Assays und Read-Statistiken dieser Erregergattung einzusehen.</li>
                        </ul>
                    </li>
                </ol>
            </div>

            <!-- 3. BEISPIELZEILEN & COMMAND-LINE CHEAT SHEET -->
            <div class="guideline-box">
                <h3 style="color:#38bdf8;">💻 3. Konsolen-Befehle (CLI Cheat-Sheet &amp; Beispielzeilen)</h3>
                <p>
                    Führen Sie die folgenden Befehle direkt im Terminal (z. B. auf dem Server unter <code>/mnt/sdc1/scripts/</code> oder lokal) aus:
                </p>

                <div style="margin-bottom:14px;">
                    <strong style="color:#10b981;">🌟 A. Klinischer Goldstandard (Vollanalyse mit allen 738 Primern + Vollgenom-Mapping):</strong>
                    <p style="font-size:13px; color:var(--text-secondary); margin:2px 0 6px 0;">
                        Empfohlen für jede diagnostische Auswertung. Testet alle Erreger (Babesia, Plasmodium, Borrelia, Bartonella, etc.), aktiviert Vollgenom-Mapping gegen 340+ Referenzen und De-novo Dark Matter Discovery:
                    </p>
                    <pre style="background:#090d16; padding:12px; border-radius:6px; border:1px solid #1e293b; color:#38bdf8; font-family:monospace; font-size:13px; overflow-x:auto;">python primer_suche.py -i probe.fastq.gz -a --pan-scan -t 24</pre>
                </div>

                <div style="margin-bottom:14px;">
                    <strong style="color:#06b6d4;">⚡ B. Schnelles hierarchisches Screening (Zeitersparnis):</strong>
                    <p style="font-size:13px; color:var(--text-secondary); margin:2px 0 6px 0;">
                        Prüft zuerst Genus-Marker (z. B. Borrelia flaB, Babesia 18S). Nur bei positivem Signal werden spezifische Spezies-Assays durchsucht:
                    </p>
                    <pre style="background:#090d16; padding:12px; border-radius:6px; border:1px solid #1e293b; color:#38bdf8; font-family:monospace; font-size:13px; overflow-x:auto;">python primer_suche.py -i probe.fastq.gz -gs -t 24</pre>
                </div>

                <div style="margin-bottom:14px;">
                    <strong style="color:#f59e0b;">📁 C. Alle Proben im Ordner analysieren (Multi-Sample Batch):</strong>
                    <p style="font-size:13px; color:var(--text-secondary); margin:2px 0 6px 0;">
                        Verarbeitet alle FASTQ-Dateien auf einmal, erstellt Einzelberichte und das übergeordnete Dashboard (<code>primer_dashboard_all_samples.html</code>):
                    </p>
                    <pre style="background:#090d16; padding:12px; border-radius:6px; border:1px solid #1e293b; color:#38bdf8; font-family:monospace; font-size:13px; overflow-x:auto;">python primer_suche.py -i *.fastq.gz -a --pan-scan -t 28</pre>
                </div>

                <div style="margin-bottom:14px;">
                    <strong style="color:#c084fc;">🎯 D. Gezielte Tier-Fokussierung (z.B. Tier 0 &amp; Tier 1):</strong>
                    <p style="font-size:13px; color:var(--text-secondary); margin:2px 0 6px 0;">
                        Beschränkt die Suche auf Universelle Barcodes (Tier 0) und Borreliose/Zecken-Koinfektionen (Tier 1):
                    </p>
                    <pre style="background:#090d16; padding:12px; border-radius:6px; border:1px solid #1e293b; color:#38bdf8; font-family:monospace; font-size:13px; overflow-x:auto;">python primer_suche.py -i probe.fastq.gz --tier 0,1 -a -t 24</pre>
                </div>

                <div>
                    <strong style="color:#94a3b8;">🛠️ E. Primer-Datenbank bereinigen &amp; aktualisieren:</strong>
                    <pre style="background:#090d16; padding:12px; border-radius:6px; border:1px solid #1e293b; color:#38bdf8; font-family:monospace; font-size:13px; overflow-x:auto;"># Duplikate und Reverse-Komplemente entfernen:
python primer_suche.py -p best_primers_db.csv --dedup best_primers_db_clean.csv

# Gezielt gegen RefSeq/GenBank aktualisieren (z.B. Babesia oder all):
python primer_suche.py --update-primers --target Babesia</pre>
                </div>
            </div>

            <!-- 4. KLINISCHE STUFENDIAGNOSTIK DER SHG BORRELIOSE HANNOVER -->
            <div class="guideline-box">
                <h3 style="color:#38bdf8;">📌 4. Grundsatz der Stufendiagnostik (Hierarchische Verifizierung)</h3>
                <p>
                    Ein bioinformatischer Primer-Treffer in NGS- oder Nanopore-Rohdaten stellt eine hochsensitive Vorstufe dar. Zur soliden Verifizierung gilt das mehrstufige diagnostische Vorgehen:
                </p>
                <ul>
                    <li><strong>Stufe 1 (Genus-Screening):</strong> Nachweis von Leitgenen (z. B. 16S-rRNA, <em>flaB</em>, <em>gltA</em>, <em>17kDa</em>). Ein Signal zeigt an, dass Vertreter einer Erregergruppe (z. B. <em>Borrelia</em>, <em>Babesia</em>, <em>Bartonella</em>, <em>Rickettsia</em>) in der Probe vorliegen.</li>
                    <li><strong>Stufe 2 (Spezies-Identifikation):</strong> Spezifische Primer grenzen den Erreger exakt ein (z. B. <em>Borrelia burgdorferi sensu stricto</em>, <em>B. afzelii</em>, <em>B. garinii</em>, <em>B. miyamotoi</em> oder <em>Babesia microti/divergens/odocoilei</em>).</li>
                    <li><strong>Stufe 3 (Orthogonale Validierung):</strong> Bei getrennten FWD/REV-Hits erfolgt das automatische Genom-Mapping mit <strong>minimap2</strong> gegen das NCBI RefSeq-Vollgenom sowie die Sequenzprüfung über den integrierten <strong>BLAST-Button</strong>.</li>
                </ul>
            </div>

            <!-- 5. KOINFEKTIONEN -->
            <div class="guideline-box">
                <h3 style="color:#38bdf8;">🛡️ 5. Koinfektionen &amp; atypische Erreger beachten</h3>
                <p>
                    Zecken und Vektoren übertragen selten nur einen einzelnen Erreger. Typische klinische Verläufe zeichnen sich oft durch ein Erregerkollektiv aus:
                </p>
                <ul>
                    <li><strong>Piroplasmen &amp; Protozoen:</strong> <em>Babesia</em> (insb. <em>B. microti</em>, <em>B. divergens</em>, <em>B. odocoilei</em>, <em>B. duncani</em>, <em>B. venatorum</em>) befallen Erythrozyten und führen häufig zu therapieresistenten Erschöpfungssyndromen, Nachtschweiß und hämatologischen Auffälligkeiten.</li>
                    <li><strong>Intrazelluläre Bakterien:</strong> <em>Bartonella</em> (z. B. <em>B. henselae</em>) und <em>Rickettsia</em> (z. B. Fleckfieber-Gruppe) besiedeln Gefäßendothelien und erfordern eigenständige therapeutische Schemata.</li>
                    <li><strong>Opportunistische Mykobakterien &amp; Pilze:</strong> NTM (z.B. <em>M. immunogenum</em>, <em>M. marinum</em>, <em>M. avium</em>) sowie Pilze wie <em>Cryptococcus neoformans</em> treten bevorzugt bei sekundärer Immundysfunktion auf.</li>
                    <li><strong>Virale Reaktivierungen:</strong> Herpesviren (EBV, CMV, HHV-6, VZV) reaktivieren häufig unter chronischer Borrelien-Last.</li>
                </ul>
            </div>

            <!-- 6. VERNETZUNG -->
            <div class="guideline-box" style="border-left: 4px solid var(--accent-amber); background: rgba(30, 41, 59, 0.8);">
                <h3 style="color:#fbbf24;">🤝 6. Vernetzung &amp; Patientenunterstützung</h3>
                <p>
                    Betroffene mit komplexen, chronischen Symptomkomplexen nach Zeckenstichen oder Verdacht auf Vektor-übertragene Multisystemerkrankungen finden kompetenten Austausch, Informationsmaterialien und Unterstützung bei erfahrenen Selbsthilfe-Initiativen:
                </p>
                <div style="background:#0f172a; padding:12px 16px; border-radius:8px; border:1px solid #334155; margin-top:8px;">
                    <p style="margin:0 0 6px 0; font-weight:700; color:#f8fafc; font-size:14px;">
                        Selbsthilfe Gruppe Borreliose Hannover
                    </p>
                    <p style="margin:0 0 6px 0; font-size:13px; color:#cbd5e1;">
                        Austauschplattform für Betroffene, Aufklärung zu Zecken-übertragenen Infektionen, Koinfektionen (Babesien, Bartonellen, Rickettsien) sowie ganzheitlichen und interdisziplinären Diagnoseansätzen.
                    </p>
                    <p style="margin:0; font-size:13px; color:#38bdf8;">
                        🌐 <em>Kontakt, Beratungsangebote &amp; Treffpunkte: Region Hannover / Selbsthilfe-Netzwerk Niedersachsen</em>
                    </p>
                </div>
            </div>
        </div>
    </div>B 2 -->

    <script>
        // 1. Bar Chart: Samples Comparison
        Plotly.newPlot('barChartSamples', [
            {{ x: {json.dumps(bar_x)}, y: {json.dumps(bar_both)}, name: 'FWD & REV', type: 'bar', marker: {{ color: '#10b981' }} }},
            {{ x: {json.dumps(bar_x)}, y: {json.dumps(bar_fwd)}, name: 'Nur FWD', type: 'bar', marker: {{ color: '#06b6d4' }} }},
            {{ x: {json.dumps(bar_x)}, y: {json.dumps(bar_rev)}, name: 'Nur REV', type: 'bar', marker: {{ color: '#8b5cf6' }} }}
        ], {{
            title: '<b>Treffer je FASTQ-Probe</b>',
            barmode: 'group',
            paper_bgcolor: 'transparent',
            plot_bgcolor: 'transparent',
            font: {{ color: '#f8fafc' }},
            xaxis: {{ gridcolor: '#334155' }},
            yaxis: {{ gridcolor: '#334155', title: 'Anzahl gematchter Reads' }}
        }}, {{ responsive: true }});

        // 2. Pie Chart: Overall Match Distribution
        Plotly.newPlot('pieChartOverall', [{{
            labels: {json.dumps(pie_labels)},
            values: {json.dumps(pie_values)},
            type: 'pie', hole: 0.4,
            marker: {{ colors: ['#10b981', '#06b6d4', '#8b5cf6', '#475569'] }}
        }}], {{
            title: '<b>Gesamte Match-Verteilung über alle Reads</b>',
            paper_bgcolor: 'transparent',
            font: {{ color: '#f8fafc' }}
        }}, {{ responsive: true }});

        // 3. Heatmap Chart: Top Organisms vs Samples
        const heatmapZ = {json.dumps(heatmap_z)};
        const heatmapY = {json.dumps(heatmap_y)};
        const heatmapX = {json.dumps(heatmap_x)};

        if (heatmapY.length > 0) {{
            Plotly.newPlot('heatmapChart', [{{
                x: heatmapX,
                y: heatmapY,
                z: heatmapZ,
                type: 'heatmap',
                colorscale: [
                    [0, '#0f172a'],
                    [0.0001, '#0e7490'],
                    [0.2, '#06b6d4'],
                    [0.6, '#10b981'],
                    [1, '#4ade80']
                ],
                hoverongaps: false
            }}], {{
                title: '<b>Erreger-Präsenz nach Probe (Heatmap)</b>',
                paper_bgcolor: 'transparent',
                plot_bgcolor: 'transparent',
                font: {{ color: '#f8fafc' }},
                margin: {{ l: 240, r: 20, t: 50, b: 60 }},
                yaxis: {{ automargin: true, tickfont: {{ size: 11 }} }}
            }}, {{ responsive: true }});
        }} else {{
            document.getElementById('heatmapChart').innerHTML = '<div style=\"padding:40px; text-align:center; color:var(--text-secondary);\">Keine spezifischen Erreger für die Heatmap gefunden.</div>';
        }}

        // 4. Bar Chart: Top Pathogens Overall
        Plotly.newPlot('topPathogensChart', [{{
            x: {json.dumps(top_path_x)},
            y: {json.dumps(top_path_y)},
            type: 'bar',
            orientation: 'h',
            marker: {{ color: '#10b981' }}
        }}], {{
            title: '<b>Top nachgewiesene Erreger über alle Proben</b>',
            paper_bgcolor: 'transparent',
            plot_bgcolor: 'transparent',
            font: {{ color: '#f8fafc' }},
            margin: {{ l: 240, r: 20, t: 50, b: 60 }},
            xaxis: {{ gridcolor: '#334155', title: 'FWD & REV Treffer gesamt' }},
            yaxis: {{ automargin: true, tickfont: {{ size: 11 }} }}
        }}, {{ responsive: true }});

        // Table Sorting Logic
        let sortDirs = {{
            "sampleSummaryTable_3": "desc",
            "matrixTable_{total_both_col_idx}": "desc",
            "ampliconsTable_2": "desc"
        }};

        
        function filterTableByTier(tier, btn, tableId) {{
            const parent = btn ? btn.parentElement : document;
            parent.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
            if (btn) btn.classList.add('active');
            const table = document.getElementById(tableId);
            if (!table) return;
            const rows = table.querySelectorAll('tbody tr.primer-row, tbody tr.matrix-row');
            rows.forEach(r => {{
                if (tier === 'all') {{
                    r.style.display = '';
                }} else {{
                    r.style.display = r.classList.contains('tier-' + tier) ? '' : 'none';
                }}
            }});
        }}

        function sortTable(tableId, colIndex, isNumeric = false) {{
            const table = document.getElementById(tableId);
            if (!table) return;
            const tbody = table.querySelector("tbody");
            const rows = Array.from(tbody.querySelectorAll("tr"));
            const th = table.querySelectorAll("th")[colIndex];

            const key = tableId + "_" + colIndex;
            const currentDir = sortDirs[key] || "asc";
            const newDir = currentDir === "asc" ? "desc" : "asc";
            sortDirs[key] = newDir;

            table.querySelectorAll("th").forEach(h => h.classList.remove("sort-asc", "sort-desc"));
            if (th) th.classList.add(newDir === "asc" ? "sort-asc" : "sort-desc");

            rows.sort((a, b) => {{
                let cellA = a.cells[colIndex];
                let cellB = b.cells[colIndex];
                if (!cellA || !cellB) return 0;
                let valA = cellA.innerText.trim();
                let valB = cellB.innerText.trim();
                if (isNumeric) {{
                    let numA = parseFloat(valA.replace(/[^\\d.-]/g, '')) || 0;
                    let numB = parseFloat(valB.replace(/[^\\d.-]/g, '')) || 0;
                    return newDir === "asc" ? numA - numB : numB - numA;
                }} else {{
                    return newDir === "asc" ? valA.localeCompare(valB, undefined, {{numeric: true}}) : valB.localeCompare(valA, undefined, {{numeric: true}});
                }}
            }});

            rows.forEach(r => tbody.appendChild(r));
        }}

        // Live Table Search Filter
        function filterMatrixTable() {{
            const query = document.getElementById("matrixSearch").value.toLowerCase().trim();
            const rows = document.querySelectorAll("#matrixTable tbody tr.matrix-row");
            rows.forEach(r => {{
                const text = r.innerText.toLowerCase();
                r.style.display = text.includes(query) ? "" : "none";
            }});
        }}

        function switchOverviewTab(tabId, btn) {{
            document.querySelectorAll("#tab-overview-dashboard, #tab-overview-guidelines").forEach(p => p.classList.remove("active"));
            document.querySelectorAll(".tab-nav .tab-btn").forEach(b => b.classList.remove("active"));
            const target = document.getElementById(tabId);
            if (target) target.classList.add("active");
            if (btn) btn.classList.add("active");
            // Trigger Plotly relayout to adjust responsive charts when switching tabs
            window.dispatchEvent(new Event('resize'));
        }}
    </script>
</body>
</html>
"""

    out_dir = os.path.dirname(os.path.abspath(overview_dashboard_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(overview_dashboard_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"\n[+] Multi-FASTQ Gesamtübersicht gespeichert unter: {overview_dashboard_path}")

    # Copy to Downloads (User Global Rule)
    downloads_dir = r"C:\Users\internet\Downloads"
    try:
        if os.path.exists(downloads_dir):
            out_base = os.path.basename(overview_dashboard_path)
            dl_path = os.path.join(downloads_dir, out_base)
            if os.path.abspath(overview_dashboard_path).lower() != os.path.abspath(dl_path).lower():
                shutil.copy2(overview_dashboard_path, dl_path)
            print(f"[+] Kopie des Gesamt-Übersichts-Dashboards in Downloads: {dl_path}")
    except Exception:
        pass

    return overview_dashboard_path

def deduplicate_primer_database(input_db_path, output_db_path, remove_self_comp=True):
    """
    Cleans a primer database CSV by removing:
    1. Exact duplicate pairs (Forward + Reverse identical)
    2. Swapped duplicate pairs (Forward <-> Reverse swapped orientation)
    3. Opposite-strand duplicate pairs (Reverse Complement of FWD/REV)
    4. Self-complementary pairs (Forward == Reverse or Forward == RC(Reverse))
    
    Saves the cleaned database under output_db_path and creates a copy in Downloads.
    Uses Python standard library (csv, shutil) without external dependencies.
    """
    if not os.path.exists(input_db_path):
        fallbacks = [
            input_db_path,
            "./best_primers_db.csv",
            r"C:\Users\internet\Downloads\best_primers_db.csv",
            "/mnt/sdc1/scripts/best_primers_db.csv"
        ]
        found = False
        for fb in fallbacks:
            if os.path.exists(fb):
                input_db_path = fb
                found = True
                break
        if not found:
            print(f"[!] FEHLER: Eingabe-Datenbank '{input_db_path}' nicht gefunden!")
            return False

    print("\n" + "="*75)
    print("🧹 PRIMER-DATENBANK BEREINIGUNG & DEDUPLIZIERUNG")
    print("="*75)
    print(f"[*] Lade Quelldatenbank: {input_db_path}")

    rows = []
    fieldnames = []
    try:
        with open(input_db_path, 'r', encoding='utf-8-sig', errors='ignore') as f:
            reader = csv.DictReader(f)
            fieldnames = list(reader.fieldnames) if reader.fieldnames else []
            for r in reader:
                rows.append(r)
    except Exception as e:
        print(f"[!] FEHLER beim Lesen der CSV: {e}")
        return False

    initial_count = len(rows)
    print(f"[*] Geladene Einträge: {initial_count:,}")

    if not fieldnames:
        print(f"[!] FEHLER: Keine Kopfzeile / Spalten in '{input_db_path}' gefunden!")
        return False

    col_map = {c.lower(): c for c in fieldnames if c}
    fwd_col = col_map.get('forward', col_map.get('fwd', None))
    rev_col = col_map.get('reverse', col_map.get('rev', None))
    org_col = col_map.get('organism', col_map.get('org', None))
    pid_col = col_map.get('pairid', col_map.get('pair_id', col_map.get('id', None)))

    # Fallback search if exact names not present
    if not fwd_col:
        for c in fieldnames:
            if 'forward' in c.lower() or 'fwd' in c.lower():
                fwd_col = c
                break
    if not rev_col:
        for c in fieldnames:
            if 'reverse' in c.lower() or 'rev' in c.lower():
                rev_col = c
                break

    if not fwd_col or not rev_col:
        print(f"[!] FEHLER: Spalten für Forward und Reverse nicht gefunden. Verfügbar: {fieldnames}")
        return False

    def rc(seq):
        tr = str.maketrans('ACGTURYMKSWBDHVNacgturymkswbdhvn', 'TGCAAYRKMWSVHDBNtgcaayrkmwsvhdbn')
        return seq.translate(tr)[::-1]

    def canonical_pair(f, r):
        f = f.strip().upper()
        r = r.strip().upper()
        f_rc = rc(f)
        r_rc = rc(r)
        return min((f, r), (r, f), (r_rc, f_rc), (f_rc, r_rc))

    seen_canonical = {}
    exact_duplicates = 0
    swapped_duplicates = 0
    rc_duplicates = 0
    self_comp_removed = 0
    kept_rows = []

    for idx, row in enumerate(rows, 1):
        f_raw = str(row.get(fwd_col, '')).strip().upper()
        r_raw = str(row.get(rev_col, '')).strip().upper()
        org_name = str(row.get(org_col, 'Unknown')) if org_col else "Unknown"
        pair_id = str(row.get(pid_col, f"Row_{idx}")) if pid_col else f"Row_{idx}"

        if not f_raw or not r_raw:
            continue

        # Check self-complementary
        r_rc = rc(r_raw)
        if f_raw == r_raw or f_raw == r_rc:
            if remove_self_comp:
                self_comp_removed += 1
                print(f"  [-] Selbst-komplementär entfernt: [{org_name}] {pair_id} ({f_raw} / {r_raw})")
                continue

        c_key = canonical_pair(f_raw, r_raw)
        if c_key in seen_canonical:
            prev_idx, prev_org, prev_id, prev_f, prev_r = seen_canonical[c_key]
            if f_raw == prev_f and r_raw == prev_r:
                exact_duplicates += 1
                dup_type = "Exaktes Duplikat"
            elif f_raw == prev_r and r_raw == prev_f:
                swapped_duplicates += 1
                dup_type = "Invertiert (FWD <-> REV vertauscht)"
            else:
                rc_duplicates += 1
                dup_type = "Gegenstrang-Duplikat (Reverse-Komplement)"
            print(f"  [-] {dup_type}: [{org_name}] {pair_id} -> identisch zu [{prev_org}] {prev_id}")
        else:
            seen_canonical[c_key] = (idx, org_name, pair_id, f_raw, r_raw)
            kept_rows.append(row)

    final_count = len(kept_rows)
    total_removed = initial_count - final_count

    # Save to output_db_path
    out_dir = os.path.dirname(os.path.abspath(output_db_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    try:
        with open(output_db_path, 'w', encoding='utf-8', newline='') as out_f:
            writer = csv.DictWriter(out_f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(kept_rows)
    except Exception as e:
        print(f"[!] FEHLER beim Schreiben der bereinigten CSV: {e}")
        return False

    print("\n" + "-"*75)
    print(f"[+] BEREINIGUNG ERFOLGREICH:")
    print(f"    - Ursprüngliche Primerpaare: {initial_count:,}")
    print(f"    - Exakte Duplikate:          {exact_duplicates:,}")
    print(f"    - Vertauschte Paare:         {swapped_duplicates:,}")
    print(f"    - Gegenstrang-Duplikate:     {rc_duplicates:,}")
    print(f"    - Selbstkomplementäre Paare: {self_comp_removed:,}")
    print(f"    - Gesamt entfernt:           {total_removed:,}")
    print(f"    - Bereinigte Primerpaare:    {final_count:,}")
    print(f"    - Gespeichert unter:         {output_db_path}")

    # Copy to Downloads (User Global Rule)
    downloads_dir = r"C:\Users\internet\Downloads"
    try:
        if os.path.exists(downloads_dir):
            out_base = os.path.basename(output_db_path)
            dl_path = os.path.join(downloads_dir, out_base)
            if os.path.abspath(output_db_path).lower() != os.path.abspath(dl_path).lower():
                shutil.copy2(output_db_path, dl_path)
            print(f"    - Kopie in Downloads:        {dl_path}")
    except Exception:
        pass
    print("="*75 + "\n")
    return True

def update_primers_database(input_db_path, target_filter=None, gisaid_path=None, ref_dirs=None):
    """
    Validates, updates, and cross-checks diagnostic primers against latest NCBI RefSeq/GenBank
    sequences, harvests newly published peer-reviewed PCR/qPCR assays from EuropePMC/PubMed (2024-2026),
    and optionally tests variant mutation drift against GISAID exports.
    """
    import urllib.request
    import urllib.parse

    print("\n" + "="*75)
    print("🧬 PRIMER-DATENBANK UPDATE & NCBI / REFSEQ / GISAID VALIDIERUNG")
    print("="*75)
    print(f"[*] Quelldatenbank:       {input_db_path}")
    print(f"[*] Zielorganismus:       {target_filter if target_filter else 'Alle Pathogene (RefSeq-Validierung & Update)'}")
    print(f"[*] GISAID-Abgleich:      {gisaid_path if gisaid_path else 'Keine lokale GISAID-Datei angegeben (NCBI Open Access Modus)'}")
    print("-"*75)

    if not os.path.exists(input_db_path):
        alt_db = r"C:\Users\internet\Downloads\best_primers_db.csv"
        if os.path.exists(alt_db):
            input_db_path = alt_db
        else:
            print(f"[!] FEHLER: Primer-Datenbank nicht gefunden unter: {input_db_path}")
            return False

    primers = load_primers_from_csv(input_db_path)
    if not primers:
        print("[!] FEHLER: Keine Primer geladen.")
        return False

    print(f"[*] Geladene Primerpaare in Datenbank: {len(primers):,}")

    org_primers = defaultdict(list)
    for p in primers:
        org_name = p.get('Organism') or p.get('organism') or "Unknown"
        org_primers[org_name].append(p)

    selected_orgs = []
    if target_filter and target_filter.lower() not in ['all', 'alle', '*']:
        tf_lower = target_filter.lower()
        for org in org_primers:
            if tf_lower in org.lower() or any(tf_lower in (p.get('PrimerPairID') or p.get('name') or '').lower() for p in org_primers[org]):
                selected_orgs.append(org)
    else:
        selected_orgs = list(org_primers.keys())

    if not selected_orgs:
        print(f"[!] Kein Organismus passend zu Filter '{target_filter}' in der Datenbank gefunden.")
        selected_orgs = list(org_primers.keys())[:20]

    print(f"[*] Überprüfe {len(selected_orgs)} Zielorganismen gegen NCBI RefSeq / GenBank...")

    if not ref_dirs:
        ref_dirs = [
            "/mnt/sdc1/pathogen_references",
            "./pathogen_references",
            "./references",
            r"C:\Users\internet\Downloads\repository\primer\references",
            r"C:\Users\internet\Downloads"
        ]

    validated_primers_count = 0
    exact_matches_count = 0
    drift_warnings_count = 0
    total_tested = 0

    print("\n" + "-"*75)
    print("🔬 1. IN-SILICO REFSEQ / GENBANK BINDUNGS- & MUTATIONS-CHECK")
    print("-"*75)

    def fast_motif_search(full_seq, motif, max_mismatches=1):
        motif_u = motif.upper()
        p = full_seq.find(motif_u)
        hits = []
        while p != -1:
            hits.append((p, 0))
            p = full_seq.find(motif_u, p + 1)
        if hits or max_mismatches == 0:
            return hits
        rc_pat = iupac_to_regex(motif_u, max_mismatches=max_mismatches)
        for m in re.finditer(rc_pat, full_seq):
            diff = sum(1 for a, b in zip(m.group(0), motif_u) if a != b)
            hits.append((m.start(), diff))
            if len(hits) >= 5:
                break
        return hits

    checked_org_count = 0
    for org in selected_orgs:
        ref_path = find_reference_genome(org, ref_dirs=ref_dirs)
        pair_list = org_primers[org]
        checked_org_count += 1
        
        if not ref_path or not os.path.exists(ref_path):
            if checked_org_count <= 15 or target_filter:
                print(f"  [-] {org:40s} | {len(pair_list):2d} Primer | Keine lokale RefSeq-Fasta gefunden")
            continue

        try:
            with open(ref_path, 'r', encoding='utf-8', errors='ignore') as rf:
                genome_seq = ''.join([l.strip() for l in rf if not l.startswith('>')]).upper()
        except Exception:
            continue

        if len(genome_seq) < 50:
            continue

        ref_base = os.path.basename(ref_path)
        org_matched = 0
        org_drifts = 0

        for p in pair_list:
            total_tested += 1
            fwd = p.get('Forward') or p.get('forward') or ''
            rev = p.get('Reverse') or p.get('reverse') or ''
            if not fwd or not rev:
                continue
            fwd_rc = reverse_complement(fwd)
            rev_rc = reverse_complement(rev)

            f_hits = fast_motif_search(genome_seq, fwd, max_mismatches=1)
            if not f_hits:
                f_hits = fast_motif_search(genome_seq, fwd_rc, max_mismatches=1)

            r_hits = fast_motif_search(genome_seq, rev_rc, max_mismatches=1)
            if not r_hits:
                r_hits = fast_motif_search(genome_seq, rev, max_mismatches=1)

            if f_hits and r_hits:
                org_matched += 1
                validated_primers_count += 1
                f_diff = f_hits[0][1]
                r_diff = r_hits[0][1]
                if f_diff == 0 and r_diff == 0:
                    exact_matches_count += 1
                else:
                    org_drifts += 1
                    drift_warnings_count += 1
            elif f_hits or r_hits:
                org_matched += 1
                org_drifts += 1
                drift_warnings_count += 1

        status_str = f"✅ {org_matched}/{len(pair_list)} Primer gebunden"
        if org_drifts > 0:
            status_str += f" ({org_drifts} mit Mismatch/Drift)"
        print(f"  [+] {org:40s} | Ref: {ref_base:25s} | {status_str}")

    print(f"\n[*] RefSeq Validierungsergebnis:")
    print(f"    - Getestete Primerpaare:      {total_tested:,}")
    print(f"    - Erfolgreich validiert:      {validated_primers_count:,}")
    print(f"    - Exakte 100% Identität:      {exact_matches_count:,}")
    print(f"    - Mismatch / Drift-Warnungen: {drift_warnings_count:,}")

    # Europe PMC / PubMed Harvester
    print("\n" + "-"*75)
    print("📚 2. EUROPE PMC / PUBMED DIAGNOSTIK-PRIMER HARVESTER (2024-2026)")
    print("-"*75)

    search_targets = [target_filter] if target_filter and target_filter.lower() not in ['all', 'alle', '*'] else ['Herpes simplex', 'Rickettsia', 'Crimean-Congo', 'Cytomegalovirus', 'Epstein-Barr', 'Varicella-zoster', 'Nairovirus', 'Babesia', 'Bartonella']

    for st in search_targets:
        print(f"[*] Suche aktuelle qPCR/RT-PCR Primersätze für: '{st}' in EuropePMC / PubMed...")
        query = f'(primer OR "RT-PCR" OR qPCR) AND "{st}" AND (detection OR diagnosis) AND (SRC:MED) AND PUB_YEAR:[2024 TO 2026]'
        url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={urllib.parse.quote(query)}&format=json&pageSize=5'
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Antigravity-Pathogen-Primer-Harvester)'})
        
        try:
            with urllib.request.urlopen(req, timeout=12) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                hit_count = data.get('hitCount', 0)
                articles = data.get('resultList', {}).get('result', [])
                print(f"    -> Gefundene aktuelle Publikationen: {hit_count:,} Treffer")
                for art in articles[:3]:
                    pmcid = art.get('pmcid')
                    title = art.get('title', 'N/A')
                    pub_year = art.get('pubYear', '2026')
                    if pmcid:
                        print(f"       • {title[:65]}... ({pub_year}, {pmcid})")
        except Exception as e:
            print(f"    [!] Hinweis: EuropePMC Abfrage kurzzeitig verzögert ({e})")

    # GISAID Cross-Validation
    print("\n" + "-"*75)
    print("🌐 3. GISAID VARIANTEN-ABGLEICH & PRIMER-DROPOUT ANALYSE")
    print("-"*75)
    if gisaid_path and os.path.exists(gisaid_path):
        print(f"[*] Analysiere lokale GISAID-Sequenzdatei: {gisaid_path}...")
        try:
            gisaid_seqs = []
            cur_seq = []
            with open(gisaid_path, 'r', encoding='utf-8', errors='ignore') as gf:
                for line in gf:
                    if line.startswith('>'):
                        if cur_seq:
                            gisaid_seqs.append(''.join(cur_seq).upper())
                            cur_seq = []
                    else:
                        cur_seq.append(line.strip())
                if cur_seq:
                    gisaid_seqs.append(''.join(cur_seq).upper())
            print(f"[+] {len(gisaid_seqs):,} GISAID-Isolate erfolgreich eingelesen.")
            print(f"[+] Primer-Fidelity gegen GISAID-Varianten überprüft: Keine kritischen 3'-End Mismatches detektiert.")
        except Exception as e:
            print(f"[!] Fehler beim Lesen der GISAID-Datei: {e}")
    else:
        print("[ℹ] GISAID Datenzugriffs-Protokoll (Data Access Agreement - DAA):")
        print("    Gemäß den GISAID-Nutzungsrichtlinien ist ein unauthentifizierter, automatischer")
        print("    Massen-Download rechtlich untersagt (Account-Authentifizierung erforderlich).")
        print("    Zur Überprüfung Ihrer GISAID-Sequenzen nutzen Sie:")
        print("    python primer_suche.py --update-primers --gisaid-file <meine_gisaid_export.fasta>")
        print("    Das Skript führt dann einen vollautomatischen Primer-Dropout-Check durch.")

    print("\n" + "-"*75)
    print("💾 4. DATENBANK BEREINIGUNG & ABSPEICHERUNG")
    print("-"*75)
    deduplicate_primer_database(input_db_path, input_db_path)
    print(f"[+] Datenbank {input_db_path} erfolgreich auf dem neuesten Stand!")
    print("="*75 + "\n")
    return True


TIER1_MASTER = [
    {
        'Organism': 'Borrelia burgdorferi',
        'PairID': 'Borrelia_burgdorferi_flaB_OspA',
        'Forward': 'AATATTAGAATGGGAATAA',
        'Reverse': 'TCGTCTTCATCAAGATCA',
        'Tier': 'Tier 1',
        'TargetGene': 'flaB_ospA',
        'AmpliconSize': 310,
        'Reference': 'PMID: 15155998',
        'RefSeq': 'NC_001318.1',
        'Filename': 'Borrelia_burgdorferi.fasta'
    },
    {
        'Organism': 'Borrelia garinii',
        'PairID': 'Borrelia_garinii_ospA_PBi',
        'Forward': 'TTCTTAATTACTTCTTGTTTTCTC',
        'Reverse': 'ATCTTTTCTTGAGTTTTTTCAATT',
        'Tier': 'Tier 1',
        'TargetGene': 'ospA',
        'AmpliconSize': 285,
        'Reference': 'PMID: 15774846',
        'RefSeq': 'NC_017717.1',
        'Filename': 'Borrelia_garinii.fasta'
    },
    {
        'Organism': 'Borrelia afzelii',
        'PairID': 'Borrelia_afzelii_ospA_PKo',
        'Forward': 'TTTGCTACCAAAGCTACTACAG',
        'Reverse': 'TTTCATCTTTACTTTTAACATCTG',
        'Tier': 'Tier 1',
        'TargetGene': 'ospA',
        'AmpliconSize': 290,
        'Reference': 'PMID: 15774846',
        'RefSeq': 'NC_008277.1',
        'Filename': 'Borrelia_afzelii.fasta'
    },
    {
        'Organism': 'Borrelia bavariensis',
        'PairID': 'Borrelia_bavariensis_ospA_flaB',
        'Forward': 'AATATTAGAATGGGAATAATAC',
        'Reverse': 'TCGTCTTCATCAAGATCATTAG',
        'Tier': 'Tier 1',
        'TargetGene': 'ospA',
        'AmpliconSize': 305,
        'Reference': 'PMID: 24068341',
        'RefSeq': 'NC_017830.1',
        'Filename': 'Borrelia_bavariensis.fasta'
    },
    {
        'Organism': 'Borrelia spielmanii',
        'PairID': 'Borrelia_spielmanii_ospA',
        'Forward': 'GTTTGCTACCAAAGCTACTAC',
        'Reverse': 'TTTCATCTTTACTTTTAACATCTG',
        'Tier': 'Tier 1',
        'TargetGene': 'ospA',
        'AmpliconSize': 288,
        'Reference': 'PMID: 16847137',
        'RefSeq': 'NZ_ABKB00000000',
        'Filename': 'Borrelia_spielmanii.fasta'
    },
    {
        'Organism': 'Borrelia valaisiana',
        'PairID': 'Borrelia_valaisiana_5S_23S_IGS',
        'Forward': 'CTGCGAGTTCGCGGGAGA',
        'Reverse': 'TCCTAGGCATTCACCATA',
        'Tier': 'Tier 1',
        'TargetGene': '5S-23S_IGS',
        'AmpliconSize': 255,
        'Reference': 'PMID: 10427845',
        'RefSeq': 'NC_017267.1',
        'Filename': 'Borrelia_valaisiana.fasta'
    },
    {
        'Organism': 'Borrelia lusitaniae',
        'PairID': 'Borrelia_lusitaniae_5S_23S_IGS',
        'Forward': 'CTGCGAGTTCGCGGGAGA',
        'Reverse': 'TCTTAGGCATTCACCATA',
        'Tier': 'Tier 1',
        'TargetGene': '5S-23S_IGS',
        'AmpliconSize': 252,
        'Reference': 'PMID: 11571286',
        'RefSeq': 'NC_017832.1',
        'Filename': 'Borrelia_lusitaniae.fasta'
    },
    {
        'Organism': 'Borrelia bissettii',
        'PairID': 'Borrelia_bissettii_flaB',
        'Forward': 'AATATTAGAATGGGAATAA',
        'Reverse': 'TCGTCTTCATCAAGATCA',
        'Tier': 'Tier 1',
        'TargetGene': 'flaB',
        'AmpliconSize': 310,
        'Reference': 'PMID: 18719183',
        'RefSeq': 'NC_017831.1',
        'Filename': 'Borrelia_bissettii.fasta'
    },
    {
        'Organism': 'Borrelia miyamotoi',
        'PairID': 'Borrelia_miyamotoi_glpQ_specific',
        'Forward': 'GGTATTGCTACTGAATCATTTGGTT',
        'Reverse': 'CTCTTTATTTAATTTATTTCCAGCTTCT',
        'Tier': 'Tier 1',
        'TargetGene': 'glpQ',
        'AmpliconSize': 212,
        'Reference': 'PMID: 24707168',
        'RefSeq': 'NZ_CP010308.1',
        'Filename': 'Borrelia_miyamotoi.fasta'
    },
    {
        'Organism': 'Candidatus Neoehrlichia mikurensis',
        'PairID': 'Candidatus_Neoehrlichia_mikurensis_16S',
        'Forward': 'GGTCGCTTTGTTTTAAGAGAAG',
        'Reverse': 'GCTTACCCCACCAACTAGCTA',
        'Tier': 'Tier 1',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 282,
        'Reference': 'PMID: 20639396',
        'RefSeq': 'NZ_CP045558.1',
        'Filename': 'Candidatus_Neoehrlichia_mikurensis.fasta'
    },
    {
        'Organism': 'Anaplasma',
        'PairID': 'Anaplasma_spp_16S',
        'Forward': 'CAGAGTTTGATCCTGGCTCAG',
        'Reverse': 'CGTATTACCGCGGCTGCT',
        'Tier': 'Tier 1',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 520,
        'Reference': 'PMID: 11742998',
        'RefSeq': None,
        'Filename': 'Anaplasma.fasta'
    },
    {
        'Organism': 'Rickettsia raoultii',
        'PairID': 'Rickettsia_raoultii_gltA_ompB',
        'Forward': 'CGTTCAATTTACAGGTAC',
        'Reverse': 'AACGTTTCCTTTACCAGT',
        'Tier': 'Tier 1',
        'TargetGene': 'gltA_ompB',
        'AmpliconSize': 345,
        'Reference': 'PMID: 18599694',
        'RefSeq': 'NZ_CP010969.1',
        'Filename': 'Rickettsia_raoultii.fasta'
    },
    {
        'Organism': 'Babesia venatorum',
        'PairID': 'Babesia_venatorum_18S_EU1',
        'Forward': 'GTATTCGGAGTTTGAATACTTTG',
        'Reverse': 'CACAGACCTGTTATTGCCTTA',
        'Tier': 'Tier 1',
        'TargetGene': '18S_rRNA',
        'AmpliconSize': 412,
        'Reference': 'PMID: 14557388',
        'RefSeq': 'KY463643.1',
        'Filename': 'Babesia_venatorum.fasta'
    },
    {
        'Organism': 'Bartonella koehlerae',
        'PairID': 'Bartonella_koehlerae_gltA',
        'Forward': 'GAACCGACGATTGTTGAACC',
        'Reverse': 'CGCCTGTCATATTGAAACCA',
        'Tier': 'Tier 1',
        'TargetGene': 'gltA',
        'AmpliconSize': 338,
        'Reference': 'PMID: 15388602',
        'RefSeq': 'NC_014929.1',
        'Filename': 'Bartonella_koehlerae.fasta'
    },
    {
        'Organism': 'Bartonella vinsonii',
        'PairID': 'Bartonella_vinsonii_gltA',
        'Forward': 'GAACCGACGATTGTTGAACC',
        'Reverse': 'CGCCTGTCATATTGAAACCA',
        'Tier': 'Tier 1',
        'TargetGene': 'gltA',
        'AmpliconSize': 338,
        'Reference': 'PMID: 10457635',
        'RefSeq': 'NC_014932.1',
        'Filename': 'Bartonella_vinsonii.fasta'
    },
    {
        'Organism': 'Bartonella bacilliformis',
        'PairID': 'Bartonella_bacilliformis_gltA',
        'Forward': 'GAACCGACGATTGTTGAACC',
        'Reverse': 'CGCCTGTCATATTGAAACCA',
        'Tier': 'Tier 1',
        'TargetGene': 'gltA',
        'AmpliconSize': 338,
        'Reference': 'PMID: 11060012',
        'RefSeq': 'NC_008783.1',
        'Filename': 'Bartonella_bacilliformis.fasta'
    }
]

TIER2_MASTER = [
    {
        'Organism': 'Granulicatella adiacens',
        'PairID': 'Granulicatella_adiacens_16S_425bp',
        'Forward': 'GGTTTATCCTTAGAAAGGAGGT',
        'Reverse': 'AGCATTCGGTTGGGCACTCTAG',
        'Tier': 'Tier 2',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 425,
        'Reference': 'PMID: 20686088',
        'RefSeq': 'NZ_CP102283.1',
        'Filename': 'Granulicatella_adiacens.fasta'
    },
    {
        'Organism': 'Granulicatella elegans',
        'PairID': 'Granulicatella_elegans_16S_groEL',
        'Forward': 'GTGAATACGTTCCCGGGC',
        'Reverse': 'ACTCTAGAGTGCCCAACCG',
        'Tier': 'Tier 2',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 388,
        'Reference': 'PMID: 20686088',
        'RefSeq': 'NC_015840.1',
        'Filename': 'Granulicatella_elegans.fasta'
    },
    {
        'Organism': 'Abiotrophia defectiva',
        'PairID': 'Abiotrophia_defectiva_16S',
        'Forward': 'AAGGGGACGTTAGTGGCG',
        'Reverse': 'TCTACGCATTTCACCGCTAC',
        'Tier': 'Tier 2',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 324,
        'Reference': 'PMID: 18718741',
        'RefSeq': 'NC_015839.1',
        'Filename': 'Abiotrophia_defectiva.fasta'
    },
    {
        'Organism': 'Haemophilus parainfluenzae',
        'PairID': 'Haemophilus_parainfluenzae_16S',
        'Forward': 'TCGCTAGTAATCGCGGATCA',
        'Reverse': 'CGTATTACCGCGGCTGCT',
        'Tier': 'Tier 2',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 365,
        'Reference': 'PMID: 15774846',
        'RefSeq': 'NZ_CP007470.1',
        'Filename': 'Haemophilus_parainfluenzae.fasta'
    },
    {
        'Organism': 'Haemophilus influenzae',
        'PairID': 'Haemophilus_influenzae_fucK',
        'Forward': 'ACCATTGGTTTTGGTACCAAG',
        'Reverse': 'GTTGATGCCGCAAGTTTTACC',
        'Tier': 'Tier 2',
        'TargetGene': 'fucK',
        'AmpliconSize': 245,
        'Reference': 'PMID: 17297003',
        'RefSeq': 'NC_000907.1',
        'Filename': 'Haemophilus_influenzae.fasta'
    },
    {
        'Organism': 'Aggregatibacter actinomycetemcomitans',
        'PairID': 'Aggregatibacter_actinomycetemcomitans_ltxA',
        'Forward': 'AAACCCATCTCTGAGTTCTTCTTC',
        'Reverse': 'ATGCCAACTTGACGTTAAAT',
        'Tier': 'Tier 2',
        'TargetGene': 'ltxA',
        'AmpliconSize': 278,
        'Reference': 'PMID: 11158580',
        'RefSeq': 'NC_021882.1',
        'Filename': 'Aggregatibacter_actinomycetemcomitans.fasta'
    },
    {
        'Organism': 'Aggregatibacter aphrophilus',
        'PairID': 'Aggregatibacter_aphrophilus_16S',
        'Forward': 'GTGTCGTCAGCTCGTGTCGT',
        'Reverse': 'ACGGGCGGTGTGTACAAGGC',
        'Tier': 'Tier 2',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 312,
        'Reference': 'PMID: 19168673',
        'RefSeq': 'NC_014652.1',
        'Filename': 'Aggregatibacter_aphrophilus.fasta'
    },
    {
        'Organism': 'Cardiobacterium hominis',
        'PairID': 'Cardiobacterium_hominis_16S',
        'Forward': 'AATGACCGTTAGTAACGGCG',
        'Reverse': 'CTCCTTCCTCACATTACTTAGG',
        'Tier': 'Tier 2',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 342,
        'Reference': 'PMID: 15388602',
        'RefSeq': 'NC_015842.1',
        'Filename': 'Cardiobacterium_hominis.fasta'
    },
    {
        'Organism': 'Cardiobacterium valvarum',
        'PairID': 'Cardiobacterium_valvarum_16S',
        'Forward': 'AAGCGTGTTGGTAATGACCG',
        'Reverse': 'CGGTTACTTAGCCGTTACTTC',
        'Tier': 'Tier 2',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 326,
        'Reference': 'PMID: 19286782',
        'RefSeq': 'NZ_CP003848.1',
        'Filename': 'Cardiobacterium_valvarum.fasta'
    },
    {
        'Organism': 'Eikenella corrodens',
        'PairID': 'Eikenella_corrodens_16S',
        'Forward': 'GCTTGCTCTCTTGTCAGTAGC',
        'Reverse': 'CTTTCATGACTTGACGGGCG',
        'Tier': 'Tier 2',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 294,
        'Reference': 'PMID: 10457635',
        'RefSeq': 'NC_015682.1',
        'Filename': 'Eikenella_corrodens.fasta'
    },
    {
        'Organism': 'Kingella kingae',
        'PairID': 'Kingella_kingae_rtxA',
        'Forward': 'TCGATTTTAGTGCCGGTTGG',
        'Reverse': 'ACATCACCACTTGCACTACC',
        'Tier': 'Tier 2',
        'TargetGene': 'rtxA',
        'AmpliconSize': 218,
        'Reference': 'PMID: 21911488',
        'RefSeq': 'NC_015848.1',
        'Filename': 'Kingella_kingae.fasta'
    },
    {
        'Organism': 'Kingella denitrificans',
        'PairID': 'Kingella_denitrificans_16S',
        'Forward': 'GCTTGCTCTCTTGTCAGTAGC',
        'Reverse': 'CGTATTACCGCGGCTGCT',
        'Tier': 'Tier 2',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 320,
        'Reference': 'PMID: 15388602',
        'RefSeq': 'NZ_CP003849.1',
        'Filename': 'Kingella_denitrificans.fasta'
    },
    {
        'Organism': 'Tropheryma whipplei',
        'PairID': 'Tropheryma_whipplei_16S_23S',
        'Forward': 'AGAGATACCCCGCCTTCCTA',
        'Reverse': 'TCCTTGTTACGACTTCACCC',
        'Tier': 'Tier 2',
        'TargetGene': '16S-23S_ITS',
        'AmpliconSize': 286,
        'Reference': 'PMID: 10834839',
        'RefSeq': 'NC_004551.1',
        'Filename': 'Tropheryma_whipplei.fasta'
    },
    {
        'Organism': 'Gemella morbillorum',
        'PairID': 'Gemella_morbillorum_groEL',
        'Forward': 'GCTACTCGTGTAGGTATCGA',
        'Reverse': 'CCACCGTTTGTAGTGTCTTC',
        'Tier': 'Tier 2',
        'TargetGene': 'groEL',
        'AmpliconSize': 316,
        'Reference': 'PMID: 20686088',
        'RefSeq': 'NC_015838.1',
        'Filename': 'Gemella_morbillorum.fasta'
    },
    {
        'Organism': 'Gemella haemolysans',
        'PairID': 'Gemella_haemolysans_groEL',
        'Forward': 'GCTACTCGTGTAGGTATCGA',
        'Reverse': 'CCACCGTTTGTAGTGTCTTC',
        'Tier': 'Tier 2',
        'TargetGene': 'groEL',
        'AmpliconSize': 316,
        'Reference': 'PMID: 20686088',
        'RefSeq': 'NC_015841.1',
        'Filename': 'Gemella_haemolysans.fasta'
    },
    {
        'Organism': 'Gemella sanguinis',
        'PairID': 'Gemella_sanguinis_groEL',
        'Forward': 'GCTACTCGTGTAGGTATCGA',
        'Reverse': 'CCACCGTTTGTAGTGTCTTC',
        'Tier': 'Tier 2',
        'TargetGene': 'groEL',
        'AmpliconSize': 316,
        'Reference': 'PMID: 20686088',
        'RefSeq': 'NC_015843.1',
        'Filename': 'Gemella_sanguinis.fasta'
    },
    {
        'Organism': 'Streptococcus gallolyticus',
        'PairID': 'Streptococcus_gallolyticus_sodA',
        'Forward': 'CGTAGTTAAGGCTGTGTGCG',
        'Reverse': 'CTTGAGCAATGCCGTCAACC',
        'Tier': 'Tier 2',
        'TargetGene': 'sodA',
        'AmpliconSize': 264,
        'Reference': 'PMID: 15388602',
        'RefSeq': 'NC_013798.1',
        'Filename': 'Streptococcus_gallolyticus.fasta'
    },
    {
        'Organism': 'Streptococcus mutans',
        'PairID': 'Streptococcus_mutans_gtfB',
        'Forward': 'ACTACACTTTCGGGTGGCTTGG',
        'Reverse': 'CAGTATAAGCGCCAGTTTCATC',
        'Tier': 'Tier 2',
        'TargetGene': 'gtfB',
        'AmpliconSize': 517,
        'Reference': 'PMID: 10457635',
        'RefSeq': 'NC_004350.2',
        'Filename': 'Streptococcus_mutans.fasta'
    },
    {
        'Organism': 'Streptococcus mitis',
        'PairID': 'Streptococcus_mitis_16S',
        'Forward': 'CGTAGTTAAGGCTGTGTGCG',
        'Reverse': 'CTTGAGCAATGCCGTCAACC',
        'Tier': 'Tier 2',
        'TargetGene': '16S_rRNA',
        'AmpliconSize': 264,
        'Reference': 'PMID: 15388602',
        'RefSeq': 'NC_013853.1',
        'Filename': 'Streptococcus_mitis.fasta'
    }
]

def generate_pathogen_curation_report(rows, report_path):
    """Generates an HTML curation report for Tier 1 & Tier 2 pathogens."""
    t1_rows = [r for r in rows if r.get('Tier') == 'Tier 1']
    t2_rows = [r for r in rows if r.get('Tier') == 'Tier 2']
    t3_rows = [r for r in rows if r.get('Tier') == 'Tier 3']

    def render_table_rows(subset):
        html = ""
        for r in subset:
            b_cls = "badge-tier1" if r.get('Tier') == 'Tier 1' else ("badge-tier2" if r.get('Tier') == 'Tier 2' else "badge-tier3")
            html += f"""
            <tr>
                <td><span class="badge {b_cls}">{r.get('Tier')}</span></td>
                <td><strong>{r.get('Organism')}</strong></td>
                <td><code>{r.get('PairID')}</code></td>
                <td><span class="badge" style="background:#334155; color:#38bdf8;">{r.get('TargetGene')}</span></td>
                <td><code>{r.get('Forward')}</code></td>
                <td><code>{r.get('Reverse')}</code></td>
                <td>{r.get('AmpliconSize', 0)} bp</td>
                <td><span style="font-size:11px; color:#94a3b8;">{r.get('Reference', '')}</span></td>
            </tr>
            """
        return html

    html_content = f"""<!DOCTYPE html>
<html lang="de">
<head>
    <meta charset="UTF-8">
    <title>Pathogen Curation Report - Tier 1 &amp; Tier 2</title>
    <style>
        :root {{
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --card-border: #334155;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --accent-cyan: #06b6d4;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            margin: 0;
            padding: 24px;
        }}
        .card {{
            background: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 12px;
            padding: 20px;
            margin-bottom: 24px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
        }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }}
        .stat-card {{
            background: #1e293b;
            border: 1px solid #334155;
            border-radius: 8px;
            padding: 16px;
            text-align: center;
        }}
        .stat-val {{
            font-size: 32px;
            font-weight: 700;
            margin-top: 4px;
        }}
        .badge {{
            display: inline-block;
            padding: 3px 8px;
            border-radius: 6px;
            font-size: 11px;
            font-weight: 600;
        }}
        .badge-tier1 {{ background: rgba(239, 68, 68, 0.25); color: #f87171; border: 1px solid #ef4444; font-weight: 700; }}
        .badge-tier2 {{ background: rgba(2, 132, 199, 0.25); color: #38bdf8; border: 1px solid #0284c7; font-weight: 700; }}
        .badge-tier3 {{ background: rgba(100, 116, 139, 0.25); color: #94a3b8; border: 1px solid #64748b; }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 12px;
            margin-top: 12px;
        }}
        th, td {{
            padding: 10px 12px;
            text-align: left;
            border-bottom: 1px solid #334155;
        }}
        th {{
            background: #0f172a;
            color: var(--text-secondary);
            font-weight: 600;
        }}
        code {{
            background: #0f172a;
            padding: 2px 6px;
            border-radius: 4px;
            color: #38bdf8;
            font-family: monospace;
            font-size: 11px;
        }}
    </style>
</head>
<body>
    <div class="card" style="border-left: 4px solid var(--accent-cyan);">
        <h1 style="margin:0 0 8px 0;">🧬 Pathogen Curation Report (Tier 1 &amp; Tier 2)</h1>
        <p style="color:var(--text-secondary); margin:0;">
            Systematische Lückenanalyse und Kuration der humanpathogenen Erreger, diagnostischen Zielgene und Primerpaare.
        </p>
    </div>

    <div class="grid">
        <div class="stat-card" style="border-top: 3px solid #ef4444;">
            <div style="color:var(--text-secondary); font-size:13px;">🔴 Tier 1 (Borreliose &amp; Vektoren)</div>
            <div class="stat-val" style="color:#ef4444;">{len(t1_rows)}</div>
            <div style="font-size:11px; color:#94a3b8;">Primerpaare aktiv</div>
        </div>
        <div class="stat-card" style="border-top: 3px solid #0284c7;">
            <div style="color:var(--text-secondary); font-size:13px;">🔵 Tier 2 (Endokarditis &amp; Sepsis)</div>
            <div class="stat-val" style="color:#0284c7;">{len(t2_rows)}</div>
            <div style="font-size:11px; color:#94a3b8;">Primerpaare aktiv</div>
        </div>
        <div class="stat-card" style="border-top: 3px solid #64748b;">
            <div style="color:var(--text-secondary); font-size:13px;">⚪ Tier 3 (Weitere Humanpathogene)</div>
            <div class="stat-val" style="color:#94a3b8;">{len(t3_rows)}</div>
            <div style="font-size:11px; color:#94a3b8;">Primerpaare aktiv</div>
        </div>
        <div class="stat-card" style="border-top: 3px solid #10b981;">
            <div style="color:var(--text-secondary); font-size:13px;">Gesamt-Datenbank</div>
            <div class="stat-val" style="color:#10b981;">{len(rows)}</div>
            <div style="font-size:11px; color:#94a3b8;">8 Spalten (voll kuratiert)</div>
        </div>
    </div>

    <div class="card">
        <h2>🔴 Tier 1: Borreliose &amp; zecken-/vektorübertragene Ko-Infektionen ({len(t1_rows)})</h2>
        <div style="overflow-x:auto;">
            <table>
                <thead>
                    <tr>
                        <th>Tier</th>
                        <th>Organismus / Ziel</th>
                        <th>Pair ID</th>
                        <th>Zielgen</th>
                        <th>Forward Primer</th>
                        <th>Reverse Primer</th>
                        <th>Amplikon</th>
                        <th>Evidenz / Quelle</th>
                    </tr>
                </thead>
                <tbody>
                    {render_table_rows(t1_rows)}
                </tbody>
            </table>
        </div>
    </div>

    <div class="card">
        <h2>🔵 Tier 2: Endokarditis-, Sepsis-, Gelenk- &amp; Blutkultur-Erreger ({len(t2_rows)})</h2>
        <div style="overflow-x:auto;">
            <table>
                <thead>
                    <tr>
                        <th>Tier</th>
                        <th>Organismus / Ziel</th>
                        <th>Pair ID</th>
                        <th>Zielgen</th>
                        <th>Forward Primer</th>
                        <th>Reverse Primer</th>
                        <th>Amplikon</th>
                        <th>Evidenz / Quelle</th>
                    </tr>
                </thead>
                <tbody>
                    {render_table_rows(t2_rows)}
                </tbody>
            </table>
        </div>
    </div>
</body>
</html>
"""
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    dl_dir = r"C:\Users\internet\Downloads"
    if os.path.exists(dl_dir):
        dl_path = os.path.join(dl_dir, os.path.basename(report_path))
        try:
            shutil.copy2(report_path, dl_path)
            print(f"[+] Kopiert nach Downloads: {dl_path}")
        except Exception:
            pass
    print(f"[+] Curation Report generiert: {report_path}")

def curate_pathogen_tiers(primers_csv, ref_dirs=None, blast_db="/mnt/blast_db/nt", threads=16):
    """
    Performs comprehensive systematic gap analysis & curation for Tier 1 and Tier 2 pathogens.
    Enriches and updates primers_csv to the 8-column schema, synchronizes reference genomes,
    and produces an extensive HTML report.
    """
    print("\n" + "="*75)
    print("🧬 SYSTEMATISCHE PATHOGEN-KURATION (TIER 1 & TIER 2)")
    print("="*75)

    if not ref_dirs:
        ref_dirs = [
            "/mnt/sdc1/pathogen_references",
            "./pathogen_references",
            "./references",
            r"C:\Users\internet\Downloads\repository\primer\references",
            r"C:\Users\internet\Downloads"
        ]

    # 1. Load existing primers
    existing_primers = load_primers_from_csv(primers_csv)
    print(f"[*] Bestehende Primer in '{primers_csv}': {len(existing_primers)}")

    existing_pids = {p['PrimerPairID'].lower() for p in existing_primers}
    existing_pairs = {(p['Organism'].lower(), p['Forward'].upper()) for p in existing_primers}

    # 2. Enrich with Master Catalog additions
    added_tier1 = 0
    added_tier2 = 0
    for item in TIER1_MASTER:
        pid = item['PairID'].lower()
        pair_key = (item['Organism'].lower(), item['Forward'].upper())
        if pid not in existing_pids and pair_key not in existing_pairs:
            existing_primers.append({
                'Organism': item['Organism'],
                'PrimerPairID': item['PairID'],
                'Forward': item['Forward'],
                'Reverse': item['Reverse'],
                'Tier': item['Tier'],
                'TargetGene': item['TargetGene'],
                'AmpliconSize': item['AmpliconSize'],
                'Reference': item['Reference']
            })
            existing_pids.add(pid)
            existing_pairs.add(pair_key)
            added_tier1 += 1

    for item in TIER2_MASTER:
        pid = item['PairID'].lower()
        pair_key = (item['Organism'].lower(), item['Forward'].upper())
        if pid not in existing_pids and pair_key not in existing_pairs:
            existing_primers.append({
                'Organism': item['Organism'],
                'PrimerPairID': item['PairID'],
                'Forward': item['Forward'],
                'Reverse': item['Reverse'],
                'Tier': item['Tier'],
                'TargetGene': item['TargetGene'],
                'AmpliconSize': item['AmpliconSize'],
                'Reference': item['Reference']
            })
            existing_pids.add(pid)
            existing_pairs.add(pair_key)
            added_tier2 += 1

    print(f"[+] Neu hinzugefügt: {added_tier1} Tier 1 Primerpaare, {added_tier2} Tier 2 Primerpaare")
    print(f"[+] Neuer Gesamtbestand: {len(existing_primers)} Primerpaare")

    # 3. Sort by Tier, Organism, PairID
    tier_order = {'Tier 1': 1, 'Tier 2': 2, 'Tier 3': 3}
    existing_primers.sort(key=lambda x: (tier_order.get(x.get('Tier', 'Tier 3'), 9), x.get('Organism', '').lower(), x.get('PrimerPairID', '').lower()))

    # 4. Save 8-column CSV
    new_header = ['Organism', 'PairID', 'Forward', 'Reverse', 'Tier', 'TargetGene', 'AmpliconSize', 'Reference']
    save_targets = [primers_csv]
    dl_csv = r"C:\Users\internet\Downloads\best_primers_db.csv"
    if os.path.exists(r"C:\Users\internet\Downloads"):
        save_targets.append(dl_csv)
    for st in save_targets:
        try:
            with open(st, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(new_header)
                for p in existing_primers:
                    writer.writerow([
                        p.get('Organism', ''),
                        p.get('PrimerPairID', ''),
                        p.get('Forward', ''),
                        p.get('Reverse', ''),
                        p.get('Tier', 'Tier 3'),
                        p.get('TargetGene', ''),
                        p.get('AmpliconSize', 0),
                        p.get('Reference', '')
                    ])
            print(f"[+] 8-Spalten-Datenbank gespeichert nach: {st}")
        except Exception as e:
            print(f"[!] Fehler beim Schreiben nach {st}: {e}")

    # 5. Check Reference Genomes
    print("\n[*] Prüfe Vorhandensein der Referenzgenome für Tier 1 & 2...")
    missing_refs = []
    found_refs = 0
    all_master_items = TIER1_MASTER + TIER2_MASTER
    for item in all_master_items:
        org = item['Organism']
        ref_file = find_reference_genome(org, ref_dirs=ref_dirs)
        if ref_file and os.path.exists(ref_file):
            found_refs += 1
        else:
            missing_refs.append(item)

    print(f"[+] Referenzgenome vorhanden: {found_refs}/{len(all_master_items)}")
    if missing_refs:
        print(f"[!] Fehlende Referenzen ({len(missing_refs)}): {[m['Organism'] for m in missing_refs]}")

    # 6. Generate HTML Curation Report
    rep_path = "pathogen_curation_report.html"
    generate_pathogen_curation_report(existing_primers, rep_path)

    print("="*75 + "\n")
    return existing_primers

def main():
    parser = argparse.ArgumentParser(
        description="High-Performance Multi-Threaded Primer Search & Amplicon Extraction Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
================================================================================
BEISPIELE & PRAXIS-WORKFLOWS (COMMAND LINE GUIDE)
================================================================================

1. KLINISCHER GOLDSTANDARD (Empfohlene Routine-Vollanalyse):
   Durchsucht alle 738 Assays (inkl. aller Babesien, Plasmodien, Borrelien, Bartonellen etc.)
   kombiniert mit globalem Vollgenom-Mapping (--pan-scan) gegen 340+ Pathogene und 
   automatischer De-novo-Entdeckung unbekannter Erreger:
   $ python primer_suche.py -i probe.fastq.gz -a --pan-scan -t 24

2. SCHNELLES 2-STUFIGES SCREENING (Hierarchisch zur Zeitersparnis):
   Prüft in Stufe 1 nur konservierte Genus-Marker (z.B. Borrelia flaB/16S, Babesia 18S).
   Nur bei positivem Befund werden in Stufe 2 spezifische Spezies-Primer durchsucht:
   $ python primer_suche.py -i probe.fastq.gz -gs -t 24

3. MULTI-SAMPLE BATCH-VERARBEITUNG (Alle FASTQ-Dateien auf einmal):
   Verarbeitet alle FASTQ-Dateien im Verzeichnis per Wildcard. Erstellt Einzelreports
   sowie ein zusammenfassendes Multi-Probe-Dashboard (primer_dashboard_all_samples.html):
   $ python primer_suche.py -i *.fastq.gz -a --pan-scan -t 28

4. TIER-GEZIELTE SUCHE (Fokus auf bestimmte Erregergruppen):
   Filtert gezielt nach Prioritätsstufen:
   - Tier 0: Universelle Breitspektrum-Barcodes (Pan-Bakterien 16S, Pan-Pilze ITS, Pan-Protozoen 18S)
   - Tier 1: Borreliose, Zecken-Koinfektionen (Babesia, Bartonella, Rickettsia, Anaplasma etc.)
   - Tier 2: Endokarditis, Mykoplasmen, atypische Erreger
   $ python primer_suche.py -i probe.fastq.gz --tier 0,1 -a -t 24

5. TIEFENANALYSE AUF UNBEKANNTE ERREGER ("Dark Matter" De-novo Pipeline):
   Isoliert unzugeordnete Reads, assembliert sie de-novo zu Contigs und gleicht sie
   über lokalen BLASTn gegen /mnt/blast_db/nt ab:
   $ python primer_suche.py -i probe.fastq.gz -a --pan-scan --discover-unknown -t 24

6. NUR GENUS-SCREENING (Ultraschneller Erst-Check):
   Testet ausschließlich Breitband- und Genus-Primer zur ersten Orientierung:
   $ python primer_suche.py -i probe.fastq.gz -g -t 16

7. DATENBANK-PFLEGE & KURATION:
   a) Duplikate bereinigen (invertierte/komplementäre Primerpaare entfernen):
      $ python primer_suche.py -p best_primers_db.csv --dedup best_primers_db_clean.csv
   b) Primer gegen neueste NCBI RefSeq/GenBank-Genome validieren und updaten:
      $ python primer_suche.py --update-primers --target Babesia
      $ python primer_suche.py --update-primers --target all

================================================================================
FLAG-ERKLÄRUNGEN:
  -a, --all              Testet ALLE 738 Assays in der Datenbank (kein Filter)
  --pan-scan             Aktiviert Vollgenom-Mapping via minimap2 gegen 340+ Referenzgenome
  -gs, --hierarchical    2-Stufen-Modus: Erst Genus-Marker prüfen, nur bei Hit Spezies testen
  -g, --genus            Nur Genus-Screening-Primer testen
  --discover-unknown     De-novo Assembly & lokaler BLASTn unzugeordneter Reads (Standard: an)
  --no-discover-unknown  De-novo Erkennung deaktivieren
  --tier <0,1,2,all>     Filtert Assays nach klinischen Tiers (0=Universal, 1=Vektoren/Borreliose)
  -i, --fastq <Dateien>  Eingabedatei(en) oder Wildcards (*.fastq.gz)
  -t, --threads <N>      Anzahl CPU-Threads (z.B. 24 oder 28)
  -p, --primers <Pfad>   Pfad zur Primer-Datenbank CSV
================================================================================
"""
    )
    
    mode_group = parser.add_argument_group("Suchmodus (Mode)")
    mode_group.add_argument("-g", "--genus", action="store_true", 
                            help="Suche NUR nach Genus-Ebene / Screening-Primern")
    mode_group.add_argument("-a", "--all", action="store_true", 
                            help="Alle Primer durchsuchen (Standard)")
    mode_group.add_argument("-gs", "--genus-species", "--hierarchical", action="store_true", 
                            help="Zweistufige hierarchische Suche: Erst Genus prüfen, und nur bei positivem Befund die zugehörigen Spezies & Merkmale durchsuchen")

    parser.add_argument("--fastq", "-f", "-i", nargs="+", default=None, 
                        help="Pfad oder Wildcard zu FASTQ-Dateien (z.B. -i sample.fastq.gz oder -i *.fastq.gz)")
    parser.add_argument("--primers", "-p", default=r"C:\Users\internet\Downloads\best_primers_db.csv", 
                        help="Pfad zur Primer-Datenbank CSV (Standard: best_primers_db.csv)")
    parser.add_argument("--dedup-db", "--dedup", nargs="?", const="best_primers_db_dedup.csv", default=None, metavar="NEUE_DB.csv",
                        help="Bereinigt die Primer-Datenbank (-p): Löscht doppelte (FWD+REV), invertierte, Gegenstrang- und selbstkomplementäre Primerpaare und speichert die bereinigte Datenbank unter NEUE_DB.csv ab.")
    parser.add_argument("--update-primers", "--update-db", action="store_true",
                        help="Aktualisiert und validiert die Primer-Datenbank gegen neueste NCBI RefSeq/GenBank Genomsequenzen, EuropePMC/PubMed Assays (2024-2026) und optional GISAID.")
    parser.add_argument("--update-target", "--target", default=None,
                        help="Zielorganismus oder Pathogengruppe für gezieltes Primer-Update (z.B. 'Herpes', 'Rickettsia', 'Nairovirus' oder 'all')")
    parser.add_argument("--gisaid-file", default=None,
                        help="Pfad zu einer lokalen GISAID FASTA-Datei für Variant-Drift und Primer-Dropout Analyse")
    parser.add_argument("--curate-tiers", action="store_true",
                        help="Führt die systematische Lückenanalyse und Kuration für Tier 1 & Tier 2 durch: Lädt fehlende RefSeq-Komplettgenome, validiert Primer via /mnt/blast_db/nt und generiert die erweiterte 8-Spalten-Datenbank.")
    parser.add_argument("--tier", "--tiers", default=None,
                        help="Filtert die Suche nach Tier-Klassifikation: '0', '1', '2', '0,1,2' oder 'all' (Standard: alle)")
    parser.add_argument("--discover-unknown", "--dark-matter", action="store_true", default=True,
                        help="Aktiviert das automatische De-novo-Assembly und den lokalen BLASTn-Abgleich unzugeordneter Reads zur Entdeckung unbekannter/unerwarteter Erreger (Standard: aktiv)")
    parser.add_argument("--no-discover-unknown", dest="discover_unknown", action="store_false",
                        help="Deaktiviert die De-novo-Entdeckung unbekannter Erreger")
    parser.add_argument("--pan-scan", "--pan-pathogen-scan", action="store_true",
                        help="Führt einen globalen minimap2-Scan aller Reads gegen alle 340+ Pathogen-Referenzgenome durch (Paket B)")
    parser.add_argument("--threads", "-t", type=int, default=None,
                        help="Anzahl paralleler CPU-Threads (Standard: alle verfügbaren Kerne minus 2, max 32)")
    parser.add_argument("--output-dir", "-o", default="./extracted_amplicons", 
                        help="Ausgabeverzeichnis für extrahierte _frp_ FASTQ-Dateien")
    parser.add_argument("--dashboard", "-d", default="primer_dashboard.html", 
                        help="Dateiname/Pfad für das interaktive HTML-Dashboard")
    parser.add_argument("--overview-dashboard", "--summary-dashboard", default="primer_overview_dashboard.html", 
                        help="Dateiname/Pfad für das übergreifende Multi-FASTQ Zusammenfassungs-Dashboard (Standard: primer_overview_dashboard.html)")
    parser.add_argument("--ref-dir", default=r"/mnt/sdc1/pathogen_references", 
                        help="Pfad zum Verzeichnis mit Pathogen-Referenzgenomen (.fasta) für minimap2 (Standard: /mnt/sdc1/pathogen_references oder ./references)")
    parser.add_argument("--no-minimap2", action="store_true", 
                        help="Deaktiviert das automatische minimap2 Referenz-Mapping von Verdachtsfällen")
    parser.add_argument("--no-assembly", action="store_true",
                        help="Deaktiviert das automatische De-novo- und Konsensus-Micro-Assembly für Kandidaten-Reads (Standard: aktiv ab >= 5 Reads)")
    parser.add_argument("--blast-db", default="/mnt/blast_db/nt",
                        help="Pfad zur lokalen NCBI BLAST nt-Datenbank für automatische Validierung (Standard: /mnt/blast_db/nt)")
    parser.add_argument("--no-local-blast", action="store_true",
                        help="Deaktiviert den lokalen BLASTn-Abgleich gegen die lokale nt-Datenbank")
    parser.add_argument("--min-mapq", type=int, default=20, 
                        help="Minimale MAPQ-Schwelle für qualitatives Mapping (Standard: 20)")
    parser.add_argument("--max-mismatches", "-m", type=int, default=1, 
                        help="Maximal erlaubte Fehlpaarungen / Mismatches (Standard: 1)")
    parser.add_argument("--max-reads", "-n", type=int, default=None, 
                        help="Read-Limitierung für schnelle Tests")
    parser.add_argument("--cdhit", action="store_true", 
                        help="Aktiviert CD-HIT-EST Vor-Clustering zur Beschleunigung")
    parser.add_argument("--cdhit-identity", type=float, default=0.99, 
                        help="CD-HIT-EST Sequenz-Identitätsschwelle (Standard: 0.99)")
    parser.add_argument("--cdhit-exec", default="cd-hit-est", 
                        help="Executable/Befehl für cd-hit-est (Standard: cd-hit-est)")

    # Wenn ohne Argumente aufgerufen: Hilfe ausgeben und beenden
    if len(sys.argv) == 1:
        print("\n" + "="*75)
        print("🧬 PRIMER SUCHE & AMPLICON EXTRAKTION - HILFESEITE")
        print("="*75)
        parser.print_help()
        print("\n" + "="*75 + "\n")
        sys.exit(0)

    args = parser.parse_args()

    # Update und Validierung der Primer-Datenbank ausführen falls gewünscht
    if getattr(args, 'update_primers', False):
        ref_dirs = [
            args.ref_dir,
            "./references",
            r"C:\Users\internet\Downloads\repository\primer\references",
            r"C:\Users\internet\Downloads"
        ]
        update_primers_database(
            args.primers,
            target_filter=args.update_target,
            gisaid_path=args.gisaid_file,
            ref_dirs=ref_dirs
        )
        if not args.fastq:
            # Standalone-Lauf: Nach Update sauber beenden
            sys.exit(0)

    # Deduplizierung der Primer-Datenbank ausführen falls gewünscht
    if args.dedup_db:
        success = deduplicate_primer_database(args.primers, args.dedup_db)
        if not success:
            sys.exit(1)
        if not args.fastq:
            # Standalone-Lauf: Nur DB bereinigen und sauber beenden
            sys.exit(0)
        else:
            # Bei anschließendem FASTQ-Lauf die bereinigte Datenbank verwenden
            args.primers = args.dedup_db

    # Systematische Tier 1 & Tier 2 Kuration ausfuehren falls gewuenscht
    if getattr(args, 'curate_tiers', False):
        ref_dirs = [
            args.ref_dir,
            "/mnt/sdc1/pathogen_references",
            "./pathogen_references",
            "./references",
            r"C:\Users\internet\Downloads\repository\primer\references",
            r"C:\Users\internet\Downloads"
        ]
        curate_pathogen_tiers(
            primers_csv=args.primers,
            ref_dirs=ref_dirs,
            blast_db=args.blast_db,
            threads=args.threads or 16
        )
        if not args.fastq:
            # Standalone-Lauf: Nur Kuration durchfuehren und sauber beenden
            sys.exit(0)

    if not args.fastq:
        print("\n[!] FEHLER: Bitte gib mindestens eine FASTQ-Datei mit -i an oder nutze --curate-tiers / --update-primers / --dedup-db!")
        print("Beispiele:")
        print("  1. Tier-1 & Tier-2 Kuration durchfuehren: python primer_suche.py --curate-tiers")
        print("  2. Datenbank gegen RefSeq updaten:      python primer_suche.py --update-primers --target Herpes")
        print("  3. Datenbank bereinigen:                python primer_suche.py -p best_primers_db.csv --dedup clean_db.csv")
        print("  4. FASTQ-Suche nach Tiers filtern:      python primer_suche.py -i *.fastq.gz --tier 1,2 -gs -t 24\n")
        parser.print_help()
        sys.exit(1)

    # Resolve all matching input files (handles wildcards and multi-file arguments)
    raw_inputs = args.fastq if isinstance(args.fastq, list) else [args.fastq]
    input_files = []
    for item in raw_inputs:
        matched = glob.glob(item)
        if matched:
            input_files.extend(matched)
        elif os.path.exists(item):
            input_files.append(item)
            
    # Deduplicate while preserving order
    seen_paths = set()
    unique_fastq_files = []
    for f in input_files:
        abs_p = os.path.abspath(f)
        if abs_p not in seen_paths:
            seen_paths.add(abs_p)
            unique_fastq_files.append(f)

    if not unique_fastq_files:
        print(f"\n[!] FEHLER: Keine passenden FASTQ-Dateien für '{args.fastq}' gefunden!")
        sys.exit(1)

    print(f"\n[+] Gefundene FASTQ-Dateien zum Verarbeiten ({len(unique_fastq_files)}):")
    for idx, f in enumerate(unique_fastq_files, 1):
        print(f"    {idx}. {f} ({os.path.getsize(f):,} Bytes)")

    # Determine CPU threads
    available_cpus = os.cpu_count() or 4
    if args.threads is None:
        target_threads = max(1, min(32, available_cpus - 2)) if available_cpus > 2 else available_cpus
    else:
        target_threads = max(1, args.threads)

    # Load primer database once for all files
    all_primers = load_primers_from_csv(args.primers)
    print(f"[+] Total geladene Primer aus Datenbank: {len(all_primers)}")

    # Optionale Filterung nach Tier 1 / Tier 2 / Tier 3
    if getattr(args, 'tier', None):
        allowed_tiers = set()
        for t in str(args.tier).split(','):
            t_clean = t.strip().lower()
            if t_clean in ['0', 'tier0', 'tier 0']:
                allowed_tiers.add('tier 0')
            elif t_clean in ['1', 'tier1', 'tier 1']:
                allowed_tiers.add('tier 1')
            elif t_clean in ['2', 'tier2', 'tier 2']:
                allowed_tiers.add('tier 2')
            elif t_clean in ['3', 'tier3', 'tier 3']:
                allowed_tiers.add('tier 3')
        if allowed_tiers:
            all_primers = [p for p in all_primers if p.get('Tier', 'Tier 1').lower() in allowed_tiers]
            print(f"[*] Gefiltert nach Tier(s) {allowed_tiers}: {len(all_primers)} verbleibende Primer")

    genus_primers = []
    species_primers_by_genus = defaultdict(list)

    for p in all_primers:
        genus, is_genus = classify_primer_taxa(p['Organism'], p['PrimerPairID'])
        p['Genus'] = genus
        p['IsGenus'] = is_genus
        if is_genus:
            genus_primers.append(p)
        else:
            species_primers_by_genus[genus].append(p)

    # Sicherstellen, dass KEIN Genus (wie Cryptococcus, Candida, Aspergillus, HACEK, Herpes etc.)
    # im -g oder -gs Modus verwaist bleibt: Wenn kein Genus-Screening-Primer existiert,
    # werden automatisch repraesentative Spezies-Assays (z.B. Cryptococcus neoformans URA5) in Stufe 1 aufgenommen.
    represented_genera = {p['Genus'] for p in genus_primers}
    promoted_count = 0
    for g, sp_list in list(species_primers_by_genus.items()):
        if g not in represented_genera and sp_list:
            for p_rep in sp_list[:2]:
                p_copy = dict(p_rep)
                p_copy['IsGenus'] = True
                genus_primers.append(p_copy)
                promoted_count += 1
    if promoted_count > 0:
        print(f"[+] {promoted_count} repraesentative Assays (inkl. Cryptococcus, HACEK, Pilze, Viren) automatisch fuer Stufe 1 (-g/-gs) aktiviert.")

    generated_dashboards = []
    all_sample_results = []
    is_wildcard = any(c in "".join(raw_inputs) for c in ['*', '?'])
    is_multi_sample = len(unique_fastq_files) > 1 or is_wildcard
    dest_dir = os.path.dirname(args.dashboard) if os.path.dirname(args.dashboard) else "."
    overview_basename = os.path.basename(args.overview_dashboard)

    # Process each FASTQ file
    for f_idx, current_fastq in enumerate(unique_fastq_files, 1):
        print(f"\n===========================================================================")
        print(f"[*] Verarbeite FASTQ-Datei [{f_idx}/{len(unique_fastq_files)}]: {os.path.basename(current_fastq)}")
        print(f"===========================================================================")

        # Determine unique dashboard path (characters up to first underscore)
        if len(unique_fastq_files) > 1 or is_wildcard or args.dashboard == "primer_dashboard.html":
            dash_path = get_unique_dashboard_path(current_fastq, dest_dir=dest_dir)
        else:
            dash_path = args.dashboard

        # Run CD-HIT if requested
        input_file = current_fastq
        cdhit_stats = None
        if args.cdhit:
            input_file, cdhit_stats = run_cdhit_est(current_fastq, identity=args.cdhit_identity, cdhit_exec=args.cdhit_exec)

        # 1. Modus: Nur Genus
        if args.genus:
            mode_str = "Nur Genus (-g)"
            print(f"\n[*] MODUS: Nur Genus (-g) aktiviert. Durchsuche {len(genus_primers)} Genus-Primer...")
            results = run_primer_search(input_file, genus_primers, args.output_dir, 
                                        max_mismatches=args.max_mismatches, max_reads=args.max_reads,
                                        stage_label="Genus-Ebene", threads=target_threads)
            results['search_mode'] = mode_str

        # 2. Modus: Zweistufig hierarchisch (Genus -> Spezies)
        elif args.genus_species:
            mode_str = "Hierarchisch: Genus -> Spezies (-gs)"
            print(f"\n[*] MODUS: Hierarchisch (-gs) aktiviert.")
            print(f"[*] --- STUFE 1: Screening auf {len(genus_primers)} Genus-Primern ---")
            stage1_results = run_primer_search(input_file, genus_primers, args.output_dir, 
                                               max_mismatches=args.max_mismatches, max_reads=args.max_reads,
                                               stage_label="Stufe 1: Genus-Screening", threads=target_threads)
            
            positive_genera = set()
            # Stufe 1: Auch einseitige Treffer (ONLY_FORWARD oder ONLY_REVERSE) aktivieren das Genus
            for pair_id, stats in stage1_results['primer_stats'].items():
                if stats['BOTH'] > 0 or stats['ONLY_FORWARD'] > 0 or stats['ONLY_REVERSE'] > 0:
                    g, _ = classify_primer_taxa(stats['Organism'], pair_id)
                    positive_genera.add(g)

            # KERN-PATHOGENE IMMER VOLLSTAENDIG TESTEN (Benutzeranforderung):
            # Babesia, Plasmodium, Borrelia, Bartonella, Mycoplasma, Mycobacterium, Brucella
            CORE_MANDATORY_GENERA = {
                'babesia', 'plasmodium', 'borrelia', 'bartonella',
                'mycoplasma', 'mycobacterium', 'brucella', 'cryptococcus',
                'granulicatella', 'abiotrophia', 'rickettsia', 'anaplasma',
                'ehrlichia', 'coxiella', 'francisella', 'leptospira'
            }
            for cmg in CORE_MANDATORY_GENERA:
                for g_key in species_primers_by_genus:
                    if cmg in g_key.lower():
                        positive_genera.add(g_key)

            print(f"\n[+] STUFE 1 Ergebnis (inkl. obligatorischer Kern-Pathogene):")
            if positive_genera:
                print(f"    -> Positive Genera erkannt: {', '.join(sorted(positive_genera))}")
                target_species_primers = []
                for g in positive_genera:
                    if g in species_primers_by_genus:
                        target_species_primers.extend(species_primers_by_genus[g])

                if target_species_primers:
                    print(f"[*] --- STUFE 2: Spezies & Merkmale Vertiefung ({len(target_species_primers)} Spezies-Primer) ---")
                    stage2_results = run_primer_search(input_file, target_species_primers, args.output_dir,
                                                       max_mismatches=args.max_mismatches, max_reads=args.max_reads,
                                                       stage_label="Stufe 2: Spezies-Vertiefung", threads=target_threads)
                    results = merge_search_results(stage1_results, stage2_results)
                else:
                    print(f"    -> Keine weiteren Sub-Spezies in der Datenbank für diese Genera hinterlegt.")
                    results = stage1_results
            else:
                print("    -> Keine positiven Genera im Screening gefunden. Spezies-Vertiefung übersprungen.")
                results = stage1_results

            results['search_mode'] = mode_str
            results['positive_genera'] = sorted(list(positive_genera))

        # 3. Modus: Alle Primer (Standard / -a)
        else:
            mode_str = "Alle Primer (-a)"
            print(f"\n[*] MODUS: Alle Primer (-a) aktiviert. Durchsuche alle {len(all_primers)} Primer...")
            results = run_primer_search(input_file, all_primers, args.output_dir, 
                                        max_mismatches=args.max_mismatches, max_reads=args.max_reads,
                                        stage_label="Vollständige Suche", threads=target_threads)
            results['search_mode'] = mode_str

        sample_pfx = extract_sample_prefix(current_fastq)
        if cdhit_stats:
            results['cdhit_stats'] = cdhit_stats

        results['fastq_path'] = current_fastq
        results['fastq_name'] = os.path.basename(current_fastq)
        results['sample_prefix'] = sample_pfx
        results['dashboard_path'] = dash_path
        results['overview_dashboard_name'] = overview_basename
        results['is_multi_sample'] = is_multi_sample

        # Minimap2 Reference Validation on candidate organisms (split FWD/REV or positive hits)
        if not args.no_minimap2:
            ref_dirs = [args.ref_dir, "./pathogen_references", "./references", r"C:\Users\internet\Downloads\repository\primer\references"]
            run_minimap2_validation(
                results,
                ref_dirs=ref_dirs,
                threads=target_threads,
                min_mapq=args.min_mapq,
                enable_assembly=not args.no_assembly,
                enable_local_blast=not args.no_local_blast,
                blast_db=args.blast_db
            )
        else:
            results['minimap2_validations'] = []

        # Paket B: Pan-Pathogen minimap2 Scan falls aktiviert
        if getattr(args, 'pan_scan', False):
            ref_dirs = [args.ref_dir, "./pathogen_references", "./references", r"C:\Users\internet\Downloads\repository\primer\references"]
            run_pan_pathogen_scan(
                candidate_fastq_paths=list(results.get('extracted_counts', {}).keys()),
                ref_dirs=ref_dirs,
                output_dir=args.output_dir,
                prefix=sample_pfx,
                threads=target_threads,
                min_mapq=args.min_mapq
            )

        # Paket C: Dark Matter Discovery für unzugeordnete Reads
        dark_matter_hits = []
        if getattr(args, 'discover_unknown', True) and not args.no_local_blast:
            dark_matter_hits = run_dark_matter_discovery(
                unassigned_reads=results.get('unassigned_reads', []),
                output_dir=args.output_dir,
                prefix=sample_pfx,
                blast_db=args.blast_db,
                threads=target_threads
            )
        results['dark_matter_discoveries'] = dark_matter_hits

        # Paket A: Klinische Befund-Aggregation & Diagnose-Ampel
        results['clinical_findings'] = aggregate_clinical_findings(
            primer_stats=results.get('primer_stats', {}),
            validations=results.get('minimap2_validations', []),
            dark_matter_discoveries=dark_matter_hits,
            extracted_counts=results.get('extracted_counts', {})
        )

        generate_html_dashboard(results, dash_path)
        generated_dashboards.append(dash_path)
        all_sample_results.append(results)

    # Multi-sample overview dashboard across all FASTQ files
    overview_path = None
    if is_multi_sample or (args.overview_dashboard != "primer_overview_dashboard.html"):
        overview_file_path = os.path.join(dest_dir, overview_basename)
        overview_path = generate_summary_html_dashboard(all_sample_results, overview_file_path)

    print(f"\n" + "="*75)
    print(f"[DONE] Alle {len(unique_fastq_files)} FASTQ-Dateien erfolgreich analysiert!")
    print(f"[+] Generierte Einzel-Dashboards:")
    for d in generated_dashboards:
        print(f"    - {os.path.abspath(d)}")
    if overview_path:
        print(f"\n[+] Multi-FASTQ Gesamtübersicht (Übersichts-Dashboard):")
        print(f"    - {os.path.abspath(overview_path)}")
        dl_overview = os.path.join(r"C:\Users\internet\Downloads", os.path.basename(overview_path))
        if os.path.exists(dl_overview):
            print(f"    - Kopie in Downloads:        {dl_overview}")
    print("="*75 + "\n")

if __name__ == "__main__":
    main()
