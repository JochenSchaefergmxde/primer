# Pathogen In-Silico Primer Search & Metagenomics Discovery Pipeline

High-Performance Multi-Threaded Primer Search, Amplicon Extraction, and Metagenomic Pathogen Discovery Pipeline.

Developed in collaboration with diagnostic experience values from **Selbsthilfe Gruppe Borreliose Hannover**.

---

## 🔬 Core Capabilities

1. **738 Curated Diagnostic Assays (`best_primers_db.csv`)**:
   - **Tier 0:** Universal broad-range barcoding (Pan-Bacteria 16S / rpoB, Pan-Fungi ITS / 28S, Pan-Protozoa 18S, Pan-Flavivirus, Pan-Herpes).
   - **Tier 1:** Borreliosis and vector-borne coinfections (*Borrelia burgdorferi sensu lato*, *Borrelia miyamotoi*, *Babesia microti/divergens/odocoilei/venatorum*, *Plasmodium*, *Bartonella*, *Rickettsia*, *Anaplasma*, *Ehrlichia*).
   - **Tier 2:** Endocarditis, atypical bacterial pathogens (*Mycoplasma*, *Mycobacterium*, *Brucella*, *Coxiella*, *Francisella*, *Tropheryma whipplei*, NVS).
   - **Tier 3:** Extended bacterial, viral, and fungal pathogen panels (*Cryptococcus neoformans*, *Candida*, *Aspergillus*, etc.).

2. **Single-Primer-Anchor & Dual 8-mer Extension (Paket B)**:
   - Recovers fragmented metagenomic reads with 1 mismatch tolerance.
   - Dual-seed alignment with reverse complement recognition.

3. **Pan-Pathogen Genome Scan (`--pan-scan`)**:
   - High-throughput `minimap2` alignment against 340+ reference complete genomes.

4. **Dark Matter De-Novo Discovery (`--discover-unknown`)**:
   - Clusters and micro-assembles unassigned reads into contigs.
   - Automated local validation via BLASTn (`/mnt/blast_db/nt`).

5. **Clinical Traffic Light Dashboard**:
   - 🟢 **Bestätigt (Confirmed):** Multi-locus amplicon detection or orthogonal BLAST-confirmed genus/species match.
   - 🟡 **Verdacht (Suspect):** Single-primer signal or partial mapping requiring review.
   - ⚪ **Gering (Low):** Isolated low-abundance reads.
   - 🚫 **Kitom / Kontaminant:** Automatic identification of environmental/kit flora (e.g. *Ralstonia*, *Cutibacterium*).

---

## 🚀 Quick Start

### 1. Clinical Gold Standard (Full Panel + Pan-Scan)
```bash
python primer_suche.py -i sample.fastq.gz -a --pan-scan -t 24
```

### 2. Fast 2-Stage Hierarchical Screening
```bash
python primer_suche.py -i sample.fastq.gz -gs -t 24
```

### 3. Multi-Sample Batch Processing
```bash
python primer_suche.py -i *.fastq.gz -a --pan-scan -t 28
```

### 4. Targeted Priority Tiers (e.g. Tiers 0 and 1)
```bash
python primer_suche.py -i sample.fastq.gz --tier 0,1 -a -t 24
```

---

## 📊 Interactive Dashboards

- **Single-Sample Dashboard (`primer_dashboard_<prefix>.html`):**
  - Tab 1: 🏥 Clinical Findings & Traffic Light Status
  - Tab 2: 🎯 BLAST-Confirmed Matches (Genus & Species)
  - Tab 3: 📊 Detailed Primer Matrix (738 Assays with 1-Click NCBI Web-BLAST)
  - Tab 4: 🔬 De-novo Discoveries (Dark Matter Contigs)
  - Tab 5: 📁 Extracted Amplicons
  - Tab 6: 📖 Clinical Guidelines & CLI Documentation (SHG Borreliose Hannover)

- **Multi-Sample Overview (`primer_overview_dashboard.html`):**
  - Cross-sample pathogen heatmaps, comparative bar charts, and multi-locus presence matrix.

- **Example Report:**
  - An interactive sample primer dashboard is available under [`examples/js_5600180_primer_dashboard_js.html`](examples/js_5600180_primer_dashboard_js.html).

---

## 🛠️ Requirements

- Python 3.10+
- `biopython`, `pandas`, `plotly`
- Optional: `minimap2`, `blastn` (with NCBI `nt` database)
