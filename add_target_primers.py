import os
import csv
import json
import re
import zipfile
import io
import urllib.request
import sys

sys.stdout.reconfigure(encoding='utf-8')

# Target genera requested by user
TARGET_GENERA = ["babesia", "plasmodium", "mycobacterium", "mycoplasma", "bartonella", "borrelia"]

# 1. Standard diagnostic PCR primer database for species & subspecies of requested genera
DIAGNOSTIC_PRIMERS = [
    # --- Babesia species & subspecies ---
    {"Organism": "Babesia microti", "PairID": "Babesia_microti_18S_Bab1_Bab4", "Forward": "CTTAGTATAAGCTTTTATACAGC", "Reverse": "GATAGGTCAGAAACTTGAATGATACAC"},
    {"Organism": "Babesia microti", "PairID": "Babesia_microti_18S_nss-1_nss-2", "Forward": "GTTATCCTTAAGCCATGCATGTCT", "Reverse": "TGATCCTTCTGCAGGTTCACCTAC"},
    {"Organism": "Babesia divergens", "PairID": "Babesia_divergens_18S_Bdiv_F_R", "Forward": "ACTATGGTTAATAGGAGCAG", "Reverse": "CAATTCCTTTAAGTTTCAGC"},
    {"Organism": "Babesia bovis", "PairID": "Babesia_bovis_18S_Bbov1_Bbov2", "Forward": "TGTATGGCTAACCGGTAGTA", "Reverse": "ACCAGTACTTCTCACTACCT"},
    {"Organism": "Babesia duncani", "PairID": "Babesia_duncani_18S_Bdunc_F_R", "Forward": "CGTTAACTACTAGCGGATTT", "Reverse": "CCTCTAAGAAGCTCGTACTA"},
    {"Organism": "Babesia canis", "PairID": "Babesia_canis_18S_Bcan_F_R", "Forward": "TGCGTTGAATACTAGCCTTT", "Reverse": "CAATTCCTTTAAGTTTCAGC"},
    {"Organism": "Babesia caballi", "PairID": "Babesia_caballi_BC48_F_R", "Forward": "TCGAACCTTGATTTCTCTGG", "Reverse": "GCTTCCTTCTTAACGTTGGC"},

    # --- Plasmodium species & subspecies ---
    {"Organism": "Plasmodium falciparum", "PairID": "Plasmodium_falciparum_18S_rFAL1_rFAL2", "Forward": "TTAAACTGGTTTGGGAAAACCAAATATATT", "Reverse": "ACACAATGAACTCAATCATGACTACCCGTC"},
    {"Organism": "Plasmodium vivax", "PairID": "Plasmodium_vivax_18S_rVIV1_rVIV2", "Forward": "CGCTTCTAGCTTAATCCACATAACTGATAC", "Reverse": "ACTTCCAAGCCGAAGCAAAGAAAGTCCTTA"},
    {"Organism": "Plasmodium malariae", "PairID": "Plasmodium_malariae_18S_rMAL1_rMAL2", "Forward": "ATAACATAGTTGTACGTTAAGAATAACCGC", "Reverse": "AAAATTCCCATGCAAAAGTTCCCATTTAAA"},
    {"Organism": "Plasmodium ovale", "PairID": "Plasmodium_ovale_18S_rOVA1_rOVA2", "Forward": "ATCTCTTTTGCTATTTTTTAGTATTGGAGA", "Reverse": "GGAAAAGGACACATTAATTGTATCCTAGTG"},
    {"Organism": "Plasmodium knowlesi", "PairID": "Plasmodium_knowlesi_18S_Pkn_F_R", "Forward": "TCCTTCCTATTATTATTCCATGCT", "Reverse": "ATCTCTTTTCTCTCCCTTTTCC"},
    {"Organism": "Plasmodium berghei", "PairID": "Plasmodium_berghei_18S_Pberg_F_R", "Forward": "AACTGGTTTGGATGACCAAATATATT", "Reverse": "ACACAATGAATTCAATCATGACTACC"},

    # --- Mycobacterium species & subspecies ---
    {"Organism": "Mycobacterium tuberculosis", "PairID": "Mycobacterium_tuberculosis_IS6110_INS1_INS2", "Forward": "CGTGAGGGCATCGAGGTGGC", "Reverse": "GCGTAGGCGTCGGTGACAAA"},
    {"Organism": "Mycobacterium tuberculosis", "PairID": "Mycobacterium_tuberculosis_rpoB_TbF_TbR", "Forward": "GGAGCGGATGACCACCCAGGAC", "Reverse": "GCCGCGATCAAGGAGTTCTTCG"},
    {"Organism": "Mycobacterium leprae", "PairID": "Mycobacterium_leprae_RLEP_RLEP1_RLEP2", "Forward": "ATATCGATGCAGGCGTGAG", "Reverse": "GGTAGCCGAAATCACCGT"},
    {"Organism": "Mycobacterium avium", "PairID": "Mycobacterium_avium_IS1245_F_R", "Forward": "GCCGCCGAAACGATCTAC", "Reverse": "AGGTGGCGTCGAGGAAGAC"},
    {"Organism": "Mycobacterium avium subsp. paratuberculosis", "PairID": "Mycobacterium_avium_paratuberculosis_IS900_P90_P91", "Forward": "GAAGGGTGTTCGGGGCCGTCGCTTAGG", "Reverse": "GGCGTCGAGGTCGATCGCCCACGTGAC"},
    {"Organism": "Mycobacterium ulcerans", "PairID": "Mycobacterium_ulcerans_IS2404_F_R", "Forward": "AGATCAAGCGCATCAAGCAC", "Reverse": "ACCAGTAACGCGACATCG"},
    {"Organism": "Mycobacterium abscessus", "PairID": "Mycobacterium_abscessus_rpoB_F_R", "Forward": "TCAAGGAGAAGCGCTACGA", "Reverse": "GGATGTTGATCAGGGTCTGC"},

    # --- Mycoplasma species & subspecies ---
    {"Organism": "Mycoplasma pneumoniae", "PairID": "Mycoplasma_pneumoniae_P1_adhesin_MP5_MP7", "Forward": "GAAGCTTATGGTACAGGTCGG", "Reverse": "ATTACCATCCTTGTTGTAAGGG"},
    {"Organism": "Mycoplasma pneumoniae", "PairID": "Mycoplasma_pneumoniae_16S_MpnF_MpnR", "Forward": "TGCTTGACAGTACCATTAG", "Reverse": "CACTGTAGTATCTAATCCTG"},
    {"Organism": "Mycoplasma genitalium", "PairID": "Mycoplasma_genitalium_mgpB_MgF_MgR", "Forward": "AGTTGATGAAACCTTAACCCCTTGG", "Reverse": "CCGTTGAGGGTTTTCCATTTTTG"},
    {"Organism": "Mycoplasma hominis", "PairID": "Mycoplasma_hominis_16S_MhF_MhR", "Forward": "CAATGGCTAATGCCGGATAC", "Reverse": "GGTGTTCTTCCTTATATCTACGC"},
    {"Organism": "Mycoplasma hyopneumoniae", "PairID": "Mycoplasma_hyopneumoniae_16S_Mhyo_F_R", "Forward": "GAGCCTTCAAGCTTCACCAAG", "Reverse": "ACTAGATCAAGGTTCCTGAAG"},
    {"Organism": "Mycoplasma gallisepticum", "PairID": "Mycoplasma_gallisepticum_16S_Mgall_F_R", "Forward": "GAGCTAATCTGTAAAGTTGGTC", "Reverse": "GCTTCCTTGCGGTTAGTCAAC"},

    # --- Bartonella species & subspecies ---
    {"Organism": "Bartonella henselae", "PairID": "Bartonella_henselae_gltA_BhF_BhR", "Forward": "TTTCTACAGTTATGGCTACC", "Reverse": "GATCYTCAATCATTTCAAA"},
    {"Organism": "Bartonella quintana", "PairID": "Bartonella_quintana_gltA_BqF_BqR", "Forward": "CATGTCAATACTAATGGTTTTA", "Reverse": "TCAGATTTTTCAATTTCATC"},
    {"Organism": "Bartonella bacilliformis", "PairID": "Bartonella_bacilliformis_gltA_BbF_BbR", "Forward": "CCTTCATCAGGTTATGCTA", "Reverse": "GATCTTCAATCATTTCAAA"},
    {"Organism": "Bartonella alsatica", "PairID": "Bartonella_alsatica_ITS_BaF_BaR", "Forward": "CCTTCATCAGGTTATGCTA", "Reverse": "GATCTTCAATCATTTCAAA"},

    # --- Borrelia species & subspecies ---
    {"Organism": "Borrelia burgdorferi", "PairID": "Borrelia_burgdorferi_flaB_FlaF_FlaR", "Forward": "TGCTAATTCAGCTTGAGAA", "Reverse": "TCAGACTCAACTACCGAAC"},
    {"Organism": "Borrelia burgdorferi sensu stricto", "PairID": "Borrelia_burgdorferi_ss_ospA_F_R", "Forward": "AATATTTCAAACGATCTAAGTCA", "Reverse": "GTAGGTTCTACAGCTACATGA"},
    {"Organism": "Borrelia afzelii", "PairID": "Borrelia_afzelii_5S_23S_IGS_F_R", "Forward": "CTGCGAGTTCGCGGGAGA", "Reverse": "TCCTAGGCATTCACCATA"},
    {"Organism": "Borrelia garinii", "PairID": "Borrelia_garinii_5S_23S_IGS_F_R", "Forward": "CTGCGAGTTCGCGGGAGA", "Reverse": "TCCTAGGCATTCACCATA"},
    {"Organism": "Borrelia miyamotoi", "PairID": "Borrelia_miyamotoi_glpQ_BmF_BmR", "Forward": "GATGCTTGAAACTTTAGCAGG", "Reverse": "AGAAATACTAAGAGGGCTAAT"},
    {"Organism": "Borrelia recurrentis", "PairID": "Borrelia_recurrentis_flaB_BrF_BrR", "Forward": "ATGAAAATCTTAGCAATGAT", "Reverse": "TTAGTTGAACATTTACTTGA"},
    {"Organism": "Borrelia duttonii", "PairID": "Borrelia_duttonii_16S_BdF_BdR", "Forward": "CATGCAAGTCGAACGGAACG", "Reverse": "CCGTCATTCCCTTACGCTAC"}
]

def load_existing_db(csv_path):
    existing = []
    if os.path.exists(csv_path):
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                existing.append(row)
    return existing

def main():
    # 1. Fetch 1,280 qPrimerDB best primers
    import build_1315_direct_fast
    with open('organisms.json', 'r', encoding='utf-8') as f:
        org_data = json.load(f)

    org_list = []
    for cat, grps in org_data.items():
        if isinstance(grps, dict):
            for grp, items in grps.items():
                if isinstance(items, list):
                    for item in items:
                        if item.get('bestPrimers', 0) > 0:
                            clean_key = item['orgKey'].replace('___', '__')
                            url = f"https://qprimerdb.biodb.org/tmp/qprimerdb/downloads/{clean_key}.best.primers.xls.zip"
                            org_list.append({
                                'name': item['orgName'],
                                'key': item['orgKey'],
                                'url': url
                            })

    print(f"[*] Fetching 1,280 qPrimerDB Best Primers...")
    qprimer_rows = []
    import concurrent.futures
    with concurrent.futures.ThreadPoolExecutor(max_workers=30) as executor:
        future_to_org = {executor.submit(build_1315_direct_fast.fetch_organism_best_primer, org): org for org in org_list}
        for future in concurrent.futures.as_completed(future_to_org):
            res = future.result()
            if res:
                qprimer_rows.append(res)

    print(f"[+] Fetched {len(qprimer_rows)} qPrimerDB primer pairs.")

    existing_pair_ids = set(r['PairID'] for r in qprimer_rows)
    for p in DIAGNOSTIC_PRIMERS:
        if p['PairID'] not in existing_pair_ids:
            qprimer_rows.append(p)
            existing_pair_ids.add(p['PairID'])

    print(f"[+] Combined Total Primer Pairs in Database: {len(qprimer_rows)}")

    csv_path = r"C:\Users\internet\Downloads\best_primers_db.csv"
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
        writer.writeheader()
        writer.writerows(qprimer_rows)

    print(f"[+] Saved complete database ({len(qprimer_rows)} primers) to: {csv_path}")

if __name__ == "__main__":
    main()
