import os
import csv
import sys

sys.stdout.reconfigure(encoding='utf-8')

VIRAL_PRIMERS = [
    # --- Herpesviruses ---
    {"Organism": "Epstein-Barr virus (EBV / HHV-4)", "PairID": "EBV_BamHI_W_F_R", "Forward": "CCCAACACTCCACCACACC", "Reverse": "TCTTAGGAGCTGTCACCAAC"},
    {"Organism": "Epstein-Barr virus (EBV / HHV-4)", "PairID": "EBV_EBNA1_F_R", "Forward": "TACAGGACCTACCCAGCTCA", "Reverse": "AAGACACCATCGACCCAGAC"},
    {"Organism": "Cytomegalovirus (CMV / HHV-5)", "PairID": "CMV_UL83_pp65_F_R", "Forward": "GTCAGCGTTCGTGTTTCCCA", "Reverse": "GGGACACAACACCGTAAAGC"},
    {"Organism": "Cytomegalovirus (CMV / HHV-5)", "PairID": "CMV_UL54_pol_F_R", "Forward": "CCTAGTGTGAGATCCATACTC", "Reverse": "ATGTTTCTCTTCCTGAATCTG"},
    {"Organism": "Herpes simplex virus 1 (HSV-1)", "PairID": "HSV1_gB_F_R", "Forward": "GCAGTTTACGTACAACCACATAC", "Reverse": "AGCTTGCGGGCCTCGTT"},
    {"Organism": "Herpes simplex virus 2 (HSV-2)", "PairID": "HSV2_gD_F_R", "Forward": "GTACATCGACTGGACCGAGC", "Reverse": "GCTCCAGATACCGGAAGTGC"},
    {"Organism": "Varicella-zoster virus (VZV / HHV-3)", "PairID": "VZV_ORF29_F_R", "Forward": "TGTACTGTAGAAAAGCCACCG", "Reverse": "CCTGGTAGGAATACCCCACAC"},
    {"Organism": "Human herpesvirus 6 (HHV-6A/B)", "PairID": "HHV6_U57_F_R", "Forward": "GCTAGAACGTATTTGCTGTC", "Reverse": "ATCCTCGAGTTCATCAACACC"},
    {"Organism": "Human herpesvirus 7 (HHV-7)", "PairID": "HHV7_U10_F_R", "Forward": "CGGAAGTCACTGGAGTAGTG", "Reverse": "GCAGTATGCCATCGTAGGAC"},
    {"Organism": "Kaposi sarcoma-associated herpesvirus (KSHV / HHV-8)", "PairID": "HHV8_ORF26_F_R", "Forward": "AGCCGAAAGGATTCCACCAT", "Reverse": "TCCGTGTTGTCTACGTCCAG"},

    # --- Respiratory Viruses ---
    {"Organism": "SARS-CoV-2 (COVID-19)", "PairID": "SARS_CoV_2_CDC_N1", "Forward": "GACCCCAAAATCAGCGAAAT", "Reverse": "TCTGGTTACTGCCAGTTGAATCTG"},
    {"Organism": "SARS-CoV-2 (COVID-19)", "PairID": "SARS_CoV_2_CDC_N2", "Forward": "TTACAAACATTGGCCGCAAA", "Reverse": "GCGCGACATTCCGAAGAA"},
    {"Organism": "SARS-CoV-2 (COVID-19)", "PairID": "SARS_CoV_2_E_gene_Charite", "Forward": "ACAGGTACGTTAATAGTTAATAGCGT", "Reverse": "ATATTGCAGCAGTACGCACACA"},
    {"Organism": "Influenza A virus", "PairID": "FluA_M1_universal", "Forward": "AAGACCAATCCTGTCACCTCTGA", "Reverse": "CAAAGCGTCTACGCTGCAGTCC"},
    {"Organism": "Influenza B virus", "PairID": "FluB_NS1_universal", "Forward": "TCCTCAACTCACTCTTCGAGCG", "Reverse": "CGGTGCTCTTGACCAAATTGG"},
    {"Organism": "Respiratory syncytial virus (RSV A)", "PairID": "RSV_A_N_gene", "Forward": "AGATCAACTTCTGTCATCCAGCAA", "Reverse": "GCACATCATAATTAGGAGTATCAAT"},
    {"Organism": "Respiratory syncytial virus (RSV B)", "PairID": "RSV_B_N_gene", "Forward": "AAGATGCAAATCATAAATTCACAGGA", "Reverse": "TGATATCCAGCATCTTTAAGTATCTCA"},
    {"Organism": "Human adenovirus (HAdV)", "PairID": "Adenovirus_Hexon_univ", "Forward": "GCCACNGAGACNTACTTCAGC", "Reverse": "TTGTACAGNGCCACNACNGG"},
    {"Organism": "Human rhinovirus / Enterovirus", "PairID": "Enterovirus_5UTR_univ", "Forward": "TCCGGCCCCTGAATGCGGCTAAT", "Reverse": "ACACGGACACCCAAAGTAGTCGGT"},

    # --- Hepatitis & Bloodborne Viruses ---
    {"Organism": "Hepatitis B virus (HBV)", "PairID": "HBV_S_gene_F_R", "Forward": "ATCCTGCTGCTATGCCTCATC", "Reverse": "ATTGAGAAGAAGAACTCCCTCGC"},
    {"Organism": "Hepatitis C virus (HCV)", "PairID": "HCV_5UTR_F_R", "Forward": "TCTGCGGAACCGGTGAGTAC", "Reverse": "GCAAGCACCCTATCAGGCAGT"},
    {"Organism": "Human immunodeficiency virus 1 (HIV-1)", "PairID": "HIV1_gag_F_R", "Forward": "ACATCAAGCAGCCATGCAAAT", "Reverse": "ATCTGGCCTGGTGCAATAGG"},
    {"Organism": "Human T-lymphotropic virus 1 (HTLV-1)", "PairID": "HTLV1_tax_F_R", "Forward": "CGGATACCCAGTCTACGTGT", "Reverse": "GAGCCGATAACGCGTCCATC"},

    # --- Polyomaviruses & Parvovirus ---
    {"Organism": "BK polyomavirus (BKV)", "PairID": "BKV_VP1_F_R", "Forward": "CTAGTGCTGCCAACAGCCTA", "Reverse": "TGGATATCTTAGAGTGCTCAA"},
    {"Organism": "JC polyomavirus (JCV)", "PairID": "JCV_VP1_F_R", "Forward": "CTGAGTGTTTAGAAGTTG", "Reverse": "GTTCTCACTATCAATTCC"},
    {"Organism": "Primate erythroparvovirus 1 (Parvovirus B19)", "PairID": "ParvoB19_NS1_F_R", "Forward": "GGGGCAGCATCAGTAGACAA", "Reverse": "TGCTTCTACTGGGCTTTTACC"},

    # --- Gastrointestinal & Papillomaviruses ---
    {"Organism": "Norovirus GI", "PairID": "Norovirus_GI_ORF1_ORF2", "Forward": "CGYTGGATGCGNTTYCATGA", "Reverse": "CTTAGACGCCATCATCATTYAC"},
    {"Organism": "Norovirus GII", "PairID": "Norovirus_GII_ORF1_ORF2", "Forward": "CARGARBCNATGTTYAGRTGGATGAG", "Reverse": "TCGACGCCATCTTCATTCACA"},
    {"Organism": "Human papillomavirus 16 (HPV-16)", "PairID": "HPV16_E6_F_R", "Forward": "CAGCAATACAACAAACCGTTGG", "Reverse": "CGCAACAAGACATACATCGAC"},
    {"Organism": "Human papillomavirus 18 (HPV-18)", "PairID": "HPV18_E6_F_R", "Forward": "TGCAAGACATAGAAATAACCT", "Reverse": "ACCTCTGTAAGTTCCAAC"},
    {"Organism": "Human papillomavirus (Broad Spectrum)", "PairID": "HPV_MY09_MY11", "Forward": "CGTCCMARRGGAWACTGATC", "Reverse": "GCMCAGGGWCATAAYAATGG"}
]

def main():
    csv_path = r"C:\Users\internet\Downloads\best_primers_db.csv"
    existing_rows = []
    if os.path.exists(csv_path):
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                existing_rows.append(row)

    print(f"[+] Loaded {len(existing_rows)} existing primer pairs.")
    existing_ids = set(r['PairID'] for r in existing_rows)

    added = 0
    for vp in VIRAL_PRIMERS:
        if vp['PairID'] not in existing_ids:
            existing_rows.append(vp)
            existing_ids.add(vp['PairID'])
            added += 1

    print(f"[+] Added {added} diagnostic viral primer pairs!")
    print(f"[+] Total primer pairs in database: {len(existing_rows)}")

    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Organism', 'PairID', 'Forward', 'Reverse'])
        writer.writeheader()
        writer.writerows(existing_rows)

    print(f"[+] Successfully saved updated database to: {csv_path}")

if __name__ == "__main__":
    main()
