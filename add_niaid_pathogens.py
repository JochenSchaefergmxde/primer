import pandas as pd
import re
import os

IUPAC = re.compile(r'^[ACGTNRYMKSwbdhvnACGTURYMKSWBDHVN]+$')

niaid_priority_additions = [
    # =========================================================================
    # NIAID CATEGORY A PRIORITY PATHOGENS
    # =========================================================================
    # Viral Hemorrhagic Fevers - Filoviruses
    ('Ebola virus (Zaire)', 'Ebola_Zaire_NP_Trombley', 'CAGTCTCAGTGACAGGCATTGA', 'ACATCTGTGTTGATGTTCCATTGA'),
    ('Ebola virus (Zaire)', 'Ebola_Zaire_GP_Trombley', 'TGGGCTGAAAACTGCTACAATC', 'CTTTGTGCACATACCGGCAC'),
    ('Ebola virus (Sudan)', 'Ebola_Sudan_GP_Drosten', 'AAGTGGTATCCAGCACCTGG', 'GCTTTCCATCCATTGCATCC'),
    ('Marburg virus', 'Marburg_virus_NP_Towner', 'ATTTTCACTGAAGGTGGTCAGTG', 'TCGGCACTTACTGAACCTAAGA'),
    ('Marburg virus', 'Marburg_virus_VP40_Drosten', 'TCCTCGTCAACAACAACCTC', 'CCTTGTCGTCCTCTTCGATG'),

    # Viral Hemorrhagic Fevers - Arenaviruses
    ('Lassa virus', 'Lassa_virus_GPC_Vieth', 'ACCGGGGATCCTAGGCATTT', 'CCCAAAGTACATGGCATCATCA'),
    ('Lassa virus', 'Lassa_virus_L_gene_Nikisins', 'ACCGACAACTTGATTGGTTG', 'CAGTATCCATCGGTTTCTTC'),
    ('Junin virus (Argentine HF)', 'Junin_virus_NP_Palacios', 'TCACTGTCCTTTCCTCCCTC', 'TCCATGATGTCCTCGCTGTC'),
    ('Machupo virus (Bolivian HF)', 'Machupo_virus_NP_Palacios', 'CAATGCTGAGCTGTCCTGTG', 'AGATGTTGTCCTCGCCCTTG'),
    ('Guanarito virus (Venezuelan HF)', 'Guanarito_virus_NP_Palacios', 'AAGCTGACTTTGAGTGTTGG', 'TCCATTTCCTCAGTGCCATC'),
    ('Chapare virus', 'Chapare_virus_NP_F_R', 'CTCTTGATGCTGTCATCCAG', 'CCATAGTAGTCTTCGCCCTC'),
    ('Lujo virus', 'Lujo_virus_GPC_F_R', 'ACACTGACTGTCACCATCAA', 'GGTCGTTGAACTCCTTGTTG'),

    # Viral Hemorrhagic Fevers - Bunyavirales
    ('Crimean-Congo hemorrhagic fever virus (CCHFV)', 'CCHFV_S_segment_Duh', 'TGGACACCTTCTCACAGACTC', 'TCAGCAAGCTCCAAATATCC'),
    ('Rift Valley fever virus (RVFV)', 'RVFV_L_segment_Bird', 'TGAAAATTCCTGAGACACATGG', 'ACTTCCTTGCATCATCTGATG'),
    ('Hantavirus (Pan-Hanta)', 'Pan_Hantavirus_S_segment_Klempa', 'ATGAACAACAACAGYATGTC', 'CCWGTBCCRTCATCATTGAT'),
    ('Sin Nombre hantavirus', 'Sin_Nombre_virus_S_segment', 'GGAAATGGTCCTCGTCATTC', 'GATGTTGTTGCTTGTCCTG'),
    ('Puumala hantavirus', 'Puumala_virus_S_segment', 'CAATGCTGTTTCTGCTGATG', 'TGCTTCATCACTGTCTTCATC'),

    # Orthopoxviruses
    ('Variola virus (Smallpox)', 'Variola_virus_B9R_Kulesh', 'TGTCTACATCACTACTCCCA', 'CCGTATTTCTTTGCCAATGCT'),
    ('Monkeypox virus (Mpox)', 'Monkeypox_virus_G2R_G_CDC', 'GGAAAATGTAATTCTCGAGGCTG', 'AAGTGCTTGATTTGTTGATTCC'),
    ('Monkeypox virus (Mpox)', 'Monkeypox_virus_C19L_generic', 'ATGACACCAGAAAACGACGA', 'CTCGGCTTGAGTGTAATCTA'),
    ('Orthopoxvirus (Pan-Orthopox)', 'Orthopoxvirus_E9L_pol_CDC', 'CGATGCAATCTCAGCATCTC', 'TGACCAAACTACCGACACCA'),

    # High-Consequence Bacterial Toxins
    ('Clostridium botulinum', 'Clostridium_botulinum_bontA_neurotoxin', 'GTGATACAAACACAGAAATC', 'CTTGTTGCTTGTAATAACTG'),
    ('Clostridium botulinum', 'Clostridium_botulinum_bontB_neurotoxin', 'TGGAATGTTTCGTTTAATCC', 'ACTACTCCTCCAGCATCAAA'),
    ('Clostridium botulinum', 'Clostridium_botulinum_bontE_neurotoxin', 'CAGATACACCGCTTAGTGAT', 'TGAACCATTAGTAGCTGCAT'),

    # =========================================================================
    # NIAID CATEGORY B PRIORITY PATHOGENS
    # =========================================================================
    ('Chlamydia psittaci (Psittacosis)', 'Chlamydia_psittaci_cps_ompA_Sachse', 'TGCAATCATAGCCGCTATTC', 'ACGATTTGACCTGTCGTCCT'),
    ('Clostridium perfringens', 'Clostridium_perfringens_etx_epsilon', 'GCGGTGTTTCAACATCTACA', 'TCCAAACTCACTCCTTCATC'),
    ('Clostridium perfringens', 'Clostridium_perfringens_cpe_enterotoxin', 'GGAGATGGTTGGATATTAGG', 'CCAACTCTACTTTTCCTTCC'),
    ('Yersinia enterocolitica', 'Yersinia_enterocolitica_ail_ISO18867', 'ACTCGATGATAACTGGGGAG', 'CCCCCAGTAATCCATAAAGG'),
    ('Cyclospora cayetanensis', 'Cyclospora_cayetanensis_18S_F_R', 'ATGTCCTCCCAAGCCCTC', 'TTAAGCGTCGAACCGTGC'),

    # Mosquito-borne Encephalitis / Alphaviruses & Flaviviruses
    ('Zika virus', 'Zika_virus_prME_CDC_Lanciotti', 'CCGCTGCCCAACACAAG', 'CCACTAACGTTCTTTTGCAGACAT'),
    ('Japanese encephalitis virus', 'Japanese_encephalitis_virus_NS3', 'TGTGTCTACAACATCATGGG', 'GTTGTTGCTTGTCCTGCTTCTT'),
    ('Venezuelan equine encephalitis virus (VEEV)', 'VEEV_nsP1_Linssen', 'CCGTCACCCTTGATGTTTCA', 'TGCTCCCGTCATCTTTACTC'),
    ('Eastern equine encephalitis virus (EEEV)', 'EEEV_nsP1_Lambert', 'CCACTAGCCAGGACGTCATT', 'TGACAGCGTCCACATTCCTC'),
    ('Western equine encephalitis virus (WEEV)', 'WEEV_nsP1_Lambert', 'TCGTCCTCGTCATCATCTTC', 'ACATCATCAGCGTCGTAGTC'),
    ('St. Louis encephalitis virus', 'SLEV_envelope_F_R', 'TCAGTGTCACAGCCCTTGTC', 'CCTTGCTTGTCCTCCATCTC'),
    ('La Crosse virus (California serogroup)', 'La_Crosse_virus_S_segment', 'ATGGACAACTTGATTGGTTG', 'CAGTATCCATCGGTTTCTTC'),

    # =========================================================================
    # NIAID CATEGORY C PRIORITY PATHOGENS & EMERGING THREATS
    # =========================================================================
    ('Nipah virus', 'Nipah_virus_N_gene_CDC', 'TCAGCAGGAAGGCAAGAGAGTAA', 'CCCCTTCATCGATATCTTGATCA'),
    ('Hendra virus', 'Hendra_virus_M_gene_Smith', 'AATCGGTCTCGCTCTTTGGT', 'CGCTTCCTCGTAAGTGTTGA'),
    ('Oropouche virus (OROV)', 'Oropouche_virus_S_segment_CDC', 'GACAAGTGCTCAATGCTKGTGT', 'CGTTGTCCGGSACTGGATT'),
    ('Heartland virus (Bandavirus)', 'Heartland_virus_S_segment_CDC', 'AAGATGGTCCTCGTCATCTC', 'TCCATCAGCGTCTTCATCAT'),
    ('Bourbon virus (Thogotovirus)', 'Bourbon_virus_PB1_CDC', 'TGCAACACAGGACTCATTCC', 'ACATCCAGCGTCTTGTCATC'),
    ('Avian influenza A (H5N1)', 'Influenza_A_H5_HA_CDC', 'ACACATGCYCARGACATACT', 'CTYTGRTTYAGTGTTGATGT'),
    ('Avian influenza A (H7N9)', 'Influenza_A_H7_HA_CDC', 'AATTCCTCGTGCCTTGTTTG', 'TCCTTCTTTGACACCATCGA'),

    # Emerging Fungi (NIAID Research Priority)
    ('Coccidioides immitis/posadasii (Valley fever)', 'Coccidioides_Ag2_PRA_F_R', 'CTTGTCGCCCACTACGTCAT', 'TCGAGCGCCTTGATGTCCTT'),
    ('Histoplasma capsulatum', 'Histoplasma_capsulatum_100kDa_F_R', 'GCCGAGTTCAATACCATCGT', 'TGATCGGACGGTCGATGTC'),
    ('Blastomyces dermatitidis', 'Blastomyces_dermatitidis_BAD1_F_R', 'TCATCGTCCTCATCATCTTC', 'CCATGATGTCCTCGCTGTC'),
    ('Rhizopus arrhizus (Mucormycosis)', 'Rhizopus_arrhizus_18S_F_R', 'AGAGTTTGATCCTGGCTCAG', 'CCATCATCAGAAAGAATCATCA'),
    ('Mucor circinelloides (Mucormycosis)', 'Mucor_circinelloides_ITS_F_R', 'GGTGTATGTCTTTAAGTCATAC', 'ACGGCTACCTTGTTACGACT'),

    # Emerging Rickettsiales & Tick Pathogens
    ('Ehrlichia chaffeensis (HME)', 'Ehrlichia_chaffeensis_16S_HE1_HE3', 'CAATTGCTTATAACCTTTTGGTTATAAAT', 'TATAGGTACCGTCATTATCTTCCCTAT'),
    ('Ehrlichia ewingii', 'Ehrlichia_ewingii_16S_EE72_EE379', 'CAATTCCTAAATAGTCTCTGACTATTTAAC', 'TAAATCGGTTTCACCAAAATTC')
]

print(f"Total NIAID additions to process: {len(niaid_priority_additions)}")

# Validation
for org, pid, f, r in niaid_priority_additions:
    f = f.strip().upper()
    r = r.strip().upper()
    assert 15 <= len(f) <= 45, f"Length F error: {pid} ({len(f)})"
    assert 15 <= len(r) <= 45, f"Length R error: {pid} ({len(r)})"
    assert IUPAC.match(f), f"IUPAC F error: {pid} ({f})"
    assert IUPAC.match(r), f"IUPAC R error: {pid} ({r})"

db = pd.read_csv('best_primers_db.csv')
print(f"Current DB size: {len(db)}")

existing = set(zip(db['Forward'].str.upper(), db['Reverse'].str.upper()))

new_rows = []
for org, pid, f, r in niaid_priority_additions:
    pair = (f.upper(), r.upper())
    if pair not in existing:
        new_rows.append({
            'Organism': org,
            'PairID': pid,
            'Forward': f.upper(),
            'Reverse': r.upper()
        })
        existing.add(pair)

print(f"Appending {len(new_rows)} new unique NIAID priority pathogen pairs.")

updated = pd.concat([db, pd.DataFrame(new_rows)], ignore_index=True)
updated = updated.drop_duplicates(subset=['Forward', 'Reverse']).reset_index(drop=True)
updated = updated.sort_values(by=['Organism', 'PairID']).reset_index(drop=True)

# Write to disk and to Downloads
updated.to_csv('best_primers_db.csv', index=False)
downloads_path = os.path.join(r"C:\Users\internet\Downloads", "best_primers_db.csv")
updated.to_csv(downloads_path, index=False)

print(f"Updated DB size: {len(updated)} primer pairs.")
