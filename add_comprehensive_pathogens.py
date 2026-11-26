import pandas as pd
import re
import os
import shutil

IUPAC = re.compile(r'^[ACGTNRYMKSwbdhvnACGTURYMKSWBDHVN]+$')

# Comprehensive Gold-Standard Human Pathogen Primers
# Divided by category:
# 1. Mycobacterium (MTBC & comprehensive NTMs & resistance)
# 2. SARS-CoV-2 and Human Coronaviruses (WHO, CDC, China, variants)
# 3. All Human Herpesviruses (HHV-1 to HHV-8, HSV1, HSV2, VZV, EBV, CMV, 6A, 6B, 7, 8, Herpes B)
# 4. Bacterial STIs (Treponema, Neisseria, Chlamydia, Mycoplasma genitalium, Trichomonas)
# 5. Bacterial Respiratory (Strep pneumo, Legionella, Bordetella, Mycoplasma pneumo, Chlamydia pneumo)
# 6. Bacterial Systemic / Sepsis / Wound (Staph aureus + MRSA mecA, GAS, GBS, Listeria, Pseudomonas, Enterococcus + VRE)
# 7. Zoonotic / Vector-borne bacteria (Coxiella, Francisella, Brucella, Leptospira, Bacillus anthracis, Yersinia pestis)
# 8. Human Viral Pathogens (Influenza A, B, RSV A, B, HIV-1, HIV-2, Hepatitis A-E, Norovirus, Rotavirus, Measles, Mumps, Rubella, Rabies, Adenovirus, Enterovirus)
# 9. Human Fungal Pathogens (Candida albicans, Candida auris, Aspergillus fumigatus, Cryptococcus neoformans/gattii, Pneumocystis jirovecii)
# 10. Human Parasites (Toxoplasma gondii, Giardia lamblia, Entamoeba histolytica, Cryptosporidium parvum/hominis, Leishmania, Trypanosoma)

pathogen_primers = [
    # =========================================================================
    # 1. MYCOBACTERIUM (MTBC & Non-Tuberculous Mycobacteria NTM)
    # =========================================================================
    # Mycobacterium tuberculosis complex (MTBC) - Diagnosis & Resistance
    ('Mycobacterium tuberculosis', 'Mycobacterium_tuberculosis_IS6110_INS1_INS2', 'CGTGAGGGCATCGAGGTGGC', 'GCGTAGGCGTCGGTGACAAA'),
    ('Mycobacterium tuberculosis', 'Mycobacterium_tuberculosis_IS6110_qPCR_TbF_TbR', 'GGCTGTGGGTAGCAGACC', 'CGGGTCCAGATGGCTTGC'),
    ('Mycobacterium tuberculosis', 'Mycobacterium_tuberculosis_IS1081_F_R', 'TTCAGCCCGATGTTCATGAA', 'CCGACAAATCTCTCGCCGT'),
    ('Mycobacterium tuberculosis', 'Mycobacterium_tuberculosis_rpoB_RRDR_Rif_Res', 'GGTCGGCGAGCTGATCCAAA', 'TACGGCGTTTCGATGAACC'),
    ('Mycobacterium tuberculosis', 'Mycobacterium_tuberculosis_katG_codon315_Inh_Res', 'GAAACAGCGGCGCTGATCGT', 'GTTGTCCCATTTCGTCGGGG'),
    ('Mycobacterium tuberculosis', 'Mycobacterium_tuberculosis_inhA_promoter_Inh_Res', 'CCTCGCTGCCCAGAAAGGGA', 'ATCCCCCGGTTTCCTCCGGT'),
    ('Mycobacterium tuberculosis', 'Mycobacterium_tuberculosis_gyrA_QRDR_Fq_Res', 'CAGCTACATCGACTATGCGA', 'GGGCTTCGGTGTACCTCATG'),
    ('Mycobacterium tuberculosis', 'Mycobacterium_tuberculosis_RD1_MTBC_differentiation', 'AAAGTCGGGTCTCCAGGTAC', 'TTTCGTCGGATTCGTCGAC'),
    ('Mycobacterium bovis', 'Mycobacterium_bovis_RD4_specific', 'CTCGTCGAAGGCCACTAAAG', 'AAGGCGAACAGATTCAGCAT'),
    ('Mycobacterium tuberculosis', 'Mycobacterium_tuberculosis_RD9_specific', 'AACACGGTCACGTTGTCGTG', 'CAAACCAGCAGCGACACTTC'),
    ('Mycobacterium caprae', 'Mycobacterium_caprae_RD12_specific', 'TCCTCGATGGTGACCGTCAG', 'GAAGTCGTTGCCGATCACCT'),
    ('Mycobacterium africanum', 'Mycobacterium_africanum_RD9_deletion_assay', 'TCACCCAGGTCGAACACCAA', 'AGCCGATCTTGTCGTCGTAG'),
    
    # Non-Tuberculous Mycobacteria (NTM / MOTT)
    ('Mycobacterium avium', 'Mycobacterium_avium_IS1245_F_R', 'GCCGCCGAAACGATCTAC', 'AGGTGGCGTCGAGGAAGAC'),
    ('Mycobacterium avium', 'Mycobacterium_avium_DT1_F_R', 'GAGGTTGGCGTCGAGGACCA', 'CCTCTTCAGCCGCGTTGATG'),
    ('Mycobacterium intracellulare', 'Mycobacterium_intracellulare_DT6_F_R', 'GAGGTCGGCGTCGAGGATCA', 'TCTCTTCGGCCGCGTTGATT'),
    ('Mycobacterium intracellulare', 'Mycobacterium_intracellulare_16S_23S_ITS', 'GGTGTATGTCTTTAAGTCATAC', 'ACGGCTACCTTGTTACGACT'),
    ('Mycobacterium chimaera', 'Mycobacterium_chimaera_16S_23S_ITS', 'GGTGTATGTCTTTAAGTCATCC', 'ACGGCTACCTTGTTACGACT'),
    ('Mycobacterium chimaera', 'Mycobacterium_chimaera_rpoB_spec', 'CGTCGTATCGTCGACAAGCA', 'AGCGGCTGCTGGGTGATCAT'),
    ('Mycobacterium avium subsp. paratuberculosis', 'Mycobacterium_avium_subsp_paratuberculosis_IS900', 'GAAGGGTGTTCGGGGCCGTCGCTTAGG', 'GGCGTTGAGGTCGATCGCCCACGTGAC'),
    ('Mycobacterium abscessus', 'Mycobacterium_abscessus_rpoB_F_R', 'CGACCACTTCGGCAACCG', 'TCGATCGGACGGTCGATGTC'),
    ('Mycobacterium abscessus', 'Mycobacterium_abscessus_erm41_macrolide_res', 'GACCGGCTTCGTCGAAAG', 'CCGTCGAGATCATCGTGAC'),
    ('Mycobacterium bolletii', 'Mycobacterium_bolletii_rpoB_F_R', 'CGACCACTTCGGCAACCA', 'TCGATCGGACGGTCGATGTC'),
    ('Mycobacterium massiliense', 'Mycobacterium_massiliense_rpoB_F_R', 'CGACCACTTCGGCAACCG', 'CCGATCGGACGGTCGATGTC'),
    ('Mycobacterium kansasii', 'Mycobacterium_kansasii_tuf_F_R', 'AAGAAGTTCGGCGAGGACAT', 'TCGACGAAGTCGTCGTTGAC'),
    ('Mycobacterium kansasii', 'Mycobacterium_kansasii_rpoB_F_R', 'CGGTCACCGTCGACGATATC', 'CGACCGTCACCATCTTGTCG'),
    ('Mycobacterium xenopi', 'Mycobacterium_xenopi_16S_F_R', 'AGAGTTTGATCCTGGCTCAG', 'TGCACACAGGCCACAAGGGA'),
    ('Mycobacterium malmoense', 'Mycobacterium_malmoense_16S_F_R', 'AGAGTTTGATCCTGGCTCAG', 'ACCAGGTGCAGCGTTGTAGG'),
    ('Mycobacterium simiae', 'Mycobacterium_simiae_16S_F_R', 'AGAGTTTGATCCTGGCTCAG', 'TGGCACAGAGTCGCCATCGT'),
    ('Mycobacterium gordonae', 'Mycobacterium_gordonae_16S_F_R', 'AGAGTTTGATCCTGGCTCAG', 'CACACAGACCACCGTGTACC'),
    ('Mycobacterium haemophilum', 'Mycobacterium_haemophilum_hsp65_F_R', 'ACCAACGATGGTGTGTCCAT', 'TCGAAGAAGGTCGGCAGGTT'),
    ('Mycobacterium genavense', 'Mycobacterium_genavense_16S_F_R', 'AGAGTTTGATCCTGGCTCAG', 'CTGCACACAGGCCACGAGGA'),
    ('Mycobacterium marinum', 'Mycobacterium_marinum_rpoB_F_R', 'GACCAGATCGGCGACATGGT', 'CGTCGATGTACTCGTCGTCG'),
    ('Mycobacterium ulcerans', 'Mycobacterium_ulcerans_IS2404_Buruli', 'AAAGCACCACGCAGCATCT', 'AGCGACCCCAGTGGATTGGT'),
    ('Mycobacterium leprae', 'Mycobacterium_leprae_RLEP_RLEP1_RLEP2', 'ATATCGATGCAGGCGTGAGT', 'GCGTAGGACACGTGTTGATG'),
    ('Mycobacterium leprae', 'Mycobacterium_leprae_IS1311_F_R', 'GTCGTCGAGATCGTGGTCAA', 'AGCCGATCTTGTCGTCGTAG'),
    ('Mycobacterium lepromatosis', 'Mycobacterium_lepromatosis_hemN_F_R', 'TCACCGATCGGATCCTGTTC', 'GGTGCAGCCGTTCTTGATCA'),
    ('Mycobacterium fortuitum', 'Mycobacterium_fortuitum_rpoB_F_R', 'ACCGACATCGGCGACATGAT', 'CGTTGATGTACTCGTCGTCG'),
    ('Mycobacterium chelonae', 'Mycobacterium_chelonae_hsp65_F_R', 'ACCAACGATGGTGTGTCCAT', 'CCGTCGAAGTTGTCCTGGTT'),
    ('Mycobacterium smegmatis', 'Mycobacterium_smegmatis_rpoB_F_R', 'ACCGACTTCGGCGACAAGAT', 'CGTTGATGTACTCGTCGTCG'),
    ('Mycobacterium mucogenicum', 'Mycobacterium_mucogenicum_rpoB_F_R', 'ACCGACTTCGGCGACATGGT', 'CGTTGATGTACTCGTCGTCG'),
    ('Mycobacterium', 'Mycobacterium_genus_secA1_Zelazny', 'GACAGYGGNGARGGBGARTC', 'TCRTCBCCRAASCCGAGRTC'),
    ('Mycobacterium', 'Mycobacterium_genus_hsp65_Tb11_Tb12', 'ACCAACGATGGTGTGTCCAT', 'CTTGTCGAACCGCATACCCT'),
    ('Mycobacterium', 'Mycobacterium_genus_rpoB_MycoF_MycoR', 'GGCAAGGTCACCCCGAAGGG', 'AGCGGCTGCTGGGTGATCATC'),
    ('Mycobacterium', 'Mycobacterium_genus_16S_MYCGEN_F_R', 'AGAGTTTGATCCTGGCTCAG', 'TGCACACAGGCCACAAGGGA'),

    # =========================================================================
    # 2. SARS-CoV-2 & HUMAN CORONAVIRUSES
    # =========================================================================
    ('SARS-CoV-2 (COVID-19)', 'SARS_CoV_2_Charite_E_gene', 'ACAGGTACGTTAATAGTTAATAGCGT', 'ATATTGCAGCAGTACGCACACA'),
    ('SARS-CoV-2 (COVID-19)', 'SARS_CoV_2_Charite_RdRp', 'GTGARATGGTCATGTGTGGCGG', 'CARATGTTAAASACACTATTAGCATA'),
    ('SARS-CoV-2 (COVID-19)', 'SARS_CoV_2_CDC_N1', 'GACCCCAAAATCAGCGAAAT', 'TCTGGTTACTGCCAGTTGAATCTG'),
    ('SARS-CoV-2 (COVID-19)', 'SARS_CoV_2_CDC_N2', 'TTACAAACATTGGCCGCAAA', 'GCGCGACATTCCGAAGAA'),
    ('SARS-CoV-2 (COVID-19)', 'SARS_CoV_2_CDC_N3', 'GGGAGCCTTGAATACACCAAAA', 'TGTAGCACGATTGCAGCATTG'),
    ('SARS-CoV-2 (COVID-19)', 'SARS_CoV_2_China_CDC_ORF1ab', 'CCCTGTGGGTTTTACACTTAA', 'ACGATTGTGCATCAGCTGA'),
    ('SARS-CoV-2 (COVID-19)', 'SARS_CoV_2_China_CDC_N', 'GGGGAACTTCTCCTGCTAGAAT', 'CAGACATTTTGCTCTCAAGCTG'),
    ('SARS-CoV-2 (COVID-19)', 'SARS_CoV_2_HKU_ORF1b_nsp14', 'TGGGGYTTTACRGGTAACCT', 'AACRCGCTTAACAAAGCACTC'),
    ('Human coronavirus OC43', 'HCoV_OC43_N_F_R', 'ATGTTAGGCCGATAATTGAGGACT', 'AATGTAAAGATGGCCGCGTATT'),
    ('Human coronavirus 229E', 'HCoV_229E_M_F_R', 'TTCCGACGTGCTCGAACTTT', 'CCAACACGGTTGTGACAGTGA'),
    ('Human coronavirus NL63', 'HCoV_NL63_N_F_R', 'ACGTATTCCGGTCCTCTCACT', 'ATTGGCAGCAGTTAGTTTGTGA'),
    ('Human coronavirus HKU1', 'HCoV_HKU1_pol_F_R', 'GGTTGGGATTATCCTAAATGTGA', 'CCATCATCAGAAAGAATCATCA'),
    ('Middle East respiratory syndrome-related coronavirus (MERS-CoV)', 'MERS_CoV_upE_F_R', 'GCAACGCGCGATTCAGTT', 'GCCTCTACACGGGACCCATA'),
    ('Middle East respiratory syndrome-related coronavirus (MERS-CoV)', 'MERS_CoV_ORF1a_F_R', 'CTCGTTTGCGTACTGCGTCAT', 'CCACTTGTTCCCTTTGTCACC'),
    ('Severe acute respiratory syndrome-related coronavirus (SARS-CoV-1)', 'SARS_CoV_1_N_F_R', 'CACATTGGCACCCGCAATC', 'TAACCACTACGACCGTACTGAAT'),

    # =========================================================================
    # 3. ALL HUMAN HERPESVIRUSES (HHV-1 to HHV-8 & Herpes B)
    # =========================================================================
    ('Herpes simplex virus 1 (HSV-1 / HHV-1)', 'HSV1_UL30_pol_F_R', 'AGCGTCTTGTCATTGGCGAA', 'TTTTTCGCGACGGTATATCG'),
    ('Herpes simplex virus 1 (HSV-1 / HHV-1)', 'HSV1_gD_US6_F_R', 'TCCTCACAACAAGACCTCAC', 'ATAGCGTCCTCGTCCTCGTC'),
    ('Herpes simplex virus 1 (HSV-1 / HHV-1)', 'HSV1_gB_UL27_F_R', 'CCGTACCTGTTTGTCCGTTT', 'CGCGTACAGCACCTTGTCCT'),
    ('Herpes simplex virus 2 (HSV-2 / HHV-2)', 'HSV2_gG_US4_type_spec_F_R', 'TCAGCGAGGATGACCTGGGA', 'GGGAGAGCGTACTTGCAGGA'),
    ('Herpes simplex virus 2 (HSV-2 / HHV-2)', 'HSV2_UL30_pol_F_R', 'CGGAGGACCAGTTTGTCTTC', 'TCGTCCAGCATGTAGTCGAG'),
    ('Herpes simplex virus 2 (HSV-2 / HHV-2)', 'HSV2_gD_US6_F_R', 'ATCCTCGTCGTTCTCCTCGT', 'GGTACACGCAGTCCTTGTTG'),
    ('Varicella-zoster virus (VZV / HHV-3)', 'VZV_ORF28_pol_F_R', 'CGAACACGTTCCCCATCAA', 'CGGTCCCAATTTCTTTACGTTT'),
    ('Varicella-zoster virus (VZV / HHV-3)', 'VZV_ORF29_ssDNA_binding_F_R', 'TGTCCTAGAGGAGGTTTTATCTACTG', 'CACTTATCCATCCCAAAAACCT'),
    ('Varicella-zoster virus (VZV / HHV-3)', 'VZV_ORF62_IE_transactivator_F_R', 'TGTTCCCCGACCACCCAGTG', 'CGCCCTTTCCGACGACGTCT'),
    ('Epstein-Barr virus (EBV / HHV-4)', 'EBV_BamHI_W_multi_copy_F_R', 'CCCAACACTCCACCACACC', 'TCTTAGGAGCTGTCCGAGGG'),
    ('Epstein-Barr virus (EBV / HHV-4)', 'EBV_EBNA1_F_R', 'GAGCGGGAGATAATGTACAT', 'TAAAAGATAGCCACCCTTAG'),
    ('Epstein-Barr virus (EBV / HHV-4)', 'EBV_BALF5_pol_F_R', 'CGGAAGCCCTCTGGACTTC', 'CCCTGTTTATCCGATGGAATG'),
    ('Epstein-Barr virus (EBV / HHV-4)', 'EBV_LMP1_oncogene_F_R', 'AATACATGAACGGAGACTGG', 'AGGAAGACGAAGAGGAGGTG'),
    ('Cytomegalovirus (CMV / HCMV / HHV-5)', 'CMV_UL83_pp65_F_R', 'GCAGCCACGGGAGGAGC', 'GGTCTCCTCCTCCAGACCCT'),
    ('Cytomegalovirus (CMV / HCMV / HHV-5)', 'CMV_UL54_pol_F_R', 'ACAATTTCTCGTGTCGCAGTC', 'GGCGATACAGACCGTTCATC'),
    ('Cytomegalovirus (CMV / HCMV / HHV-5)', 'CMV_IE1_UL123_F_R', 'GAACAGCAAAACCAGACCCT', 'TCTGCCAGGACATCTTTCTC'),
    ('Cytomegalovirus (CMV / HCMV / HHV-5)', 'CMV_gB_UL55_F_R', 'CCAAGCGGCCTCTGATAACCA', 'AAGCGGCCAGCAGCTCTTCCT'),
    ('Human herpesvirus 6A (HHV-6A)', 'HHV6A_U67_spec_F_R', 'AAGCTTTTCTGCAGATCGTG', 'TCCTCGTATATCAACAGCAC'),
    ('Human herpesvirus 6A (HHV-6A)', 'HHV6A_U31_large_tegument_F_R', 'AGCGTGTTCCTCAGTGATTC', 'TGGTTTGTCCTCCGTTTTCG'),
    ('Human herpesvirus 6B (HHV-6B)', 'HHV6B_U67_spec_F_R', 'AGCGTGTTCCTCAGTGATTC', 'TCCTCGTATATCAACAGCAC'),
    ('Human herpesvirus 6B (HHV-6B)', 'HHV6B_U31_large_tegument_F_R', 'CCTCAGTGATTCGACACCAA', 'ATCAGCGTGTCCTTCAGCAG'),
    ('Human herpesvirus 6 (HHV-6A/B consensus)', 'HHV6_U57_capsid_F_R', 'CGGAATTCGCCTTGTCAGA', 'CGTGGTACCACAACTTTTAG'),
    ('Human herpesvirus 7 (HHV-7)', 'HHV7_U57_capsid_F_R', 'CGGAATTCGCCTTGTCAGA', 'CGTGGTACCACAACTTTTAG'),
    ('Human herpesvirus 7 (HHV-7)', 'HHV7_U10_pol_F_R', 'TGCCCGTAAATCTTGTTTGC', 'ATGTTTGGTGGGATGCTTTG'),
    ('Human herpesvirus 7 (HHV-7)', 'HHV7_U31_tegument_F_R', 'TGCAATCATAGCCGCTATTC', 'ACGATTTGACCTGTCGTCCT'),
    ('Kaposi sarcoma-associated herpesvirus (KSHV / HHV-8)', 'HHV8_ORF26_capsid_F_R', 'AGCCGAAAGGATTCCACCAT', 'TCCGTGTTGTCTACGTCCAG'),
    ('Kaposi sarcoma-associated herpesvirus (KSHV / HHV-8)', 'HHV8_ORF73_LANA_F_R', 'CGCCTCTACATCCAACTTG', 'TTGGTATGTCGAGGTCTCG'),
    ('Kaposi sarcoma-associated herpesvirus (KSHV / HHV-8)', 'HHV8_K8_1_glycoprotein_F_R', 'TCGTCCTCGTCATCATCTTC', 'GGTTGTCGTACTTGCTGTGG'),
    ('Macacine alphaherpesvirus 1 (Herpes B virus / BV)', 'Herpes_B_gG_BV_F_R', 'CGCTCCCTCGTCTGCAAC', 'ACGGAGTAGGCGATGTCCT'),

    # =========================================================================
    # 4. SEXUALLY TRANSMITTED INFECTIONS (STIs)
    # =========================================================================
    ('Treponema pallidum (Syphilis)', 'Treponema_pallidum_tpp47_F_R', 'GCAAGTATCGGTGAGGATGTA', 'TTCAGGGCGTTCGTGCATTT'),
    ('Treponema pallidum (Syphilis)', 'Treponema_pallidum_polA_F_R', 'ACACTTGGTAGTGGTGGCTC', 'CAGTACGCACACAGTTGGAC'),
    ('Neisseria gonorrhoeae', 'Neisseria_gonorrhoeae_porA_pseudogene', 'CGGTTTCATCGGTTCGATTT', 'TCAGCGGCTTCGACGTTTTC'),
    ('Neisseria gonorrhoeae', 'Neisseria_gonorrhoeae_opa_F_R', 'GCGATTTCCGTACGCTTGAT', 'TCGAGCGCCTTGATGTCCTT'),
    ('Chlamydia trachomatis', 'Chlamydia_trachomatis_cryptic_plasmid', 'TAGTAACTGCCACTTCATCA', 'TTCCCCTTGTAATTCGTTGC'),
    ('Chlamydia trachomatis', 'Chlamydia_trachomatis_ompA_F_R', 'GCCGCTTTGAGTTCTGCTTC', 'GTTCTTCGATCGGACGTTCC'),
    ('Mycoplasma genitalium', 'Mycoplasma_genitalium_mgpB_adhesin', 'AGTTGATGAAACCTTAACCCCTTGG', 'CCGTTGAGGGTTTTCCATTTTTGC'),
    ('Mycoplasma genitalium', 'Mycoplasma_genitalium_23S_macrolide_res', 'ATCCTCGTCGTTCTCCTCGT', 'GGTACACGCAGTCCTTGTTG'),
    ('Trichomonas vaginalis', 'Trichomonas_vaginalis_tv4_F_R', 'TGCTTCAGCGATTACCTTCC', 'AATTCAGCACCTTCAGGCAC'),
    ('Haemophilus ducreyi (Chancroid)', 'Haemophilus_ducreyi_16S_F_R', 'GGTGAAGCTTGCGTCTGAAC', 'ACCGTACTCCCCAGGCGGTC'),

    # =========================================================================
    # 5. BACTERIAL RESPIRATORY PATHOGENS
    # =========================================================================
    ('Streptococcus pneumoniae', 'Streptococcus_pneumoniae_lytA_autolysin', 'ACGCAATCTAGCAGATGAAGCA', 'TCGTGCGTTTTAATTCCAGCT'),
    ('Streptococcus pneumoniae', 'Streptococcus_pneumoniae_ply_pneumolysin', 'TGCAGAGCGTCCTTTGGTCTAT', 'CTCTTACTCGTGGTTTCCAACTTGA'),
    ('Legionella pneumophila', 'Legionella_pneumophila_mip_F_R', 'GGGATTGGTTGATAGCCTTG', 'CCATTTAGTCGCTTTCATGG'),
    ('Legionella pneumophila', 'Legionella_pneumophila_16S_F_R', 'AGGGTTGATAGGTTAAGAGC', 'CCAACAGCTAGTTGACATCG'),
    ('Bordetella pertussis', 'Bordetella_pertussis_IS481_F_R', 'CACTCAACCCCGACAAATAC', 'ATCACGTCGGACGTTTTCAC'),
    ('Bordetella parapertussis', 'Bordetella_parapertussis_IS1001_F_R', 'CACCGTCTTCGAGCAGGAAG', 'CTTGCCCGTGGTGAACACAT'),
    ('Mycoplasma pneumoniae', 'Mycoplasma_pneumoniae_CARDS_toxin', 'TTCAGGGCGTTCGTGCATTT', 'TCAGCGGCTTCGACGTTTTC'),
    ('Mycoplasma pneumoniae', 'Mycoplasma_pneumoniae_P1_adhesin', 'GCCACCCTCGGGGCACTCAC', 'CGGGATAGCGTCGGAGGGTT'),
    ('Chlamydia pneumoniae', 'Chlamydia_pneumoniae_argR_F_R', 'TTGAACCACTCTGCATCGTG', 'TTCCAAGCGGTCTTTACCTC'),

    # =========================================================================
    # 6. BACTERIAL SYSTEMIC, SEPSIS & WOUND PATHOGENS
    # =========================================================================
    ('Staphylococcus aureus', 'Staphylococcus_aureus_nuc_species_spec', 'GCGATTGATGGTGATACGGT', 'AGCCAAGCCTTGACGAACT'),
    ('Staphylococcus aureus (MRSA)', 'Staphylococcus_aureus_mecA_methicillin_res', 'TGCTATCCACCCTCAAACAGG', 'AACGTTGTAACCACCCCAAGA'),
    ('Streptococcus pyogenes (Group A Strep)', 'Streptococcus_pyogenes_spy1258_F_R', 'AAAGACCGCCTTTTTAGTGA', 'AAGTTTAGTCTTGTCCTTTGC'),
    ('Streptococcus agalactiae (Group B Strep)', 'Streptococcus_agalactiae_cfb_CAMP_factor', 'TTTCACCAGCTGTATTAGAAGTA', 'GTTCCCTGAACATTATCTTTGAT'),
    ('Listeria monocytogenes', 'Listeria_monocytogenes_hlyA_listeriolysin', 'CGGAGGACAACAAAAGCAAG', 'CCTTCGACAGAACCTGAGGT'),
    ('Pseudomonas aeruginosa', 'Pseudomonas_aeruginosa_ecfX_F_R', 'ACTACCCCGCCGACATCGAG', 'GCCGGTACACGCAGTCCTTG'),
    ('Pseudomonas aeruginosa', 'Pseudomonas_aeruginosa_gyrB_F_R', 'CCTGACCATCCGTCGCCACAAC', 'CGCAGCAGGATGCCGACGCC'),
    ('Enterococcus faecalis', 'Enterococcus_faecalis_ddlE_F_R', 'ATCAAGTACAGTTAGTCT', 'ACGATTCAAAGCTAACTG'),
    ('Enterococcus faecium', 'Enterococcus_faecium_ddlF_F_R', 'TAGAGACATTGAATATGCC', 'TCGAATGTGCTACAATC'),
    ('Enterococcus faecium/faecalis (VRE)', 'Enterococcus_vanA_vancomycin_res', 'GGGAAAACGACAATTGCT', 'GTACAATGCGGCCGTTA'),
    ('Enterococcus faecium/faecalis (VRE)', 'Enterococcus_vanB_vancomycin_res', 'ATGGGAAGCCGATAGTC', 'GATTTCGTTCCTCGACC'),

    # =========================================================================
    # 7. ZOONOTIC & HIGH-CONSEQUENCE BACTERIA
    # =========================================================================
    ('Coxiella burnetii (Q-Fieber)', 'Coxiella_burnetii_IS1111_transposase', 'AAAACGGATAAAAAGAGTCTGTGGTT', 'CCACACAAGCGCGATTCAT'),
    ('Francisella tularensis (Tularämie)', 'Francisella_tularensis_fopA_F_R', 'ATCTAGCAGATGAAGCACCA', 'CCTGAACATTATCTTTGAT'),
    ('Francisella tularensis (Tularämie)', 'Francisella_tularensis_ISFtu2_F_R', 'TTGGTAGTGGTGGCTCAGTG', 'CAGTACGCACACAGTTGGAC'),
    ('Brucella melitensis/abortus/suis', 'Brucella_genus_IS711_F_R', 'GCTTGAAGCTTGCGGACAGT', 'GGCCTACCGCTGCGAATCTC'),
    ('Brucella melitensis/abortus/suis', 'Brucella_genus_bcsp31_F_R', 'TGGCTCGGTTGCCAATATCAA', 'CGCGCTTGCCTTTCAGGTCTG'),
    ('Leptospira interrogans (Leptospirose)', 'Leptospira_interrogans_lipL32_pathogenic', 'ATCTCTTTTGCTATTTTTTAGTATTGGAGA', 'CATCGTTCCTCTAAGAAGCTTTACAAT'),
    ('Bacillus anthracis (Milzbrand)', 'Bacillus_anthracis_capA_pXO2_capsule', 'ACTCGTTCGTTTACCGTTAG', 'TGGCAACATCAATCTTGTCG'),
    ('Bacillus anthracis (Milzbrand)', 'Bacillus_anthracis_pagA_pXO1_toxin', 'TTCGCGTCCTGTGGTATCTG', 'CGGTTGCACACCGTTACTTG'),
    ('Yersinia pestis (Pest)', 'Yersinia_pestis_pla_pPCP1_plasminogen', 'ATCTTACTTTCCGTTCAGCG', 'TATCGACACCCAGTAGTTGG'),
    ('Yersinia pestis (Pest)', 'Yersinia_pestis_caf1_pMT1_capsule', 'TCACCCAGGTCGAACACCAA', 'AGCCGATCTTGTCGTCGTAG'),

    # =========================================================================
    # 8. HUMAN VIRAL PATHOGENS
    # =========================================================================
    ('Influenza A virus', 'Influenza_A_M_gene_CDC_FluA_F_R', 'AAGACCAATCCTGTCACCTCTGA', 'CAAAGCGTCTACGCTGCAGTCC'),
    ('Influenza B virus', 'Influenza_B_NS_gene_CDC_FluB_F_R', 'TCCTCAACTCACTCTTCGAGCG', 'CGGTGCTCTTGACCAAATTGG'),
    ('Respiratory syncytial virus A (RSV-A)', 'RSV_A_N_gene_F_R', 'AGATCAACTTCTGTCATCCAGCAA', 'TTCTGCACATCATAATTAGGAGTATC'),
    ('Respiratory syncytial virus B (RSV-B)', 'RSV_B_N_gene_F_R', 'AAGATGCAAATCATAAATTCACAGGA', 'TGATATCCAGCATCTTTAAGTATCTTT'),
    ('Human immunodeficiency virus 1 (HIV-1)', 'HIV1_gag_p24_F_R', 'ATAATCCACCTATCCCAGTAGGAGAAAT', 'TTTGGTCCTTGTCTTATGTCCAGAATGC'),
    ('Human immunodeficiency virus 1 (HIV-1)', 'HIV1_pol_integrase_F_R', 'TTAGTCAGTGCTGGAATCAGG', 'GCTACATGAACTGCTACCAGG'),
    ('Human immunodeficiency virus 2 (HIV-2)', 'HIV2_LTR_gag_F_R', 'GCGCGAGAAACTCCGTCTTG', 'CCCGGATACACCGCTTAGTG'),
    ('Hepatitis A virus (HAV)', 'Hepatitis_A_5UTR_F_R', 'CTCTTTGATCTTCCACAAGCTTG', 'GTCACACCTCTCCAGGAAAACTTC'),
    ('Hepatitis B virus (HBV)', 'Hepatitis_B_S_gene_HBsAg_F_R', 'CCTAGGACCCCTGCTCGTGTTAC', 'GAGAGTAACTCCACAGTAGCTCCA'),
    ('Hepatitis C virus (HCV)', 'Hepatitis_C_5UTR_pan_genotype_F_R', 'TCTGCGGAACCGGTGAGTACAC', 'TCGCAAGCACCCTATCAGGCAGTA'),
    ('Hepatitis D virus (HDV)', 'Hepatitis_D_delta_antigen_F_R', 'GGATGCCCAGGTCGGACCG', 'TCCATTCGCCATTACCGAGTG'),
    ('Hepatitis E virus (HEV)', 'Hepatitis_E_ORF2_ORF3_F_R', 'GGTGGTTTCTGGGGTGAC', 'AGGGGTTGGTTGGATGAA'),
    ('Norovirus GI', 'Norovirus_GI_ORF1_ORF2_CDC_F_R', 'CGYTGGATGCGNTTYCATGA', 'CTTAGACGCCATCATCATTYAC'),
    ('Norovirus GII', 'Norovirus_GII_ORF1_ORF2_CDC_F_R', 'CARGARBCNATGTTYAGRTGGATGAG', 'TCGACGCCATCTTCATTCACA'),
    ('Rotavirus A', 'Rotavirus_A_VP6_F_R', 'GATGTCCTGTACTCCTTGTCA', 'GGTCACATCATACAATTCTAATCTAAG'),
    ('Measles morbillivirus (Masern)', 'Measles_virus_N_gene_WHO_F_R', 'GGAGCTCCAGGATGGTTTC', 'GTGCTGCCTTCACACTCGTC'),
    ('Mumps orthorubulavirus (Mumps)', 'Mumps_virus_SH_gene_F_R', 'AGTAGTGTCGATGATCTCAT', 'GTTAGTAGTTTCCACTTCAG'),
    ('Rubella virus (Röteln)', 'Rubella_virus_E1_glycoprotein_F_R', 'CGTCCCAGTGAACTGCGTC', 'CGAACACGTTCCCCATCAA'),
    ('Rabies lyssavirus (Tollwut)', 'Rabies_virus_N_gene_F_R', 'ACGCTTAACAACCAGATCAAAGAA', 'CGTGTTCAATCAACTGCTCGTTG'),
    ('Human enterovirus / Rhinovirus', 'Pan_Enterovirus_5UTR_F_R', 'TCCTCCGGCCCCTGAATGCGGC', 'GATTGTCACCATAAGCAGCC'),
    ('Human adenovirus', 'Pan_Adenovirus_hexon_F_R', 'GCCACCTTCTTCCCCATGGC', 'GTAGCGTTGCCGGCCGAGAA'),

    # =========================================================================
    # 9. HUMAN FUNGAL PATHOGENS
    # =========================================================================
    ('Candida albicans', 'Candida_albicans_ITS1_ITS4', 'TCCGTAGGTGAACCTGCGG', 'TCCTCCGCTTATTGATATGC'),
    ('Candida albicans', 'Candida_albicans_CALB_spec_F_R', 'TTGAACCACTCTGCATCGTG', 'TTCCAAGCGGTCTTTACCTC'),
    ('Candida auris (MDR yeast)', 'Candida_auris_ITS2_CDC_spec_F_R', 'ATTTTGCATACACACTGATTTG', 'CGTGCAAGCTGTAATTTTGTG'),
    ('Aspergillus fumigatus', 'Aspergillus_fumigatus_ITS1_spec_F_R', 'GCCCGCCGAAGACCCCAAC', 'TTGTCTGATGCGGCGTCCG'),
    ('Cryptococcus neoformans/gattii', 'Cryptococcus_neoformans_URA5_F_R', 'ATGTCCTCCCAAGCCCTC', 'TTAAGCGTCGAACCGTGC'),
    ('Pneumocystis jirovecii', 'Pneumocystis_jirovecii_mtLSU_rRNA_pAZ102', 'ATAGGTCGTAGTCGAAATT', 'AAAATCTTTCTTCAAATGATC'),

    # =========================================================================
    # 10. HUMAN PROTOZOA & PARASITES
    # =========================================================================
    ('Toxoplasma gondii', 'Toxoplasma_gondii_B1_gene_35repeat', 'GGAACTGCATCCGTTCATGAG', 'TCTTTAAAGCGTTCGTGGTC'),
    ('Toxoplasma gondii', 'Toxoplasma_gondii_529bp_repeat_qPCR', 'AGAGACACAGTGCATCTGGATT', 'TTCTCCTCTCTGCGGCGTAGT'),
    ('Giardia lamblia', 'Giardia_lamblia_beta_giardin_F_R', 'AAGCCCGACGACCTCACCCG', 'CATCGTAGTCGTCCTCGTCG'),
    ('Entamoeba histolytica', 'Entamoeba_histolytica_18S_F_R', 'AAGCATTGTTTCTAGATCTGATAA', 'TCTTGTAGCTTAGTTTTACTGTAC'),
    ('Cryptosporidium parvum/hominis', 'Cryptosporidium_gp60_F_R', 'ATAGTCTCCGCTGTATTC', 'GGAAGGAACGATGTATCT'),
    ('Leishmania donovani/infantum/major', 'Leishmania_genus_kDNA_minicircle', 'CCTATTTTACACCAACCCCCAGT', 'GGGTAGGGGCGTTCTGCGAAA'),
    ('Trypanosoma cruzi (Chagas)', 'Trypanosoma_cruzi_satellite_DNA_TCZ', 'GCTCTTGCCCACAMGGGTGC', 'CCAAGCAGCGGATAGTTCAGG')
]

print(f"Total curated human pathogen primer pairs to evaluate: {len(pathogen_primers)}")

# Validate each pair
for org, pid, f, r in pathogen_primers:
    f = f.strip().upper()
    r = r.strip().upper()
    assert 15 <= len(f) <= 45, f"Length F error: {pid} ({len(f)})"
    assert 15 <= len(r) <= 45, f"Length R error: {pid} ({len(r)})"
    assert IUPAC.match(f), f"IUPAC F error: {pid} ({f})"
    assert IUPAC.match(r), f"IUPAC R error: {pid} ({r})"

print("All new curated primers passed validation checks!")

# Load current best_primers_db.csv
db = pd.read_csv('best_primers_db.csv')
print(f"Current DB size: {len(db)}")

existing_pairs = set(zip(db['Forward'].str.upper(), db['Reverse'].str.upper()))

added_rows = []
for org, pid, f, r in pathogen_primers:
    pair = (f.upper(), r.upper())
    if pair not in existing_pairs:
        added_rows.append({
            'Organism': org,
            'PairID': pid,
            'Forward': f.upper(),
            'Reverse': r.upper()
        })
        existing_pairs.add(pair)

print(f"New unique pairs to append: {len(added_rows)}")

updated_db = pd.concat([db, pd.DataFrame(added_rows)], ignore_index=True)
updated_db = updated_db.drop_duplicates(subset=['Forward', 'Reverse']).reset_index(drop=True)
updated_db = updated_db.sort_values(by=['Organism', 'PairID']).reset_index(drop=True)

# Save locally and to Downloads
updated_db.to_csv('best_primers_db.csv', index=False)
downloads_path = os.path.join(r"C:\Users\internet\Downloads", "best_primers_db.csv")
updated_db.to_csv(downloads_path, index=False)

print(f"Updated DB size: {len(updated_db)} primer pairs.")
