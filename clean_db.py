import pandas as pd
import re
import os
import shutil

df = pd.read_csv('best_primers_db.csv')
initial_len = len(df)
print(f"Initial length: {initial_len}")

# 1. Remove non-pathogen / control / off-target contamination from PMC tables:
drop_patterns = [
    'PMC13630449', # Mpox, SARS2, Candida albicans in Mycobacterium paper
    'PMC13414585', # Tick host cox genes (Hlong, Ixodes)
    'PMC13516315_Borrelia_dsb_beta-tubulin',
    'beta-tubulin',
    'Bactin',
    'Borrelia_18s_NS2A',
    'Borrelia_18s_S_segment',
    'Borrelia_18s_VP12',
    'Borrelia_dsb_S_segment',
    'Borrelia_ospA_NS5'
]

mask_drop = df['PairID'].str.contains('|'.join(drop_patterns), case=False, na=False)
dropped = df[mask_drop]
print(f"Dropping {len(dropped)} off-target / host / virus contaminant rows:")
for idx, row in dropped.iterrows():
    print(f"  - [{row['Organism']}] {row['PairID']}")

df = df[~mask_drop].copy()

# 2. Correct misassigned organisms:
# Brucella primers misassigned to Bartonella
brucella_mask = df['PairID'].str.contains('BMEII|BruAb2|IS711|IS1111|BR0952|TWSP', case=False, na=False)
df.loc[brucella_mask, 'Organism'] = 'Brucella'

# Mycoplasma primers misassigned to Mycobacterium
mycoplasma_mask = df['PairID'].str.contains('VlhA|mgc2', case=False, na=False)
df.loc[mycoplasma_mask, 'Organism'] = 'Mycoplasma'

# Anaplasma msp2 misassigned to Borrelia
anaplasma_mask = df['PairID'].str.contains('msp2', case=False, na=False)
df.loc[anaplasma_mask, 'Organism'] = 'Anaplasma phagocytophilum'

# Ehrlichia dsb misassigned to Borrelia
ehrlichia_mask = df['PairID'].str.contains('dsb_dsb', case=False, na=False)
df.loc[ehrlichia_mask, 'Organism'] = 'Ehrlichia'

# Babesia / Theileria piroplasmid misassigned to Borrelia
piroplasmid_mask = df['PairID'].str.contains('Piroplasmid|18s_rRNA_GGGACTTTGCGTTCATAAAACGC', case=False, na=False)
df.loc[piroplasmid_mask, 'Organism'] = 'Babesia'

# Rickettsia gltA & ompB misassigned to Borrelia
rick_mask = df['PairID'].str.contains('ompB_CAAATGTTGCAGTTCCTCTAAA|gltA_GAGAGAAAATTATATATCCAAATGTTGAT', case=False, na=False)
df.loc[rick_mask, 'Organism'] = 'Rickettsia'

# Bartonella pap13 misassigned to Borrelia
bartonella_mask = df['PairID'].str.contains('pap13', case=False, na=False)
df.loc[bartonella_mask, 'Organism'] = 'Bartonella'

# 3. Add legitimate tick-borne viruses (e.g. FSME / TBEV)
tbev_primers = [
    {'Organism': 'Tick-borne encephalitis virus (FSME / TBEV)', 'PairID': 'TBEV_FSME_3UTR_qPCR_Schwaiger', 'Forward': 'GGGCGGTTCTTGTTCTCC', 'Reverse': 'ACACATCACCTCCTTGTCAG'},
    {'Organism': 'Tick-borne encephalitis virus (FSME / TBEV)', 'PairID': 'TBEV_FSME_NS5_Gaeumann', 'Forward': 'AATGGAAGATGAAATGAAC', 'Reverse': 'CCGAGCCAAAGTGAGGATGT'}
]
df = pd.concat([df, pd.DataFrame(tbev_primers)], ignore_index=True)

# Deduplicate
df = df.drop_duplicates(subset=['Forward', 'Reverse']).reset_index(drop=True)
df = df.sort_values(by=['Organism', 'PairID']).reset_index(drop=True)

print(f"Cleaned length: {len(df)}")
df.to_csv('best_primers_db.csv', index=False)
downloads_path = os.path.join(r"C:\Users\internet\Downloads", "best_primers_db.csv")
df.to_csv(downloads_path, index=False)
print("Saved to best_primers_db.csv and Downloads!")
