import pandas as pd
import shutil
import os
import paramiko

original_path = r"C:\Users\internet\Downloads\best_primers_db.csv"
backup_path = r"C:\Users\internet\Downloads\best_primers_db_backup.csv"
harvested_path = r"C:\Users\internet\Downloads\harvested_pathogen_primers.csv"
repo_path = r"C:\Users\internet\Downloads\repository\primer\best_primers_db.csv"

# 1. Create backup
if not os.path.exists(backup_path):
    shutil.copy2(original_path, backup_path)
    print(f"Created backup at: {backup_path}")

df_orig = pd.read_csv(original_path)
print(f"Original entries in best_primers_db: {len(df_orig)}")

# 2. Fix specific known errors & deduplicate:
# Remove known bogus duplicate: Mycobacterium marinum 16S alt (which is actually M. fortuitum rpoB)
df_orig = df_orig[~(
    (df_orig['Organism'] == 'Mycobacterium marinum') & 
    (df_orig['PairID'] == 'Mycobacterium_marinum_16S_F_R_alt')
)]

# Remove redundant M. kansasii row that copied the genus-level Telenti Tb11/Tb12 primer
df_orig = df_orig[~(
    (df_orig['Organism'] == 'Mycobacterium kansasii') & 
    (df_orig['PairID'] == 'Mycobacterium_kansasii_hsp65_F_R')
)]

# Fix Bartonella alsatica mislabeled as ITS (it is actually gltA)
mask_ba = (df_orig['Organism'] == 'Bartonella alsatica') & (df_orig['PairID'] == 'Bartonella_alsatica_ITS_BaF_BaR')
df_orig.loc[mask_ba, 'PairID'] = 'Bartonella_alsatica_gltA_BaF_BaR'

# Remove taxonomic duplicate: Propionibacterium acnes (synonym of Cutibacterium acnes)
df_orig = df_orig[~(
    (df_orig['Organism'] == 'Propionibacterium acnes') & 
    (df_orig['PairID'] == 'Propionibacterium_acnes_16S_univ')
)]

# Deduplicate any remaining exact Forward & Reverse duplicates, keeping the first
initial_len = len(df_orig)
df_clean = df_orig.drop_duplicates(subset=['Forward', 'Reverse'], keep='first').copy()
print(f"Removed {initial_len - len(df_clean)} duplicate rows from original database.")

# 3. Merge harvested literature primers
df_harvested = pd.read_csv(harvested_path)
print(f"Loaded {len(df_harvested)} harvested literature primers.")

new_rows = []
existing_pairs = set(zip(df_clean['Forward'].str.upper(), df_clean['Reverse'].str.upper()))

added_count = 0
for _, row in df_harvested.iterrows():
    f_seq = str(row['Forward']).strip().upper()
    r_seq = str(row['Reverse']).strip().upper()
    
    if (f_seq, r_seq) not in existing_pairs:
        existing_pairs.add((f_seq, r_seq))
        
        # Build clean PairID
        gene = row.get('Target_Gene', 'PCR')
        pmcid = row.get('PMCID', '')
        pid = row.get('PairID', 'F_R')
        clean_pair_id = f"{row['Organism']}_{gene}_{pid}_{pmcid}".replace(" ", "_")
        
        new_rows.append({
            'Organism': row['Organism'],
            'PairID': clean_pair_id,
            'Forward': f_seq,
            'Reverse': r_seq
        })
        added_count += 1

if new_rows:
    df_new = pd.DataFrame(new_rows)
    df_combined = pd.concat([df_clean, df_new], ignore_index=True)
else:
    df_combined = df_clean

print(f"Added {added_count} newly harvested primer pairs!")
print(f"Final combined database size: {len(df_combined)} pairs.")

# Save locally to Downloads and repo
df_combined.to_csv(original_path, index=False)
df_combined.to_csv(repo_path, index=False)
print(f"Saved updated database to {original_path}")

# 4. Upload to remote server so pipeline & search scripts use it immediately
try:
    print("Uploading updated best_primers_db.csv to server (192.168.178.90)...")
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect("192.168.178.90", username="internet", password="fever123!")
    sftp = ssh.open_sftp()
    sftp.put(original_path, "/mnt/sdc1/scripts/best_primers_db.csv")
    sftp.close()
    ssh.close()
    print("Successfully updated /mnt/sdc1/scripts/best_primers_db.csv on server!")
except Exception as e:
    print(f"Server upload error: {e}")
