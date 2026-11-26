import paramiko
import os

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

print(f"[*] Connecting to {username}@{host} to upload full 1,280 primer database...")

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=15)

sftp = ssh.open_sftp()
local_csv = r"C:\Users\internet\Downloads\best_primers_db.csv"
remote_csv = "/mnt/sdc1/scripts/best_primers_db.csv"
print(f"[*] SFTP uploading full database ({os.path.getsize(local_csv)} bytes)...")
sftp.put(local_csv, remote_csv)
sftp.close()

# 1. Run 16S Long Read Fastq Analysis with 1,280 Primers
fastq_16s = "/mnt/fastq/js_16S_Long_Read_pathogen.fastq.gz"
out_dir_16s = "/mnt/sdc1/scripts/extracted_amplicons_16S_LongRead"
remote_dash_16s = "/mnt/sdc1/scripts/primer_dashboard_16S_LongRead.html"
local_dash_16s = r"C:\Users\internet\Downloads\primer_dashboard_16S_LongRead.html"

cmd_16s = f"python3 /mnt/sdc1/scripts/primer_suche.py --fastq {fastq_16s} --primers {remote_csv} --output-dir {out_dir_16s} --dashboard {remote_dash_16s}"
print(f"[*] Running 16S Long Read analysis with 1,280 primers...")
stdin, stdout, stderr = ssh.exec_command(cmd_16s)
out1 = stdout.read().decode('utf-8', errors='ignore')
err1 = stderr.read().decode('utf-8', errors='ignore')
print("[16S STDOUT]:\n", out1)
if err1:
    print("[16S STDERR]:\n", err1)

# Download 16S Dashboard
sftp = ssh.open_sftp()
try:
    sftp.get(remote_dash_16s, local_dash_16s)
    print(f"[+] Downloaded 16S Dashboard to: {local_dash_16s}")
except Exception as e:
    print(f"[!] 16S Dashboard download error: {e}")

# 2. Run TSAA1868 Fastq Analysis with 1,280 Primers
fastq_tsaa = "/mnt/fastq/js_u_t2t_2022_11_TSAA1868_100.fastq.gz"
out_dir_tsaa = "/mnt/sdc1/scripts/extracted_amplicons"
remote_dash_tsaa = "/mnt/sdc1/scripts/primer_dashboard_TSAA1868.html"
local_dash_tsaa = r"C:\Users\internet\Downloads\primer_dashboard_TSAA1868.html"

cmd_tsaa = f"python3 /mnt/sdc1/scripts/primer_suche.py --fastq {fastq_tsaa} --primers {remote_csv} --output-dir {out_dir_tsaa} --dashboard {remote_dash_tsaa}"
print(f"[*] Running TSAA1868 analysis with 1,280 primers...")
stdin, stdout, stderr = ssh.exec_command(cmd_tsaa)
out2 = stdout.read().decode('utf-8', errors='ignore')
err2 = stderr.read().decode('utf-8', errors='ignore')
print("[TSAA1868 STDOUT]:\n", out2)

try:
    sftp.get(remote_dash_tsaa, local_dash_tsaa)
    print(f"[+] Downloaded TSAA1868 Dashboard to: {local_dash_tsaa}")
except Exception as e:
    print(f"[!] TSAA1868 Dashboard download error: {e}")

sftp.close()
ssh.close()
