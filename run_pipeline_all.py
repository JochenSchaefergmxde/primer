import paramiko
import os
import sys

sys.stdout.reconfigure(encoding='utf-8')

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

# 1. Generate all_primers_db.csv if not present or small
all_csv_path = r"C:\Users\internet\Downloads\all_primers_db.csv"
best_csv_path = r"C:\Users\internet\Downloads\best_primers_db.csv"

if not os.path.exists(all_csv_path) or os.path.getsize(all_csv_path) < 1000:
    print("[*] Generating all_primers_db.csv using build_all_primers_db.py...")
    import build_all_primers_db
    build_all_primers_db.main()

target_csv = all_csv_path if (os.path.exists(all_csv_path) and os.path.getsize(all_csv_path) > 1000) else best_csv_path

print(f"[*] Connecting to {username}@{host}...")
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=15)

sftp = ssh.open_sftp()

remote_csv = "/mnt/sdc1/scripts/all_primers_db.csv"
print(f"[*] SFTP uploading database ({os.path.getsize(target_csv):,} bytes) to {remote_csv}...")
sftp.put(target_csv, remote_csv)

local_script = r"C:\Users\internet\Downloads\repository\primer\primer_suche.py"
remote_script = "/mnt/sdc1/scripts/primer_suche.py"
sftp.put(local_script, remote_script)
sftp.close()

ssh.exec_command("chmod +x /mnt/sdc1/scripts/primer_suche.py")

# 2. Run FASTQ analysis
fastq_tsaa = "/mnt/fastq/js_u_t2t_2022_11_TSAA1868_100.fastq.gz"
out_dir_tsaa = "/mnt/sdc1/scripts/extracted_amplicons"
remote_dash_tsaa = "/mnt/sdc1/scripts/primer_dashboard_TSAA1868.html"
local_dash_tsaa = r"C:\Users\internet\Downloads\primer_dashboard_TSAA1868.html"

cmd = f"python3 /mnt/sdc1/scripts/primer_suche.py --fastq {fastq_tsaa} --primers {remote_csv} --output-dir {out_dir_tsaa} --dashboard {remote_dash_tsaa}"
print(f"[*] Running primer search on {fastq_tsaa}...")

stdin, stdout, stderr = ssh.exec_command(cmd)

for line in stdout:
    print(line, end="")

err = stderr.read().decode('utf-8', errors='ignore')
if err:
    print("[STDERR]:", err)

# Download resulting dashboard
sftp = ssh.open_sftp()
try:
    sftp.get(remote_dash_tsaa, local_dash_tsaa)
    print(f"[+] Downloaded Dashboard to: {local_dash_tsaa}")
except Exception as e:
    print(f"[!] Dashboard download error: {e}")

sftp.close()
ssh.close()
