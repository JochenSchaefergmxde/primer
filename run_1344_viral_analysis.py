import paramiko
import os

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

print(f"[*] Connecting to {username}@{host} to run analysis with 1,344 primers (including viruses)...")

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=15)

sftp = ssh.open_sftp()
local_csv = r"C:\Users\internet\Downloads\best_primers_db.csv"
remote_csv = "/mnt/sdc1/scripts/best_primers_db.csv"
print(f"[*] Verifying SFTP upload of 1,344 primers database ({os.path.getsize(local_csv)} bytes)...")
sftp.put(local_csv, remote_csv)

fastq_tsaa = "/mnt/fastq/js_u_t2t_2022_11_TSAA1868_100.fastq.gz"
out_dir_tsaa = "/mnt/sdc1/scripts/extracted_amplicons"
remote_dash_tsaa = "/mnt/sdc1/scripts/primer_dashboard_TSAA1868.html"
local_dash_tsaa = r"C:\Users\internet\Downloads\primer_dashboard_TSAA1868.html"

cmd_tsaa = f"python3 /mnt/sdc1/scripts/primer_suche.py --fastq {fastq_tsaa} --primers {remote_csv} --output-dir {out_dir_tsaa} --dashboard {remote_dash_tsaa}"
print(f"[*] Running TSAA1868 analysis with 1,344 primers...")
stdin, stdout, stderr = ssh.exec_command(cmd_tsaa)

for line in stdout:
    print(line, end="")

err = stderr.read().decode('utf-8', errors='ignore')
if err:
    print("[STDERR]:", err)

try:
    sftp.get(remote_dash_tsaa, local_dash_tsaa)
    print(f"[+] Downloaded updated Dashboard to: {local_dash_tsaa}")
except Exception as e:
    print(f"[!] Dashboard download error: {e}")

sftp.close()
ssh.close()
