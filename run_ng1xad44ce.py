import paramiko
import os

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=10)

fastq_file = "/mnt/fastq/js_2022_05_NG1XAD44CE_t2t_u_s.fastq.gz"
out_dir = "/mnt/sdc1/scripts/extracted_amplicons_NG1XAD44CE"
remote_dash = "/mnt/sdc1/scripts/primer_dashboard_NG1XAD44CE.html"
local_dash = r"C:\Users\internet\Downloads\primer_dashboard_NG1XAD44CE.html"

cmd = f"python3 /mnt/sdc1/scripts/primer_suche.py --fastq {fastq_file} --primers /mnt/sdc1/scripts/best_primers_db.csv --output-dir {out_dir} --dashboard {remote_dash}"

print(f"[*] Executing analysis for {fastq_file}...")
stdin, stdout, stderr = ssh.exec_command(cmd)

out = stdout.read().decode('utf-8', errors='ignore')
err = stderr.read().decode('utf-8', errors='ignore')

print("[STDOUT]:\n", out)
if err:
    print("[STDERR]:\n", err)

# Download dashboard if created
sftp = ssh.open_sftp()
try:
    sftp.get(remote_dash, local_dash)
    print(f"[+] Downloaded Dashboard to: {local_dash}")
except Exception as e:
    print(f"[!] Dashboard download error: {e}")

sftp.close()
ssh.close()
