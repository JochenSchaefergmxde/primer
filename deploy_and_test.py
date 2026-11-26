import paramiko
import os

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=10)

sftp = ssh.open_sftp()
local_script = r"C:\Users\internet\Downloads\repository\primer\primer_suche.py"
remote_script = "/mnt/sdc1/scripts/primer_suche.py"
print(f"[*] SFTP uploading fixed {local_script} -> {remote_script}...")
sftp.put(local_script, remote_script)
sftp.close()

ssh.exec_command("chmod +x /mnt/sdc1/scripts/primer_suche.py")

cmd = "python3 /mnt/sdc1/scripts/primer_suche.py --fastq /mnt/fastq/js_u_t2t_2022_11_TSAA1868_100.fastq.gz --primers /mnt/sdc1/scripts/best_primers_db.csv --output-dir /mnt/sdc1/scripts/extracted_amplicons --dashboard /mnt/sdc1/scripts/primer_dashboard_TSAA1868.html"

print(f"[*] Running analysis on server: {cmd}")
stdin, stdout, stderr = ssh.exec_command(cmd)

out = stdout.read().decode('utf-8', errors='ignore')
err = stderr.read().decode('utf-8', errors='ignore')

print("[STDOUT]:\n", out)
if err:
    print("[STDERR]:\n", err)

# Launch background screen for full run if needed
screen_cmd = f"screen -dmS primer_run bash -c '{cmd}'"
ssh.exec_command(screen_cmd)

# Download dashboard if generated
sftp = ssh.open_sftp()
remote_dash = "/mnt/sdc1/scripts/primer_dashboard_TSAA1868.html"
local_dash = r"C:\Users\internet\Downloads\primer_dashboard_TSAA1868.html"

try:
    sftp.get(remote_dash, local_dash)
    print(f"[+] Downloaded Dashboard to: {local_dash}")
except Exception as e:
    print(f"[!] Could not download dashboard yet: {e}")

sftp.close()
ssh.close()
