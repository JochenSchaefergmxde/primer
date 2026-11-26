import paramiko

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

print(f"[*] Connecting to {username}@{host} to start primer search in screen session...")

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=15)

cmd = 'screen -dmS primer_run bash -c "python3 /mnt/sdc1/scripts/primer_suche.py --fastq /mnt/fastq/js_u_t2t_2022_11_TSAA1868_100.fastq.gz --primers /mnt/sdc1/scripts/best_primers_db.csv --output-dir /mnt/sdc1/scripts/extracted_amplicons --dashboard /mnt/sdc1/scripts/primer_dashboard_TSAA1868.html > /mnt/sdc1/scripts/run.log 2>&1"'
ssh.exec_command(cmd)

print("[+] Successfully launched screen session 'primer_run' on server!")
ssh.close()
