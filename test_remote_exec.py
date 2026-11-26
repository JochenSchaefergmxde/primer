import paramiko

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=10)

cmd = "python3 /mnt/sdc1/scripts/primer_suche.py --fastq /mnt/fastq/js_u_t2t_2022_11_TSAA1868_100.fastq.gz --primers /mnt/sdc1/scripts/best_primers_db.csv --output-dir /mnt/sdc1/scripts/extracted_amplicons --dashboard /mnt/sdc1/scripts/primer_dashboard_TSAA1868.html --max-reads 50000"

print(f"[*] Running command directly to test output: {cmd}")
stdin, stdout, stderr = ssh.exec_command(cmd)

out = stdout.read().decode('utf-8', errors='ignore')
err = stderr.read().decode('utf-8', errors='ignore')

print("[STDOUT]:")
print(out)
print("[STDERR]:")
print(err)

ssh.close()
