import paramiko

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
print(f"[*] SFTP uploading updated {local_script} -> {remote_script}...")
sftp.put(local_script, remote_script)
sftp.close()

ssh.exec_command("chmod +x /mnt/sdc1/scripts/primer_suche.py")

fastq_file = "/mnt/fastq/js_16S_Long_Read_pathogen.fastq.gz"
out_dir = "/mnt/sdc1/scripts/extracted_amplicons_16S_LongRead"
remote_dash = "/mnt/sdc1/scripts/primer_dashboard_16S_LongRead.html"
local_dash = r"C:\Users\internet\Downloads\primer_dashboard_16S_LongRead.html"

cmd = f"python3 /mnt/sdc1/scripts/primer_suche.py --fastq {fastq_file} --primers /mnt/sdc1/scripts/best_primers_db.csv --output-dir {out_dir} --dashboard {remote_dash}"

print(f"[*] Running analysis for {fastq_file}...")
stdin, stdout, stderr = ssh.exec_command(cmd)

out = stdout.read().decode('utf-8', errors='ignore')
err = stderr.read().decode('utf-8', errors='ignore')

print("[STDOUT]:\n", out)
if err:
    print("[STDERR]:\n", err)

sftp = ssh.open_sftp()
try:
    sftp.get(remote_dash, local_dash)
    print(f"[+] Downloaded Dashboard to: {local_dash}")
except Exception as e:
    print(f"[!] Dashboard download error: {e}")

sftp.close()
ssh.close()
