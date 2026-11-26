import paramiko
import os

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

print(f"[*] Uploading expanded 1,344 primer database (including viruses) to {username}@{host}...")

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=15)

sftp = ssh.open_sftp()
local_csv = r"C:\Users\internet\Downloads\best_primers_db.csv"
remote_csv = "/mnt/sdc1/scripts/best_primers_db.csv"
sftp.put(local_csv, remote_csv)
sftp.close()

print(f"[+] Successfully uploaded updated database ({os.path.getsize(local_csv)} bytes, 1,344 primer pairs) to /mnt/sdc1/scripts/best_primers_db.csv!")

ssh.close()
