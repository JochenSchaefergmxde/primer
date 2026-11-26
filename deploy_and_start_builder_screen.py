import paramiko
import os

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

print(f"[*] Connecting to {username}@{host} to deploy and launch build_all_primers_db.py in screen...")

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=15)

sftp = ssh.open_sftp()

local_builder = r"C:\Users\internet\Downloads\repository\primer\build_all_primers_db.py"
remote_builder = "/mnt/sdc1/scripts/build_all_primers_db.py"
print(f"[*] SFTP uploading {local_builder}...")
sftp.put(local_builder, remote_builder)

local_orgs = r"C:\Users\internet\Downloads\repository\primer\organisms.json"
remote_orgs = "/mnt/sdc1/scripts/organisms.json"
print(f"[*] SFTP uploading {local_orgs}...")
sftp.put(local_orgs, remote_orgs)

sftp.close()

# Start screen session
cmd = 'cd /mnt/sdc1/scripts && screen -dmS build_db bash -c "python3 /mnt/sdc1/scripts/build_all_primers_db.py > /mnt/sdc1/scripts/build_db.log 2>&1"'
ssh.exec_command(cmd)

print("[+] Successfully started 'build_db' screen session on server!")
ssh.close()
