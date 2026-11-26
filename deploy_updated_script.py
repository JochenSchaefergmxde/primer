import paramiko

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

print(f"[*] Uploading updated primer_suche.py (with CD-HIT-EST support) to {username}@{host}...")

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=15)

sftp = ssh.open_sftp()
local_script = r"C:\Users\internet\Downloads\repository\primer\primer_suche.py"
remote_script = "/mnt/sdc1/scripts/primer_suche.py"
sftp.put(local_script, remote_script)
sftp.close()

ssh.exec_command("chmod +x /mnt/sdc1/scripts/primer_suche.py")
print("[+] Successfully uploaded updated primer_suche.py to /mnt/sdc1/scripts/primer_suche.py!")

ssh.close()
