import paramiko
import time
import os

host = "192.168.178.90"
user = "internet"
password = "fever123!"

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, username=user, password=password)

print("Monitoring primer_suche.py on server...")

while True:
    stdin, stdout, stderr = ssh.exec_command('ps aux | grep "[p]rimer_suche.py --fastq /mnt/fastq/tm_NG1YETBR8Q"')
    out = stdout.read().decode().strip()
    
    if not out:
        print("Process finished! Attempting to download the dashboard...")
        break
    else:
        time.sleep(30)

# Download
sftp = ssh.open_sftp()
remote_file = "/mnt/sdc1/scripts/primer_dashboard_tm_NG1YETBR8Q.html"
local_file = r"C:\Users\internet\Downloads\primer_dashboard_tm_NG1YETBR8Q.html"
try:
    sftp.get(remote_file, local_file)
    print(f"File successfully downloaded to {local_file}")
except Exception as e:
    print(f"Error downloading file: {e}")
finally:
    sftp.close()

ssh.close()
