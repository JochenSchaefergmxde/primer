import paramiko
import time

host = "192.168.178.90"
user = "internet"
password = "fever123!"

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, username=user, password=password)

for _ in range(5):
    stdin, stdout, stderr = ssh.exec_command("ls -lh /mnt/sdc1/scripts/all_primers_db.csv && ps aux | grep build_all_primers_db")
    print(stdout.read().decode())
    time.sleep(3)

ssh.close()
