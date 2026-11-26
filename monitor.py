import paramiko
import time

host = "192.168.178.90"
user = "internet"
password = "fever123!"

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, username=user, password=password)

running = True
while running:
    stdin, stdout, stderr = ssh.exec_command("ps aux | grep build_all_primers_db | grep -v grep | wc -l")
    count = int(stdout.read().decode().strip())
    
    stdin, stdout, stderr = ssh.exec_command("ls -lh /mnt/sdc1/scripts/all_primers_db.csv")
    print(stdout.read().decode().strip())
    
    if count == 0:
        print("Process finished!")
        break
        
    time.sleep(5)

ssh.close()
