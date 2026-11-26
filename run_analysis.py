import paramiko
import os
import sys

host = "192.168.178.90"
user = "internet"
password = "fever123!"

try:
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host, username=user, password=password)

    # Note: Using `--cdhit` for performance as discussed previously
    cmd = "cd /mnt/sdc1/scripts && python3 primer_suche.py --fastq /mnt/fastq/tm_NG1YETBR8Q_t2t_u.fastq.gz --primers best_primers_db.csv --cdhit --dashboard primer_dashboard_tm_NG1YETBR8Q.html"
    print(f"Running command on server: {cmd}")
    
    stdin, stdout, stderr = ssh.exec_command(cmd)
    
    # Stream the output directly so we can see progress in the log
    for line in iter(stdout.readline, ""):
        print(line, end="")
    
    exit_status = stdout.channel.recv_exit_status()
    
    err = stderr.read().decode()
    if err:
        print("STDERR:")
        print(err)

    if exit_status == 0:
        print("Analysis finished successfully. Downloading the dashboard...")
        sftp = ssh.open_sftp()
        remote_file = "/mnt/sdc1/scripts/primer_dashboard_tm_NG1YETBR8Q.html"
        local_file = r"C:\Users\internet\Downloads\primer_dashboard_tm_NG1YETBR8Q.html"
        sftp.get(remote_file, local_file)
        sftp.close()
        print(f"File successfully downloaded to {local_file}")
    else:
        print(f"Analysis failed with exit code {exit_status}")
        
finally:
    ssh.close()
