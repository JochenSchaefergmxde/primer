import paramiko
import os
import sys

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

print(f"[*] Connecting to {username}@{host}:{port}...")

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

try:
    ssh.connect(host, port=port, username=username, password=password, timeout=10)
    print("[+] SSH connection established successfully!")

    # 1. Create target directory
    print("[*] Creating /mnt/sdc1/scripts directory if missing...")
    stdin, stdout, stderr = ssh.exec_command("mkdir -p /mnt/sdc1/scripts /mnt/sdc1/scripts/extracted_amplicons")
    stdout.channel.recv_exit_status()

    # 2. Upload files via SFTP
    sftp = ssh.open_sftp()
    
    local_script = r"C:\Users\internet\Downloads\repository\primer\primer_suche.py"
    remote_script = "/mnt/sdc1/scripts/primer_suche.py"
    print(f"[*] SFTP uploading {local_script} -> {remote_script}...")
    sftp.put(local_script, remote_script)

    local_csv = r"C:\Users\internet\Downloads\best_primers_db.csv"
    remote_csv = "/mnt/sdc1/scripts/best_primers_db.csv"
    print(f"[*] SFTP uploading {local_csv} -> {remote_csv}...")
    sftp.put(local_csv, remote_csv)
    
    sftp.close()

    # 3. Make executable
    ssh.exec_command("chmod +x /mnt/sdc1/scripts/primer_suche.py")

    # 4. Check files on server
    stdin, stdout, stderr = ssh.exec_command("ls -la /mnt/sdc1/scripts/")
    print("[+] Remote /mnt/sdc1/scripts/ files:")
    print(stdout.read().decode('utf-8', errors='ignore'))

    # 5. Launch screen session
    screen_cmd = (
        "screen -dmS primer_run bash -c "
        "'cd /mnt/sdc1/scripts && python3 primer_suche.py "
        "--fastq /mnt/fastq/js_u_t2t_2022_11_TSAA1868_100.fastq.gz "
        "--primers /mnt/sdc1/scripts/best_primers_db.csv "
        "--output-dir /mnt/sdc1/scripts/extracted_amplicons "
        "--dashboard /mnt/sdc1/scripts/primer_dashboard_TSAA1868.html'"
    )
    print(f"[*] Launching screen session...")
    stdin, stdout, stderr = ssh.exec_command(screen_cmd)
    stdout.channel.recv_exit_status()

    # 6. Verify screen session
    stdin, stdout, stderr = ssh.exec_command("screen -ls")
    screen_output = stdout.read().decode('utf-8', errors='ignore')
    print("[+] Active screen sessions:")
    print(screen_output)

    ssh.close()
    print("[+] Deployment and job launch completed successfully!")

except Exception as e:
    print(f"[!] Error connecting or deploying to {host}: {e}")
    sys.exit(1)
