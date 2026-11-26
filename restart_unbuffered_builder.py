import paramiko
import shutil

host = "192.168.178.90"
port = 22
username = "internet"
password = "fever123!"

# Copy locally to Downloads per rule
local_script = r"C:\Users\internet\Downloads\repository\primer\build_all_primers_db.py"
shutil.copy(local_script, r"C:\Users\internet\Downloads\build_all_primers_db.py")

print(f"[*] Uploading updated unbuffered script to {username}@{host}...")

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, port=port, username=username, password=password, timeout=15)

sftp = ssh.open_sftp()
sftp.put(local_script, "/mnt/sdc1/scripts/build_all_primers_db.py")
sftp.close()

# Kill previous build_db screen session
ssh.exec_command("screen -X -S build_db quit; pkill -f build_all_primers_db.py")

# Launch new unbuffered python3 -u session
cmd = 'cd /mnt/sdc1/scripts && screen -dmS build_db bash -c "python3 -u /mnt/sdc1/scripts/build_all_primers_db.py > /mnt/sdc1/scripts/build_db.log 2>&1"'
ssh.exec_command(cmd)

print("[+] Successfully restarted unbuffered build_db screen session!")
ssh.close()
