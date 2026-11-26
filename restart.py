import paramiko

host = "192.168.178.90"
user = "internet"
password = "fever123!"

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, username=user, password=password)

# Kill old processes
ssh.exec_command("pkill -f build_all_primers_db.py")
ssh.exec_command("pkill -f SCREEN")

# Upload new script
sftp = ssh.open_sftp()
local_path = r"C:\Users\internet\Downloads\repository\primer\build_all_primers_db.py"
remote_path = "/mnt/sdc1/scripts/build_all_primers_db.py"
sftp.put(local_path, remote_path)
sftp.close()

# Run it
stdin, stdout, stderr = ssh.exec_command("cd /mnt/sdc1/scripts && rm -f all_primers_db.csv && screen -dmS build_primers python3 build_all_primers_db.py")
ssh.close()

print("Restarted cleanly.")
