import paramiko, os, glob

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect('192.168.178.90', port=22, username='internet', password='fever123!', timeout=20)
sftp = ssh.open_sftp()

local_ref_dir = r"C:\Users\internet\Downloads\repository\primer\references"
remote_ref_dir = "/mnt/sdc1/pathogen_references"

patterns = ["*Theileria*", "*Cytauxzoon*", "*Babesia*", "*Piroplasmida*"]
uploaded = 0
for pat in patterns:
    for f in glob.glob(os.path.join(local_ref_dir, pat)):
        fname = os.path.basename(f)
        remote_path = f"{remote_ref_dir}/{fname}"
        sftp.put(f, remote_path)
        uploaded += 1
        print(f"[+] Hochgeladen: {fname} -> {remote_path}")

sftp.close()
print(f"Insgesamt {uploaded} Referenzdateien hochgeladen!")

stdin, stdout, stderr = ssh.exec_command("ls -la /mnt/sdc1/pathogen_references/*Theileria* /mnt/sdc1/pathogen_references/*Cytauxzoon* /mnt/sdc1/pathogen_references/*Piroplasmida*")
print("\n--- Server Listing ---")
print(stdout.read().decode())

ssh.close()
