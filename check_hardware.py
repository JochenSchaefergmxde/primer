import paramiko

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect('192.168.178.90', username='internet', password='fever123!')

stdin, stdout, stderr = ssh.exec_command('lscpu')
print("=== CPU INFO ===")
for line in stdout:
    if any(k in line for k in ['Model name', 'CPU(s):', 'Thread(s) per core', 'Core(s) per socket', 'Socket(s)']):
        print(line.strip())

stdin, stdout, stderr = ssh.exec_command('nvidia-smi')
gpu_out = stdout.read().decode()
if gpu_out:
    print("\n=== GPU INFO ===")
    print(gpu_out[:500])
else:
    print("\nNo NVIDIA GPU or nvidia-smi not available.")

stdin, stdout, stderr = ssh.exec_command('free -h')
print("\n=== MEMORY INFO ===")
print(stdout.read().decode())

stdin, stdout, stderr = ssh.exec_command('which pigz ripgrep seqkit bwa minimap2')
print("=== INSTALLED BIO TOOLS ===")
print(stdout.read().decode())

ssh.close()
