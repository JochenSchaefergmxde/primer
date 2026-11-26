import paramiko

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect('192.168.178.90', username='internet', password='fever123!', timeout=15)

test_code = """
import subprocess, time

q_fa = '/tmp/test_blast_q.fa'
# Let's test with a 300 bp 18S sequence from Cytauxzoon felis
with open(q_fa, 'w') as f:
    f.write('>test_query\\nACTTTCGTTCCTTTAATTGGTTTCAGTTTAAGTCATACATGTATCAGTCTCAGATTCCCCCAACTGTTCCCTATTAACACATTCATATTTTGCCTTGTTGAAGTCATTGCTTGCCATTACTTTAGCATGGAATAACTAAGTAGGACTTTGGTTCTATTTTGTTGGTTTAAGAGCCAAAGTAATGATATATAGGACAAGTTGGGGGCATTCGTATTTAAATTGTCAGAGGTGAAATTCTTAGATTTGTTAAAGACGAACTACTGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAAGAACGA\\n')

# Run blastn against /mnt/blast_db/nt with outfmt 6 (tabular: qseqid sseqid pident length evalue bitscore stitle)
cmd = "blastn -db /mnt/blast_db/nt -query " + q_fa + " -max_target_seqs 3 -max_hsps 1 -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' -num_threads 4"
t0 = time.time()
res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
elapsed = time.time() - t0
print(f"Blastn finished in {elapsed:.2f}s, exit code: {res.returncode}")
print("STDOUT:")
print(res.stdout)
print("STDERR:")
print(res.stderr)
"""

stdin, stdout, stderr = ssh.exec_command('python3 -')
stdin.write(test_code)
stdin.channel.shutdown_write()

print(stdout.read().decode('utf-8'))
print(stderr.read().decode('utf-8'))
ssh.close()
