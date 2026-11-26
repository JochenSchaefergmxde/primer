import paramiko

host = "192.168.178.90"
user = "internet"
password = "fever123!"

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host, username=user, password=password)

ssh.exec_command('pkill -f SCREEN')
ssh.exec_command('pkill -f build_all_primers_db.py')

cmd = 'screen -dmS pipeline bash -c "cd /mnt/sdc1/scripts && rm -f all_primers_db.csv && python3 build_all_primers_db.py && python3 primer_suche.py -i js_u_t2t_2022_11_TSAA1868_100.fastq.gz --primers all_primers_db.csv --cdhit > pipeline.log 2>&1"'
ssh.exec_command(cmd)
ssh.close()
print("Pipeline started!")
