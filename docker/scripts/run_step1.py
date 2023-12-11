#!/opt/conda/bin/python
import os, sys, subprocess
from multiprocessing import Pool

def run_cmd(sample,cmd):
    status=subprocess.run(cmd,shell=True).returncode
    if status==0:
        return sample+' succeed'
    else:
        return sample+' fail'

def batch_run(cmds,process_num):
    log=[]
    pool = Pool(processes=process_num)
    result = []
    for cmd in cmds:
        sample=cmd.split()[2]
        process = pool.apply_async(run_cmd,(sample,cmd,))
        result.append(process)
    for process in result:
        log.append(process.get())
    return log



dataset = sys.argv[1]
ncore = int(sys.argv[2])

if not os.path.exists('step1/'):
    os.system('mkdir step1')
if not os.path.exists('step1/' + dataset):
    os.system('mkdir step1/' + dataset)

fp = open(dataset + '_struc.list', 'r')
cases = []
for line in fp:
    words = line.split()
    accession = words[0]
    cases.append(accession)
#    cases.append([accession, version])
fp.close()

cmds = []
for case in cases:
    if os.path.exists('step1/' + dataset + '/' + case + '.fa'):
        fp = open('step1/' + dataset + '/' + case + '.fa', 'r')
        check_header = 0
        check_seq = 0
        check_length = 0
        for line in fp:
            check_length += 1
            if line[0] == '>':
                if line[1:-1] == case:
                    check_header = 1
            else:
                if len(line) > 10:
                    check_seq = 1
        fp.close()
        if check_header and check_seq and check_length == 2:
            pass
        else:
            os.system('rm step1/' + dataset + '/' + case + '.fa')
            cmds.append('python /opt/DPAM/scripts/step1_get_AFDB_seqs.py ' + dataset + ' ' + case)
    else:
        cmds.append('python /opt/DPAM/scripts/step1_get_AFDB_seqs.py ' + dataset + ' ' + case)


if cmds:
    logs = batch_run(cmds,ncore)
    fail = [i for i in logs if 'fail' in i]
    if fail:
        with open(dataset + '_step1.log','w') as f:
            for i in fail:
                f.write(i+'\n')
    else:
        with open(dataset + '_step1.log','w') as f:
            f.write('done\n')
else:
    with open(dataset + '_step1.log','w') as f:
        f.write('done\n')
