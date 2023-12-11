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
if not os.path.exists('step2'):
    os.system('mkdir step2/')
if not os.path.exists('step2/' + dataset):
    os.system('mkdir step2/' + dataset)

fp = open(dataset + '_struc.list', 'r')
cases = []
for line in fp:
    words = line.split()
    cases.append(words[0])
fp.close()

cmds = []
for case in cases:
    fasta_length = 0
    if os.path.exists('step1/' + dataset + '/' + case + '.fa'):
        fp = open('step1/' + dataset + '/' + case + '.fa', 'r')
        for line in fp:
            if line[0] != '>':
                fasta_length = len(line[:-1])
        fp.close()

    pdb_resids = set([])
    if os.path.exists('step2/' + dataset + '/' + case + '.pdb'):
        fp = open('step2/' + dataset + '/' + case  + '.pdb', 'r')
        for line in fp:
            if len(line) >= 50:
                if line[:4] == 'ATOM':
                    resid = int(line[22:26])
                    pdb_resids.add(resid)
        fp.close()
    pdb_length = len(pdb_resids)

    if fasta_length == pdb_length:
        if fasta_length:
            pass
        else:
            if os.path.exists('step2/' + dataset + '/' + case + '.pdb'):
                os.system('rm step2/' + dataset + '/' + case + '.pdb')
            cmds.append('python /opt/DPAM/scripts/step2_get_AFDB_pdbs.py ' + dataset + ' ' + case + ' ' + case + '\n')
    else:
        if os.path.exists('step2/' + dataset + '/' + case + '.pdb'):
            os.system('rm step2/' + dataset + '/' + case + '.pdb')
        cmds.append('python /opt/DPAM/scripts/step2_get_AFDB_pdbs.py ' + dataset + ' ' + case + ' ' + case + '\n')

if cmds:
    logs=batch_run(cmds, ncore)
    fail = [i for i in logs if 'fail' in i]
    if fail:
        with open(dataset + '_step2.log','w') as f:
            for i in fail:
                f.write(i+'\n')
    else:
        with open(dataset + '_step2.log','w') as f:
            f.write('done\n')
else:
    with open(dataset + '_step2.log','w') as f:
        f.write('done\n')
