#!/opt/conda/bin/python
import os, sys, subprocess
from multiprocessing import Pool

def run_cmd(cmd):
    status=subprocess.run(cmd,shell=True).returncode
    if status==0:
        return cmd+' succeed'
    else:
        return cmd+' fail'

def batch_run(cmds,process_num):
    log=[]
    pool = Pool(processes=process_num)
    result = []
    for cmd in cmds:
        process = pool.apply_async(run_cmd,(cmd,))
        result.append(process)
    for process in result:
        log.append(process.get())
    return log


dataset = sys.argv[1]
ncore = int(sys.argv[2])
os.system('rm step21_' + dataset + '_*.list')

fp = os.popen('ls -1 step19/' + dataset + '/*.result')
prots = []
for line in fp:
    prot = line.split('/')[2].split('.')[0]
    prots.append(prot)
fp.close()

cases = []
all_cases = set([])
for prot in prots:
    fp = open('step19/' + dataset + '/' + prot + '.result', 'r')
    for line in fp:
        words = line.split()
        domain1 = words[0]
        resids1 = words[1]
        domain2 = words[2]
        resids2 = words[3]
        cases.append([prot, domain1, resids1, domain2, resids2])
        all_cases.add(prot + '_' + domain1 + '_' + domain2)
    fp.close()

get_cases = set([])
if os.path.exists('step21_' + dataset + '.result'):
    fp = open('step21_' + dataset + '.result','r')
    for line in fp:
        words = line.split()
        get_cases.add(words[0] + '_' + words[1] + '_' + words[2])
    fp.close()

if all_cases == get_cases:
    with open(dataset + '_step21.log','w') as f:
        f.write('done\n')
else:
    total = len(cases)
    batchsize = total // ncore + 1
    cmds = []
    for i in range(ncore):
        rp = open('step21_' + dataset + '_' + str(i) + '.list', 'w')
        for case in cases[batchsize * i : batchsize * i + batchsize]:
            rp.write(case[0] + '\t' + case[1] + '\t' + case[2] + '\t' + case[3] + '\t' + case[4] + '\n')
        rp.close()
        cmds.append('step21_compare_domains.py ' + dataset + ' ' + str(i))
    logs = batch_run(cmds, ncore)
    fail = [i for i in logs if 'fail' in i]
    if fail:
        with open(dataset + '_step21.log','w') as f:
            for i in fail:
                f.write(i+'\n')
    else:
        with open(dataset + '_step21.log','w') as f:
            f.write('done\n')
    status=subprocess.run('cat step21_' + dataset + '_*.result >> step21_' + dataset + '.result',shell=True).returncode
    os.system('rm step21_' + dataset + '_*.list')
    os.system('rm step21_' + dataset + '_*.result')
