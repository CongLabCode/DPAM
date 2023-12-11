#!/opt/conda/bin/python
import os, sys, subprocess
from multiprocessing import Pool
def run_cmd(cmd):
    status=subprocess.run(cmd,shell=True).returncode
    if status==0:
        return cmd + ' succeed'
    else:
        return cmd + ' fail'

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

if not os.path.exists('step14'):
    os.system('mkdir step14/')
if not os.path.exists('step14/' + dataset):
    os.system('mkdir step14/' + dataset)

fp = open(dataset + '_struc.list', 'r')
cases = []
for line in fp:
    words = line.split()
    cases.append(words[0])
fp.close()


need_cases = []
for case in cases:
    prot = case
    if os.path.exists('step14/' + dataset + '/' + prot + '.domains'):
        word_counts = set([])
        fp = open('step14/' + dataset + '/' + prot + '.domains', 'r')
        for line in fp:
            words = line.split()
            word_counts.add(len(words))
        fp.close()
        if len(word_counts) == 1 and 2 in word_counts:
            pass
        else:
            os.system('rm step14/' + dataset + '/' + prot + '.domains')
            need_cases.append(case)
    elif os.path.exists('step14/' + dataset + '/' + prot + '.done'):
        pass
    else:
        need_cases.append(case)

if need_cases:
    cmds = []
    for case in need_cases:
        cmds.append('step14_parse_domains.py ' + dataset + ' ' + case)
    logs = batch_run(cmds, ncore)
    fail = [i for i in logs if 'fail' in i]
    if fail:
        with open(dataset + '_step14.log','w') as f:
            for i in fail:
                f.write(i+'\n')
    else:
        with open(dataset + '_step14.log','w') as f:
            f.write('done\n')
else:
    with open(dataset + '_step14.log','w') as f:
        f.write('done\n')
