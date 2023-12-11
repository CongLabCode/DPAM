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

if not os.path.exists('step12'):
    os.system('mkdir step12/')
if not os.path.exists('step12/' + dataset):
    os.system('mkdir step12/' + dataset)

fp = open(dataset + '_struc.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

need_prots = []
for prot in prots:
    if os.path.exists('step12/' + dataset + '/' + prot + '.sse'):
        fp = open('step12/' + dataset + '/' + prot + '.sse', 'r')
        word_counts = set([])
        for line in fp:
            words = line.split()
            word_counts.add(len(words))
        fp.close()
        if len(word_counts) == 1 and 4 in word_counts:
            pass
        else:
            os.system('rm step12/' + dataset + '/' + prot + '.sse')
            need_prots.append(prot)
    else:
        need_prots.append(prot)

if need_prots:
    cmds = []
    for prot in need_prots:
        cmds.append('step12_get_sse.py ' + dataset + ' ' + prot + '\n')
    logs = batch_run(cmds, ncore)
    fail = [i for i in logs if 'fail' in i]
    if fail:
        with open(dataset + '_step12.log','w') as f:
            for i in fail:
                f.write(i+'\n')
    else:
        with open(dataset + '_step12.log','w') as f:
            f.write('done\n')
else:
    with open(dataset + '_step12.log','w') as f:
        f.write('done\n')
