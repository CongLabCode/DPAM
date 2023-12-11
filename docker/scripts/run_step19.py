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
if not os.path.exists('step19'):
    os.system('mkdir step19/')
if not os.path.exists('step19/' + dataset):
    os.system('mkdir step19/' + dataset)

fp = open(dataset + '_struc.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

need_prots = set([])
for prot in prots:
    check_info = 0
    if os.path.exists('step19/' + dataset + '/' + prot + '.info'):
        word_counts = set([])
        fp = open('step19/' + dataset + '/' + prot + '.info', 'r')
        for line in fp:
            words = line.split()
            word_counts.add(len(words))
        fp.close()
        if len(word_counts) == 1 and 2 in word_counts:
            check_info = 1
        else:
            os.system('rm step19/' + dataset + '/' + prot + '.info')
            need_prots.add(prot)

    check_result = 0
    if os.path.exists('step19/' + dataset + '/' + prot + '.result'):
        word_counts = set([])
        fp = open('step19/' + dataset + '/' + prot + '.result', 'r')
        for line in fp:
            words = line.split()
            word_counts.add(len(words))
        fp.close()
        if len(word_counts) == 1 and 4 in word_counts:
            check_result = 1
        else:
            os.system('rm step19/' + dataset + '/' + prot + '.result')
            need_prots.add(prot)

    if check_info and check_result:
        pass
    elif os.path.exists('step19/' + dataset + '/' + prot + '.done'):
        pass
    else:
        need_prots.add(prot)

if need_prots:
    cmds = []
    for prot in need_prots:
        cmds.append('step19_get_merge_candidates.py ' + dataset + ' ' + prot)
    logs = batch_run(cmds, ncore)
    fail = [i for i in logs if 'fail' in i]
    if fail:
        with open(dataset + '_step19.log','w') as f:
            for i in fail:
                f.write(i+'\n')
    else:
        with open(dataset + '_step19.log','w') as f:
            f.write('done\n')
else:
    with open(dataset + '_step19.log','w') as f:
        f.write('done\n')
