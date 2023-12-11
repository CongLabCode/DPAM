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

if not os.path.exists('step10'):
    os.system('mkdir step10/')
if not os.path.exists('step10/' + dataset):
    os.system('mkdir step10/' + dataset)

fp = open(dataset + '_struc.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

need_prots = set([])
for prot in prots:
    get_seq = 0
    if os.path.exists('step10/' + dataset + '/' + prot + '_sequence.result'):
        fp = open('step10/' + dataset + '/' + prot + '_sequence.result','r')
        word_counts = set([])
        for line in fp:
            words = line.split()
            word_counts.add(len(words))
        fp.close()
        if len(word_counts) == 1 and 8 in word_counts:
            get_seq = 1
        else:
            os.system('rm step10/' + dataset + '/' + prot + '_sequence.result')
            need_prots.add(prot)

    get_str = 0
    if os.path.exists('step10/' + dataset + '/' + prot + '_structure.result'):
        fp = open('step10/' + dataset + '/' + prot + '_structure.result','r')
        word_counts = set([])
        for line in fp:
            words = line.split()
            word_counts.add(len(words))
        fp.close()
        if len(word_counts) == 1 and 12 in word_counts:
            get_str = 1
        else:
            os.system('rm step10/' + dataset + '/' + prot + '_structure.result')
            need_prots.add(prot)

    if get_seq and get_str:
        pass
    elif os.path.exists('step10/' + dataset + '/' + prot + '.done'):
        pass
    else:
        need_prots.add(prot)


if need_prots:
    cmds = []
    for prot in need_prots:
        cmds.append('step10_get_support.py ' + dataset + ' ' + prot + '\n')
    logs = batch_run(cmds, ncore)
    fail = [i for i in logs if 'fail' in i]
    if fail:
        with open(dataset + '_step10.log','w') as f:
            for i in fail:
                f.write(i+'\n')
    else:
        with open(dataset + '_step10.log','w') as f:
            f.write('done\n')
else:
    with open(dataset + '_step10.log','w') as f:
        f.write('done\n')
