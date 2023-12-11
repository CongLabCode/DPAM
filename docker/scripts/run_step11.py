#!/opt/conda/bin/python
import os, sys, subprocess
from  multiprocessing import Pool

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

if not os.path.exists('step11'):
    os.system('mkdir step11/')
if not os.path.exists('step11/' + dataset):
    os.system('mkdir step11/' + dataset)

fp = open(dataset + '_struc.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

need_prots = []
for prot in prots:
    if os.path.exists('step11/' + dataset + '/' + prot + '.goodDomains'):
        fp = open('step11/' + dataset + '/' + prot + '.goodDomains','r')
        seq_word_counts = set([])
        str_word_counts = set([])
        for line in fp:
            words = line.split()
            if words[0] == 'sequence':
                seq_word_counts.add(len(words))
            elif words[0] == 'structure':
                str_word_counts.add(len(words))
        fp.close()

        bad_seq = 0
        bad_str = 0
        if seq_word_counts:
            if len(seq_word_counts) == 1 and 10 in seq_word_counts:
                pass
            else:
                bad_seq = 1
        if str_word_counts:
            if len(str_word_counts) == 1 and 16 in str_word_counts:
                pass
            else:
                bad_str = 1

        if bad_seq or bad_str:
            os.system('rm step11/' + dataset + '/' + prot + '.goodDomains')
            need_prots.append(prot)
    elif os.path.exists('step11/' + dataset + '/' + prot + '.done'):
        pass
    else:
        need_prots.append(prot)


if need_prots:
    cmds = []
    for prot in need_prots:
        cmds.append('step11_get_good_domains.py ' + dataset + ' ' + prot)
    logs = batch_run(cmds, ncore)
    fail = [i for i in logs if 'fail' in i]
    if fail:
        with open(dataset + '_step11.log','w') as f:
            for i in fail:
                f.write(i+'\n')
    else:
        with open(dataset + '_step11.log','w') as f:
            f.write('done\n')
else:
    with open(dataset + '_step11.log','w') as f:
        f.write('done\n')
