#!/opt/conda/bin/python

import os, sys, subprocess
def run_cmd(cmd):
    status = subprocess.run(cmd,shell = True).returncode
    if status == 0:
        return cmd + ' succeed'
    else:
        return cmd + ' fail'

dataset = sys.argv[1]
ncore = sys.argv[2]

if not os.path.exists('step4'):
    os.system('mkdir step4/')

if not os.path.exists('step4/' + dataset):
    os.system('mkdir step4/' + dataset)

fp = open(dataset + '_struc.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

need_prots = []
for prot in prots:
    if os.path.exists('step4/' + dataset + '/' + prot + '.foldseek'):
        fp = open('step4/' + dataset + '/' + prot + '.foldseek','r')
        word_counts = set([])
        for line in fp:
            words = line.split()
            word_counts.add(len(words))
        fp.close()
        if len(word_counts) == 1 and 12 in word_counts:
            pass
        elif os.path.exists('step4/' + dataset + '/' + prot + '.done'):
            pass
        else:
            os.system('rm step4/' + dataset + '/' + prot + '.foldseek')
            need_prots.append(prot)
    else:
        need_prots.append(prot)


if need_prots:
    with open('step4/' + dataset + '_step4.list','w') as f:
        for i in need_prots:
            f.write(i+'\n')
    log = run_cmd('step4_run_foldseek.py ' + dataset + ' ' + ncore + ' \n')
    if 'fail' in log:
        with open(dataset + '_step4.log','w') as f:
            f.write(dataset + ' fail\n')
    else:
        with open(dataset + '_step4.log','w') as f:
            f.write('done\n')
        os.system('rm step4/' + dataset + '_step4.list')
else:
    if os.path.exists('step4/' + dataset + '_step4.list'):
        os.system('rm step4/' + dataset + '_step4.list')
    with open(dataset + '_step4.log','w') as f:
        f.write('done\n')
