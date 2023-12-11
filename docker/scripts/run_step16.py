#!/opt/conda/bin/python
import os, sys, subprocess

dataset = sys.argv[1]

if not os.path.exists('step16'):
    os.system('mkdir step16/')
if not os.path.exists('step16/' + dataset):
    os.system('mkdir step16/' + dataset)

fp = open(dataset + '_struc.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

need_prots = []
for prot in prots:
    if os.path.exists('step16/' + dataset + '/' + prot + '.result'):
        word_counts = set([])
        fp = open('step16/' + dataset + '/' + prot + '.result', 'r')
        for line in fp:
            words = line.split()
            word_counts.add(len(words))
        fp.close()
        if len(word_counts) == 1 and 21 in word_counts:
            pass
        else:
            os.system('rm step16/' + dataset + '/' + prot + '.result')
            need_prots.append(prot)
    elif os.path.exists('step16/' + dataset + '/' + prot + '.done'):
        pass
    else:
        need_prots.append(prot)


if need_prots:
    rp = open('step16_' + dataset + '.list', 'w')
    for prot in need_prots:
        rp.write(prot + '\n')
    rp.close()
    rcode=subprocess.run('step16_run_domass.py ' + dataset,shell=True).returncode
    if rcode!=0:
        with open(dataset + '_step16.log','w')as f:
            f.write(' '.join(need_prots)+' fail\n')
    else:
        with open(dataset + '_step16.log','w')as f:
            f.write('done\n')
    os.system('rm step16_' + dataset + '*.list\n')
else:
    with open(dataset + '_step16.log','w')as f:
        f.write('done\n')
