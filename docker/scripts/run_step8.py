#!/opt/conda/bin/python
import os, sys, subprocess
def run_cmd(cmd):
    status=subprocess.run(cmd,shell=True).returncode
    if status==0:
        return cmd +' succeed'
    else:
        return cmd +' fail'

dataset = sys.argv[1]
ncore = sys.argv[2]
if not os.path.exists('step8'):
    os.system('mkdir step8/')

if not os.path.exists('step8/' + dataset):
    os.system('mkdir step8/' + dataset)

fp = open(dataset + '_struc.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

need_prots = []
for prot in prots:
    if os.path.exists('step8/' + dataset + '/' + prot + '_hits'):
        hit_count = 0
        fp = open('step8/' + dataset + '/' + prot + '_hits', 'r')
        hit_lines = []
        hit_line_count = 0
        bad = 0
        for line in fp:
            if line[0] == '>':
                hit_count += 1
                if hit_line_count:
                    if hit_line_count + 4 != len(hit_lines):
                        bad = 1
                        break
                words = line.split()
                hit_line_count = int(words[2])
                hit_lines = []
            else:
                hit_lines.append(line)
        fp.close()
        if hit_line_count:
            if hit_line_count + 4 != len(hit_lines):
                bad = 1
        if bad:
            os.system('rm step8/' + dataset + '/' + prot + '_hits')
            need_prots.append(prot)
        elif hit_count:
            pass
        else:
            if os.path.exists('step8/' + dataset + '/' + prot + '.done'):
                pass
            else:
                os.system('rm step8/' + dataset + '/' + prot + '_hits')
                need_prots.append(prot)
    else:
        if os.path.exists('step8/' + dataset + '/' + prot + '.done'):
            pass
        else:
            need_prots.append(prot)

if need_prots:
    print(need_prots)
    fail = []
    for prot in need_prots:
        log = run_cmd ('step8_iterative_dali.py ' + dataset + ' ' + prot + ' ' + ncore )
        if 'fail' in log:
            fail.append(log)
    if fail:
        with open(dataset + '_step8.log','w') as f:
            for i in fail:
                f.write(i+'\n')
    else:
        with open(dataset + '_step8.log','w') as f:
            f.write('done\n')
else:
    with open(dataset + '_step8.log','w') as f:
        f.write('done\n')
