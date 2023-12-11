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
        sample=cmd.split()[2]
        process = pool.apply_async(run_cmd,(cmd,))
        result.append(process)
    for process in result:
        log.append(process.get())
    return log


dataset = sys.argv[1]
ncore = int(sys.argv[2])
if not os.path.exists('step7'):
    os.system('mkdir step7/')
if not os.path.exists('step7/' + dataset):
    os.system('mkdir step7/' + dataset)

fp = open(dataset + '_struc.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

need_prots = []
for prot in prots:
    if os.path.exists('step7/' + dataset + '/' + prot + '_hits'):
        pass
    elif os.path.exists('step7/' + dataset + '/' + prot + '.done'):
        pass
    else:
        need_prots.append(prot)

if need_prots:
    cmds = []
    for prot in need_prots:
        cmds.append('step7_prepare_dali.py ' + dataset + ' ' + prot)
    logs = batch_run(cmds, ncore)
    fail = [i for i in logs if 'fail' in i]
    if fail:
        with open(dataset + '_step7.log','w') as f :
            for i in fail:
                f.write(i+'\n')
    else:
        with open(dataset + '_step7.log','w') as f:
            f.write('done\n')
else:
    with open(dataset + '_step7.log','w') as f:
        f.write('done\n')
