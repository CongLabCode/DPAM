#!/opt/conda/bin/python
import os, sys, subprocess,time
def run_cmd(cmd):
    status = subprocess.run(cmd,shell = True).returncode
    if status == 0:
        return cmd + ' succeed'
    else:
        return cmd + ' fail'


dataset = sys.argv[1]
ncore = sys.argv[2]
if not os.path.exists('step3'):
    os.system('mkdir step3/')
if not os.path.exists('step3/' + dataset):
    os.system('mkdir step3/' + dataset)

fp = open(dataset + '_struc.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

cmds = []
for prot in prots:
    if os.path.exists('step3/' + dataset + '/' + prot + '.hmm') and os.path.exists('step3/' + dataset + '/' + prot + '.hhsearch'):
        fp = open('step3/' + dataset + '/' + prot + '.hmm', 'r')
        get_sspred = 0
        get_ssconf = 0
        for line in fp:
            if len(line) >= 10:
                if line[0] == '>' and line[1:8] == 'ss_pred':
                    get_sspred = 1
                elif line[0] == '>' and line[1:8] == 'ss_conf':
                    get_ssconf = 1
            if get_sspred and get_ssconf:
                break
        fp.close()

        if get_sspred and get_ssconf:
            pass
        elif os.path.exists('step3/' + dataset + '/' + prot + '.a3m'):
            fp = open('step3/' + dataset + '/' + prot + '.a3m', 'r')
            count_line = 0
            for line in fp:
                count_line += 1
            fp.close()
            if count_line == 2:
                get_sspred = 1
                get_ssconf = 1

        fp = open('step3/' + dataset + '/' + prot + '.hhsearch', 'r')
        start = 0
        end = 0
        hitsA = set([])
        hitsB = set([])
        for line in fp:
            words = line.split()
            if len(words) >= 2:
                if words[0] == 'No' and words[1] == 'Hit':
                    start = 1
                elif words[0] == 'No' and words[1] == '1':
                    hitsB.add(int(words[1]))
                    end = 1
                elif start and not end:
                    hitsA.add(int(words[0]))
                elif end:
                    if words[0] == 'No':
                        hitsB.add(int(words[1]))
        fp.close()
        last_words = line.split()
        
        if get_sspred and get_ssconf and hitsA == hitsB and not last_words:
            pass
        else:
            os.system('rm step1/' + dataset + '/' + prot + '.hhr')
            os.system('rm step3/' + dataset + '/' + prot + '.a3m')
            os.system('rm step3/' + dataset + '/' + prot + '.hmm')
            os.system('rm step3/' + dataset + '/' + prot + '.hhsearch')
            cmds.append('step3_run_hhsearch.py ' + dataset + ' ' + prot + ' ' + ncore)
    else:
        if os.path.exists('step1/' + dataset + '/' + prot + '.hhr'):
            os.system('rm step1/' + dataset + '/' + prot + '.hhr')
        if os.path.exists('step3/' + dataset + '/' + prot + '.a3m'):
            os.system('rm step3/' + dataset + '/' + prot + '.a3m')
        if os.path.exists('step3/' + dataset + '/' + prot + '.hmm'):
            os.system('rm step3/' + dataset + '/' + prot + '.hmm')
        if os.path.exists('step3/' + dataset + '/' + prot + '.hhsearch'):
            os.system('rm step3/' + dataset + '/' + prot + '.hhsearch')
        cmds.append('step3_run_hhsearch.py ' + dataset + ' ' + prot + ' ' + ncore)

if cmds:
    fail=[]
    for cmd in cmds:
        for i in range(5):
            log=run_cmd(cmd)
            if 'succeed' in log:
                break
            time.sleep(1)
        else:
            fail.append(log)
    if fail:
        with open(dataset + '_step3.log','w') as f:
            for i in fail:
                f.write(i + '\n')
        sys.exit(1)
    else:
        with open(dataset + '_step3.log','w') as f:
            f.write('done\n')
