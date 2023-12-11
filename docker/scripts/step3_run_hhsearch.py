#!/opt/conda/bin/python
import os, sys, subprocess
def run_cmd(cmd):
    status = subprocess.run(cmd,shell = True).returncode
    return status

dataset = sys.argv[1]
prot = sys.argv[2]
cpu = sys.argv[3]

if os.path.exists('step3/' + dataset + '/' + prot + '.hhsearch'):
    pass
else:
    if os.path.exists('step3/' + dataset + '/' + prot + '.hmm'):
        status = run_cmd('hhsearch -cpu ' + cpu + ' -Z 100000 -B 100000 -i step3/' + dataset + '/' + prot + '.hmm -d /mnt/databases/pdb70/pdb70 -o step3/' + dataset + '/' + prot + '.hhsearch')
        if status != 0:
            sys.exit(1)
    elif os.path.exists('step3/' + dataset + '/' + prot + '.hhm'):
        os.system('mv step3/' + dataset + '/' + prot + '.hhm step3/' + dataset + '/' + prot + '.hmm')
        status = run_cmd('hhsearch -cpu ' + cpu + ' -Z 100000 -B 100000 -i step3/' + dataset + '/' + prot + '.hmm -d /mnt/databases/pdb70/pdb70 -o step3/' + dataset + '/' + prot + '.hhsearch')
        if status != 0:
            sys.exit(1)
    else:
        cmds= ['hhblits -cpu ' + cpu + ' -i step1/' + dataset + '/' + prot + '.fa -d /mnt/databases/UniRef30_2022_02/UniRef30_2022_02 -oa3m step3/' + dataset + '/' + prot + '.a3m','addss.pl step3/' + dataset + '/' + prot + '.a3m step3/' + dataset + '/' + prot + '.a3m.ss -a3m','mv step3/' + dataset + '/' + prot + '.a3m.ss step3/' + dataset + '/' + prot + '.a3m','hhmake -i step3/' + dataset + '/' + prot + '.a3m -o step3/' + dataset + '/' + prot + '.hmm','hhsearch -cpu ' + cpu + ' -Z 100000 -B 100000 -i step3/' + dataset + '/' + prot + '.hmm -d /mnt/databases/pdb70/pdb70 -o step3/' + dataset + '/' + prot + '.hhsearch']
        for cmd in cmds:
            status = run_cmd(cmd)
            if status != 0:
                sys.exit(1)
