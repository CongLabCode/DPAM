#!/usr/bin/env python
import os, sys, time, subprocess
dataset = sys.argv[1]
ncore = sys.argv[2]
wd = os.getcwd()

for step in range(1,25):
    if 1 <= step <= 19 or step == 21 or step == 23: 
        if os.path.exists(dataset + '_step' + str(step) + '.log'):
            with open(dataset + '_step' + str(step) + '.log') as f:
                step_logs = f.read()
            if 'done\n' != step_logs:
                rcode = subprocess.run('run_step' + str(step) + '.py ' + dataset + ' ' + ncore,shell=True).returncode
                if rcode != 0:
                    print(f'Error in step{step}')
                    sys.exit()
        else:
            for s in range(step,25):
                os.system('rm ' + dataset + '_step' + str(s) + '.log')
                os.system('rm -rf step_' + str(s) + '/' + dataset + '/*')
            rcode = subprocess.run('run_step' + str(step) + '.py ' + dataset + ' ' + ncore,shell=True).returncode
            if rcode != 0:
                print(f'Error in step{step}')
                sys.exit()
    elif step == 20:
        run_flag = 0
        if os.path.exists(dataset + '_step' + str(step) + '.log'):
            with open(dataset + '_step' + str(step) + '.log') as f:
                step_logs=f.read()
            if 'done\n' != step_logs:
                run_flag = 1
        else:
            run_flag = 1
        if run_flag == 1:
            for s in range(step,25):
                os.system('rm ' + dataset + '_step' + str(s) + '.log')
                os.system('rm -rf step_' + str(s) + '/' + dataset + '/*')
            status_code = subprocess.run('step20_extract_domains.py ' + dataset,shell=True).returncode
            if status_code == 0:
                with open(dataset + '_step20.log','w') as f:
                    f.write('done\n')
            else:
                with open(dataset + '_step20.log','w') as f:
                    f.write('fail\n')
                print(f'Error in step{step}')
                sys.exit()
    elif step == 22:
        run_flag = 0
        if os.path.exists(dataset + '_step' + str(step) + '.log'):
            with open(dataset + '_step' + str(step) + '.log') as f:
                step_logs=f.read()
            if 'done\n' != step_logs:
                run_flag = 1
        else:
            run_flag = 1
        if run_flag == 1:
            for s in range(step,25):
                os.system('rm ' + dataset + '_step' + str(s) + '.log')
                os.system('rm -rf step_' + str(s) + '/' + dataset + '/*')
            status_code = subprocess.run('step22_merge_domains.py ' + dataset,shell=True).returncode
            if status_code == 0:
                with open(dataset + '_step22.log','w') as f:
                    f.write('done\n')
            else:
                with open(dataset + '_step22.log','w') as f:
                    f.write('fail\n')
                print(f'Error in step{step}')
                sys.exit()
    elif step == 24:
        run_flag = 0
        if os.path.exists(dataset + '_step' + str(step) + '.log'):
            with open(dataset + '_step' + str(step) + '.log') as f:
                step_logs = f.read()
            if 'done\n' != step_logs:
                run_flag = 1
        else:
            run_flag = 1
        if run_flag == 1:
            for s in range(step,25):
                os.system('rm ' + dataset + '_step' + str(s) + '.log')
                os.system('rm -rf step_' + str(s) + '/' + dataset + '/*')
            status_code = subprocess.run('step24_integrate_results.py ' + dataset,shell=True).returncode
            if status_code == 0:
                with open(dataset + '_step24.log','w') as f:
                    f.write('done\n')
            else:
                with open(dataset + '_step24.log','w') as f:
                    f.write('fail\n')
                print(f'Error in step{step}')
                sys.exit()
filelist=[wd + '/' + dataset + '_step' + str(k)+'.log' for k in range(1,25)]
undone = 24
for name in filelist:
    with open(name) as f:
        info = f.read()
    if info.strip() == 'done':
        undone = undone - 1
    else:
        print(dataset + ' ' + name.split('/')[-1].split(dataset + '_')[1]+' has errors..Fail')
        break
if undone == 0:
    print(dataset + ' done')
