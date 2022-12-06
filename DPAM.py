import sys,os,time
from datetime import datetime
import subprocess
script_dir=os.path.dirname(os.path.realpath(__file__))

def print_usage ():
    print("usage: DPAM.py <input_cif/pdb> <input_pae> <accession> <output_dir> <threads> <datadir>")

def check_progress(basedir, basename):
    full_progress=range(1,13)
    if os.path.exists(basedir):
        if os.path.exists(basedir + '/' + basename + '_progress_logs'):
            with open(basedir + '/' + basename + '_progress_logs') as f:
                logs=f.readlines()
            logs=[i.strip() for  i in logs if i.strip()!='']
            if logs:
                logs=[int(i.split()[0]) for i in logs]
                full_progress=set(full_progress)-set(logs)
                full_progress=sorted(full_progress)
    if full_progress:
        progress=full_progress[0]
    else:
        progress=[]
    return progress

if len(sys.argv) != 7:
    print_usage()
else:
    logs=[]
    input_struc = sys.argv[1]
    input_pae = sys.argv[2]
    basename = sys.argv[3]
    basedir = sys.argv[4]
    threads = sys.argv[5]
    datadir = sys.argv[6]
    if basedir[0] != '/':
        basedir = os.getcwd() + basedir
    if not os.path.exists(basedir):
        os.system('mkdir ' + basedir)
    if '.cif' == input_struc[-4:]:
        os.system('cp ' + input_struc + ' ' + basedir + '/' + basename + '.cif')
    elif '.pdb' == input_struc[-4:]:
        os.system('cp ' + input_struc + ' ' + basedir + '/' + basename + '.pdb')
    else:
        print("Cannot recognize the structure file.Please use either mmcif or PDB as input. Exiting...")
        sys.exit()
    os.system('cp ' + input_pae + ' ' + basedir + '/' + basename + '.json')
    print('start input processing', datetime.now())
    status = subprocess.call(f'python {script_dir}/step1_get_AFDB_seqs.py {basename} {basedir}',shell=True)
    if status != 0:
        print('Cannot get protein sequence. Exiting...')
        sys.exit()
    status = subprocess.call(f'python {script_dir}/step1_get_AFDB_pdbs.py {basename} {basedir}',shell=True)
    if status != 0:
        print('Cannot process structure file.Exiting...')
        sys.exit()
    logs.append('0')
    progress=check_progress(basedir, basename)
    if progress!=[]:
        cmds=[]
        cmds.append(f'python {script_dir}/step2_run_hhsearch.py {basename} {threads} {basedir} {datadir}')
        cmds.append(f'python {script_dir}/step3_run_foldseek.py {basename} {threads} {basedir} {datadir}')
        cmds.append(f'python {script_dir}/step4_filter_foldseek.py {basename} {basedir}')
        cmds.append(f'python {script_dir}/step5_map_to_ecod.py {basename} {basedir} {datadir}')
        cmds.append(f'python {script_dir}/step6_get_dali_candidates.py {basename} {basedir}')
        cmds.append(f'python {script_dir}/step7_iterative_dali_aug_multi.py {basename} {threads} {basedir} {datadir}')
        cmds.append(f'python {script_dir}/step8_analyze_dali.py {basename} {basedir} {datadir}')
        cmds.append(f'python {script_dir}/step9_get_support.py {basename} {basedir} {datadir}')
        cmds.append(f'python {script_dir}/step10_get_good_domains.py {basename} {basedir} {datadir}')
        cmds.append(f'python {script_dir}/step11_get_sse.py {basename} {basedir}')
        cmds.append(f'python {script_dir}/step12_get_diso.py {basename} {basedir}')
        cmds.append(f'python {script_dir}/step13_parse_domains.py {basename} {basedir}')
        step=1
        for cmd in cmds[progress-1:]:
            print(f'start {cmd}', datetime.now())
            status = subprocess.call(cmd,shell=True)
            if status != 0:
                print(f"Error in {cmd}.Exiting...")
                with open(f'{basedir}/{basename}_progress_logs','w') as f:
                    for i in logs:
                        f.write(i+'\n')
                sys.exit()
            else:
                logs.append(str(step))
                step = step + 1
            print(f'end {cmd}', datetime.now()) 
        print(f'Domain Parsing for {basename} done')
        with open(f'{basedir}/{basename}_progress_logs','w') as f:
            for i in logs:
                f.write(i+'\n')
    else:
        print(f'Previous domain parsing result for {basename} is complete')
