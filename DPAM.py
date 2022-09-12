import sys,os,time
script_dir=os.path.dirname(os.path.realpath(__file__))


def check_progress(mName,output_dir):
    full_progress=range(1,12)
    if os.path.exists(output_dir):
        if os.path.exists(output_dir+'/log'):
            with open(output_dir+'/log') as f:
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

model=sys.argv[1]
output_dir=sys.argv[2]
pae=model.split('-model')[0]+'-predicted_aligned_error_v2.json'

if model[0]!='/':
    cwd=os.getcwd()
    model=cwd+'/'+model
    pae=cwd+'/'+pae

if not os.path.exists(output_dir):
    os.system('mkdir '+output_dir)

os.chdir(output_dir)
if not os.path.exists(model.split('/')[-1]):
    os.system(f'cp {model} {output_dir}')
if not os.path.exists(pae.split('/')[-1]):
    os.system(f'cp {pae} {output_dir}')


if model[-4:]=='.pdb':
    mName=model.split('/')[-1].split('.pdb')[0]
    mFormat='pdb'
    os.system(f'python {script_dir}/get_seq.py {mName} {mFormat}')
elif model[-4:]=='.cif': 
    mName=model.split('/')[-1].split('.cif')[0]
    mFormat='cif'
    os.system(f'python {script_dir}/get_seq.py {mName} {mFormat}')
    if not os.path.exists(output_dir+'/'+mName+'.pdb'):
        print(f'python {script_dir}/get_pdbs.py {mName}')
        os.system(f'python {script_dir}/get_pdbs.py {mName}')
else:
    print('Unknow format of model detected, please add extension to indicate the format of the model')

progress=check_progress(model,output_dir)
if progress!=[]:
    cmds=[f'python {script_dir}/run_hhsearch.py {mName}',f'python {script_dir}/run_dali.py {mName}',f'python {script_dir}/map_hhsearch_to_ecod.py {mName}',f'python {script_dir}/filter_dali.py {mName}',f'python {script_dir}/iterative_dali.py {mName}',f'python {script_dir}/analyze_PDB.py {mName}',f'python {script_dir}/get_support.py {mName}',f'python {script_dir}/get_good_domains.py {mName}',f'python {script_dir}/get_sse.py {mName}',f'python {script_dir}/get_diso.py {mName}',f'python {script_dir}/parse_domains.py {mName}']
    for cmd in cmds[progress-1:]:
        os.system(cmd)
    print(f'Domain Parsing for {mName} done')
else:
    print(f'Previous domain parsing result for {mName} is complete.Exit')

