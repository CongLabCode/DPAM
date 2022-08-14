import sys,os,time
script_dir=os.path.dirname(os.path.realpath(__file__))
with open(script_dir+'/config_file') as f:
    configs=f.readlines()
configs=[i.strip() for i in configs if i.strip()!='']
configs={i.split(':')[0]:i.split(':')[1] for i in configs}


mName=sys.argv[1]
cwd=os.getcwd()
print(cwd)
ecod_database_dir=configs['ecod_database_dir']
if not os.path.exists('ECOD_database'):
    os.system(f'ln -s {ecod_database_dir} ./ECOD_database')
if os.path.exists('dali_tmp'):
    os.system('rm -rf dali_tmp')
os.system('mkdir dali_tmp')
os.chdir('dali_tmp')
alldomains=os.listdir(ecod_database_dir)
logp=open(f'../{mName}_dali.failed','w')
outputp=open(f'../{mName}_dali.result','w')
for domain in alldomains:
    cmd=f"dali.pl --pdbfile1 ../{mName}.pdb --pdbfile2 ../ECOD_database/{domain} --dat1 ./ --dat2 ./ --outfmt summary,alignments,transrot"
    os.system(cmd)
    fp = os.popen("ls -1 mol*.txt")
    info = fp.readlines()
    fp.close()
    filenames = []
    for line in info:
        filenames.append(line[:-1])
    if filenames:
        outputp.write(">" + mName + "_" + domain + "\n")
        for filename in filenames:
            fp = open(filename, "r")
            newlines = fp.readlines()
            fp.close()
            for newline in newlines:
                outputp.write(newline)
    else:
        logp.write(mName + "_" + domain + "\n")
    if os.getcwd()==cwd+'/dali_tmp':
        os.system('rm -rf *')
        time.sleep(0.02)
outputp.close()
logp.close()
with open(f'{cwd}/log','a') as f:
    f.write('2\n')
