import sys,os
script_dir=os.path.dirname(os.path.realpath(__file__))
with open(script_dir+'/config_file') as f:
    configs=f.readlines()
configs=[i.strip() for i in configs if i.strip()!='']
configs={i.split(':')[0]:i.split(':')[1] for i in configs}

mName=sys.argv[1]
hhuniprotDB=configs['hhuniprot']
hhpdb70=configs['hhpdb70']
cmd=f'hhblits -cpu 4 -i {mName}.fa -d {hhuniprotDB} -oa3m {mName}.a3m'
os.system(cmd)
cmd=f'addss.pl {mName}.a3m {mName}.a3m.ss -a3m'
os.system(cmd)
os.system(f'mv  {mName}.a3m.ss {mName}.a3m')
os.system(f'hhmake -i {mName}.a3m -o {mName}.hmm')
os.system(f'hhsearch -cpu 4 -i {mName}.hmm -d {hhpdb70} -o {mName}.hhsearch')
with open('log','a') as f:
    f.write('1\n')
