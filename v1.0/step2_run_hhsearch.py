import os, sys

prefix = sys.argv[1]
CPUs = sys.argv[2]
wd = sys.argv[3]
data_dir=sys.argv[4]
if os.getcwd() != wd:
    os.chdir(wd)

print (f'hhblits -cpu {CPUs} -i {prefix}.fa -d {data_dir}/UniRef30_2022_02/UniRef30_2022_02 -oa3m {prefix}.a3m')
os.system(f'hhblits -cpu {CPUs} -i {prefix}.fa -d {data_dir}/UniRef30_2022_02/UniRef30_2022_02 -oa3m {prefix}.a3m')
print (f'addss.pl {prefix}.a3m {prefix}.a3m.ss -a3m')
os.system(f'addss.pl {prefix}.a3m  {prefix}.a3m.ss -a3m')
os.system(f'mv {prefix}.a3m.ss {prefix}.a3m')
print (f'hhmake -i {prefix}.a3m -o {prefix}.hmm')
os.system(f'hhmake -i {prefix}.a3m -o {prefix}.hmm')
print (f'hhsearch -cpu {CPUs} -Z 100000 -B 100000 -i {prefix}.hmm -d {data_dir}/pdb70/pdb70 -o {prefix}.hhsearch')
os.system(f'hhsearch -cpu {CPUs} -Z 100000 -B 100000 -i {prefix}.hmm -d {data_dir}/pdb70/pdb70 -o {prefix}.hhsearch')
