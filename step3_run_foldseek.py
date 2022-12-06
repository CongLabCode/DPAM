import os, sys,time


prefix = sys.argv[1]
threads = sys.argv[2]
wd = sys.argv[3]
data_dir = sys.argv[4]

if os.getcwd() != wd:
    os.chdir(wd)
if not os.path.exists('foldseek_tmp'):
    os.system('mkdir foldseek_tmp')

os.system(f'foldseek easy-search {prefix}.pdb {data_dir}/ECOD_foldseek_DB/ECOD_foldseek_DB {prefix}.foldseek foldseek_tmp -e 1000000 --max-seqs 1000000 --threads {threads} > {prefix}_foldseek.log')
os.system('rm -rf foldseek_tmp')
