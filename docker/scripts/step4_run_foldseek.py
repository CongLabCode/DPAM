#!/opt/conda/bin/python
import os, sys, random, string
def generate_random_directory_name():
    characters = string.ascii_lowercase + string.digits
    random_string = ''.join(random.choice(characters) for _ in range(8))
    return random_string


dataset = sys.argv[1]
ncore = sys.argv[2]
tmp_dir_name = generate_random_directory_name()

if not os.path.exists('/tmp/' + dataset + '_' + tmp_dir_name):
    os.system('mkdir /tmp/' + dataset + '_' + tmp_dir_name)

fp = open('step4/' + dataset + '_step4.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

for prot in prots:
    os.system('foldseek easy-search step2/' + dataset + '/' + prot + '.pdb /mnt/databases/ECOD_foldseek_DB/ECOD_foldseek_DB step4/' + dataset + '/' + prot + '.foldseek /tmp/' + dataset + '_' + tmp_dir_name + ' -e 1000 --max-seqs 1000000 --threads ' + ncore + ' >> /tmp/step4_' + dataset + '_' + tmp_dir_name + '.log')
    fp = open('step4/' + dataset + '/' + prot + '.foldseek', 'r')
    countline = 0
    for line in fp:
        countline += 1
    fp.close()
    if not countline:
        os.system('echo \'done\' > step4/' + dataset + '/' + prot + '.done')

os.system('rm -rf /tmp/' + dataset + '_' + tmp_dir_name)
os.system('mv /tmp/step4_' + dataset + '_' + tmp_dir_name + '.log ./')
