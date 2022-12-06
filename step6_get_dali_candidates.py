import sys,os

prefix = sys.argv[1]
wd = sys.argv[2]

if os.getcwd() != wd:
    os.chdir(wd)

domains = set([])
fp = open(prefix + '.map2ecod.result', 'r')
for countl, line in enumerate(fp):
    if countl:
        words = line.split()
        domains.add(words[0])
fp.close()

fp = open(prefix + '.foldseek.flt.result','r')
for countl, line in enumerate(fp):
    if countl:
        words = line.split()
        domains.add(words[0])
fp.close()

rp = open(prefix + '_hits4Dali', 'w')
for domain in domains:
    rp.write(domain + '\n')
rp.close()
