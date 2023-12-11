#!/opt/conda/bin/python
import os, sys

dataset = sys.argv[1]
prot = sys.argv[2]


domains = set([])
if os.path.exists('step5/' + dataset + '/' + prot + '.result'):
    fp = open('step5/' + dataset + '/' + prot + '.result', 'r')
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            domains.add(words[1])
    fp.close()

if os.path.exists('step6/' + dataset + '/' + prot + '.result'):
    fp = open('step6/' + dataset + '/' + prot + '.result','r')
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            domains.add(words[0])
    fp.close()

if domains:
    rp = open('step7/' + dataset + '/' + prot + '_hits', 'w')
    for domain in domains:
        rp.write(domain + '\n')
    rp.close()
else:
    os.system('echo \'done\' > step7/' + dataset + '/' + prot + '.done')
