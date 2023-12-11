#!/opt/conda/bin/python
import os, sys

dataset = sys.argv[1]
if not os.path.exists('step20'):
    os.system('mkdir step20')
os.system('mkdir step20/' + dataset)

fp = os.popen('ls -1 step19/' + dataset + '/*.result')
prots = []
for line in fp:
    prot = line.split('/')[2].split('.result')[0]
    prots.append(prot)
fp.close()

domains = []
for prot in prots:
    get_domains = set([])
    fp = open('step19/' + dataset + '/' + prot + '.result', 'r')
    for line in fp:
        words = line.split()
        domain1 = words[0]
        if not domain1 in get_domains:
            get_domains.add(domain1)
            resids1 = set([])
            for seg in words[1].split(','):
                if '-' in seg:
                    start = int(seg.split('-')[0])
                    end = int(seg.split('-')[1])
                    for res in range(start, end + 1):
                        resids1.add(res)
                else:
                    resids1.add(int(seg))
            domains.append([prot, domain1, resids1])
        
        domain2 = words[2]
        if not domain2 in get_domains:
            get_domains.add(domain2)
            resids2 = set([])
            for seg in words[3].split(','):
                if '-' in seg:
                    start = int(seg.split('-')[0])
                    end = int(seg.split('-')[1])
                    for res in range(start, end + 1):
                        resids2.add(res)
                else:
                    resids2.add(int(seg))
            domains.append([prot, domain2, resids2])
fp.close()

for item in domains:
    prot = item[0]
    dname = item[1]
    resids = item[2]
    fp = open('step2/' + dataset + '/' + prot + '.pdb', 'r')
    rp = open('step20/' + dataset + '/' + prot + '_' + dname + '.pdb', 'w')
    for line in fp:
        resid = int(line[22:26])
        if resid in resids:
            rp.write(line)
    fp.close()
    rp.close()
