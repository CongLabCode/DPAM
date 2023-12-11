#!/opt/conda/bin/python
import sys

dataset = sys.argv[1]
prot = sys.argv[2]

fp = open('step1/' + dataset + '/' + prot + '.fa','r')
query_seq = ''
for line in fp:
    if line[0] != '>':
        query_seq += line[:-1]
fp.close()
qlen = len(query_seq)

fp = open('step4/' + dataset + '/' + prot + '.foldseek', 'r')
hits = []
for line in fp:
    words = line.split()
    dnum = words[1].split('.')[0]
    qstart = int(words[6])
    qend = int(words[7])
    qresids = set([])
    for qres in range(qstart, qend + 1):
        qresids.add(qres)
    evalue = float(words[10])
    hits.append([dnum, evalue, qstart, qend, qresids])
fp.close()
hits.sort(key = lambda x:x[1])

qres2count = {}
for res in range(1, qlen + 1):
    qres2count[res] = 0

rp = open('step6/' + dataset + '/' + prot + '.result', 'w')
rp.write('ecodnum\tevalue\trange\n')
for hit in hits:
    dnum = hit[0]
    evalue = hit[1]
    qstart = hit[2]
    qend = hit[3]
    qresids = hit[4]
    for res in qresids:
        qres2count[res] += 1
    good_res = 0
    for res in qresids:
        if qres2count[res] <= 100:
            good_res += 1
    if good_res >= 10:
        rp.write(dnum + '\t' + str(evalue) + '\t' + str(qstart) + '-' + str(qend) + '\n')
rp.close()
