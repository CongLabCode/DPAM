#!/opt/conda/bin/python
import os, sys
import numpy as np

dataset = sys.argv[1]
prots = []
if os.path.exists(dataset + '_struc.list'):
    fp = open(dataset + '_struc.list','r')
    for line in fp:
        words = line.split()
        prots.append(words[0])
    fp.close()
else:
    fp = open(dataset + '.list','r')
    for line in fp:
        words = line.split()
        prots.append(words[0])
    fp.close()

full_domains = []
part_domains = []
miss_domains = []
for prot in prots:
    if os.path.exists('step23/' + dataset + '/' + prot + '.assign'):
        fp = open('step23/' + dataset + '/' + prot + '.assign', 'r')
        for line in fp:
            words = line.split()
            if words[0] == 'full':
                full_domains.append([prot] + words[1:-1])
            elif words[0] == 'part':
                part_domains.append([prot] + words[1:-1])
            elif words[0] == 'miss':
                miss_domains.append([prot] + words[1:-1])
        fp.close()


prot2goodres = {}
protres2sse = {}
prot2Hsses = {}
prot2Esses = {}
for prot in prots:
    prot2goodres[prot] = set([])
    protres2sse[prot] = {}
    prot2Hsses[prot] = set([])
    prot2Esses[prot] = set([])

    fp = open('step12/' + dataset + '/' + prot + '.sse', 'r')
    for line in fp:
        words = line.split()
        if words[2] != 'na':
            res = int(words[0])
            sse = int(words[2])
            prot2goodres[prot].add(res)
            protres2sse[prot][res] = sse
            if words[3] == 'H':
                prot2Hsses[prot].add(sse)
            if words[3] == 'E':
                prot2Esses[prot].add(sse)
    fp.close()


all_results = []
for words in miss_domains:
    prot = words[0]
    dom = words[1]
    resids = []
    for seg in words[2].split(','):
        if '-' in seg:
            start = int(seg.split('-')[0])
            end = int(seg.split('-')[1])
            for res in range(start, end + 1):
                resids.append(res)
        else:
            resids.append(int(seg))

    sses = []
    sse2count = {}
    for resid in resids:
        if resid in prot2goodres[prot]:
            sse = protres2sse[prot][resid]
            try:
                sse2count[sse] += 1
            except KeyError:
                sses.append(sse)
                sse2count[sse] = 1

    Hcount = 0
    Scount = 0
    for sse in sses:
        if sse in prot2Hsses[prot] and sse2count[sse] >= 6:
            Hcount += 1
        if sse in prot2Esses[prot] and sse2count[sse] >= 3:
            Scount += 1
    if Hcount + Scount < 3:
        judge = 'simple_topology'
    else:
        judge = 'low_confidence'
    all_results.append(words + [judge, str(Hcount), str(Scount)])


for words in part_domains:
    prot = words[0]
    dom = words[1]
    HHprob = float(words[6])
    ratio1 = float(words[8])
    ratio2 = float(words[9])
    resids = []
    for seg in words[2].split(','):
        if '-' in seg:
            start = int(seg.split('-')[0])
            end = int(seg.split('-')[1])
            for res in range(start, end + 1):
                resids.append(res)
        else:
            resids.append(int(seg))

    sses = []
    sse2count = {}
    for resid in resids:
        if resid in prot2goodres[prot]:
            sse = protres2sse[prot][resid]
            try:
                sse2count[sse] += 1
            except KeyError:
                sses.append(sse)
                sse2count[sse] = 1

    Hcount = 0
    Scount = 0
    for sse in sses:
        if sse in prot2Hsses[prot] and sse2count[sse] >= 6:
            Hcount += 1
        if sse in prot2Esses[prot] and sse2count[sse] >= 3:
            Scount += 1
    if Hcount + Scount < 3:
        if HHprob >= 0.95 and ratio1 >= 0.8 and ratio2 >= 0.8:
            judge = 'partial_domain'
        else:
            judge = 'simple_topology'
    else:
        judge = 'partial_domain'
    all_results.append(words + [judge, str(Hcount), str(Scount)])


for words in full_domains:
    prot = words[0]
    dom = words[1]
    HHprob = float(words[6])
    ratio1 = float(words[8])
    ratio2 = float(words[9])
    resids = []
    for seg in words[2].split(','):
        if '-' in seg:
            start = int(seg.split('-')[0])
            end = int(seg.split('-')[1])
            for res in range(start, end + 1):
                resids.append(res)
        else:
            resids.append(int(seg))

    sses = []
    sse2count = {}
    for resid in resids:
        if resid in prot2goodres[prot]:
            sse = protres2sse[prot][resid]
            try:
                sse2count[sse] += 1
            except KeyError:
                sses.append(sse)
                sse2count[sse] = 1

    Hcount = 0
    Scount = 0
    for sse in sses:
        if sse in prot2Hsses[prot] and sse2count[sse] >= 6:
            Hcount += 1
        if sse in prot2Esses[prot] and sse2count[sse] >= 3:
            Scount += 1
    if Hcount + Scount < 3:
        if HHprob >= 0.95 and ratio1 >= 0.8 and ratio2 >= 0.8:
            judge = 'good_domain'
        else:
            judge = 'simple_topology'
    else:
        judge = 'good_domain'
    all_results.append(words + [judge, str(Hcount), str(Scount)])


prot2results = {}
for item in all_results:
    prot = item[0]
    resids = []
    for seg in item[2].split(','):
        if '-' in seg:
            start = int(seg.split('-')[0])
            end = int(seg.split('-')[1])
            for res in range(start, end + 1):
                resids.append(res)
        else:
            resids.append(int(seg))
    mean_resid = np.mean(resids)
    if item[7] == 'na':
        try:
            prot2results[prot].append([mean_resid, item[2], item[3], item[4], item[5], item[6], item[7], item[8], item[9], item[10], item[11], item[12]])
        except KeyError:
            prot2results[prot] = [[mean_resid, item[2], item[3], item[4], item[5], item[6], item[7], item[8], item[9], item[10], item[11], item[12]]]
    else:
        try:
            prot2results[prot].append([mean_resid, item[2], item[3], item[4], item[5], item[6], str(round(float(item[7]) * 10, 1)), item[8], item[9], item[10], item[11], item[12]])
        except KeyError:
            prot2results[prot] = [[mean_resid, item[2], item[3], item[4], item[5], item[6], str(round(float(item[7]) * 10, 1)), item[8], item[9], item[10], item[11], item[12]]]


fp = open('/mnt/databases/ecod.latest.domains','r')
domain2keyword = {}
for line in fp:
    if line[0] != '#':
        words = line.split()
        domain2keyword[words[0]] = words[1]
fp.close()

os.system('mkdir step24/')
os.system('mkdir step24/' + dataset)
sump = open(dataset + '_domains', 'w')
sump.write('Protein\tDomain\tRange\tECOD_num\tECOD_key\tT-group\tDPAM_prob\tHH_prob\tDALI_zscore\tHit_cov\tTgroup_cov\tJudge\tHcount\tScount\n')
for prot in prots:
    try:
        results = prot2results[prot]
    except KeyError:
        results = []
    if results:
        results.sort(key = lambda x:x[0])
        rp = open('step24/' + dataset + '/' + prot + '_domains','w')
        for counti, item in enumerate(results):
            domain = item[2]
            try:
                keyword = domain2keyword[domain]
            except KeyError:
                keyword = 'na'
            rp.write('nD' + str(counti + 1) + '\t' + item[1] + '\t' + item[2] + '\t' + keyword + '\t' + item[3] + '\t' + item[4] + '\t' + item[5] + '\t' + item[6] + '\t' + item[7] + '\t' + item[8] + '\t' + item[9] + '\t' + item[10] + '\t' + item[11] + '\n')
            sump.write(prot + '\tnD' + str(counti + 1) + '\t' + item[1] + '\t' + item[2] + '\t' + keyword + '\t' + item[3] + '\t' + item[4] + '\t' + item[5] + '\t' + item[6] + '\t' + item[7] + '\t' + item[8] + '\t' + item[9] + '\t' + item[10] + '\t' + item[11] + '\n')
        rp.close()
sump.close()
