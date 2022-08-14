import os, sys

def get_domain_range(resids):
    segs = []
    resids.sort()
    cutoff1 = 5
    cutoff2 = len(resids) * 0.05
    cutoff = max(cutoff1, cutoff2)
    for resid in resids:
        if not segs:
            segs.append([resid])
        else:
            if resid > segs[-1][-1] + cutoff:
                segs.append([resid])
            else:
                segs[-1].append(resid)
    seg_string = []
    for seg in segs:
        start = str(seg[0])
        end = str(seg[-1])
        seg_string.append(start + '-' + end)
    return ','.join(seg_string)


prot = sys.argv[1]
fp = open(prot + '.pdb', 'r')
resids_set = set([])
for line in fp:
    if line[:4]=="ATOM":
        resid = int(line[22:26])
        resids_set.add(resid)
fp.close()
resids = list(resids_set)
resids.sort()


fp = open( prot + '_dali.result', 'r')
dnum = ''
qali = ''
sali = ''
zscore = 0
alignments = []
getit = 0
for line in fp:
    if line[0] == '>':
        if dnum and zscore and qali and sali:
            alignments.append([dnum, zscore, qali, sali])
#        myprot = line[1:-1].split('_')[0]
        myprot=prot
        dnum=line[1:-1].split(myprot+'_')[1].split('.pdb')[0]
#        dnum = line[1:-1].split('_')[1]
        if prot != myprot:
            print ('error1\t' + prot + '\t' + dnum)
        qali = ''
        sali = ''
        getit = 1
    elif getit:
        words = line.split()
        if len(words) >= 2:
            if words[0] == 'Query':
                qali += words[1]
            elif words[0] == 'Sbjct':
                sali += words[1]
            elif words[0] == 'No' and words[1] == '1:':
                for word in words:
                    if '=' in word:
                        subwords = word.split('=')
                        if subwords[0] == 'Z-score':
                            zinfo = subwords[1].split('.')
                            zscore = float(zinfo[0] + '.' + zinfo[1])
            elif words[0] == 'No' and words[1] == '2:':
                getit = 0
fp.close()
if dnum and zscore and qali and sali:
    alignments.append([dnum, zscore, qali, sali])


hits = []
for item in alignments:
    dnum = item[0]
    zscore = item[1]
    qali = item[2]
    sali = item[3]
    qinds = set([])
    length = len(qali)
    qposi = 0
    sposi = 0
    match = 0
    for i in range(length):
        if qali[i] != '-':
            qposi += 1
        if sali[i] != '-':
            sposi += 1
        if qali[i] != '-' and sali[i] != '-':
            if qali[i].isupper() and sali[i].isupper():
                match += 1
                qinds.add(qposi)
    hits.append([dnum, zscore, qposi, sposi, match, qinds])


hits.sort(key = lambda x:x[1], reverse = True)
all_qres2cov = {}
good_qres2cov = {}
good_hits = []
for hit in hits:
    dnum = hit[0]
    zscore = hit[1]
    qlen = hit[2]
    if qlen != len(resids):
        print ('error2\t' + prot + '\t' + dnum)
    else:
        slen = hit[3]
        match = hit[4]
        qinds = hit[5]
        
        if match >= 30:
            for qres in qinds:
                try:
                    all_qres2cov[qres] += 1
                except KeyError:
                    all_qres2cov[qres] = 1

            getit = 0
            if match / slen >= 0.5:
                for qres in qinds:
                    try:
                        good_qres2cov[qres] += 1
                        all_qres2cov[qres] += 1
                    except KeyError:
                        good_qres2cov[qres] = 1
                        all_qres2cov[qres] = 1
                    if good_qres2cov[qres] <= 5:
                        getit += 1

            else:
                for qres in qinds:
                    try:
                        all_qres2cov[qres] += 1
                    except KeyError:
                        all_qres2cov[qres] = 1
                    if all_qres2cov[qres] <= 5:
                        getit += 1

            if getit >= 10:
                good_hits.append(hit)


rp = open(prot + '_hits','w')
for hit in good_hits:
    dnum = hit[0]
    zscore = hit[1]
    qlen = hit[2]
    slen = hit[3]
    match = hit[4]
    raw_qinds = list(hit[5])
    raw_qinds.sort()
    raw_qresids = []
    for qind in raw_qinds:
        raw_qresids.append(resids[qind - 1])
    qrange = get_domain_range(raw_qresids)
    qresids = set([])
    qsegs = qrange.split(',')
    for qseg in qsegs:
        qedges = qseg.split('-')
        qstart = int(qedges[0])
        qend = int(qedges[1])
        for qres in range(qstart, qend + 1):
            qresids.add(qres)
    remain_resids = resids_set.difference(qresids)
    rp.write(dnum + '\t' + str(zscore) + '\t' + str(match) + '\t' + str(qlen) + '\t' + str(slen) + '\t' + qrange + '\t' + str(len(remain_resids)) + '\n')
    
    if not os.path.exists(prot):
        os.system('mkdir '+prot)
rp.close()
with open('log','a') as f:
    f.write('4\n')
