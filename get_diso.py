#!/usr/bin/python
import numpy as np
import os, sys, json, math, string


prot = sys.argv[1]
insses = set([])
res2sse = {}
fp = open(prot + '.sse','r')
for line in fp:
    words = line.split()
    if words[2] != 'na':
        sseid = int(words[2])
        resid = int(words[0])
        insses.add(resid)
        res2sse[resid] = sseid
fp.close()


hit_resids = set([])
if os.path.exists(prot + '.goodDomains.result'):
    fp = open(prot + '.goodDomains.result','r')
    for line in fp:
        words = line.split()
        if words[0] == 'sequence':
            segs = words[8].split(',')
        elif words[0] == 'structure':
            segs = words[14].split(',')
        for seg in segs:
            if '-' in seg:
                start = int(seg.split('-')[0])
                end = int(seg.split('-')[1])
                for resid in range(start, end + 1):
                    hit_resids.add(resid)
            else:
                resid = int(seg)
                hit_resids.add(resid)
    fp.close()


fp = open(prot.split('-model')[0]+'-predicted_aligned_error_v2.json','r')
text = fp.read()[1:-1]
fp.close()
json_dict = json.loads(text)
resid1s = json_dict['residue1']
resid2s = json_dict['residue2']
prot_len1 = max(resid1s)
prot_len2 = max(resid2s)
if prot_len1 != prot_len2:
    print ('error, matrix is not a square with shape (' + str(prot_len1) + ', ' + str(prot_len2) + ')')
else:
    length = prot_len1

allerrors = json_dict['distance']
mtx_size = len(allerrors)

rpair2error = {}
for i in range(mtx_size):
    res1 = resid1s[i]
    res2 = resid2s[i]
    err = allerrors[i]
    try:
        rpair2error[res1]
    except KeyError:
        rpair2error[res1] = {}
    rpair2error[res1][res2] = err


res2contacts = {}
for i in range(length):
    for j in range(length):
        res1 = i + 1
        res2 = j + 1
        error = rpair2error[res1][res2]
        if res1 + 10 <= res2 and error < 12:
            if res2 in insses:
                if res1 in insses and res2sse[res1] == res2sse[res2]:
                    pass
                else:
                    try:
                        res2contacts[res1].append(res2)
                    except KeyError:
                        res2contacts[res1] = [res2]
            if res1 in insses:
                if res2 in insses and res2sse[res2] == res2sse[res1]:
                    pass
                else:
                    try:
                        res2contacts[res2].append(res1)
                    except KeyError:
                        res2contacts[res2] = [res1]


diso_resids = set([])
for start in range(1, length - 9):
    total_contact = 0
    hitres_count = 0
    for res in range(start, start + 10):
        if res in hit_resids:
            hitres_count += 1
        if res in insses:
            try:
                total_contact += len(res2contacts[res])
            except KeyError:
                pass
    if total_contact <= 30 and hitres_count <= 5:
        for res in range(start, start + 10):
            diso_resids.add(res)

diso_resids_list = list(diso_resids)
diso_resids_list.sort()
rp = open(prot + '.diso','w')
for resid in diso_resids_list:
    rp.write(str(resid) + '\n')
rp.close()
with open('log','a') as f:
    f.write('10\n')
