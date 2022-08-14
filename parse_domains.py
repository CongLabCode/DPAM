#!/usr/bin/python
import numpy as np
import os, sys, json, math, string


def get_domain_range(resids):
    segs = []
    resids.sort()
    for resid in resids:
        if not segs:
            segs.append([resid])
        else:
            if resid > segs[-1][-1] + 1:
                segs.append([resid])
            else:
                segs[-1].append(resid)
    seg_string = []
    for seg in segs:
        start = str(seg[0])
        end = str(seg[-1])
        seg_string.append(start + '-' + end)
    return ','.join(seg_string)


def get_PDB_prob(dist):
    if dist <= 3:
        prob = 0.95
    elif dist <= 6:
        prob = 0.94
    elif dist <= 9:
        prob = 0.93
    elif dist <= 12:
        prob = 0.91
    elif dist <= 15:
        prob = 0.89
    elif dist <= 18:
        prob = 0.85
    elif dist <= 21:
        prob = 0.81
    elif dist <= 24:
        prob = 0.77
    elif dist <= 27:
        prob = 0.71
    elif dist <= 30:
        prob = 0.66
    elif dist <= 35:
        prob = 0.58
    elif dist <= 40:
        prob = 0.48
    elif dist <= 45:
        prob = 0.40
    elif dist <= 50:
        prob = 0.33
    elif dist <= 55:
        prob = 0.28
    elif dist <= 60:
        prob = 0.24
    elif dist <= 70:
        prob = 0.22
    elif dist <= 80:
        prob = 0.20
    elif dist <= 120:
        prob = 0.19
    elif dist <= 160:
        prob = 0.15
    elif dist <= 200:
        prob = 0.1
    else:
        prob = 0.06
    return prob


def get_PAE_prob(error):
    if error <= 1:
        prob = 0.97
    elif error <= 2:
        prob = 0.89
    elif error <= 3:
        prob = 0.77
    elif error <= 4:
        prob = 0.67
    elif error <= 5:
        prob = 0.61
    elif error <= 6:
        prob = 0.57
    elif error <= 7:
        prob = 0.54
    elif error <= 8:
        prob = 0.52
    elif error <= 9:
        prob = 0.50
    elif error <= 10:
        prob = 0.48
    elif error <= 11:
        prob = 0.47
    elif error <= 12:
        prob = 0.45
    elif error <= 14:
        prob = 0.44
    elif error <= 16:
        prob = 0.42
    elif error <= 18:
        prob = 0.41
    elif error <= 20:
        prob = 0.39
    elif error <= 22:
        prob = 0.37
    elif error <= 24:
        prob = 0.32
    elif error <= 26:
        prob = 0.25
    elif error <= 28:
        prob = 0.16
    else:
        prob = 0.11
    return prob


def get_HHS_prob(hhpro):
    if hhpro >= 180:
        prob = 0.98
    elif hhpro >= 160:
        prob = 0.94
    elif hhpro >= 140:
        prob = 0.92
    elif hhpro >= 120:
        prob = 0.88
    elif hhpro >= 110:
        prob = 0.87
    elif hhpro >= 100:
        prob = 0.81
    elif hhpro >= 50:
        prob = 0.76
    else:
        prob = 0.5
    return prob


def get_DALI_prob(daliz):
    if daliz >= 35:
        prob = 0.95
    elif daliz >= 25:
        prob = 0.94
    elif daliz >= 20:
        prob = 0.93
    elif daliz >= 18:
        prob = 0.9
    elif daliz >= 16:
        prob = 0.87
    elif daliz >= 14:
        prob = 0.85
    elif daliz >= 12:
        prob = 0.8
    elif daliz >= 11:
        prob = 0.77
    elif daliz >= 10:
        prob = 0.74
    elif daliz >= 9:
        prob = 0.71
    elif daliz >= 8:
        prob = 0.68
    elif daliz >= 7:
        prob = 0.63
    elif daliz >= 6:
        prob = 0.60
    elif daliz >= 5:
        prob = 0.57
    elif daliz >= 4:
        prob = 0.54
    elif daliz >= 3:
        prob = 0.53
    elif daliz >= 2:
        prob = 0.52
    else:
        prob = 0.5
    return prob


prot = sys.argv[1]
#param1 = float(sys.argv[2])
#param2 = float(sys.argv[3])
param1 = 0.64
param2 = 1.1

fp = open(prot + '.fa', 'r')
protseq = ''
for line in fp:
    if line[0] != '>':
        protseq += line[:-1]
fp.close()
protlen = len(protseq)


diso_resids = set([])
fp = open(prot + '.diso', 'r')
for line in fp:
    words = line.split()
    diso_resids.add(int(words[0]))
fp.close()


fp = open(prot + '.pdb','r')
resid2coors = {}
resids = []
for line in fp:
    resid = int(line[22:26])
    coorx = float(line[30:38])
    coory = float(line[38:46])
    coorz = float(line[46:54])
    atomtype = line[13:17].strip()
    try:
        resid2coors[resid].append([coorx, coory, coorz])
    except KeyError:
        resid2coors[resid] = [[coorx, coory, coorz]]
        resids.append(resid)
fp.close()

rpair2dist = {}
for res1 in resids:
    for res2 in resids:
        if res1 < res2:
            dists = []
            for coor1 in resid2coors[res1]:
                for coor2 in resid2coors[res2]:
                    dist = ((coor1[0] - coor2[0]) ** 2 + (coor1[1] - coor2[1]) ** 2 + (coor1[2] - coor2[2]) ** 2) ** 0.5
                    dists.append(dist)
            mindist = min(dists)
            try:
                rpair2dist[res1]
            except KeyError:
                rpair2dist[res1] = {}
            try:
                rpair2dist[res2]
            except KeyError:
                rpair2dist[res2] = {}
            rpair2dist[res1][res2] = mindist
            rpair2dist[res2][res1] = mindist


fp = open(prot.split('-model')[0] + '-predicted_aligned_error_v2.json','r')
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


rpair2hhpros = {}
rpair2dalizs = {}
if os.path.exists(prot + '.goodDomains.result'):
    fp = open(prot + '.goodDomains.result','r')
    for line in fp:
        words = line.split()
        if words[0] == 'sequence':
            resids = []
            for seg in words[8].split(','):
                if '-' in seg:
                    start = int(seg.split('-')[0])
                    end = int(seg.split('-')[1])
                    for res in range(start, end + 1):
                        resids.append(res)
                else:
                    resids.append(int(seg))
            for c1, res1 in enumerate(resids):
                for c2, res2 in enumerate(resids):
                    if c1 < c2:
                        try:
                            rpair2hhpros[res1]
                        except KeyError:
                            rpair2hhpros[res1] = {}
                        try:
                            rpair2hhpros[res2]
                        except KeyError:
                            rpair2hhpros[res2] = {}
                        try:
                            rpair2hhpros[res1][res2].append(float(words[5]))
                        except KeyError:
                            rpair2hhpros[res1][res2] = [float(words[5])]
                        try:
                            rpair2hhpros[res2][res1].append(float(words[5]))
                        except KeyError:
                            rpair2hhpros[res2][res1] = [float(words[5])]
        
        elif words[0] == 'structure':
            resids = []
            for seg in words[14].split(','):
                if '-' in seg:
                    start = int(seg.split('-')[0])
                    end = int(seg.split('-')[1])
                    for res in range(start, end + 1):
                        resids.append(res)
                else:
                    resids.append(int(seg))
            for c1, res1 in enumerate(resids):
                for c2, res2 in enumerate(resids):
                    if c1 < c2:
                        try:
                            rpair2dalizs[res1]
                        except KeyError:
                            rpair2dalizs[res1] = {}
                        try:
                            rpair2dalizs[res2]
                        except KeyError:
                            rpair2dalizs[res2] = {}
                        try:
                            rpair2dalizs[res1][res2].append(float(words[7]))
                        except KeyError:
                            rpair2dalizs[res1][res2] = [float(words[7])]
                        try:
                            rpair2dalizs[res2][res1].append(float(words[7]))
                        except KeyError:
                            rpair2dalizs[res2][res1] = [float(words[7])]
    fp.close()


rpair2hhpro = {}
for res1 in rpair2hhpros.keys():
    rpair2hhpro[res1] = {}
    for res2 in rpair2hhpros[res1].keys():
        if len(rpair2hhpros[res1][res2]) > 10:
            rpair2hhpro[res1][res2] = max(rpair2hhpros[res1][res2]) + 100
        else:
            rpair2hhpro[res1][res2] = max(rpair2hhpros[res1][res2]) + len(rpair2hhpros[res1][res2]) * 10 - 10

rpair2daliz = {}
for res1 in rpair2dalizs.keys():
    rpair2daliz[res1] = {}
    for res2 in rpair2dalizs[res1].keys():
        if len(rpair2dalizs[res1][res2]) > 5:
            rpair2daliz[res1][res2] = max(rpair2dalizs[res1][res2]) + 5
        else:
            rpair2daliz[res1][res2] = max(rpair2dalizs[res1][res2]) + len(rpair2dalizs[res1][res2]) - 1


RESIDS = list(range(1, protlen + 1))
for res1 in RESIDS:
    for res2 in RESIDS:
        if res1 < res2:
            try:
                rpair2hhpro[res1]
            except KeyError:
                rpair2hhpro[res1] = {}
            try:
                rpair2hhpro[res1][res2]
            except KeyError:
                rpair2hhpro[res1][res2] = 20
            try:
                rpair2hhpro[res2]
            except KeyError:
                rpair2hhpro[res2] = {}
            try:
                rpair2hhpro[res2][res1]
            except KeyError:
                rpair2hhpro[res2][res1] = 20
            try:
                rpair2daliz[res1]
            except KeyError:
                rpair2daliz[res1] = {}
            try:
                rpair2daliz[res1][res2]
            except KeyError:
                rpair2daliz[res1][res2] = 1
            try:
                rpair2daliz[res2]
            except KeyError:
                rpair2daliz[res2] = {}
            try:
                rpair2daliz[res2][res1]
            except KeyError:
                rpair2daliz[res2][res1] = 1
rpair2prob = {}
for res1 in RESIDS:
    for res2 in RESIDS:
        if res1 < res2:
            try:
                mydist = 0.5 * (rpair2dist[res1][res2] + rpair2dist[res2][res1])
                check_dist = 1
            except KeyError:
                check_dist = 0
            try:
                myerror = 0.5 * (rpair2error[res1][res2] + rpair2error[res2][res1])
                check_error = 1
            except KeyError:
                check_error = 0
            try:
                myhhpro = 0.5 * (rpair2hhpro[res1][res2] + rpair2hhpro[res2][res1])
                check_hhpro = 1
            except KeyError:
                check_hhpro = 0
            try:
                mydaliz = 0.5 * (rpair2daliz[res1][res2] + rpair2daliz[res2][res1])
                check_daliz = 1
            except KeyError:
                check_daliz = 0
        
            if check_dist and check_error and check_hhpro and check_daliz:
                dist_prob = get_PDB_prob(mydist)
                error_prob = get_PAE_prob(myerror)
                hhpro_prob = get_HHS_prob(myhhpro)
                daliz_prob = get_DALI_prob(mydaliz)
                total_prob = (dist_prob * error_prob * hhpro_prob * daliz_prob) ** 0.25
                try:
                    rpair2prob[res1]
                except KeyError:
                    rpair2prob[res1] = {}
                try:
                    rpair2prob[res2]
                except KeyError:
                    rpair2prob[res2] = {}
                rpair2prob[res1][res2] = total_prob
                rpair2prob[res2][res1] = total_prob
            else:
                print ('error\t' + prot + '\t' + str(res1) + '\t' + str(res2))
START = min(RESIDS)
END = max(RESIDS)
mylist = list(range(START, END + 1))
chunks = []
for x in range(0, len(mylist), 5):
    chunks.append(mylist[x : x + 5])
segments_V0 = []
for chunk in chunks:
    segment = []
    for res in chunk:
        if not res in diso_resids:
            segment.append(res)
    if len(segment) >= 3:
        segments_V0.append(segment)
segment_probs = []
spair2count = {}
spair2total = {}
for i, segi in enumerate(segments_V0):
    for j, segj in enumerate(segments_V0):
        if i < j:
            try:
                spair2count[i]
            except KeyError:
                spair2count[i] = {}
            try:
                spair2total[i]
            except KeyError:
                spair2total[i] = {}
            try:
                spair2count[j]
            except KeyError:
                spair2count[j] = {}
            try:
                spair2total[j]
            except KeyError:
                spair2total[j] = {}
            spair2count[i][j] = 0
            spair2total[i][j] = 0
            spair2count[j][i] = 0
            spair2total[j][i] = 0
        
            good = 0
            total = 0
            probs = []
            for resi in segi:
                for resj in segj:
                    if resi + 5 < resj:
                        spair2count[i][j] += 1
                        spair2total[i][j] += rpair2prob[resi][resj]
                        spair2count[j][i] += 1
                        spair2total[j][i] += rpair2prob[resi][resj]
                        probs.append(rpair2prob[resi][resj])
            if probs:
                meanprob = np.mean(probs)
                if meanprob > param1:
                    segment_probs.append([i, j, meanprob])
segment_probs.sort(key = lambda x:x[2], reverse = True)
segments_V1 = []
for item in segment_probs:
    segi = item[0]
    segj = item[1]
    if not segments_V1:
        segments_V1.append(set([segi, segj]))
    else:
        isdone = 0
        candidates = []
        for counts, segment in enumerate(segments_V1):
            if segi in segment and segj in segment:
                isdone = 1
            elif segi in segment:
                candidates.append(counts)
            elif segj in segment:
                candidates.append(counts)
        if not isdone:
            if len(candidates) == 2:
                c1 = candidates[0]
                c2 = candidates[1]
                intra_count1 = 0
                intra_total1 = 0
                intra_count2 = 0
                intra_total2 = 0
                inter_count = 0
                inter_total = 0
                for i in segments_V1[c1]:
                    for j in segments_V1[c1]:
                        if i < j:
                            intra_count1 += spair2count[i][j]
                            intra_total1 += spair2total[i][j]
                for i in segments_V1[c2]:
                    for j in segments_V1[c2]:
                        if i < j:
                            intra_count2 += spair2count[i][j]
                            intra_total2 += spair2total[i][j]
                for i in segments_V1[c1]:
                    for j in segments_V1[c2]:
                        inter_count += spair2count[i][j]
                        inter_total += spair2total[i][j]
                merge = 0
                if intra_count1 <= 20 or intra_count2 <= 20:
                    merge = 1
                else:
                    intra_prob1 = intra_total1 / intra_count1
                    intra_prob2 = intra_total2 / intra_count2
                    inter_prob = inter_total / inter_count
                    if inter_prob * param2 >= intra_prob1 or inter_prob * param2 >= intra_prob2:
                        merge = 1
                if merge:
                    new_segments = []
                    new_segment = set([])
                    for counts, segment in enumerate(segments_V1):
                        if counts in candidates:
                            new_segment = new_segment.union(segment)
                        else:
                            new_segments.append(segment)
                    new_segments.append(new_segment)
                    segments_V1 = new_segments
                else:
                    pass
            elif len(candidates) == 1:
                c0 = candidates[0]
                intra_count = 0
                intra_total = 0
                inter_count = 0
                inter_total = 0
                for i in segments_V1[c0]:
                    for j in segments_V1[c0]:
                        if i < j:
                            intra_count += spair2count[i][j]
                            intra_total += spair2total[i][j]    
                if segi in segments_V1[c0]:
                    for k in segments_V1[c0]:
                        if segj != k:
                            inter_count += spair2count[k][segj]
                            inter_total += spair2total[k][segj]
                elif segj in segments_V1[c0]:
                    for k in segments_V1[c0]:
                        if segi != k:
                            inter_count += spair2count[k][segi]
                            inter_total += spair2total[k][segi]
                else:
                    print ('error0')
                merge = 0
                if intra_total <= 20:
                    merge = 1
                else:
                    intra_prob = intra_total / intra_count
                    inter_prob = inter_total / inter_count
                    if inter_prob * param2 >= intra_prob:
                        merge = 1
                if merge:
                    segments_V1[c0].add(segi)
                    segments_V1[c0].add(segj)
                else:
                    pass
            elif len(candidates) == 0:
                segments_V1.append(set([segi, segj]))
            else:
                print ('error1')
sorted_segments = []
for item in segments_V1:
    resids = []
    for segind in item:
        for res in segments_V0[segind]:
            resids.append(res)
    resids.sort()
    if resids:
        sorted_segments.append([resids, np.mean(resids)])
sorted_segments.sort(key = lambda x:x[1])
domains_v0 = []
domain_resids = set([])
for item in sorted_segments:
    if len(item[0]) >= 20:
        domain_range = get_domain_range(item[0])
        for res in item[0]:
            domain_resids.add(res)
        domains_v0.append([item[0], domain_range])
#print ("****** v0 *** " + str(len(domains_v0)) + " ******")
#for item in domains_v0:
#    print (item[0], item[1])
domains_v1 = []
for item in domains_v0:
    domain = item[0]
    if len(domain) >= 20:
        segs = []
        for res in domain:
            if not segs:
                segs.append([res])
            elif res == segs[-1][-1] + 1:
                segs[-1].append(res)
            else:
                segs.append([res])
        if len(segs) > 1:
            newdomain = []
            for counts, seg in enumerate(segs):
                if not counts:
                    for residue in seg:
                        newdomain.append(residue)
                else:
                    lastseg = segs[counts - 1]
                    count_other = 0
                    count_good = 0
                    count_all = 0
                    interseg = range(lastseg[-1] + 1, seg[0])
                    for residue in interseg:
                        if residue in domain_resids:
                            count_other += 1
                        else:
                            count_good += 1
                        count_all += 1
                    getit = 0
                    if count_all <= 10:
                        getit = 1
                    elif count_all <= 20 and count_other <= 10:
                        getit = 1
                    if getit:
                        for residue in interseg:
                            newdomain.append(residue)
                    for residue in seg:
                        newdomain.append(residue)
            domain_range = get_domain_range(newdomain)
            domains_v1.append([newdomain, domain_range])
        else:
            domain_range = get_domain_range(segs[0])
            domains_v1.append([segs[0], domain_range])
#print ("****** v1 *** " + str(len(domains_v1)) + " ******")
#for item in domains_v1:
#    print (item[0], item[1])
domains_v2 = []
for counti, itemi in enumerate(domains_v1):
    domain = itemi[0]
    other_resids = set([])
    for countj, itemj in enumerate(domains_v1):
        if counti != countj:
            for res in itemj[0]:
                other_resids.add(res)
    segs = []
    for res in domain:
        if not segs:
            segs.append([res])
        elif res == segs[-1][-1] + 1:
            segs[-1].append(res)
        else:
            segs.append([res])
    newdomain = []
    for counts, seg in enumerate(segs):
        count_good = 0
        for resid in seg:
            if not resid in other_resids:
                count_good += 1
        if count_good >= 10:
            for resid in seg:
                newdomain.append(resid)
    if len(newdomain) >= 25:
        domain_range = get_domain_range(newdomain)
        domains_v2.append([newdomain, domain_range])
#print ("****** v2 *** " + str(len(domains_v2)) + " ******")
#for item in domains_v2:
#    print (item[0], item[1])
rp = open(prot + '.domains', 'w')
for countd, domain in enumerate(domains_v2):
    rp.write('D' + str(countd + 1) + '\t' + domain[1] + '\n')
rp.close()
with open('log','a') as f:
    f.write('11\n')
