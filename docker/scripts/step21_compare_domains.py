#!/opt/conda/bin/python
import os, sys

def get_seq_dist(residsA, residsB, good_resids):
    indsA = []
    for ind, resid in enumerate(good_resids):
        if resid in residsA:
            indsA.append(ind)
    indsB = []
    for ind, resid in enumerate(good_resids):
        if resid in residsB:
            indsB.append(ind)

    connected = 0
    for indA in indsA:
        for indB in indsB:
            if abs(indA - indB) <= 5:
                connected = 1
                break
        if connected:
            break
    return connected

dataset = sys.argv[1]
part = sys.argv[2]
fp = open('step21_' + dataset + '_' + part + '.list','r')
cases = []
for line in fp:
    words = line.split()
    cases.append(words)
fp.close()

rp = open('step21_' + dataset + '_' + part + '.result', 'w')
for case in cases:
    prot = case[0]
    good_resids = []
    fp = open('step14/' + dataset + '/' + prot + '.domains','r')
    for line in fp:
        words = line.split()
        for seg in words[1].split(','):
            if '-' in seg:
                start = int(seg.split('-')[0])
                end = int(seg.split('-')[1])
                for res in range(start, end + 1):
                    good_resids.append(res)
            else:
                good_resids.append(int(seg))
    fp.close()
    good_resids.sort()

    dom1 = case[1]
    segs1 = case[2]
    residsA = []
    for seg in segs1.split(','):
        if '-' in seg:
            start = int(seg.split('-')[0])
            end = int(seg.split('-')[1])
            for res in range(start, end + 1):
                residsA.append(res)
        else:
            residsA.append(int(seg))

    dom2 = case[3]
    segs2 = case[4]
    residsB = []
    for seg in segs2.split(','):
        if '-' in seg:
            start = int(seg.split('-')[0])
            end = int(seg.split('-')[1])
            for res in range(start, end + 1):
                residsB.append(res)
        else:
            residsB.append(int(seg))

    if get_seq_dist(set(residsA), set(residsB), good_resids):
        judge = 1
    else:
        resid2coors = {}
        fp = open('step20/' + dataset + '/' + prot + '_' + dom1 + '.pdb', 'r')
        for line in fp:
            resid = int(line[22:26])
            coorx = float(line[30:38])
            coory = float(line[38:46])
            coorz = float(line[46:54])
            try:
                resid2coors[resid].append([coorx, coory, coorz])
            except KeyError:
                resid2coors[resid] = [[coorx, coory, coorz]]
        fp.close()

        fp = open('step20/' + dataset + '/' + prot + '_' + dom2 + '.pdb', 'r')
        for line in fp:
            resid = int(line[22:26])
            coorx = float(line[30:38])
            coory = float(line[38:46])
            coorz = float(line[46:54])
            try:
                resid2coors[resid].append([coorx, coory, coorz])
            except KeyError:
                resid2coors[resid] = [[coorx, coory, coorz]]
        fp.close()

        interface_count = 0
        for residA in residsA:
            for residB in residsB:
                dists = []
                coorsA = resid2coors[residA]
                coorsB = resid2coors[residB]
                for coorA in coorsA:
                    for coorB in coorsB:
                        dist = ((coorA[0] - coorB[0]) ** 2 + (coorA[1] - coorB[1]) ** 2 + (coorA[2] - coorB[2]) ** 2) ** 0.5
                        dists.append(dist)
                min_dist = min(dists)
                if min_dist <= 8:
                    interface_count += 1
        if interface_count >= 9:
            judge = 2
        else:
            judge = 0
    rp.write(prot + '\t' + dom1 + '\t' + dom2 + '\t' + str(judge) + '\t' + segs1 + '\t' + segs2 + '\n')
rp.close()
