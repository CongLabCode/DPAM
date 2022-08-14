import os, sys
import numpy as np

prot = sys.argv[1]
os.system('mkdssp '+prot + '.pdb > '+prot + '.dssp')

fp = open(prot + '.fa', 'r')
for line in fp:
    if line[0] != '>':
        seq = line[:-1]
fp.close()

fp = open(prot + '.dssp', 'r')
start = 0
dssp_result = ""
resids = []
for line in fp:
    words = line.split()
    if len(words) > 3:
        if words[0] == "#" and words[1] == "RESIDUE":
            start = 1
        elif start:
            try:
                resid = int(line[5:10])
                getit = 1
            except ValueError:
                getit = 0

            if getit:
                pred = line[16]
                resids.append(resid)
                pred = line[16]
                if pred == "E" or pred == "B":
                    newpred = "E"
                elif pred == "G" or pred == "H" or pred == "I":
                    newpred = "H"
                else:
                    newpred = "-"
                dssp_result += newpred
fp.close()

res2sse = {}
dssp_segs = dssp_result.split("---")
posi = 0
Nsse = 0
for dssp_seg in dssp_segs:
    if dssp_seg.count("E") >= 3 or dssp_seg.count("H") >= 6:
        Nsse += 1
    for char in dssp_seg:
        resid = resids[posi]
        if char != "-":
            res2sse[resid] = [Nsse, char]
        posi += 1
    posi += 3

if len(resids) != len(seq):
    print ('error\t' + prot + '\t' + str(len(resids)) + '\t' + str(len(seq)))
else:
    rp = open(prot + '.sse', 'w')
    for resid in resids:
        try:
            rp.write(str(resid) + '\t' + seq[resid - 1] + '\t' + str(res2sse[resid][0]) + '\t' + res2sse[resid][1] + '\n')
        except KeyError:
            rp.write(str(resid) + '\t' + seq[resid - 1] + '\tna\tC\n')
    rp.close()
with open('log','a') as f:
    f.write('9\n')
