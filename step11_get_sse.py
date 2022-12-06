import os, sys
import numpy as np

prefix = sys.argv[1]
wd = sys.argv[2]
if os.getcwd() != wd:
    os.chdir(wd)

os.system(f'mkdssp -i {prefix}.pdb -o {prefix}.dssp')
fp = open(f'{prefix}.fa', 'r')
for line in fp:
	if line[0] != '>':
		seq = line[:-1]
fp.close()

fp = open(f'{prefix}.dssp', 'r')
start = 0
dssp_result = ''
resids = []
for line in fp:
	words = line.split()
	if len(words) > 3:
		if words[0] == '#' and words[1] == 'RESIDUE':
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
				if pred == 'E' or pred == 'B':
					newpred = 'E'
				elif pred == 'G' or pred == 'H' or pred == 'I':
					newpred = 'H'
				else:
					newpred = '-'
				dssp_result += newpred
fp.close()

res2sse = {}
dssp_segs = dssp_result.split('--')
posi = 0
Nsse = 0
for dssp_seg in dssp_segs:
	judge = 0
	if dssp_seg.count('E') >= 3 or dssp_seg.count('H') >= 6:
		Nsse += 1
		judge = 1
	for char in dssp_seg:
		resid = resids[posi]
		if char != '-':
			if judge:
				res2sse[resid] = [Nsse, char]
		posi += 1
	posi += 2

os.system(f'rm {prefix}.dssp')
if len(resids) != len(seq):
	print (f'error\t{prefix}\t{len(resids)}\t{len(seq)}')
else:
	rp = open(f'{prefix}.sse', 'w')
	for resid in resids:
		try:
			rp.write(f'{resid}\t{seq[resid - 1]}\t{res2sse[resid][0]}\t{res2sse[resid][1]}\n')
		except KeyError:
			rp.write(f'{resid}\t{seq[resid - 1]}\tna\tC\n')
	rp.close()
