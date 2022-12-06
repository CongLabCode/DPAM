import os, sys

def get_range(resids):
	resids.sort()
	segs = []
	for resid in resids:
		if not segs:
			segs.append([resid])
		else:
			if resid > segs[-1][-1] + 1:
				segs.append([resid])
			else:
				segs[-1].append(resid)
	ranges = []
	for seg in segs:
		ranges.append(f'{seg[0]}-{seg[-1]}')
	return ','.join(ranges)

prefix = sys.argv[1]
wd = sys.argv[2]
data_dir = sys.argv[3]
if os.getcwd() != wd:
    os.chdir(wd)

fp = open(f'{data_dir}/ECOD_length','r')
ecod2len = {}
for line in fp:
	words = line.split()
	ecod2len[words[0]] = int(words[2])
fp.close()

fp = open(f'{data_dir}/ecod.latest.domains','r')
ecod2id = {}
ecod2fam = {}
for line in fp:
	if line[0] != '#':
		words = line[:-1].split('\t')
		ecodnum = words[0]
		ecodid = words[1]
		ecodfam = '.'.join(words[3].split('.')[:2])
		ecod2id[ecodnum] = ecodid
		ecod2fam[ecodnum] = ecodfam
fp.close()


seqhits = []
fp = open(f'{prefix}.map2ecod.result', 'r')
for countl, line in enumerate(fp):
	if countl:
		words = line.split()
		ecodnum = words[0]
		ecodlen = ecod2len[ecodnum]
		ecodfam = ecod2fam[ecodnum]
		prob = float(words[2])
		Qsegs = words[11].split(',')
		Tsegs = words[12].split(',')
		Qresids = []
		for seg in Qsegs:
			start = int(seg.split('-')[0])
			end = int(seg.split('-')[1])
			for res in range(start, end + 1):
				Qresids.append(res)
		Tresids = []
		for seg in Tsegs:
			start = int(seg.split('-')[0])
			end = int(seg.split('-')[1])
			for res in range(start, end + 1):
				Tresids.append(res)
		seqhits.append([ecodnum, ecodlen, ecodfam, prob, Qresids, Tresids])
fp.close()


fam2hits = {}
fams = set([])
for hit in seqhits:
	fam = hit[2]
	fams.add(fam)
	try:
		fam2hits[fam].append([hit[3], hit[1], hit[4], hit[5]])
	except KeyError:
		fam2hits[fam] = [[hit[3], hit[1], hit[4], hit[5]]]

ecods = []
ecod2hits = {}
for hit in seqhits:
	ecodnum = hit[0]
	ecodlen = hit[1]
	ecodfam = hit[2]
	prob = hit[3]
	Qresids = hit[4]
	Tresids = hit[5]
	Qset = set(Qresids)
	try:
		ecod2hits[ecodnum].append([prob, ecodfam, ecodlen, Qresids, Tresids, Qset])
	except KeyError:
		ecod2hits[ecodnum] = [[prob, ecodfam, ecodlen, Qresids, Tresids, Qset]]
		ecods.append(ecodnum)


rp = open(f'{prefix}_sequence.result', 'w')
for ecodnum in ecods:
	ecodid = ecod2id[ecodnum]
	ecod2hits[ecodnum].sort(key = lambda x:x[0], reverse = True)
	get_resids = set([])
	mycount = 0
	for hit in ecod2hits[ecodnum]:
		hit_prob = hit[0]
		hit_fam = hit[1]
		hit_ecodlen = hit[2]
		query_resids = hit[3]
		query_range = get_range(query_resids)
		hit_resids = hit[4]
		hit_range = get_range(hit_resids)
		hit_resids_set = hit[5]
		hit_coverage = round(len(hit_resids_set) / hit_ecodlen, 2)
		if hit_coverage >= 0.4 and hit_prob >= 50:
			new_resids = hit_resids_set.difference(get_resids)
			if len(new_resids) >= len(hit_resids_set) * 0.5:
				mycount += 1
				get_resids = get_resids.union(hit_resids_set)
				rp.write(f'{ecodnum}_{str(mycount)}\t{ecodid}\t{hit_fam}\t{hit_prob}\t{hit_coverage}\t{hit_ecodlen}\t{query_range}\t{hit_range}\n')
rp.close()


if os.path.exists(f'{prefix}_good_hits'):
	fp = open(f'{prefix}_good_hits', 'r')
	rp = open(f'{prefix}_structure.result', 'w')
	for countl, line in enumerate(fp):
		if countl:
			words = line.split()
			hitname = words[0]
			ecodnum = words[1]
			ecodid = words[2]
			ecodfam = words[3]
			zscore = words[4]
			qscore = words[5]
			ztile = words[6]
			qtile = words[7]
			rank = words[8]
			qsegments = words[9]
			ssegments = words[10]
			segs = []
			for seg in qsegments.split(','):
				start = int(seg.split('-')[0])
				end = int(seg.split('-')[1])
				for res in range(start, end + 1):
					if not segs:
						segs.append([res])
					else:
						if res > segs[-1][-1] + 10:
							segs.append([res])
						else:
							segs[-1].append(res)
			resids = set([])
			for seg in segs:
				start = seg[0]
				end = seg[-1]
				for res in range(start, end + 1):
					resids.add(res)

			good_hits = []
			try:
				for hit in fam2hits[ecodfam]:
					prob = float(hit[0])
					Tlen = hit[1]
					Qresids = hit[2]
					Tresids = hit[3]
					get_Tresids = set([])
					for i in range(len(Qresids)):
						if Qresids[i] in resids:
							get_Tresids.add(Tresids[i])
					Tcov = len(get_Tresids) / Tlen
					good_hits.append([prob, Tcov])
			except KeyError:
				pass

			bestprob = 0
			bestcov = 0
			if good_hits:
				for item in good_hits:
					if item[0] > bestprob:
						bestprob = item[0]
				bestcovs = []
				for item in good_hits:
					if item[0] >= bestprob - 0.1:
						bestcovs.append(item[1])
				bestcov = round(max(bestcovs), 2)
			rp.write(f'{hitname}\t{ecodid}\t{ecodfam}\t{zscore}\t{qscore}\t{ztile}\t{qtile}\t{rank}\t{bestprob}\t{bestcov}\t{qsegments}\t{ssegments}\n')
	fp.close()
	rp.close()
