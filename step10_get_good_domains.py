import os, sys

prefix = sys.argv[1]
wd = sys.argv[2]
data_dir = sys.argv[3]
if os.getcwd() != wd:
    os.chdir(wd)

fp = open(f'{data_dir}/ECOD_norms', 'r')
ecod2norm = {}
for line in fp:
    words = line.split()
    ecod2norm[words[0]] = float(words[1])
fp.close()

results = []
fp = open(f'{prefix}_sequence.result', 'r')
for line in fp:
    words = line.split()
    filt_segs = []
    for seg in words[6].split(','):
        start = int(seg.split('-')[0])
        end = int(seg.split('-')[1])
        for res in range(start, end + 1):
            if not filt_segs:
                filt_segs.append([res])
            else:
                if res > filt_segs[-1][-1] + 10:
                    filt_segs.append([res])
                else:
                    filt_segs[-1].append(res)

    filt_seg_strings = []
    total_good_count = 0
    for seg in filt_segs:
        start = seg[0]
        end = seg[-1]
        good_count = 0
        for res in range(start, end + 1):
            good_count += 1
        if good_count >= 5:
            total_good_count += good_count
            filt_seg_strings.append(f'{start}-{end}')
    if total_good_count >= 25:
        results.append('sequence\t' + prefix + '\t' + '\t'.join(words[:7]) + '\t' + ','.join(filt_seg_strings) + '\n')
fp.close()

if os.path.exists(f'{prefix}_structure.result'):
    fp = open(f'{prefix}_structure.result', 'r')
    for line in fp:
        words = line.split()
        ecodnum = words[0].split('_')[0]
        edomain = words[1]
        zscore = float(words[3])
        try:
            znorm = round(zscore / ecod2norm[ecodnum], 2)
        except KeyError:
            znorm = 0.0
        qscore = float(words[4])
        ztile = float(words[5])
        qtile = float(words[6])
        rank = float(words[7])
        bestprob = float(words[8])
        bestcov = float(words[9])

        judge = 0
        if rank < 1.5:
            judge += 1
        if qscore > 0.5:
            judge += 1
        if ztile < 0.75 and ztile >= 0:
            judge += 1
        if qtile < 0.75 and qtile >= 0:
            judge += 1
        if znorm > 0.225:
            judge += 1

        seqjudge = 'no'
        if bestprob >= 20 and bestcov >= 0.2:
            judge += 1
            seqjudge = 'low'
        if bestprob >= 50 and bestcov >= 0.3:
            judge += 1
            seqjudge = 'medium'
        if bestprob >= 80 and bestcov >= 0.4:
            judge += 1
            seqjudge = 'high'
        if bestprob >= 95 and bestcov >= 0.6:
            judge += 1
            seqjudge = 'superb'

        if judge:
            seg_strings = words[10].split(',')
            filt_segs = []
            for seg in words[10].split(','):
                start = int(seg.split('-')[0])
                end = int(seg.split('-')[1])
                for res in range(start, end + 1):
                    if not filt_segs:
                        filt_segs.append([res])
                    else:
                        if res > filt_segs[-1][-1] + 10:
                            filt_segs.append([res])
                        else:
                            filt_segs[-1].append(res)
            
            filt_seg_strings = []
            total_good_count = 0
            for seg in filt_segs:
                start = seg[0]
                end = seg[-1]
                good_count = 0
                for res in range(start, end + 1):
                    good_count += 1
                if good_count >= 5:
                    total_good_count += good_count
                    filt_seg_strings.append(f'{str(start)}-{str(end)}')
            if total_good_count >= 25:
                results.append('structure\t' + seqjudge + '\t' + prefix + '\t' + str(znorm) + '\t' + '\t'.join(words[:10]) + '\t' + ','.join(seg_strings) + '\t' + ','.join(filt_seg_strings) + '\n')
    fp.close()

if results:
    rp = open(f'{prefix}.goodDomains', 'w')
    for line in results:
        rp.write(line)
    rp.close()
