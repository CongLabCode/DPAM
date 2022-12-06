import os, sys
import numpy as np


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

if os.path.exists(f'{prefix}_iterativdDali_hits'):
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

    fp = open(f'{prefix}_iterativdDali_hits','r')
    ecodnum = ''
    ecodid = ''
    ecodfam = ''
    hitname = ''
    maps = []
    hits = []
    for line in fp:
        if line[0] == '>':
            if ecodnum and ecodid and ecodfam and hitname and zscore and maps:
                hits.append([hitname, ecodnum, ecodid, ecodfam, zscore, maps])
            words = line[1:].split()
            zscore = float(words[1])
            hitname = words[0]
            ecodnum = hitname.split('_')[0]
            ecodid = ecod2id[ecodnum]
            ecodfam = ecod2fam[ecodnum]
            maps = []
        else:
            words = line.split()
            pres = int(words[0])
            eres = int(words[1])
            maps.append([pres, eres])
    fp.close()
    if ecodnum and ecodid and ecodfam and hitname and zscore and maps:
        hits.append([hitname, ecodnum, ecodid, ecodfam, zscore, maps])


    newhits = []
    for hit in hits:
        hitname = hit[0]
        ecodnum = hit[1]
        total_weight = 0
        posi2weight = {}
        zscores = []
        qscores = []
        if os.path.exists(f'{data_dir}/ecod_weights/{ecodnum}.weight'):
            fp = open(f'{data_dir}/ecod_weights/{ecodnum}.weight','r')
            posi2weight = {}
            for line in fp:
                words = line.split()
                total_weight += float(words[3])
                posi2weight[int(words[0])] = float(words[3])
            fp.close()
        if os.path.exists(f'{data_dir}/ecod_domain_info/{ecodnum}.info'):
            fp = open(f'{data_dir}/ecod_domain_info/{ecodnum}.info','r')
            for line in fp:
                words = line.split()
                zscores.append(float(words[1]))
                qscores.append(float(words[2]))
            fp.close()
        ecodid = hit[2]
        ecodfam = hit[3]
        zscore = hit[4]
        maps = hit[5]

        if zscores and qscores:
            qscore = 0
            for item in maps:
                try:
                    qscore += posi2weight[item[1]]
                except KeyError:
                    pass

            better = 0
            worse = 0
            for other_qscore in qscores:
                if other_qscore > qscore:
                    better += 1
                else:
                    worse += 1
            qtile = better / (better + worse)

            better = 0
            worse = 0
            for other_zscore in zscores:
                if other_zscore > zscore:
                    better += 1
                else:
                    worse += 1
            ztile = better / (better + worse)   
            newhits.append([hitname, ecodnum, ecodid, ecodfam, zscore, qscore / total_weight, ztile, qtile, maps])
        else:
            newhits.append([hitname, ecodnum, ecodid, ecodfam, zscore, -1, -1, -1, maps])


    newhits.sort(key = lambda x:x[4], reverse = True)
    finalhits = []
    posi2fams = {}
    for hit in newhits:
        ecodfam = hit[3]
        maps = hit[8]
        qposis = []
        eposis = []
        ranks = []
        for item in maps:
            qposis.append(item[0])
            eposis.append(item[1])
            try:
                posi2fams[item[0]].add(ecodfam)
            except KeyError:
                posi2fams[item[0]] = set([ecodfam])
            ranks.append(len(posi2fams[item[0]]))
        ave_rank = round(np.mean(ranks), 2)
        qrange = get_range(qposis)
        erange = get_range(eposis)
        finalhits.append([hit[0], hit[1], hit[2], hit[3], round(hit[4], 2), round(hit[5], 2), round(hit[6], 2), round(hit[7], 2), ave_rank, qrange, erange])

    rp = open(f'{prefix}_good_hits', 'w')
    rp.write('hitname\tecodnum\tecodkey\thgroup\tzscore\tqscore\tztile\tqtile\trank\tqrange\terange\n')
    for hit in finalhits:
        rp.write(f'{hit[0]}\t{hit[1]}\t{hit[2]}\t{hit[3]}\t{hit[4]}\t{hit[5]}\t{hit[6]}\t{hit[7]}\t{hit[8]}\t{hit[9]}\t{hit[10]}\n')
    rp.close()
