import os, sys
import numpy as np
script_dir=os.path.dirname(os.path.realpath(__file__))
with open(script_dir+'/config_file') as f:
    configs=f.readlines()
configs=[i.strip() for i in configs if i.strip()!='']
configs={i.split(':')[0]:i.split(':')[1] for i in configs}
ECOD_resWeight=configs['ECOD_resWeight']
ECOD_domain_info=configs['ECOD_domain_info']

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
        ranges.append(str(seg[0]) + '-' + str(seg[-1]))
    return ','.join(ranges)


fp = open(script_dir+'/ecod.latest.domains','r')
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


prot = sys.argv[1]
fp = open(prot + '_iterative_hits','r')
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
    posi2weight = {}
    zscores = []
    qscores = []
    if os.path.exists(f'{ECOD_resWeight}/{ecodnum}.weight'):   #ECOD step18 weight
        fp = open(f'{ECOD_resWeight}/{ecodnum}.weight','r')
        posi2weight = {}
        for line in fp:
            words = line.split()
            posi2weight[int(words[0])] = float(words[3])
        fp.close()
    if os.path.exists(f'{ECOD_domain_info}/{ecodnum}.info'):
        fp = open(f'{ECOD_doamin_info}/{ecodnum}.info','r')   #ECOD step19 info
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
        newhits.append([hitname, ecodnum, ecodid, ecodfam, zscore, qscore, ztile, qtile, maps])
    else:
        newhits.append([hitname, ecodnum, ecodid, ecodfam, zscore, 0, 1, 1, maps])


newhits.sort(key = lambda x:x[4], reverse = True)
finalhits = []
posi2fams = {}
for hit in newhits:
    ecodfam = hit[3]
    maps = hit[8]
    qposis = []
    ranks = []
    for item in maps:
        qposis.append(item[0])
        try:
            posi2fams[item[0]].add(ecodfam)
        except KeyError:
            posi2fams[item[0]] = set([ecodfam])
        ranks.append(len(posi2fams[item[0]]))
    ave_rank = round(np.mean(ranks), 2)
    qrange = get_range(qposis)
    finalhits.append([hit[0], hit[1], hit[2], hit[3], round(hit[4], 2), round(hit[5], 2), round(hit[6], 2), round(hit[7], 2), ave_rank, qrange])

rp = open(prot + '.ecod_analysis.result', 'w')
for hit in finalhits:
    rp.write(hit[0] + '\t' + hit[1] + '\t' + hit[2] + '\t' + hit[3] + '\t' + str(hit[4]) + '\t' + str(hit[5]) + '\t' + str(hit[6]) + '\t' + str(hit[7]) + '\t' + str(hit[8]) + '\t' + hit[9] + '\n')
rp.close()
with open('log','a') as f:
    f.write('6\n')
