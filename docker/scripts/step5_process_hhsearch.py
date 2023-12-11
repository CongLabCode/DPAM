#!/opt/conda/bin/python
import os, sys

def get_range(resids, chainid):
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
        if chainid:
            ranges.append(chainid + ':' + str(seg[0]) + '-' + str(seg[-1]))
        else:
            ranges.append(str(seg[0]) + '-' + str(seg[-1]))
    return ','.join(ranges)

dataset = sys.argv[1]
prot = sys.argv[2]
fp = open('step3/' + dataset + '/' + prot + '.hhsearch', 'r')
info = fp.read().split('\n>')
fp.close()
allhits = []
need_pdbchains = set([])
need_pdbs = set([])
for hit in info[1:]:
    lines = hit.split('\n')
    qstart = 0
    qend = 0
    qseq = ''
    hstart = 0
    hend = 0
    hseq = ''
    for line in lines:
        if len(line) >= 6:
            if line[:6] == 'Probab':
                words = line.split()
                for word in words:
                    subwords = word.split('=')
                    if subwords[0] == 'Probab':
                        hh_prob = subwords[1]
                    elif subwords[0] == 'E-value':
                        hh_eval = subwords[1]
                    elif subwords[0] == 'Score':
                        hh_score = subwords[1]
                    elif subwords[0] == 'Aligned_cols':
                        aligned_cols = subwords[1]
                    elif subwords[0] == 'Identities':
                        idents = subwords[1]
                    elif subwords[0] == 'Similarity':
                        similarities = subwords[1]
                    elif subwords[0] == 'Sum_probs':
                        sum_probs = subwords[1]

            elif line[:2] == 'Q ':
                words = line.split()
                if words[1] != 'ss_pred' and words[1] != 'Consensus':
                    qseq += words[3]
                    if not qstart:
                        qstart = int(words[2])
                    qend = int(words[4])

            elif line[:2] == 'T ':
                words = line.split()
                if words[1] != 'Consensus' and words[1] != 'ss_dssp' and words[1] != 'ss_pred':
                    hid = words[1]
                    hseq += words[3]
                    if not hstart:
                        hstart = int(words[2])
                    hend = int(words[4])
    allhits.append([hid, hh_prob, hh_eval, hh_score, aligned_cols, idents, similarities, sum_probs, qstart, qend, qseq, hstart, hend, hseq])
    need_pdbchains.add(hid)
    need_pdbs.add(hid.split('_')[0].lower())


fp = open('/mnt/databases/ECOD_pdbmap','r')
pdb2ecod = {}
good_hids = set([])
for line in fp:
    words = line.split()
    pdbid = words[1]
    segments = words[2].split(',')
    chainids = set([])
    resids = []
    for segment in segments:
        chainids.add(segment.split(':')[0])
        if '-' in segment:
            start = int(segment.split(':')[1].split('-')[0])
            end = int(segment.split(':')[1].split('-')[1])
            for res in range(start, end + 1):
                resids.append(res)
        else:
            resid = int(segment.split(':')[1])
            resids.append(resid)
    if len(chainids) == 1:
        chainid = list(chainids)[0]
        pdbchain = pdbid.upper() + '_' + chainid
        if pdbchain in need_pdbchains:
            good_hids.add(pdbchain)
            try:
                pdb2ecod[pdbchain]
            except KeyError:
                pdb2ecod[pdbchain] = {}
            for i, resid in enumerate(resids):
                pdb2ecod[pdbchain][resid] = words[0] + ':' + str(i + 1)
    else:
        print (line[:-1])
fp.close()

ecod2key = {}
ecod2len = {}
fp = open('/mnt/databases/ECOD_length','r')
for line in fp:
    words = line.split()
    ecod2key[words[0]] = words[1]
    ecod2len[words[0]] = int(words[2])
fp.close()

if allhits:
    results = []
    hid2count = {}
    for hit in allhits:
        hid = hit[0]
        try:
            hid2count[hid] += 1
        except KeyError:
            hid2count[hid] = 1
        hitname = hid + '_' + str(hid2count[hid])
        pdbid = hid.split('_')[0]
        chainid = hid.split('_')[1]
        ecods = []
        ecod2hres = {}
        ecod2hresmap = {}
        if hid in good_hids:
            for pdbres in pdb2ecod[hid].keys():
                for item in pdb2ecod[hid][pdbres].split(','):
                    ecod = item.split(':')[0]
                    ecodres = int(item.split(':')[1])
                    try:
                        ecod2hres[ecod]
                        ecod2hresmap[ecod]
                    except KeyError:
                        ecods.append(ecod)
                        ecod2hres[ecod] = set([])
                        ecod2hresmap[ecod] = {}
                    ecod2hres[ecod].add(pdbres)
                    ecod2hresmap[ecod][pdbres] = ecodres

            hh_prob = hit[1]
            hh_eval = hit[2]
            hh_score = hit[3]
            aligned_cols = hit[4]
            idents = hit[5]
            similarities = hit[6]
            sum_probs = hit[7]
            qstart = hit[8]
            qseq = hit[10]
            hstart = hit[11]
            hseq = hit[13]

            for ecod in ecods:
                ecodkey = ecod2key[ecod]
                ecodlen = ecod2len[ecod]
                qposi = qstart - 1
                hposi = hstart - 1
                qresids = []
                hresids = []
                eresids = []
                if len(qseq) == len(hseq):
                    for i in range(len(hseq)):
                        if qseq[i] != '-':
                            qposi += 1
                        if hseq[i] != '-':
                            hposi += 1
                        if qseq[i] != '-' and hseq[i] != '-':
                            if hposi in ecod2hres[ecod]:
                                eposi = ecod2hresmap[ecod][hposi]
                                qresids.append(qposi)
                                hresids.append(hposi)
                                eresids.append(eposi)
                    if len(qresids) >= 10 and len(eresids) >= 10:
                        qrange = get_range(qresids,'')
                        hrange = get_range(hresids, chainid)
                        erange = get_range(eresids,'')
                        coverage = round(len(eresids) / ecodlen, 3)
                        ungapped_coverage = round((max(eresids) - min(eresids) + 1) / ecodlen, 3)
                        results.append(hitname + '\t' + ecod + '\t' + ecodkey + '\t' + hh_prob + '\t' + hh_eval + '\t' + hh_score + '\t' + aligned_cols + '\t' + idents + '\t' + similarities + '\t' + sum_probs + '\t' + str(coverage) + '\t' + str(ungapped_coverage) + '\t' + qrange + '\t' + erange + '\t' + hrange + '\n')
                else:
                    print ('error\t' + prot + '\t' + ecod)

    if results:
        rp = open('step5/' + dataset + '/' + prot + '.result', 'w')
        rp.write('hitname\tecodnum\tecodkey\thh_prob\thh_eval\thh_score\taligned_cols\tidents\tsimilarities\tsum_probs\tcoverage\tungapped_coverage\tquery_range\ttemplate_range\ttemplate_seqid_range\n')
        for line in results:
            rp.write(line)
        rp.close()
    else:
        os.system('echo \'done\' > step5/' + dataset + '/' + prot + '.done')
else:
    os.system('echo \'done\' > step5/' + dataset + '/' + prot + '.done')
