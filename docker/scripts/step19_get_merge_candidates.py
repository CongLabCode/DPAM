#!/opt/conda/bin/python
import os, sys

def get_resids(domain_range):
    domain_resids = []
    for seg in domain_range.split(','):
        if '-' in seg:
            start = int(seg.split('-')[0])
            end = int(seg.split('-')[1])
            for res in range(start, end + 1):
                domain_resids.append(res)
        else:
            domain_resids.append(int(seg))
    return domain_resids


spname = sys.argv[1]
prot = sys.argv[2]
if os.path.exists('step18/' + spname + '/' + prot + '.data'):
    fp = open('step18/' + spname + '/' + prot + '.data','r')
    need_ecods = set([])
    for line in fp:
        words = line.split()
        need_ecods.add(words[2])
    fp.close()

    fp = open('/mnt/databases/ECOD_length','r')
    ecod2length = {}
    for line in fp:
        words = line.split()
        if words[0] in need_ecods:
            ecod2length[words[0]] = int(words[2])
    fp.close()

    ecod2totW = {}
    ecod2posW = {}
    for ecod in need_ecods:
        ecod2totW[ecod] = 0
        ecod2posW[ecod] = {}
        if os.path.exists('/mnt/databases/posi_weights/' + ecod + '.weight'):
            fp = open('/mnt/databases/posi_weights/' + ecod + '.weight','r')
            for line in fp:
                words = line.split()
                resid = int(words[0])
                weight = float(words[3])
                ecod2totW[ecod] += weight
                ecod2posW[ecod][resid] = weight
            fp.close()
        else:
            ecod2totW[ecod] = ecod2length[ecod]
            for i in range(ecod2length[ecod]):
                ecod2posW[ecod][i + 1] = 1


    fp = open('step18/' + spname + '/' + prot + '.data','r')
    domains = []
    ecods = []
    domain2def = {}
    domain2prob = {}
    domain2hits = {}
    ecod2hits = {}
    for line in fp:
        words = line.split()
        domain = words[0]
        domdef = words[1]
        ecod = words[2]
        tgroup = words[3]
        prob = float(words[4])
        try:
            if prob > domain2prob[domain]:
                domain2prob[domain] = prob
        except KeyError:
            domain2def[domain] = domdef
            domain2prob[domain] = prob

        if words[6] == 'na':
            Hresids = []
        else:
            Hresids = get_resids(words[6])
        if words[7] == 'na':
            Dresids = []
        else:
            Dresids = get_resids(words[7])
        if len(Dresids) > len(Hresids) * 0.5:
            HDresids = set(Dresids)
        else:
            HDresids = set(Hresids)
    
        total_weight = ecod2totW[ecod]
        get_weight = 0
        for resid in HDresids:
            try:
                get_weight += ecod2posW[ecod][resid]
            except KeyError:
                print (prot, ecod, resid)
        try:
            domain2hits[domain].append([ecod, tgroup, prob, get_weight / total_weight])
        except KeyError:
            domains.append(domain)
            domain2hits[domain] = [[ecod, tgroup, prob, get_weight / total_weight]]
        try:
            ecod2hits[ecod].append([domain, tgroup, prob, HDresids])
        except KeyError:
            ecods.append(ecod)
            ecod2hits[ecod] = [[domain, tgroup, prob, HDresids]]
    fp.close()


    domain_pairs = []
    dpair2supports = {}
    for ecod in ecods:
        if len(ecod2hits[ecod]) > 1:
            for c1, hit1 in enumerate(ecod2hits[ecod]):
                for c2, hit2 in enumerate(ecod2hits[ecod]):
                    if c1 < c2:
                        domain1 = hit1[0]
                        tgroup1 = hit1[1]
                        prob1 = hit1[2]
                        get_resids1 = hit1[3]
                        domain2 = hit2[0]
                        tgroup2 = hit2[1]
                        prob2 = hit2[2]
                        get_resids2 = hit2[3]
                        if prob1 + 0.1 > domain2prob[domain1] and prob2 + 0.1 > domain2prob[domain2]:
                            common_resids = get_resids1.intersection(get_resids2)
                            if len(common_resids) < 0.25 * len(get_resids1) or len(common_resids) < 0.25 * len(get_resids2):
                                domain_pair = domain1 + '_' + domain2
                                try:
                                    dpair2supports[domain_pair].append([ecod, tgroup1, prob1, prob2])
                                except KeyError:
                                    domain_pairs.append(domain_pair)
                                    dpair2supports[domain_pair] = [[ecod, tgroup1, prob1, prob2]]

    
    merge_pairs = []
    merge_info = []
    for domain_pair in domain_pairs:
        domain1 = domain_pair.split('_')[0]
        domain2 = domain_pair.split('_')[1]
        support_ecods = set([])
        for item in dpair2supports[domain_pair]:
            ecod = item[0]
            support_ecods.add(ecod)
        against_ecods1 = set([])
        against_ecods2 = set([])
        merge_info.append(domain1 + ',' + domain2 + '\t' + ','.join(support_ecods))

        for item in domain2hits[domain1]:
            ecod = item[0]
            tgroup = item[1]
            prob = item[2]
            ratio = item[3]
            if prob + 0.1 > domain2prob[domain1]:
                if ratio > 0.5:
                    if not ecod in support_ecods:
                        against_ecods1.add(ecod)

        for item in domain2hits[domain2]:
            ecod = item[0]
            tgroup = item[1]
            prob = item[2]
            ratio = item[3]
            if prob + 0.1 > domain2prob[domain2]:
                if ratio > 0.5:
                    if not ecod in support_ecods:
                        against_ecods2.add(ecod)

        if len(support_ecods) > len(against_ecods1) or len(support_ecods) > len(against_ecods2):
            merge_pairs.append([domain1, domain2])
 
    if merge_info:
        rp = open('step19/' + spname + '/' + prot + '.info','w')
        for merge_line in merge_info:
            rp.write(merge_line + '\n')
        rp.close()

    if merge_pairs:
        rp = open('step19/' + spname + '/' + prot + '.result','w')
        for merge_pair in merge_pairs:
            domain1 = merge_pair[0]
            domain2 = merge_pair[1]
            rp.write(domain1 + '\t' + domain2def[domain1] + '\t' + domain2 + '\t' + domain2def[domain2] + '\n')
        rp.close()

    if merge_info and merge_pairs:
        pass
    else:
        os.system('echo \'done\' > step19/' + spname + '/' + prot + '.done')
else:
    os.system('echo \'done\' > step19/' + spname + '/' + prot + '.done')
