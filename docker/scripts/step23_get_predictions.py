#!/opt/conda/bin/python
import os, sys

def get_resids(segs):
    resids = []
    for seg in segs.split(','):
        if '-' in seg:
            start = int(seg.split('-')[0])
            end = int(seg.split('-')[1])
            for res in range(start, end + 1):
                resids.append(res)
        else:
            resids.append(int(seg))
    return resids


dataset = sys.argv[1]
prot = sys.argv[2]
fp = open('/mnt/databases/tgroup_length','r')
tgroup2length = {}
for line in fp:
    words = line.split()
    tgroup2length[words[0]] = float(words[1])
fp.close()

domain2hits = {}
if os.path.exists('step16/' + dataset + '/' + prot + '.result'):
    fp = open('step16/' + dataset + '/' + prot + '.result', 'r')
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            domain = words[0]
            tgroup = words[2]
            ecodid = words[3]
            DPAMprob = float(words[4])
            HHprob = float(words[5])
            DALIzscore = float(words[8])
            try:
                domain2hits[domain].append([ecodid, tgroup, DPAMprob, HHprob, DALIzscore])
            except KeyError:
                domain2hits[domain] = [[ecodid, tgroup, DPAMprob, HHprob, DALIzscore]]
    fp.close()

mergdoms = set([])
merges = []
if os.path.exists('step22_' + dataset + '.result'):
    fp = open('step22_' + dataset + '.result', 'r')
    for line in fp:
        words = line.split()
        if words[0] == prot:
            merges.append([words[1], words[2]])
            for dom in words[1].split(','):
                mergdoms.add(dom)
    fp.close()

need_ecods = set([])
if os.path.exists('step18/' + dataset + '/' + prot + '.data'):
    fp = open('step18/' + dataset + '/' + prot + '.data','r')
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
    if os.path.exists('posi_weights/' + ecod + '.weight'):
        fp = open('posi_weights/' + ecod + '.weight','r')
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

all_domains = set([])
for item in merges:
    all_domains.add(item[0])
domain2range = {}
if os.path.exists('step14/' + dataset + '/' + prot + '.domains'):
    fp = open('step14/' + dataset + '/' + prot + '.domains','r')
    for line in fp:
        words = line.split()
        if not words[0] in mergdoms:
            all_domains.add(words[0])
            domain2range[words[0]] = words[1]
    fp.close()
full_domains = set([])
part_domains = set([])
domain2fullhit = {}
domain2parthit = {}

single_domains = []
domain2info = {}
if os.path.exists('step18/' + dataset + '/' + prot + '.data'):
    fp = open('step18/' + dataset + '/' + prot + '.data', 'r')
    for line in fp:
        words = line.split()
        domain = words[0]
        domlen = len(get_resids(words[1]))
        if words[6] == 'na':
            Hresids = []
        else:
            Hresids = get_resids(words[6])
        if words[7] == 'na':
            Dresids = []
        else:
            Dresids = get_resids(words[7])
        if domain in mergdoms:
            pass
        else:
            single_domains.append(domain)
        try:
            domain2info[domain].append([words[2], words[3], float(words[4]), words[5], Hresids, Dresids, float(words[4]) * domlen])
        except KeyError:
            domain2info[domain] = [[words[2], words[3], float(words[4]), words[5], Hresids, Dresids, float(words[4]) * domlen]]
    fp.close()

multi_domains = []
for item in merges:
    new_domain = item[0]
    multi_domains.append(new_domain)
    domain2range[new_domain] = item[1]
    for domain in new_domain.split(','):
        for hit in domain2info[domain]:
            try:
                for ohit in domain2info[new_domain]:
                    if hit[0] == ohit[0]:
                        if ohit[2] < hit[2]:
                            ohit[2] = hit[2]
                        if ohit[3] == 'good' or hit[3] == 'good':
                            ohit[3] = 'good'
                        elif ohit[3] == 'ok' or hit[3] == 'ok':
                            ohit[3] = 'ok'
                        else:
                            ohit[3] == 'bad'
                        for Hresid in hit[4]:
                            ohit[4].append(Hresid)
                        for Dresid in hit[5]:
                            ohit[5].append(Dresid)
                        break
                else:
                    domain2info[new_domain].append(hit)
            except KeyError:
                domain2info[new_domain] = [hit]

for domain in single_domains:
    domdef = domain2range[domain]
    domresids = get_resids(domdef)
    hitinfo = domain2info[domain]
    hitinfo.sort(key = lambda x:x[2], reverse = True)
    for item in hitinfo:
        ecodid = item[0]
        tgroup = item[1]
        prob = item[2]
        judge = item[3]
        Hresids = item[4]
        Dresids = item[5]
        if len(Dresids) > len(Hresids) * 0.5:
            HDresids = set(Dresids)
        else:
            HDresids = set(Hresids)

        total_weight = ecod2totW[ecodid]
        get_weight = 0
        for resid in HDresids:
            try:
                get_weight += ecod2posW[ecodid][resid]
            except KeyError:
                print (prot, ecodid, resid)
        ratio1 = get_weight / total_weight
        ratio2 = len(domresids) / tgroup2length[tgroup]
        besthit = [ecodid, tgroup, prob, ratio1, ratio2, judge]
        try:
            domain2parthit[domain]
        except KeyError:
            domain2parthit[domain] = besthit
        if ratio1 >= 0.66 or ratio2 >= 0.66:
            if ratio1 >= 0.33 and ratio2 >= 0.33:
                domain2fullhit[domain] = besthit
                break


for domain in multi_domains:
    domdef = domain2range[domain]
    domresids = get_resids(domdef)
    hitinfo = domain2info[domain]
    hitinfo.sort(key = lambda x:x[6], reverse = True)
    for item in hitinfo:
        ecodid = item[0]
        tgroup = item[1]
        prob = item[2]
        judge = item[3]
        Hresids = item[4]
        Dresids = item[5]
        if len(Dresids) > len(Hresids) * 0.5:
            HDresids = set(Dresids)
        else:
            HDresids = set(Hresids)

        total_weight = ecod2totW[ecodid]
        get_weight = 0
        for resid in HDresids:
            try:
                get_weight += ecod2posW[ecodid][resid]
            except KeyError:
                print (prot, ecodid, resid)
        ratio1 = get_weight / total_weight
        ratio2 = len(domresids) / tgroup2length[tgroup]
        besthit = [ecodid, tgroup, prob, ratio1, ratio2, judge]
        try:
            domain2parthit[domain]
        except KeyError:
            domain2parthit[domain] = besthit
        if ratio1 >= 0.66 or ratio2 >= 0.66:
            if ratio1 >= 0.33 and ratio2 >= 0.33:
                domain2fullhit[domain] = besthit
                break

for domain in domain2fullhit.keys():
    besthit = domain2fullhit[domain]
    prob = besthit[2]
    if prob >= 0.85:
        full_domains.add(domain)
for domain in domain2parthit.keys():
    if not domain in full_domains:
        besthit = domain2parthit[domain]
        prob = besthit[2]
        if prob >= 0.85:
            part_domains.add(domain)

results = []
get_domains = set([])
for domain in full_domains:
    get_domains.add(domain)
    domdef = domain2range[domain]
    besthit = domain2fullhit[domain]
    ecodid = besthit[0]
    DPAMprobs = []
    HHprobs = []
    DALIzs = []
    for dom in domain.split(','):
        for item in domain2hits[dom]:
            if item[0] == ecodid:
                DPAMprobs.append(item[2])
                HHprobs.append(item[3])
                DALIzs.append(item[4])

    if DPAMprobs and HHprobs:
        DPAMprob = round(max(DPAMprobs), 3)
        HHprob = round(max(HHprobs), 3)
        DALIz = round(max(DALIzs), 3)
        tgroup = besthit[1]
        prob = round(besthit[2], 3)
        ratio1 = round(besthit[3], 3)
        ratio2 = round(besthit[4], 3)
        judge = besthit[5]
        if DPAMprob / prob > 0.95 and DPAMprob / prob < 1.05:
            results.append('full\t' + domain + '\t' + domdef + '\t' + ecodid + '\t' + tgroup + '\t' + str(prob) + '\t' + str(HHprob) + '\t' + str(DALIz) + '\t' + str(ratio1) + '\t' + str(ratio2) + '\t' + judge + '\n')
        else:
            print ('error1\t' + dataset + '\t' + prot + '\t' + domain)
    else:
        print ('error2\t' + dataset + '\t' + prot + '\t' + domain)


for domain in part_domains:
    get_domains.add(domain)
    domdef = domain2range[domain]
    besthit = domain2parthit[domain]
    ecodid = besthit[0]
    DPAMprobs = []
    HHprobs = []
    DALIzs = []
    for dom in domain.split(','):
        for item in domain2hits[dom]:
            if item[0] == ecodid:
                DPAMprobs.append(item[2])
                HHprobs.append(item[3])
                DALIzs.append(item[4])

    if DPAMprobs and HHprobs:
        DPAMprob = round(max(DPAMprobs), 3)
        HHprob = round(max(HHprobs), 3)
        DALIz = round(max(DALIzs), 3)
        tgroup = besthit[1]
        prob = round(besthit[2], 3)
        ratio1 = round(besthit[3], 3)
        ratio2 = round(besthit[4], 3)
        judge = besthit[5]

        if DPAMprob / prob > 0.95 and DPAMprob / prob < 1.05:
            results.append('part\t' + domain + '\t' + domdef + '\t' + ecodid + '\t' + tgroup + '\t' + str(prob) + '\t' + str(HHprob) + '\t' + str(DALIz) + '\t' + str(ratio1) + '\t' + str(ratio2) + '\t' + judge + '\n')
        else:
            print ('error1\t' + dataset + '\t' + prot + '\t' + domain)
    else:
        print ('error2\t' + dataset + '\t' + prot + '\t' + domain)


other_domains = all_domains.difference(get_domains)
for domain in other_domains:
    domdef = domain2range[domain]
    getit = 0
    try:
        domain2hits[domain]
        getit = 1
    except KeyError:
        pass
    if getit:
        domain2hits[domain].sort(key = lambda x:x[2], reverse = True)
        tophit = domain2hits[domain][0]
        ecodid = tophit[0]
        tgroup = tophit[1]
        DPAMprob = round(tophit[2], 3)
        HHprob = round(tophit[3], 3)
        DALIz = round(tophit[4], 3)
        results.append('miss\t' + domain + '\t' + domdef + '\t' + ecodid + '\t' + tgroup + '\t' + str(DPAMprob) + '\t' + str(HHprob) + '\t' + str(DALIz) + '\tna\tna\tna\n')
    else:
        results.append('miss\t' + domain + '\t' + domdef + '\tna\tna\tna\tna\tna\tna\tna\tna\n')
    
if results:
    rp = open('step23/' + dataset + '/' + prot + '.assign', 'w')
    for line in results:
        rp.write(line)
    rp.close()
else:
    os.system('echo \'done\' > step23/' + dataset + '/' + prot + '.done')
