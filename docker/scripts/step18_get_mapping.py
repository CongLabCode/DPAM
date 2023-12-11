#!/opt/conda/bin/python
import os, sys
import numpy as np

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

def check_overlap(residsA, residsB):
    overlap = set(residsA).intersection(set(residsB))
    if len(overlap) >= len(residsA) * 0.33:
        if len(overlap) >= len(residsA) * 0.5 or len(overlap) >= len(residsB) * 0.5:
            return 1
        else:
            return 0
    else:
        return 0

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
        ranges.append(f'{str(seg[0])}-{str(seg[-1])}')
    return ','.join(ranges)


spname = sys.argv[1]
prot = sys.argv[2]
HHhits = []
if os.path.exists('step5/' + spname + '/' + prot + '.result'):
    fp = open('step5/' + spname + '/' + prot + '.result', 'r')
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            ecodid = words[1]
            getres = set([])
            resmap = {}
            fp1 = open('/mnt/databases/ECOD_maps/' + ecodid + '.map', 'r')
            for line1 in fp1:
                words1 = line1.split()
                getres.add(int(words1[1]))
                resmap[int(words1[1])] = int(words1[0])
            fp1.close()
            hhprob = float(words[3]) / 100

            raw_qresids = get_resids(words[12])
            raw_tresids = get_resids(words[13])
            qresids = []
            tresids = []
            for i in range(len(raw_qresids)):
                if raw_tresids[i] in getres:
                    qresid = raw_qresids[i]
                    tresid = resmap[raw_tresids[i]]
                    qresids.append(qresid)
                    tresids.append(tresid)
            HHhits.append([ecodid, hhprob, qresids, tresids])
    fp.close()

DALIhits = []
if os.path.exists('step9/' + spname + '/' + prot + '_good_hits'):
    fp = open('step9/' + spname + '/' + prot + '_good_hits', 'r')
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            ecodid = words[1]
            zscore = float(words[4]) / 10
            qresids = get_resids(words[9])
            tresids = get_resids(words[10])
            DALIhits.append([ecodid, zscore, qresids, tresids])
    fp.close()

if os.path.exists('step17/' + spname + '/' + prot + '.result'):
    fp = open('step17/' + spname + '/' + prot + '.result', 'r')
    domains = []
    domain2def = {}
    domain2resids = {}
    domain2hits = {}
    for line in fp:
        words = line.split()
        dname = words[0]
        try:
            domain2resids[dname]
        except KeyError:
            domains.append(dname)
            domain2resids[dname] = get_resids(words[1])
            domain2def[dname] = words[1]
        tgroup = words[2]
        ecodhit = words[3]
        prob = float(words[4])
        judge = words[5]
        try:
            domain2hits[dname]
        except KeyError:
            domain2hits[dname] = {}
        domain2hits[dname][ecodhit] = [prob, tgroup, judge]
    fp.close()

    results = []
    for domain in domains:
        domain_resids = domain2resids[domain]
        domain_residset = set(domain_resids)
        hitinfo = domain2hits[domain]
        good_hits = list(hitinfo.keys())

        Hecods = set([])
        ecod2Hhit = {}
        for hit in HHhits:
            ecodid = hit[0]
            Hprob = hit[1]
            Hqresids = hit[2]
            Htresids = hit[3]
            if check_overlap(domain_resids, Hqresids):
                try:
                    if Hprob > ecod2Hhit[ecodid][0]:
                        ecod2Hhit[ecodid] = [Hprob, Hqresids, Htresids]
                except KeyError:
                    Hecods.add(ecodid)
                    ecod2Hhit[ecodid] = [Hprob, Hqresids, Htresids]

        Decods = set([])
        ecod2Dhit = {}
        for hit in DALIhits:
            ecodid = hit[0]
            Dzscore = hit[1]
            Dqresids = hit[2]
            Dtresids = hit[3]
            if check_overlap(domain_resids, Dqresids):
                try:
                    if Dzscore > ecod2Dhit[ecodid][0]:
                        ecod2Dhit[ecodid] = [Dzscore, Dqresids, Dtresids]
                except KeyError:
                    Decods.add(ecodid)
                    ecod2Dhit[ecodid] = [Dzscore, Dqresids, Dtresids]

        for hit in good_hits:
            [DPAMprob, tgroup, judge] = hitinfo[hit]
            if hit in Hecods:
                HQresids = ecod2Hhit[hit][1]
                HTresids = ecod2Hhit[hit][2]
                Hresids = []
                if len(HQresids) != len(HTresids):
                    print (spname, prot, domain, hit)
                    Hresid_string = 'na'
                elif HQresids:
                    for i in range(len(HQresids)):
                        if HQresids[i] in domain_residset:
                            Hresids.append(HTresids[i])
                    Hresid_string = get_range(Hresids)
                else:
                    Hresid_string = 'na'
            else:
                Hresid_string = 'na'

            if hit in Decods:
                DQresids = ecod2Dhit[hit][1]
                DTresids = ecod2Dhit[hit][2]
                Dresids = []
                if len(DQresids) != len(DTresids):
                    print (spname, prot, domain, hit)
                    Dresid_string = 'na'
                elif DQresids:
                    for i in range(len(DQresids)):
                        if DQresids[i] in domain_residset:
                            Dresids.append(DTresids[i])
                    Dresid_string = get_range(Dresids)
                else:
                    Dresid_string = 'na'
            else:
                Dresid_string = 'na'
            results.append(domain + '\t' + domain2def[domain] + '\t' + hit + '\t' + tgroup + '\t' + str(DPAMprob) + '\t' + judge + '\t' + Hresid_string + '\t' + Dresid_string + '\n')
    
    if results:
        rp = open('step18/' + spname + '/' + prot + '.data', 'w')
        for line in results:
            rp.write(line)
        rp.close()
    else:
        os.system('echo \'done\' > step18/' + spname + '/' + prot + '.done')
else:
    os.system('echo \'done\' > step18/' + spname + '/' + prot + '.done')
