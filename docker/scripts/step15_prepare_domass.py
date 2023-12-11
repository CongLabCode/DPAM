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
    if len(overlap) >= len(residsA) * 0.5 or len(overlap) >= len(residsB) * 0.5:
        return 1
    else:
        return 0

fp = open('/mnt/databases/ecod.latest.domains', 'r')
ecod2tgroup = {}
ecod2hgroup = {}
for line in fp:
    if line[0] != '#':
        words = line.split()
        ecodnum = words[0]
        tgroup = '.'.join(words[3].split('.')[:3])
        hgroup = '.'.join(words[3].split('.')[:2])
        ecod2tgroup[ecodnum] = tgroup
        ecod2hgroup[ecodnum] = hgroup
fp.close()

dataset = sys.argv[1]
prot = sys.argv[2]
fp = open('step12/' + dataset + '/' + prot + '.sse', 'r')
res2sse = {}
for line in fp:
    words = line.split()
    resid = int(words[0])
    res2sse[resid] = [words[2], words[3]]
fp.close()

domains = []
if os.path.exists('step14/' + dataset + '/' + prot + '.domains'):
    fp = open('step14/' + dataset + '/' + prot + '.domains', 'r')
    for line in fp:
        words = line.split()
        dname = words[0]
        resids = get_resids(words[1])
        sse2count = {}
        sse2type = {}
        sses = []
        for res in res2sse.keys():
            sse = res2sse[res][0]
            sse_type = res2sse[res][1]
            try:
                sse2count[sse] += 1
            except KeyError:
                sses.append(sse)
                sse2count[sse] = 1
                sse2type[sse] = sse_type

        helix_count = 0
        strand_count = 0
        for sse in sses:
            if sse2type[sse] == 'H':
                if sse2count[sse] >= 6:
                    helix_count += 1
            elif sse2type[sse] == 'E':
                if sse2count[sse] >= 3:
                    strand_count += 1
        domains.append([dname, set(resids), helix_count, strand_count, words[1]])
    fp.close()


HHhits = []
max_hrank = 0
if os.path.exists('step5/' + dataset + '/' + prot + '.result'):
    fp = open('step5/' + dataset + '/' + prot + '.result', 'r')
    qres2hgroups = {}
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            hitname = words[0]
            ecodid = words[1]
            hgroup = ecod2hgroup[ecodid]
            hhprob = float(words[3]) / 100
            coverage = float(words[10])
            qresids = get_resids(words[12])
            tresids = get_resids(words[13])
            for qres in qresids:
                try:
                    qres2hgroups[qres].add(hgroup)
                except KeyError:
                    qres2hgroups[qres] = set([hgroup])
            ranks = []
            for qres in qresids:
                ranks.append(len(qres2hgroups[qres]))
            hrank = np.mean(ranks)
            if hrank > max_hrank:
                max_hrank = hrank
            if len(qresids) != len(tresids):
                print (prot + '\t' + ecodid + '\thhsearch')
            else:
                HHhits.append([ecodid, hhprob, coverage, hrank / 10, qresids, tresids, hitname])
    fp.close()
if max_hrank > 100:
    pass
else:
    max_hrank = 100


DALIhits = []
max_drank = 0
if os.path.exists('step9/' + dataset + '/' + prot + '_good_hits'):
    fp = open('step9/' + dataset + '/' + prot + '_good_hits', 'r')
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            hitname = words[0]
            ecodid = words[1]
            resmap = {}
            fp1 = open('/mnt/databases/ECOD_maps/' + ecodid + '.map','r')
            for line1 in fp1:
                words1 = line1.split()
                resmap[int(words1[0])] = int(words1[1])
            fp1.close()

            zscore = float(words[4]) / 10
            qscore = float(words[5])
            ztile = float(words[6])
            qtile = float(words[7])
            drank = float(words[8])
            if drank > max_drank:
                max_drank = drank
            qresids = get_resids(words[9])
            raw_tresids = get_resids(words[10])
            tresids = []
            for resid in raw_tresids:
                tresids.append(resmap[resid])
            rot1 = words[11]
            rot2 = words[12]
            rot3 = words[13]
            trans = words[14]
            if len(qresids) != len(tresids):
                print (prot + '\t' + ecodid + '\tdali')
            else:
                DALIhits.append([ecodid, zscore, qscore, ztile, qtile, drank / 10, qresids, tresids, hitname, rot1, rot2, rot3, trans])
    fp.close()
if max_drank > 100:
    pass
else:
    max_drank = 100


if domains:
    rp = open('step15/' + dataset + '/' + prot + '.data', 'w')
    rp.write('domID\tdomRange\ttgroup\tecodid\tdomLen\tHelix_num\tStrand_num\tHHprob\tHHcov\tHHrank\tDzscore\tDqscore\tDztile\tDqtile\tDrank\tCdiff\tCcov\tHHname\tDname\tDrot1\tDrot2\tDrot3\tDtrans\n')
    for countd, domain in enumerate(domains):
        domain_name = domain[0]
        domain_resids = domain[1]
        domain_length = len(domain_resids)
        domain_helix = domain[2]
        domain_strand = domain[3]
        domain_range = domain[4]

        ecod2Hhit = {}
        Hecods = set([])
        for hit in HHhits:
            ecodid = hit[0]
            Hprob = hit[1]
            Hcov = hit[2]
            Hrank = hit[3]
            Hqresids = hit[4]
            Htresids = hit[5]
            Hname = hit[6]
            if check_overlap(domain_resids, Hqresids):
                try:
                    if Hprob > ecod2Hhit[ecodid][0]:
                        ecod2Hhit[ecodid] = [Hprob, Hcov, Hrank, Hqresids, Htresids, Hname]
                except KeyError:
                    Hecods.add(ecodid)
                    ecod2Hhit[ecodid] = [Hprob, Hcov, Hrank, Hqresids, Htresids, Hname]

        ecod2Dhit = {}
        Decods = set([])
        for hit in DALIhits:
            ecodid = hit[0]
            Dzscore = hit[1]
            Dqscore = hit[2]
            Dztile = hit[3]
            Dqtile = hit[4]
            Drank = hit[5]
            Dqresids = hit[6]
            Dtresids = hit[7]
            Dname = hit[8]
            Drot1 = hit[9]
            Drot2 = hit[10]
            Drot3 = hit[11]
            Dtrans = hit[12]
            if check_overlap(domain_resids, Dqresids):
                try:
                    if Dzscore > ecod2Dhit[ecodid][0]:
                        ecod2Dhit[ecodid] = [Dzscore, Dqscore, Dztile, Dqtile, Drank, Dqresids, Dtresids, Dname, Drot1, Drot2, Drot3, Dtrans]
                except KeyError:
                    Decods.add(ecodid)
                    ecod2Dhit[ecodid] = [Dzscore, Dqscore, Dztile, Dqtile, Drank, Dqresids, Dtresids, Dname, Drot1, Drot2, Drot3, Dtrans]

        for ecodid in Hecods.intersection(Decods):
            tgroup = ecod2tgroup[ecodid]
            Hhit = ecod2Hhit[ecodid]
            Hprob = Hhit[0]
            Hcov = Hhit[1]
            Hrank = Hhit[2]
            Hqres = Hhit[3]
            Htres = Hhit[4]
            Hname = Hhit[5]
            Hmap = {}
            for i in range(len(Hqres)):
                Hmap[Hqres[i]] = Htres[i]

            Dhit = ecod2Dhit[ecodid]
            Dzscore = Dhit[0]
            Dqscore = Dhit[1]
            Dztile = Dhit[2]
            Dqtile = Dhit[3]
            Drank = Dhit[4]
            Dqres = Dhit[5]
            Dtres = Dhit[6]
            Dname = Dhit[7]
            Drot1 = Dhit[8]
            Drot2 = Dhit[9]
            Drot3 = Dhit[10]
            Dtrans = Dhit[11]
            Dmap = {}
            for i in range(len(Dqres)):
                Dmap[Dqres[i]] = Dtres[i]

            Cqres = set(Hqres).intersection(set(Dqres))
            Ccov = len(Cqres) / domain_length
            Cdiffs = []
            for res in Cqres:
                Hres = Hmap[res]
                Dres = Dmap[res]
                if Hres > Dres:
                    Cdiffs.append(Hres - Dres)
                else:
                    Cdiffs.append(Dres - Hres)
            if Cdiffs:
                Cdiff = np.mean(Cdiffs)
            else:
                Cdiff = -1
            rp.write(domain_name + '\t' + domain_range + '\t' + tgroup + '\t' + ecodid + '\t' + str(domain_length) + '\t' + str(domain_helix) + '\t' + str(domain_strand) + '\t' + str(round(Hprob, 3)) + '\t' + str(round(Hcov, 3)) + '\t' + str(round(Hrank, 2)) + '\t' + str(round(Dzscore, 3)) + '\t' + str(round(Dqscore, 3)) + '\t' + str(round(Dztile, 3)) + '\t' + str(round(Dqtile, 3)) + '\t' + str(round(Drank, 2)) + '\t' + str(round(Cdiff, 2)) + '\t' + str(round(Ccov, 3)) + '\t' + Hname + '\t' + Dname + '\t' + Drot1 + '\t' + Drot2 + '\t' + Drot3 + '\t' + Dtrans + '\n')

        for ecodid in Hecods.difference(Decods):
            tgroup = ecod2tgroup[ecodid]
            Hhit = ecod2Hhit[ecodid]
            Hprob = Hhit[0]
            Hcov = Hhit[1]
            Hrank = Hhit[2]
            Hname = Hhit[5]
            Dzscore = 0
            Dqscore = 0
            Dztile = 10
            Dqtile = 10
            Drank = max_drank
            Ccov = 0
            Cdiff = -1
            rp.write(domain_name + '\t' + domain_range + '\t' + tgroup + '\t' + ecodid + '\t' + str(domain_length) + '\t' + str(domain_helix) + '\t' + str(domain_strand) + '\t' + str(round(Hprob, 3)) + '\t' + str(round(Hcov, 3)) + '\t' + str(round(Hrank, 2)) + '\t' + str(round(Dzscore, 3)) + '\t' + str(round(Dqscore, 3)) + '\t' + str(round(Dztile, 3)) + '\t' + str(round(Dqtile, 3)) + '\t' + str(round(Drank, 2)) + '\t' + str(round(Cdiff, 2)) + '\t' + str(round(Ccov, 3)) + '\t' + Hname + '\tna\tna\tna\tna\tna\n')

        for ecodid in Decods.difference(Hecods):
            tgroup = ecod2tgroup[ecodid]
            Dhit = ecod2Dhit[ecodid]
            Dzscore = Dhit[0]
            Dqscore = Dhit[1]
            Dztile = Dhit[2]
            Dqtile = Dhit[3]
            Drank = Dhit[4]
            Dname = Dhit[7]
            Drot1 = Dhit[8]
            Drot2 = Dhit[9]
            Drot3 = Dhit[10]
            Dtrans = Dhit[11]
            Hprob = 0
            Hcov = 0
            Hrank = max_hrank
            Ccov = 0
            Cdiff = -1
            rp.write(domain_name + '\t' + domain_range + '\t' + tgroup + '\t' + ecodid + '\t' + str(domain_length) + '\t' + str(domain_helix) + '\t' + str(domain_strand) + '\t' + str(round(Hprob, 3)) + '\t' + str(round(Hcov, 3)) + '\t' + str(round(Hrank, 2)) + '\t' + str(round(Dzscore, 3)) + '\t' + str(round(Dqscore, 3)) + '\t' + str(round(Dztile, 3)) + '\t' + str(round(Dqtile, 3)) + '\t' + str(round(Drank, 2)) + '\t' + str(round(Cdiff, 2)) + '\t' + str(round(Ccov, 3)) + '\tna\t' + Dname + '\t' + Drot1 + '\t' + Drot2 + '\t' + Drot3 + '\t' + Dtrans + '\n')
    rp.close()
else:
    os.system('echo \'done\' > step15/' + dataset + '/' + prot + '.done')
