#!/opt/conda/bin/python

import os, sys

dataset = sys.argv[1]
prot = sys.argv[2]
if os.path.exists('step16/' + dataset + '/' + prot + '.result'):
    fp = open('step16/' + dataset + '/' + prot + '.result','r')
    domains = []
    domain2range = {}
    domain2hits = {}
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            domain = words[0]
            drange = words[1]
            tgroup = words[2]
            refdom = words[3]
            prob = float(words[4])
            domain2range[domain] = drange
            try:
                domain2hits[domain].append([tgroup, refdom, prob])
            except KeyError:
                domains.append(domain)
                domain2hits[domain] = [[tgroup, refdom, prob]]
    fp.close()

    results = []
    for domain in domains:
        drange = domain2range[domain]
        tgroups = []
        tgroup2best = {}
        for hit in domain2hits[domain]:
            tgroup = hit[0]
            refdom = hit[1]
            prob = hit[2]
            try:
                if prob > tgroup2best[tgroup]:
                    tgroup2best[tgroup] = prob
            except KeyError:
                tgroups.append(tgroup)
                tgroup2best[tgroup] = prob

        domain2hits[domain].sort(key = lambda x:x[2], reverse = True)
        for hit in domain2hits[domain]:
            tgroup = hit[0]
            refdom = hit[1]
            prob = hit[2]
            if prob >= 0.6:
                similar_tgroups = set([])
                for ogroup in tgroups:
                    if prob < tgroup2best[ogroup] + 0.05:
                        similar_tgroups.add(ogroup)
                similar_hgroups = set([])
                for group in similar_tgroups:
                    hgroup = group.split('.')[0] + '.' + group.split('.')[1]
                    similar_hgroups.add(hgroup)

                if len(similar_tgroups) == 1:
                    judge = 'good'
                elif len(similar_hgroups) == 1:
                    judge = 'ok'
                else:
                    judge = 'bad'
                results.append(domain + '\t' + drange + '\t' + tgroup + '\t' + refdom + '\t' + str(prob) + '\t' + judge + '\n')
    
    if results:
        rp = open('step17/' + dataset + '/' + prot + '.result','w')
        for line in results:
            rp.write(line)
        rp.close()
    else:
        os.system('echo \'done\' > step17/' + dataset + '/' + prot + '.done')
else:
    os.system('echo \'done\' > step17/' + dataset + '/' + prot + '.done')
