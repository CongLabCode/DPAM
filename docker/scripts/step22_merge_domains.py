#!/opt/conda/bin/python
import sys

def get_range(resids):
    if resids:
        resids = list(resids)
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
            start = seg[0]
            end = seg[-1]
            ranges.append(str(start) + '-' + str(end))
        return ','.join(ranges)
    else:
        return 'na'


dataset = sys.argv[1]
fp = open('step21_' + dataset + '.result', 'r')
get_prots = set([])
prot2merges = {}
domain2resids = {}
for line in fp:
    words = line.split()
    prot = words[0]
    dom1 = words[1]
    dom2 = words[2]
    resids1 = []
    for seg in words[4].split(','):
        if '-' in seg:
            start = int(seg.split('-')[0])
            end = int(seg.split('-')[1])
            for res in range(start, end + 1):
                resids1.append(res)
        else:
            resids1.append(int(seg))
    resids2 = []
    for seg in words[5].split(','):
        if '-' in seg:
            start = int(seg.split('-')[0])
            end = int(seg.split('-')[1])
            for res in range(start, end + 1):
                resids2.append(res)
        else:
            resids2.append(int(seg))

    if int(words[3]) > 0:
        try:
            prot2merges[prot].append(set([dom1, dom2]))
        except KeyError:
            prot2merges[prot] = [set([dom1, dom2])]
            get_prots.add(prot)
        try:
            domain2resids[prot]
        except KeyError:
            domain2resids[prot] = {}
        domain2resids[prot][dom1] = resids1
        domain2resids[prot][dom2] = resids2
fp.close()


rp = open('step22_' + dataset + '.result','w')
for prot in get_prots:
    pairs = prot2merges[prot]
    groups = []
    for pair in pairs:
        groups.append(pair)
    while 1:
        newgroups = []
        for group in groups:
            if not newgroups:
                newgroups.append(group)
            else:
                for newgroup in newgroups:
                    if group.intersection(newgroup):
                        for item in group:
                            newgroup.add(item)
                        break
                else:
                    newgroups.append(group)

        if len(groups) == len(newgroups):
            break
        groups = []
        for newgroup in newgroups:
            groups.append(newgroup)
    
    for group in groups:
        group_domains = []
        group_resids = set([])
        for domain in group:
            group_domains.append(domain)
            group_resids = group_resids.union(domain2resids[prot][domain])
        group_range = get_range(group_resids)
        rp.write(prot + '\t' + ','.join(group_domains) + '\t' + group_range + '\n')
rp.close()
