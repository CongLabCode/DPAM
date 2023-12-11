import os, sys
import numpy as np

def get_good_coor(value):
    if len(str(round(value,3))) <= 7:
        newvalue = str(round(value,3)).rjust(8, ' ')
    elif len(str(round(value,2))) <= 8:
        newvalue = str(round(value,2)).rjust(8, ' ')
    else:
        newvalue = str(round(value,1)).rjust(8, ' ')
        print ('error\t' + str(value))
    return newvalue

def get_domain_range(resids):
    segs = []
    resids.sort()
    for resid in resids:
        if not segs:
            segs.append([resid])
        else:
            if resid > segs[-1][-1] + 1:
                segs.append([resid])
            else:
                segs[-1].append(resid)
    seg_string = []
    for seg in segs:
        start = str(seg[0])
        end = str(seg[-1])
        seg_string.append(f'{start}-{end}')
    return ','.join(seg_string)


fp = open("/mnt/databases/ecod.latest.domains","r")
ecod2keyword = {}
for line in fp:
    if line[0] != "#":
        words = line.split()
        ecod2keyword[words[0]] = words[1]
fp.close()

dataset = sys.argv[1]
prot = sys.argv[2]
need_ecodnums = set([])
names = []
name2resids = {}
name2ecodnum = {}
name2domrange = {}

fp = open('step24/' + dataset + '/' + prot + '_domains','r')
for line in fp:
    words = line.split()
    name = words[0]
    if words[2] != 'na':
        name2domrange[name] = '+'.join(words[1].split(','))
        resids = set([])
        for seg in words[1].split(','):
            if '-' in seg:
                start = int(seg.split('-')[0])
                end = int(seg.split('-')[1])
                for res in range(start, end + 1):
                    resids.add(res)
            else:
                resids.add(int(seg))
        names.append(name)
        name2resids[name] = resids
        name2ecodnum[name] = words[2]
        need_ecodnums.add(words[2])
fp.close()


all_hits = []
ecod2info = {}
if os.path.exists('step16/' + dataset + '/' + prot + '.result'):
    fp = open('step16/' + dataset + '/' + prot + '.result','r')
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            resids = set([])
            for seg in words[1].split(','):
                if '-' in seg:
                    start = int(seg.split('-')[0])
                    end = int(seg.split('-')[1])
                    for res in range(start, end + 1):
                        resids.add(res)
                else:
                    resids.add(int(seg))

            tgroup = words[2]
            ecod = words[3]
            keyword = ecod2keyword[ecod]
            DPAMprob = words[4]
            HHprob = str(round(float(words[5]) * 100, 1))
            HHcov = words[6]
            HHrank = words[7]
            DALIzscore = str(round(float(words[8]) * 10, 1))
            DALIqscore = words[9]
            DALIztile = words[10]
            DALIqtile = words[11]
            DALIrank = words[12]
            Cdiff = words[13]
            Ccov = words[14]
            all_hits.append([resids, tgroup, keyword, float(DPAMprob), HHprob, HHcov, HHrank, DALIzscore, DALIqscore, DALIztile, DALIqtile, DALIrank, Cdiff, Ccov])

            if ecod in need_ecodnums:
                bad_mtx = 0
                Tvec = []
                for value in words[20].split(','):
                    if value == 'na':
                        bad_mtx = 1
                    else:
                        Tvec.append(float(value))
                Rmtx = [[],[],[]]
                for value in words[17].split(','):
                    if value == 'na':
                        bad_mtx = 1
                    else:
                        Rmtx[0].append(float(value))
                for value in words[18].split(','):
                    if value == 'na':
                        bad_mtx = 1
                    else:
                        Rmtx[1].append(float(value))
                for value in words[19].split(','):
                    if value == 'na':
                        bad_mtx = 1
                    else:
                        Rmtx[2].append(float(value))
                Rmtx = np.array(Rmtx)
                Tvec = np.array(Tvec)

                if not bad_mtx:
                    resids = set([])
                    for seg in words[1].split(','):
                        if '-' in seg:
                            start = int(seg.split('-')[0])
                            end = int(seg.split('-')[1])
                            for res in range(start, end + 1):
                                resids.add(res)
                        else:
                            resids.add(int(seg))
                    try:
                        ecod2info[ecod].append([resids, Rmtx, Tvec])
                    except KeyError:
                        ecod2info[ecod] = [[resids, Rmtx, Tvec]]
    fp.close()
all_hits.sort(key = lambda x:x[3], reverse = True)


name2hits = {}
for name in names:
    name2hits[name] = []
    domain_resids = name2resids[name]
    for counti, item in enumerate(all_hits):
        if item[0].issubset(domain_resids):
            name2hits[name].append(item)
    ecodnum = name2ecodnum[name]
    ecodkey = ecod2keyword[ecodnum]
    hits = []
    try:
        for item in ecod2info[ecodnum]:
            resids = item[0]
            Rmtx = item[1]
            Tvec = item[2]
            overlap = len(resids.intersection(domain_resids))
            hits.append([overlap, Rmtx, Tvec])
    except KeyError:
        pass

    if hits:
        hits.sort(key = lambda x:x[0], reverse = True)
        best_hit = hits[0]
        Rmtx = best_hit[1]
        Tvec = best_hit[2]
    
        fp = open('ECOD_pdbs/' + ecodnum + '.pdb', 'r')
        rp = open('step25/' + dataset + '/' + prot + '_' + name + '_' + ecodnum + '.pdb', 'w')
        for line in fp:
            if len(line) >= 60:
                if line[:4] == 'ATOM':
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coor = np.array([x, y, z])
                    newcoor = Rmtx.dot(coor) + Tvec
                    newx = get_good_coor(newcoor[0])
                    newy = get_good_coor(newcoor[1])
                    newz = get_good_coor(newcoor[2])
                    rp.write(line[:30] + newx + newy + newz + line[54:])
                else:
                    rp.write(line)
            else:
                rp.write(line)
        fp.close()
        rp.close()

        domrange = name2domrange[name]
        rp = open('step25/' + dataset + '/' + prot + '_' + name + '_' + ecodnum + '.pml', 'w')
        rp.write('load step2/' + dataset + '/' + prot + '.pdb, ' + prot + '\n')
        rp.write('load step25/' + dataset + '/' + prot + '_' + name + '_' + ecodnum + '.pdb, ' + ecodkey + '\n')
        rp.write('color grey, ' + prot + '\n')
        rp.write('spectrum count, rainbow, ' + ecodkey + '\n')
        rp.write('select ' + prot + '_' + name + ', resi ' + domrange + ' & ' + prot + '\n')
        rp.write('spectrum count, rainbow, ' + prot + '_' + name + '\n')
        rp.write('delete selection ' + prot + '_' + name + '\n')
        rp.write('save ' + dataset + '_web/' + prot + '/' + prot + '_' + name + '_' + ecodnum + '.pse\n')
        rp.write('delete ' + prot + '\n')
        rp.write('delete ' + ecodkey + '\n')
        rp.close()
        os.system('pymol -c step25/' + dataset + '/' + prot + '_' + name + '_' + ecodnum + '.pml >& pml_log')


if all_hits:
    for name in names:
        rp = open(dataset + '_web/' + prot + '/' + prot + '_' + name + '.html','w')
        rp.write('<!doctype html>\n\n')
        rp.write('<html lang=\'en\'>\n')
        rp.write('<head>\n')
        rp.write('  <meta charset=\'utf-8\'>\n')
        rp.write('  <link rel=\'stylesheet\' href=\'../table.css\'>\n')
        rp.write('</head>\n')
        rp.write('<body>\n')
        rp.write(' <div id=\'wrapper\'>\n')
        rp.write('  <h1>Candidate homologous ECOD domains for ' + prot + '</h1>\n\n')
        rp.write('  <table id=\'keywords\' cellspacing=\'0\' cellpadding=\'0\'>\n')
        rp.write('    <thead>\n')
        rp.write('      <tr>\n')
        rp.write('        <th><span>Query</span></th>\n')
        rp.write('        <th><span>Domain</span></th>\n')
        rp.write('        <th><span>Range</span></th>\n')
        rp.write('        <th><span>T-group</span></th>\n')
        rp.write('        <th><span>ECOD ref</span></th>\n')
        rp.write('        <th><span>DPAM prob</span></th>\n')
        rp.write('        <th><span>HH prob</span></th>\n')
        rp.write('        <th><span>HH cov</span></th>\n')
        rp.write('        <th><span>HH rank</span></th>\n')
        rp.write('        <th><span>DALI zscore</span></th>\n')
        rp.write('        <th><span>DALI qscore</span></th>\n')
        rp.write('        <th><span>DALI ztile</span></th>\n')
        rp.write('        <th><span>DALI qtile</span></th>\n')
        rp.write('        <th><span>DALI rank</span></th>\n')
        rp.write('        <th><span>CON diff</span></th>\n')
        rp.write('        <th><span>CON cov</span></th>\n')
        rp.write('      </tr>\n')
        rp.write('    </thead>\n')
        rp.write('    <tbody>\n')
    
        for hit in name2hits[name]:
            resids = list(hit[0])
            domrange = get_domain_range(resids)
            tgroup = hit[1]
            ecodkey = hit[2]
            DPAMprob = str(hit[3])
            HHprob = hit[4]
            HHcov = hit[5]
            HHrank = hit[6]
            DALIzscore = hit[7]
            DALIqscore = hit[8]
            DALIztile = hit[9]
            DALIqtile = hit[10]
            DALIrank = hit[11]
            Cdiff = hit[12]
            Ccov = hit[13]

            rp.write('      <tr>\n')
            rp.write('        <td><a href=\'https://alphafold.ebi.ac.uk/entry/' + prot + '\'>' + prot + '</a></td>\n')
            rp.write('        <td>' + name + '</td>\n')
            rp.write('        <td>' + domrange + '</td>\n')
            rp.write('        <td><a href=\'http://prodata.swmed.edu/ecod/complete/tree?id=' + tgroup + '\'>' + tgroup + '</a></td>\n')
            rp.write('        <td><a href=\'http://prodata.swmed.edu/ecod/complete/domain/' + ecodkey + '\'>' + ecodkey + '</a></td>\n')
            rp.write('        <td>' + DPAMprob + '</td>\n')
            rp.write('        <td>' + HHprob + '</td>\n')
            rp.write('        <td>' + HHcov + '</td>\n')
            rp.write('        <td>' + HHrank + '</td>\n')
            rp.write('        <td>' + DALIzscore + '</td>\n')
            rp.write('        <td>' + DALIqscore + '</td>\n')
            rp.write('        <td>' + DALIztile + '</td>\n')
            rp.write('        <td>' + DALIqtile + '</td>\n')
            rp.write('        <td>' + DALIrank + '</td>\n')
            rp.write('        <td>' + Cdiff + '</td>\n')
            rp.write('        <td>' + Ccov + '</td>\n')
            rp.write('      </tr>\n')

        rp.write('    </tbody>\n')
        rp.write('  </table>\n')
        rp.write(' </div>\n')
        rp.write('</body>\n')
        rp.write('<script src=\'https://ajax.googleapis.com/ajax/libs/jquery/2.2.2/jquery.min.js\'></script>\n')
        rp.write('<script type=\'text/javascript\' src=\'../tablesorter-master/dist/js/jquery.tablesorter.js\'></script>\n')
        rp.write('<script>\n')
        rp.write('  $(function(){\n')
        rp.write('    $(\'#keywords\').tablesorter();\n')
        rp.write('  });\n')
        rp.write('</script>\n')
        rp.write('</html>\n')
        rp.close()
