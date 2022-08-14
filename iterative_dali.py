import os, sys
import time

script_dir=os.path.dirname(os.path.realpath(__file__))
with open(script_dir+'/config_file') as f:
    configs=f.readlines()
configs=[i.strip() for i in configs if i.strip()!='']
configs={i.split(':')[0]:i.split(':')[1] for i in configs}



def get_domain_range(resids):
    segs = []
    resids.sort()
    cutoff1 = 5
    cutoff2 = len(resids) * 0.05
    cutoff = max(cutoff1, cutoff2)
    for resid in resids:
        if not segs:
            segs.append([resid])
        else:
            if resid > segs[-1][-1] + cutoff:
                segs.append([resid])
            else:
                segs[-1].append(resid)
    seg_string = []
    for seg in segs:
        start = str(seg[0])
        end = str(seg[-1])
        seg_string.append(start + '-' + end)
    return ','.join(seg_string)


prot = sys.argv[1]
ecod_database_dir=configs['ecod_database_dir']
if not os.path.exists('ECOD_database'):
    os.system(f'ln -s {ecod_database_dir} ./ECOD_database')
cwd=os.getcwd()
if os.path.exists(prot+'_hits'):
    fp = open(prot + '_hits','r')
    edomains = []
    for line in fp:
        words = line.split()
        edomains.append(words[0])
    fp.close()
    if os.path.exists('iterative_dali_tmp'):
        os.system('rm -rf iterative_dali_tmp')
    os.system('mkdir iterative_dali_tmp')
    os.system('mkdir iterative_dali_tmp/tmp')
    os.chdir('iterative_dali_tmp/tmp')
    for edomain in edomains:
        print('cp '+cwd+'/'+prot + '.pdb '+cwd+'/iterative_dali_tmp/'+prot + '__' + edomain + '.pdb')
        os.system('cp '+cwd+'/'+prot + '.pdb '+cwd+'/iterative_dali_tmp/'+prot + '__' + edomain + '.pdb')
        alicount = 0
        while 1:
            print(f'dali.pl --pdbfile1 ../{prot}__{edomain}.pdb --pdbfile2 ../../ECOD_database/{edomain}.pdb --dat1 . --dat2 . --outfmt summary,alignments,transrot >& log')
            os.system(f'dali.pl --pdbfile1 ../{prot}__{edomain}.pdb --pdbfile2 ../../ECOD_database/{edomain}.pdb --dat1 . --dat2 . --outfmt summary,alignments,transrot >& log')
            fp = os.popen('ls -1 mol*.txt')
            info = fp.readlines()
            fp.close()
            filenames = []
            for line in info:
                filenames.append(line[:-1])

            if filenames:
                fp = open('../'+prot + '__' + edomain + '.pdb', 'r')
                Qresids_set = set([])
                for line in fp:
                    if line[:4]=='ATOM':
                        resid = int(line[22:26])
                        Qresids_set.add(resid)
                fp.close()
                Qresids = list(Qresids_set)

                info = []
                for filename in filenames:
                    fp = open(filename, "r")
                    lines = fp.readlines()
                    fp.close()
                    for line in lines:
                        info.append(line)

                qali = ''
                sali = ''
                getit = 1
                zscore = 0
                for line in info:
                    words = line.split()
                    if len(words) >= 2 and getit:
                        if words[0] == 'Query':
                            qali += words[1]
                        elif words[0] == 'Sbjct':
                            sali += words[1]
                        elif words[0] == 'No' and words[1] == '1:':
                            for word in words:
                                if '=' in word:
                                    subwords = word.split('=')
                                    if subwords[0] == 'Z-score':
                                        zinfo = subwords[1].split('.')
                                        zscore = float(zinfo[0] + '.' + zinfo[1])
                        elif words[0] == 'No' and words[1] == '2:':
                            getit = 0

                print (zscore, qali, sali)
                qinds = []
                sinds = []
                length = len(qali)
                qposi = 0
                sposi = 0
                match = 0
                for i in range(length):
                    if qali[i] != '-':
                        qposi += 1
                    if sali[i] != '-':
                        sposi += 1
                    if qali[i] != '-' and sali[i] != '-':
                        if qali[i].isupper() and sali[i].isupper():
                            match += 1
                            qinds.append(qposi)
                            sinds.append(sposi)
                qlen = qposi
                slen = sposi

                if match >= 20:
                    alicount += 1
                    rp = open(cwd+'/'+prot + '_iterative_hits', 'a')
                    rp.write('>' + edomain + '_' + str(alicount) + '\t' + str(zscore) + '\t' + str(match) + '\t' + str(qlen) + '\t' + str(slen) + '\n')
                    for i in range(len(qinds)):
                        qind = qinds[i] - 1
                        sind = sinds[i]
                        rp.write(str(Qresids[qind]) + '\t' + str(sind) + '\n')
                    rp.close()

                    raw_qresids = []
                    for qind in qinds:
                        raw_qresids.append(Qresids[qind - 1])
                    qrange = get_domain_range(raw_qresids)
                    qresids = set([])
                    qsegs = qrange.split(',')
                    for qseg in qsegs:
                        qedges = qseg.split('-')
                        qstart = int(qedges[0])
                        qend = int(qedges[1])
                        for qres in range(qstart, qend + 1):
                            qresids.add(qres)
                    remain_resids = Qresids_set.difference(qresids)
    
                    if len(remain_resids) >= 40:
                        rp = open(cwd+'/iterative_dali_tmp/'+prot + '__' + edomain + '.pdbnew', 'w')
                        fp = open(cwd+'/iterative_dali_tmp/'+prot + '__' + edomain + '.pdb', 'r')
                        for line in fp:
                            if line[:4]=='ATOM':
                                resid = int(line[22:26])
                                if resid in remain_resids:
                                    rp.write(line)
                        fp.close()
                        rp.close()
                        os.system('mv '+cwd+'/iterative_dali_tmp/'+prot + '__' + edomain + '.pdbnew '+cwd+'/iterative_dali_tmp/'+prot + '__' + edomain + '.pdb')
                        os.system('rm  '+cwd+'/iterative_dali_tmp/tmp/*')
                    else:
                        os.system('rm '+cwd+'/iterative_dali_tmp/tmp/*')
                        break
                else:
                    os.system('rm '+cwd+'/iterative_dali_tmp/tmp/*')
                    break
            else:
                os.system('rm '+cwd+'/iterative_dali_tmp/tmp/*')
                break
    os.system('rmdir '+cwd+'/iterative_dali_tmp/tmp')
#    os.system('echo \'done\' > /gscratch/scrubbed/congq/ECOD_domains/step13/' + prot + '.done')
    os.chdir(cwd)
with open('log','a') as f:
    f.write('5\n')
