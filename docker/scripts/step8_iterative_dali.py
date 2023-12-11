#!/opt/conda/bin/python
import os, sys
import time, subprocess
from multiprocessing import Pool

wdir = os.getcwd()

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

def run_dali(edomain):
    os.system('cp ' + wdir + '/step2/' + dataset + '/' + prot + '.pdb /tmp/' + prot + '_' + edomain + '.pdb')
    alicount = 0
    os.system('mkdir /tmp/tmp_' + prot + '_' + edomain)
    os.chdir('/tmp/tmp_' + prot + '_' + edomain)
    while 1:
        rcode=subprocess.run('dali.pl --pdbfile1 /tmp/' + prot + '_' + edomain + '.pdb --pdbfile2 /mnt/databases/ECOD_pdbs/' + edomain + '.pdb --dat1 ./ -dat2 ./ --outfmt summary,alignments,transrot > log',shell=True).returncode
        if str(rcode) != '0':
            sys.exit(1)
        fp = os.popen('ls -1 mol*.txt')
        info = fp.readlines()
        fp.close()
        filenames = []
        for line in info:
            filenames.append(line[:-1])

        if filenames:
            fp = open('/tmp/' + prot + '_' + edomain + '.pdb', 'r')
            Qresids_set = set([])
            for line in fp:
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
            rotdata = []
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
                    elif len(words) >= 8 and words[0] == "-matrix":
                        rotdata.append(line)
                    elif words[0] == 'No' and words[1] == '2:':
                        getit = 0

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
                rp = open('/tmp/' + prot + '_' + edomain + '_hits', 'a')
                rp.write('>' + edomain + '_' + str(alicount) + '\t' + str(zscore) + '\t' + str(match) + '\t' + str(qlen) + '\t' + str(slen) + '\n')
                rotation = ''
                translation = 'translation'
                for rotline in rotdata:
                    rotwords = rotline[:-2].split()
                    rotation += 'rotation\t' + '\t'.join(rotwords[4:7]) + '\n'
                    translation += '\t' + rotwords[7]
                rp.write(rotation)
                rp.write(translation + '\n')
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

                if len(remain_resids) >= 20:
                    rp = open('/tmp/' + prot + '_' + edomain + '.pdbnew', 'w')
                    fp = open('/tmp/' + prot + '_' + edomain + '.pdb', 'r')
                    for line in fp:
                        resid = int(line[22:26])
                        if resid in remain_resids:
                            rp.write(line)
                    fp.close()
                    rp.close()
                    os.system('mv /tmp/' + prot + '_' + edomain + '.pdbnew /tmp/' + prot + '_' + edomain + '.pdb')
                    if os.getcwd() == '/tmp/tmp_' + prot + '_' + edomain:
                        pass
                        os.system('rm *')
                else:
                    if os.getcwd() == '/tmp/tmp_' + prot + '_' + edomain:
                        pass
                        os.system('rm *')
                    break
            else:
                if os.getcwd() == '/tmp/tmp_' + prot + '_' + edomain:
                    pass
                    os.system('rm *')
                break
        else:
            if os.getcwd() == '/tmp/tmp_' + prot + '_' + edomain:
                pass
                os.system('rm *')
            break
    os.chdir('/tmp/')
    os.system('rmdir /tmp/tmp_' + prot + '_' + edomain)
    os.system('rm /tmp/' + prot + '_' + edomain + '.pdb')


dataset = sys.argv[1]
prot = sys.argv[2]
ncore = int(sys.argv[3])

if os.path.exists(wdir + '/step7/' + dataset + '/' + prot + '_hits'):
    fp = open(wdir + '/step7/' + dataset + '/' + prot + '_hits','r')
    edomains = []
    for line in fp:
        words = line.split()
        edomains.append(words[0])
    fp.close()

    inputs = []
    for edomain in edomains:
        inputs.append([edomain])
    pool = Pool(processes = ncore)
    results = []
    for item in inputs:
        process = pool.apply_async(run_dali, item)
        results.append(process)
    for process in results:
        process.get()
    fp = os.popen('ls -1 /tmp/' + prot + '_*_hits')
    get_hit_count = 0
    for line in fp:
        get_hit_count += 1
    fp.close()
    if get_hit_count:
        os.system('cat /tmp/' + prot + '_*_hits > ' + wdir + '/step8/' + dataset + '/' + prot + '_hits')
        os.system('rm /tmp/' + prot + '_*_hits')
    else:
        os.system('echo \'done\' > ' + wdir + '/step8/' + dataset + '/' + prot + '.done')
else:
    os.system('echo \'done\' > ' + wdir + '/step8/' + dataset + '/' + prot + '.done')
