import os, sys
import time
from multiprocessing import Pool

prefix = sys.argv[1]
CPUs = sys.argv[2]
wd = sys.argv[3]
data_dir=sys.argv[4]


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
    alicount = 0
    os.system(f'mkdir {wd}/iterativeDali_{prefix}/tmp_{prefix}_{edomain}')
    os.system(f'cp {wd}/{prefix}.pdb {wd}/iterativeDali_{prefix}/tmp_{prefix}_{edomain}/{prefix}_{edomain}.pdb')
    os.system(f'mkdir {wd}/iterativeDali_{prefix}/tmp_{prefix}_{edomain}/output_tmp')
    os.chdir(f'{wd}/iterativeDali_{prefix}/tmp_{prefix}_{edomain}/output_tmp')
    while True:
        os.system(f'dali.pl --pdbfile1 ../{prefix}_{edomain}.pdb --pdbfile2 {data_dir}/ECOD70/{edomain}.pdb --dat1 ./ -dat2 ./ --outfmt summary,alignments,transrot >& log')
        fp = os.popen('ls -1 mol*.txt')
        info = fp.readlines()
        fp.close()
        filenames = []
        for line in info:
            filenames.append(line[:-1])

        if filenames:
            fp = open('../' + prefix + '_' + edomain + '.pdb', 'r')
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
                rp = open(f'{wd}/iterativeDali_{prefix}/{prefix}_{edomain}_hits', 'a')
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

                if len(remain_resids) >= 20:
                    rp = open('../' + prefix + '_' + edomain + '.pdbnew', 'w')
                    fp = open('../' + prefix + '_' + edomain + '.pdb', 'r')
                    for line in fp:
                        resid = int(line[22:26])
                        if resid in remain_resids:
                            rp.write(line)
                    fp.close()
                    rp.close()
                    os.system('mv ../' + prefix + '_' + edomain + '.pdbnew ../' + prefix + '_' + edomain + '.pdb')
                    if os.getcwd() == f'{wd}/iterativeDali_{prefix}/tmp_{prefix}_{edomain}/output_tmp':
                        os.system('rm *')
                else:
                    if os.getcwd() == f'{wd}/iterativeDali_{prefix}/tmp_{prefix}_{edomain}/output_tmp':
                        os.system('rm *')
                    break
            else:
                if os.getcwd() == f'{wd}/iterativeDali_{prefix}/tmp_{prefix}_{edomain}/output_tmp':
                    os.system('rm *')
                break
        else:
            if os.getcwd() == f'{wd}/iterativeDali_{prefix}/tmp_{prefix}_{edomain}/output_tmp':
                os.system('rm *')
            break
    os.chdir(wd)
    time.sleep(1)
    os.system(f'rm -rf {wd}/iterativeDali_{prefix}/tmp_{prefix}_{edomain}')

if os.getcwd() != wd:
    os.chdir(wd)

if os.path.exists(prefix + '.iterativeDali.done'):
    pass
else:
    if not os.path.exists(f'{wd}/iterativeDali_{prefix}'):
        os.system(f'mkdir {wd}/iterativeDali_{prefix}')
    fp = open(prefix + '_hits4Dali','r')
    edomains = []
    for line in fp:
        words = line.split()
        edomains.append(words[0])
    fp.close()

    inputs = []
    for edomain in edomains:
        inputs.append([edomain])
    pool = Pool(processes = int(CPUs))
    results = []
    for item in inputs:
        process = pool.apply_async(run_dali, item)
        results.append(process)
    for process in results:
        process.get()

    os.system(f'cat {wd}/iterativeDali_{prefix}/{prefix}_*_hits > {prefix}_iterativdDali_hits')
    os.system(f'rm -rf {wd}/iterativeDal_{prefix}/tmp_*')
    os.system(f'rm -rf {wd}/iterativeDal_{prefix}')
    os.system(f'echo "done" > {prefix}.iterativeDali.done')
