import shutil
hhsuite=['hhblits','hhsearch','hhmake','addss.pl']
programs=['hhblits','hhsearch','hhmake','addss.pl','foldseek','dali.pl','mkdssp']
missing=[]
pdbx=0
for prog in programs:
    check = shutil.which(prog)
    if check == None:
        missing.append(prog)
try:
    import pdbx
    from pdbx.reader.PdbxReader import PdbxReader
except:
    pdbx = 1


if missing or pdbx == 1:
    if missing:
        text = "Please add"
        hhsuite_missing = [i for i in missing if i in hhsuite]
        if hhsuite_missing:
            if len(hhsuite_missing) >= 2:
                hhsuite_missing = ','.join(hhsuite_missing[:-1]) +' and '+hhsuite_missing[-1] + " in HH-suite3"
            else:
                hhsuite_missing = hhsuite_missing[0] + " in HH-suite3"
            text = text + " " + hhsuite_missing
        others = [i for i in missing if i not in hhsuite_missing]
        if others:
            if len(others) >= 2:
                others = ','.join(others[:-1]) + ' and ' + others[-1]
            else:
                others = others[0]
            text = text + " and " + others
        text = text + " to envirnoment path"
        print(text)
    if pdbx == 1:
        print('pdbx is not installed properly.Please refer to https://github.com/soedinglab/pdbx for installation')
else:
    print('HH-suite, Foldseek and dali.pl are found')
