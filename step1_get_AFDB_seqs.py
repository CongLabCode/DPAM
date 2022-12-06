#!/usr1/local/bin/python
import os, sys
import pdbx
from pdbx.reader.PdbxReader import PdbxReader

three2one = {}
three2one["ALA"] = 'A'
three2one["CYS"] = 'C'
three2one["ASP"] = 'D'
three2one["GLU"] = 'E'
three2one["PHE"] = 'F'
three2one["GLY"] = 'G'
three2one["HIS"] = 'H'
three2one["ILE"] = 'I'
three2one["LYS"] = 'K'
three2one["LEU"] = 'L'
three2one["MET"] = 'M'
three2one["MSE"] = 'M'
three2one["ASN"] = 'N'
three2one["PRO"] = 'P'
three2one["GLN"] = 'Q'
three2one["ARG"] = 'R'
three2one["SER"] = 'S'
three2one["THR"] = 'T'
three2one["VAL"] = 'V'
three2one["TRP"] = 'W'
three2one["TYR"] = 'Y'


prefix = sys.argv[1]
wd = sys.argv[2]
if os.getcwd() != wd:
    os.chdir(wd)


struc_fn = prefix + '.cif'
if os.path.exists(struc_fn):
    cif = open(prefix + ".cif")
else:
    if os.path.exists(prefix + ".pdb"):
        os.system(f'pdb2fasta '+ prefix + ".pdb > " + prefix + ".fa")
        with open(prefix + ".fa") as f:
            fa = f.readlines()
        fa[0] = fa[0].split(':')[0] + '\n'
        with open(prefix + ".fa",'w') as f:
            f.write(''.join(fa))
    else:
        print("No recognized structure file (*.cif or *.pdb). Existing...")
    sys.exit()

pRd = PdbxReader(cif)
data = []
pRd.read(data)
block = data[0]

modinfo = {}
mod_residues = block.getObj("pdbx_struct_mod_residue")
if mod_residues:
    chainid = mod_residues.getIndex("label_asym_id")
    posiid = mod_residues.getIndex("label_seq_id")
    parentid = mod_residues.getIndex("parent_comp_id")
    resiid = mod_residues.getIndex("label_comp_id")
    for i in range(mod_residues.getRowCount()):
        words = mod_residues.getRow(i)
        try:
            modinfo[words[chainid]]
        except KeyError:
            modinfo[words[chainid]] = {}
        modinfo[words[chainid]][words[posiid]] = [words[resiid], words[parentid]]

entity_poly = block.getObj("entity_poly")
pdbx_poly_seq_scheme = block.getObj("pdbx_poly_seq_scheme")
if pdbx_poly_seq_scheme and entity_poly:
    typeid = entity_poly.getIndex("type")
    entityid1 = entity_poly.getIndex("entity_id")
    entityid2 = pdbx_poly_seq_scheme.getIndex("entity_id")
    chainid = pdbx_poly_seq_scheme.getIndex("asym_id")
    resiid = pdbx_poly_seq_scheme.getIndex("mon_id")
    posiid = pdbx_poly_seq_scheme.getIndex("seq_id")

    good_entities = []
    for i in range(entity_poly.getRowCount()):
        words = entity_poly.getRow(i)
        entity = words[entityid1]
        type = words[typeid]
        if type == "polypeptide(L)":
             good_entities.append(entity)

    if good_entities:
        chains = []
        residues = {}
        seqs = {}
        rp = open(prefix + ".fa","w")
        for i in range(pdbx_poly_seq_scheme.getRowCount()):
            words = pdbx_poly_seq_scheme.getRow(i)
            entity = words[entityid2]
            if entity in good_entities:
                chain = words[chainid]

                try:
                    aa = three2one[words[resiid]]
                except KeyError:
                    try:
                        modinfo[chain][words[posiid]]
                        resiname = modinfo[chain][words[posiid]][0]
                        if words[resiid] == resiname:
                            new_resiname = modinfo[chain][words[posiid]][1]
                            try:
                                aa = three2one[new_resiname]
                            except KeyError:
                                aa = "X"
                                print ("error1 " + new_resiname)
                        else:
                            aa = "X"
                            print ("error2 " + words[resiid] + " " + resiname)
                    except KeyError:
                        print (modinfo)
                        print (words[resiid])
                        aa = "X"

                try:
                    seqs[chain]
                except KeyError:
                    chains.append(chain)
                    seqs[chain] = {}

                try:
                    if seqs[chain][int(words[posiid])] == "X" and aa != "X":
                        seqs[chain][int(words[posiid])] = aa
                except KeyError:
                    seqs[chain][int(words[posiid])] = aa

                try:
                    residues[chain].add(int(words[posiid]))
                except KeyError:
                    residues[chain] = set([int(words[posiid])])

        for chain in chains:
            for i in range(len(residues[chain])):
                if not i + 1 in residues[chain]:
                    print ("error3 " + prefix + " " + chain)
                    break
            else:
                rp.write(">" + prefix + "\n")
                finalseq = []
                for i in range(len(residues[chain])):
                    finalseq.append(seqs[chain][i+1])
                rp.write("".join(finalseq) + "\n")
        rp.close()
    else:
        print ("empty " + prefix)
else:
    print ("bad " + prefix)
