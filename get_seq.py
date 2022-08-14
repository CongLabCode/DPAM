#!/usr1/local/bin/python
import os, sys
import subprocess
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

prot = sys.argv[1]
mFormat=sys.argv[2]
if mFormat=='cif':
#cif = open("CIFs/" + myfile + ".cif")
    cif = open(prot+'.cif')
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
	    newchainid = pdbx_poly_seq_scheme.getIndex("asym_id")
	    oldchainid = pdbx_poly_seq_scheme.getIndex("pdb_strand_id")
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
		    mymap = {}
		    residues = {}
		    seqs = {}
		    rp = open(prot + ".fa","w")
		    for i in range(pdbx_poly_seq_scheme.getRowCount()):
			    words = pdbx_poly_seq_scheme.getRow(i)
			    entity = words[entityid2]
			    if entity in good_entities:
				    chain1 = words[newchainid]
				    chain2 = words[oldchainid]

				    try:
					    aa = three2one[words[resiid]]
				    except KeyError:
					    try:
						    modinfo[chain1][words[posiid]]
						    resiname = modinfo[chain1][words[posiid]][0]
						    if words[resiid] == resiname:
							    new_resiname = modinfo[chain1][words[posiid]][1]
							    try:
								    aa = three2one[new_resiname]
							    except KeyError:
								    aa = "X"
								    print ("error2 " + new_resiname)
						    else:
							    aa = "X"
							    print ("error1 " + words[resiid] + " " + resiname)
					    except KeyError:
						    print (modinfo)
						    print (words[resiid])
						    aa = "X"

				    try:
					    seqs[chain1]
				    except KeyError:
					    seqs[chain1] = {}

				    try:
					    if seqs[chain1][int(words[posiid])] == "X" and aa != "X":
						    seqs[chain1][int(words[posiid])] = aa
				    except KeyError:
					    seqs[chain1][int(words[posiid])] = aa

				    try:
					    residues[chain1].add(int(words[posiid]))
					    if not chain2 in mymap[chain1]:
						    mymap[chain1].append(chain2)
				    except KeyError:
					    residues[chain1] = set([int(words[posiid])])
					    mymap[chain1] = [chain2]

		    for chain1 in mymap.keys():
			    if len(mymap[chain1]) > 1:
				    print ("error1 " + prot + " " + chain1)
			    else:
				    for i in range(len(residues[chain1])):
					    if not i + 1 in residues[chain1]:
						    print ("error2 " + prot + " " + chain1)
						    break
				    else:
					    chain2 = mymap[chain1][0]
					    rp.write(">" + prot + "_" + chain1 + " " + prot + "_" + chain2 + "\n")
					    finalseq = []
					    for i in range(len(residues[chain1])):
						    finalseq.append(seqs[chain1][i+1])
					    rp.write("".join(finalseq) + "\n")
		    rp.close()
	    else:
		    print ("empty " + prot)
    else:
	    print ("bad " + prot)
else:
    status=subprocess.call(f'{script_dir}/pdb2fasta {prot}.pdb > {prot}.fa')
    if status==0:
        with open(prot+'.fa') as f:
            fasta=f.readlines()
        fasta[0]=fasta[0].split()[0].replace(':','_')+'\n'
        with open(prot+'.fa','w')as f:
            for i in fasta:
                f.write(i)
    else:
        print('cannot get sequence of the file')
        sys.exit()
