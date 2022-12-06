#!/usr1/local/bin/python
import os, sys, string
import pdbx
from pdbx.reader.PdbxReader import PdbxReader

three2one = {}
three2one["ALA"] = "A"
three2one["CYS"] = "C"
three2one["ASP"] = "D"
three2one["GLU"] = "E"
three2one["PHE"] = "F"
three2one["GLY"] = "G"
three2one["HIS"] = "H"
three2one["ILE"] = "I"
three2one["LYS"] = "K"
three2one["LEU"] = "L"
three2one["MET"] = "M"
three2one["MSE"] = "M"
three2one["ASN"] = "N"
three2one["PRO"] = "P"
three2one["GLN"] = "Q"
three2one["ARG"] = "R"
three2one["SER"] = "S"
three2one["THR"] = "T"
three2one["VAL"] = "V"
three2one["TRP"] = "W"
three2one["TYR"] = "Y"

prefix = sys.argv[1]
wd=sys.argv[2]
if os.getcwd() != wd:
    os.chdir(wd)

if os.path.exists(prefix+'.cif') and os.path.exists(prefix + ".fa"):
    fp = open(prefix + ".fa", "r")
    myseq = ""
    for line in fp:
        if line[0] == ">":
            pass
        else:
            myseq += line[:-1]
    fp.close()

    cif = open(prefix + ".cif", "r")
    pRd = PdbxReader(cif)
    data = []
    pRd.read(data)
    block = data[0]

    atom_site = block.getObj("atom_site")
    record_type_index = atom_site.getIndex("group_PDB")
    atom_type_index = atom_site.getIndex("type_symbol")
    atom_identity_index = atom_site.getIndex("label_atom_id")
    residue_type_index = atom_site.getIndex("label_comp_id")
    chain_id_index = atom_site.getIndex("label_asym_id")
    residue_id_index = atom_site.getIndex("label_seq_id")
    coor_x_index = atom_site.getIndex("Cartn_x")
    coor_y_index = atom_site.getIndex("Cartn_y")
    coor_z_index = atom_site.getIndex("Cartn_z")
    alt_id_index = atom_site.getIndex("label_alt_id")
    model_num_index = atom_site.getIndex("pdbx_PDB_model_num")

    if model_num_index == -1:
        mylines = []
        for i in range(atom_site.getRowCount()):
            words = atom_site.getRow(i)
            chain_id = words[chain_id_index]
            record_type = words[record_type_index]
            if chain_id == "A" and record_type == "ATOM":
                mylines.append(words)
    else:
        model2lines = {}
        models = []
        for i in range(atom_site.getRowCount()):
            words = atom_site.getRow(i)
            chain_id = words[chain_id_index]
            record_type = words[record_type_index]
            model_num = int(words[model_num_index])
            if chain_id == "A" and record_type == "ATOM":
                try:
                    model2lines[model_num].append(words)
                except KeyError:
                    model2lines[model_num] = [words]
                    models.append(model_num)
        best_model = min(models)
        mylines = model2lines[best_model]

    goodlines = []
    resid2altid = {}
    resid2aa = {}
    atom_count = 0
    for words in mylines:
        atom_type = words[atom_type_index]
        atom_identity = words[atom_identity_index]
        residue_type = words[residue_type_index]
        residue_id = int(words[residue_id_index])
        alt_id = words[alt_id_index]

        if atom_identity == "CA":
            try:
                resid2aa[residue_id] = three2one[residue_type]
            except KeyError:
                resid2aa[residue_id] = "X"

        get_line = 0
        if alt_id == ".":
            get_line = 1
        else:
            try:
                if resid2altid[residue_id] == alt_id:
                    get_line = 1
                else:
                    get_line = 0
            except KeyError:
                resid2altid[residue_id] = alt_id
                get_line = 1

        if get_line:
            atom_count += 1
            coor_x_info = words[coor_x_index].split(".")
            if len(coor_x_info) >= 2:
                coor_x = coor_x_info[0] + "." + coor_x_info[1][:3]
            else:
                coor_x = coor_x_info[0]
            coor_y_info = words[coor_y_index].split(".")
            if len(coor_y_info) >= 2:
                coor_y = coor_y_info[0] + "." + coor_y_info[1][:3]
            else:
                coor_y = coor_y_info[0]
            coor_z_info = words[coor_z_index].split(".")
            if len(coor_z_info) >= 2:
                coor_z = coor_z_info[0] + "." + coor_z_info[1][:3]
            else:
                coor_z = coor_z_info[0]
            if len(atom_identity) < 4:
                goodlines.append("ATOM  " + str(atom_count).rjust(5) + "  " +  atom_identity.ljust(3) + " " + residue_type.ljust(3) + " A" + str(residue_id).rjust(4) + "    " + coor_x.rjust(8) + coor_y.rjust(8) + coor_z.rjust(8) + "  1.00  0.00           " + atom_type + "\n")
            elif len(atom_identity) == 4:
                goodlines.append("ATOM  " + str(atom_count).rjust(5) + " " +  atom_identity + " " + residue_type.ljust(3) + " A" + str(residue_id).rjust(4) + "    " + coor_x.rjust(8) + coor_y.rjust(8) + coor_z.rjust(8) + "  1.00  0.00           " + atom_type + "\n")

    newseq = ""
    for i in range(len(myseq)):
        resid = i + 1
        try:
            newseq += resid2aa[resid]
            if resid2aa[resid] == "X":
                pass
            elif resid2aa[resid] == myseq[i]:
                pass
            else:
                print ("error\t" + prefix)
        except KeyError:
            newseq += "-"
    if newseq == myseq:
        rp = open(prefix + ".pdb","w")
        for goodline in goodlines:
            rp.write(goodline)
        rp.close()
    else:
        sys.exit(1)
elif os.path.exists(prefix + ".pdb") and os.path.exists(prefix + ".fa"):
    fp = open(prefix + ".fa", "r")
    myseq = ""
    for line in fp:
        if line[0] == ">":
            pass
        else:
            myseq += line[:-1]
    fp.close()

    pdb_new = []
    res = []
    f = open(prefix + ".pdb")
    for line in f:
        if line[:4] == 'ATOM':
            pdb_new.append(line)
            if line[13:16].strip()=='CA':
                try:
                    res.append(three2one[line[17:20]])
                except:
                    res.append('X')
    f.close()
    if ''.join(res) == myseq:
        rp = open(prefix + ".pdb","w")
        for i in pdb_new:
            rp.write(i)
        rp.close()
    else:
        sys.exit(1)
else:
    sys.exit(1)
