#!/usr1/local/bin/python
import os, sys, string
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

filename = sys.argv[1]
cwd=os.getcwd()
chain2lines = {}
if os.path.exists(filename + ".cif") and os.path.exists(filename + ".fa"):
    chain2seq = {}
    fp = open(filename + ".fa", "r")
    chainid = ""
    seq = ""
    for line in fp:
        if line[0] == ">":
            if chainid and seq:
                chain2seq[chainid] = seq
            chainid = line[1:].split()[0].split("_")[-1]
            seq = ""
        else:
            seq += line[:-1]
    fp.close()
    if chainid and seq:
        chain2seq[chainid] = seq

    cif = open(filename + ".cif", "r")
    pRd = PdbxReader(cif)
    data = []
    pRd.read(data)
    block = data[0]
    good_entities = []
    good_chains = []

    entity_poly = block.getObj("entity_poly")
    pdbx_poly_seq_scheme = block.getObj("pdbx_poly_seq_scheme")
    if pdbx_poly_seq_scheme and entity_poly:
        typeid = entity_poly.getIndex("type")
        entityid1 = entity_poly.getIndex("entity_id")
        entityid2 = pdbx_poly_seq_scheme.getIndex("entity_id")
        chainid = pdbx_poly_seq_scheme.getIndex("asym_id")
        for i in range(entity_poly.getRowCount()):
            words = entity_poly.getRow(i)
            entity = words[entityid1]
            chaintype = words[typeid]
            if chaintype == "polypeptide(L)":
                good_entities.append(entity)

    if good_entities:
        good_chains = []
        for i in range(pdbx_poly_seq_scheme.getRowCount()):
            words = pdbx_poly_seq_scheme.getRow(i)
            entity = words[entityid2]
            if entity in good_entities:
                chain = words[chainid]
                if not chain in good_chains:
                    good_chains.append(chain)

    if good_chains:
        chain2lines = {}
        chain2models = {}
        for chain in good_chains:
            chain2lines[chain] = {}
            chain2models[chain] = []
        chain2select = {}

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
            for chain in good_chains:
                chain2lines[chain][1] = []
            for i in range(atom_site.getRowCount()):
                words = atom_site.getRow(i)
                chain_id = words[chain_id_index]
                record_type = words[record_type_index]
                if chain_id in good_chains:
                    if record_type == "ATOM":
                        chain2lines[chain_id][1].append(words)
            for chain in good_chains:
                chain2select[chain] = 1
        else:
            for i in range(atom_site.getRowCount()):
                words = atom_site.getRow(i)
                chain_id = words[chain_id_index]
                record_type = words[record_type_index]
                model_num = int(words[model_num_index])
                if chain_id in good_chains:
                    if record_type == "ATOM":
                        try:
                            chain2lines[chain_id][model_num].append(words)
                        except KeyError:
                            chain2lines[chain_id][model_num] = [words]
                            chain2models[chain_id].append(model_num)
            for chain in good_chains:
                if chain2models[chain]:
                    chain2select[chain] = min(chain2models[chain])

        for chain in good_chains:
            if chain2models[chain]:
                goodlines = []
                rp = open(cwd+'/'+filename + ".pdb","w")
                resid2altid = {}
                resid2aa = {}
                atom_count = 0
                for words in chain2lines[chain][chain2select[chain]]:
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
                        if len(chain) < 3:
                            if len(atom_identity) < 4:
                                goodlines.append("ATOM  " + str(atom_count).rjust(5) + "  " +  atom_identity.ljust(3) + " " + residue_type.ljust(3) + chain.rjust(2) + str(residue_id).rjust(4) + "    " + coor_x.rjust(8) + coor_y.rjust(8) + coor_z.rjust(8) + "  1.00  0.00           " + atom_type + "\n")
                            elif len(atom_identity) == 4:
                                goodlines.append("ATOM  " + str(atom_count).rjust(5) + " " +  atom_identity + " " + residue_type.ljust(3) + chain.rjust(2) + str(residue_id).rjust(4) + "    " + coor_x.rjust(8) + coor_y.rjust(8) + coor_z.rjust(8) + "  1.00  0.00           " + atom_type + "\n")
                        else:
                            if len(atom_identity) < 4:
                                goodlines.append("ATOM  " + str(atom_count).rjust(5) + "  " +  atom_identity.ljust(3) + " " + residue_type.ljust(3) + " ." + str(residue_id).rjust(4) + "    " + coor_x.rjust(8) + coor_y.rjust(8) + coor_z.rjust(8) + "  1.00  0.00           " + atom_type + "\n")
                            elif len(atom_identity) == 4:
                                goodlines.append("ATOM  " + str(atom_count).rjust(5) + " " +  atom_identity + " " + residue_type.ljust(3) + " ." + str(residue_id).rjust(4) + "    " + coor_x.rjust(8) + coor_y.rjust(8) + coor_z.rjust(8) + "  1.00  0.00           " + atom_type + "\n")

                newseq = ""
                oldseq = chain2seq[chain]
                for i in range(len(oldseq)):
                    resid = i + 1
                    try:
                        newseq += resid2aa[resid]
                        if resid2aa[resid] == "X":
                            pass
                        elif resid2aa[resid] == oldseq[i]:
                            pass
                        else:
                            print ("error " + filename + " " + chain)
                    except KeyError:
                        newseq += "-"
#                rp.write("#" + oldseq + "\n")
#                rp.write("#" + newseq + "\n")
                for goodline in goodlines:
                    rp.write(goodline)
                rp.close()
