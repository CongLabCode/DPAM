# DPAM
A **D**omain **P**arser for **A**lphafold **M**odels

## Updates:
Replace Dali with Foldseek for initial hits searching. (2022-11-30)
Fix a bug in analyze_PDB.py which prevents the proper usage of Dali results. (2022-10-31)
## Prerequisites:
### Software and packages
- HH-suite3: https://github.com/soedinglab/hh-suite (enable addss.pl to add secondary structure)
- DaliLite.v5: http://ekhidna2.biocenter.helsinki.fi/dali/
- Python 3.8 
- Foldseek 
- pdbx: https://github.com/soedinglab/pdbx

Please add above software to environment path for DPAM. 
### Supporting database:
- hhsearch UniRef database (https://wwwuser.gwdg.de/~compbiol/uniclust/2022_02/)
- pdb70 (https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)
- ECOD database 
  - ECOD ID map to pdb
  - ECOD domain length
  - ECOD domain list (need to decompress)
  - ECOD norms 
  - ECOD domain quality information
  - ECOD residue weight in domains 
  - ECOD domain structures 

We provide a script scripts/download_all_data.sh that can be used to download all of these databases.

`bash scripts/download_all_data.sh <DOWNLOAD_DIR>`

After downloading the databases, please decompress files. All supporting database files should be put in the same directory and the directory should be provided to `DPAM.py` as `<datadir>`. The <datadir> should have the following structure and files. 
`<datadir>'

## Installation
git clone https://github.com/CongLabCode/DPAM.git


## Usage
`python DPAM.py <input_cif/pdb> <input_pae> <accession> <output_dir> <threads> <datadir>`

## Future improvments
- Incoperate mmseq to improve search speed
- Provide public server and incoperate with ECOD database 
