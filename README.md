# DPAM
A **D**omain **P**arser for **A**lphafold **M**odels 

DPAM: A Domain Parser for AlphaFold Models （https://www.biorxiv.org/content/10.1101/2022.09.22.509116v1, accepted by Protein Science)

## Updates:
A docker image can be dowloaded by **docker pull conglab/dpam:latest** (this is an enhanced version of current DPAM, we will soon update the repository too)

Upload domain parser results for six model organisms.  (2022-12-6)

Replace Dali with Foldseek for initial hits searching. (2022-11-30)

Fix a bug in analyze_PDB.py which prevents the proper usage of Dali results. (2022-10-31)
## Prerequisites:

### Software and packages
- HH-suite3: https://github.com/soedinglab/hh-suite (enable addss.pl to add secondary structure)
- DaliLite.v5: http://ekhidna2.biocenter.helsinki.fi/dali/
- Python 3.8 
- Foldseek 
- mkdssp
- pdbx: https://github.com/soedinglab/pdbx
- pdb2fasta (https://zhanggroup.org/pdb2fasta)

Please add above software to environment path for DPAM. We also provide a script `check_dependencies.py` to check if above programs can be found. 
### Supporting database:
- hhsearch UniRef database (https://wwwuser.gwdg.de/~compbiol/uniclust/2022_02/)
- pdb70 (https://conglab.swmed.edu/DPAM/pdb70.tgz)
- ECOD database 
  - ECOD ID map to pdb
  - ECOD domain length
  - ECOD domain list
  - ECOD norms 
  - ECOD domain quality information
  - ECOD residue weight in domains 
  - ECOD70 domain structures 
  - ECOD70 foldseek database

We provide a script download_all_data.sh that can be used to download all of these databases.

`bash download_all_data.sh <DOWNLOAD_DIR>`

After downloading the databases, please decompress files. All supporting database files should be put in the same directory and the directory should be provided to `DPAM.py` as `<datadir>`. The `<datadir>` should have the following structure and files. 
  
    <datadir>/
        ECOD70/
        ecod_domain_info/
        ECOD_foldseek_DB/
        ecod_weights/
        pdb70/
        UniRef30_2022_02/
        ecod.latest.domains
        ECOD_length
        ECOD_norms
        ECOD_pdbmap
    

## Installation
git clone https://github.com/CongLabCode/DPAM.git

conda install -c qianlabcode dpam

## Usage
`python DPAM.py <input_cif/pdb> <input_pae> <accession> <output_dir> <threads> <datadir>`

## Future improvments
- Incoperate mmseq to improve search speed
- Provide public server and incoperate with ECOD database 
