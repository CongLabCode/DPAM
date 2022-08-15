# DPAM
A domain parser for Alphafold models
## Prerequisites:
### Software
- HH-suite3: https://github.com/soedinglab/hh-suite (enable addss.pl to add secondary structure)
- DaliLite.v5: http://ekhidna2.biocenter.helsinki.fi/dali/
- Python 3.8 
- pdbx: https://github.com/soedinglab/pdbx
Above software is required to be added to $PATH for smooth execution of the DPAM. 
### Database:
- hhsearch UniRef database (https://wwwuser.gwdg.de/~compbiol/uniclust/2022_02/)
- pdb70 (https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)
- ECOD database 
### Other auxiliary data required for the software:
- ECOD map to pdb
- ECOD domain length
- ECOD domain list
- ECOD norms
- ECOD domain quality
- ECOD residue weight

## Installation
After installing required software and downloading the DPAM and necessary auxiliary data, please modify the config_file so DPAM can access the required data. 

## Usage
DPAM.py [model name] [output_dir]

## Future improvments
- Incoperate mmseq and foldseek to accelerate speed of the search
- Provide server for public usage and integrate with ECOD
