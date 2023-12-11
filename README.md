# DPAM
A **D**omain **P**arser for **A**lphafold **M**odels 

DPAM: A Domain Parser for AlphaFold Models （https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4548)

## Updates:
A docker image for DPAM v2.0 can be dowloaded by **docker pull conglab/dpam:latest** (2023-12-10). New version includes domain classification based on ECOD database and addresses over-segmentation for some proteins. 

Upload domain parser results for six model organisms.  (2022-12-6)

Replace Dali with Foldseek for initial hits searching. (2022-11-30)

Fix a bug in analyze_PDB.py which prevents the proper usage of Dali results. (2022-10-31)
## Prerequisites: 
docker image conglab/dpam:latest 
run_dpam_docker.py 

### Software and packages used by DPAM
- HH-suite3: https://github.com/soedinglab/hh-suite (enable addss.pl to add secondary structure)
- DaliLite.v5: http://ekhidna2.biocenter.helsinki.fi/dali/
- Python 3.7
- Foldseek 
- mkdssp
- pdbx: https://github.com/soedinglab/pdbx
- pdb2fasta (https://zhanggroup.org/pdb2fasta)
- tensorflow=1.14


### Supporting databases for DPAM:

The databases required for DPAM and all other supporting files can be download from our lab sever https://conglab.swmed.edu/DPAM/. 
  
It includes the following databases and files. </p>
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

The databases required for DPAM and all other supporting files can be download from our lab web https://conglab.swmed.edu/DPAM/

After downloading the databases, please decompress files. All databases and supporting files should be put in the same directory and the directory should be provided to `run_dpam_docker.py` 
    
## Installation
**docker pull conglab/dpam:latest**

## Usage
<pre>run_dpam_docker.py [-h] --databases_dir DATABASES_DIR --input_dir
                    INPUT_DIR --dataset DATASET
                    [--image_name IMAGE_NAME] [--threads THREADS]
                    [--log_file LOG_FILE]</pre>

### Arguments

- `-h`, `--help`  
  Show this help message and exit. Use this argument if you need information about different command options.

- `--databases_dir DATABASES_DIR`  
  **(Required)** Specify the path to the databases directory that needs to be mounted.

- `--input_dir INPUT_DIR`  
  **(Required)** Specify the path to the input directory that needs to be mounted.

- `--dataset DATASET`  
  **(Required)** Provide the name of the dataset for domain segmentation and classification.

- `--image_name IMAGE_NAME`  
  Specify the Docker image name. If not provided, a default image name will be used.

- `--threads THREADS`  
  Define the number of threads to be used. By default, the script is configured to utilize all available CPUs.

- `--log_file LOG_FILE`  
  Specify a file where the logs should be saved. If not provided, logs will be displayed in the standard output.

### Input organization

Before running the wrapper, the `INPUT_DIR` needs to be in the following structure:
    
    <INPUT_DIR>/
        dataset1/
        dataset1_struc.list
        dataset2/
        dataset2_struc.list
        ...


The `dataset1/` and `dataset2/` directories include PDB/mmCIF files and json file for PAE and `dataset1_struc.list` and `dataset2_struc.list` include targets (prefix of PDB/mmCIF and json), one line per one target.

`example/` is provided.

`example/` is considered `<INPUT_DIR>` and `test` under `example/` is considered as `<dataset>`

**exmaple command**:

`python ./run_dpam_docker.py --dataset test --input_dir example  --databases_dir databases --threads 32`

`databases` is the directory uncompressed fromd databases.tar.gz from our lab server. 

Final output should be <dataset>_domins under <INPUT_DIR>. For the example case, it should be test_domains under `example/`
