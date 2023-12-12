# DPAM
A **D**omain **P**arser for **A**lphafold **M**odels 

DPAM: A Domain Parser for AlphaFold Models （https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4548)

## Updates:
A docker image for DPAM v2.0 can be dowloaded by **docker pull conglab/dpam:latest** and previous version (v1.0) is moved to v1.0 directory (2023-12-10) . New version includes domain classification based on ECOD database and addresses over-segmentation for some proteins. 

Upload domain parser results for six model organisms.  (2022-12-6)

Replace Dali with Foldseek for initial hits searching. (2022-11-30)

Fix a bug in analyze_PDB.py which prevents the proper usage of Dali results. (2022-10-31)
## Prerequisites (required): 
Docker 

Python3

[Databases and supporting files](https://conglab.swmed.edu/DPAM/databases.tar.gz)

### Software and packages used by DPAM (already installed in the docker image)
- HH-suite3: https://github.com/soedinglab/hh-suite (enable addss.pl to add secondary structure)
- DaliLite.v5: http://ekhidna2.biocenter.helsinki.fi/dali/
- Python 3.7
- Foldseek 
- mkdssp
- pdbx: https://github.com/soedinglab/pdbx
- pdb2fasta (https://zhanggroup.org/pdb2fasta)
- tensorflow=1.14


### Supporting databases for DPAM:

The databases necessary for DPAM, along with all supporting files, are available for download from our lab server at [https://conglab.swmed.edu/DPAM/](https://conglab.swmed.edu/DPAM/). The compressed file size is approximately 89GB, while the size of the databases when uncompressed reaches around 400GB. It is essential to ensure that you have sufficient hard drive space to accommodate these databases. Additionally, due to their substantial size, downloading these databases might require several hours to a few days, depending on your internet connection speed.

After downloading the databases, please decompress files. All databases and supporting files should be put in the same directory and the directory should be provided to `run_dpam_docker.py` 
    
## Installation
    docker pull conglab/dpam:latest
    git clone https://github.com/CongLabCode/DPAM
    wget https://conglab.swmed.edu/DPAM/databases.tar.gz
    tar -xzf databases.tar.gz

### Quick test
`python run_dpam_docker.py --dataset test --input_dir example  --databases_dir databases --threads 32`

## Usage
<pre>python run_dpam_docker.py [-h] --databases_dir DATABASES_DIR --input_dir
                    INPUT_DIR --dataset DATASET
                    [--image_name IMAGE_NAME] [--threads THREADS]
                    [--log_file LOG_FILE]</pre>

### Arguments

- `-h`, `--help`  
  Show this help message and exit. Use this argument if you need information about different command options.

- `--databases_dir DATABASES_DIR`  
  **(Required)** Specify the path to the databases directory (downloaded before and uncompressed) that needs to be mounted to the docker. Please make sure you download the databases before running

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
        <dataset1>/
        <dataset1>_struc.list
        <dataset2>/
        <dataset2>_struc.list
        ...


The `<dataset1>/` and `<dataset2>/` directories include PDB/mmCIF files and json file for PAE and `dataset1_struc.list` and `dataset2_struc.list` include targets (prefix of PDB/mmCIF and json), one line per one target. <dataset> can be any name and postfix _struc.list has to be maintained. 

In the example test in **Quick test** above, 

`example/` is `<INPUT_DIR>` and `test` under `example/` is `<dataset>`

**exmaple command**:

`python run_dpam_docker.py --dataset test --input_dir example  --databases_dir databases --threads 32`

`databases` is the directory uncompressed fromd databases.tar.gz from our lab server. 

### Output
The pipeline will generate log files for each step for debugging. 

Final output is \<dataset\>_domains under <INPUT_DIR>. 

For the example, it should be `test_domains` under `example/`. 
