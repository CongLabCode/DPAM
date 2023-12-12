# DPAM
A **D**omain **P**arser for **A**lphafold **M**odels 

DPAM: A Domain Parser for AlphaFold Models （https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4548)

## Updates:
A docker image for DPAM v2.0 can be dowloaded by **docker pull conglab/dpam:latest** and previous version (v1.0) is moved to v1.0 directory (2023-12-10) . New version includes domain classification based on ECOD database and addresses over-segmentation for some proteins. Warning: current Docker image only works on AMD x86, not Apple M series chip. We're updating it for the compatibility. Stay tuned!
Upload domain parser results for six model organisms.  (2022-12-6)

Replace Dali with Foldseek for initial hits searching. (2022-11-30)

Fix a bug in analyze_PDB.py which prevents the proper usage of Dali results. (2022-10-31)
## Prerequisites (required): 
Docker/Singularity

Python3

[Databases and supporting files](https://conglab.swmed.edu/DPAM/databases.tar.gz)

### Supporting databases for DPAM:

The databases necessary for DPAM, along with all supporting files, are available for download from our lab server at [https://conglab.swmed.edu/DPAM/](https://conglab.swmed.edu/DPAM/). The compressed file is around 89GB, expanding to about **400GB** when uncompressed (best run on a computing cluster/workstation due to the substantial storage needs, which may surpass the capacity of typical personal computers). It is essential to ensure that you have sufficient hard drive space to accommodate these databases. Additionally, due to their substantial size, downloading these databases might require several hours to a few days, depending on your internet connection speed.

After downloading the databases.tar.gz, please decompress the file. And the directory(`[download_path]/databases`) must be provided to `run_dpam_docker.py` as `--databases_dir`
    
## Installation
For Docker:
    
    docker pull conglab/dpam:latest
    git clone https://github.com/CongLabCode/DPAM
    cd ./DPAM
    wget https://conglab.swmed.edu/DPAM/databases.tar.gz
    tar -xzf databases.tar.gz

For Singularity:

    git clone https://github.com/CongLabCode/DPAM
    cd ./DPAM
    wget https://conglab.swmed.edu/DPAM/databases.tar.gz
    tar -xzf databases.tar.gz
    singularity pull dpam.sif docker://conglab/dpam
    
    
    

### Quick test
For Docker:

    python run_dpam_docker.py --dataset test --input_dir example  --databases_dir databases --threads 32

For Singularity:

    python ./run_dpam_singularity.py --databases_dir databases --input_dir example --dataset test --threads 32 --image_name dpam.sif`

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
