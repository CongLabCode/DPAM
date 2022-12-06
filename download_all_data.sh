#!/bin/bash
# Usage: bash download_all_data.sh /path/to/download/directory
set -e

if [[ $# -eq 0 ]]; then
  echo "Error: download directory must be provided as an input argument."
  exit 1
fi

if ! command -v aria2c &> /dev/null ; then
  echo "Error: aria2c could not be found. Please install aria2c."
  exit 1
fi


DOWNLOAD_DIR="$1"
### Download ECOD70
echo "Downloading ECOD70..."
SOURCE_URL="https://conglab.swmed.edu/DPAM/ECOD70.tgz"
wget --no-check-certificate "${SOURCE_URL}"  -P "${DOWNLOAD_DIR}"

### Download pdb70
echo "Downloading pdb70..."
SOURCE_URL="https://conglab.swmed.edu/DPAM/pdb70.tgz"
wget --no-check-certificate "${SOURCE_URL}"  -P "${DOWNLOAD_DIR}"

### Download UniRef30
echo "Downloading UniRef30..."
SOURCE_URL="https://wwwuser.gwdg.de/~compbiol/uniclust/2022_02/UniRef30_2022_02_hhsuite.tar.gz"
wget --no-check-certificate "${SOURCE_URL}"  -P "${DOWNLOAD_DIR}"

### Download ECOD70 foldseek database
echo "Downloading ECOD70 foldseek database"
SOURCE_URL="https://conglab.swmed.edu/DPAM/ECOD_foldseek_DB.tgz"
wget --no-check-certificate "${SOURCE_URL}"  -P "${DOWNLOAD_DIR}"

### Download ECOD position weights
echo "Downloading ECOD position weights"
SOURCE_URL="https://conglab.swmed.edu/DPAM/ecod_weights.tgz"
wget --no-check-certificate "${SOURCE_URL}"  -P "${DOWNLOAD_DIR}"

### Download ECOD domain information
echo "Downloading ECOD domain information"
SOURCE_URL="https://conglab.swmed.edu/DPAM/ecod_domain_info.tgz"
wget --no-check-certificate "${SOURCE_URL}"  -P "${DOWNLOAD_DIR}"


### Download ECOD domain list, length, relationship to pdb and normalization
echo "Downloading other ECOD related data"
files=("ECOD_norms", "ecod.latest.domains", "ECOD_length", "ECOD_pdbmap")
for str in ${files[@]} 
do
    SOURCE_URL="https://conglab.swmed.edu/DPAM/${str}"
    wget --no-check-certificate "${SOURCE_URL}"  -P "${DOWNLOAD_DIR}"
done

### Download benchmark data
echo "Download benchmark data"
SOURCE_URL="https://conglab.swmed.edu/DPAM/ECOD_benchmark.tgz"
wget --no-check-certificate "${SOURCE_URL}"  -P "${DOWNLOAD_DIR}"

echo "Download done"
