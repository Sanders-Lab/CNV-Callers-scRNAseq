#!/bin/bash

# this can be automated for every env
# exporting each of the env serially
# will parallelize later
# copykat is in infercnv env only.
# since both rely on jags and r-jags

source ~/.bashrc
DATE=$(date -I)


conda activate numbat
printf "\nExporting currently active env: %s" $(echo $CONDA_DEFAULT_ENV)
conda env export > ../environments/${CONDA_DEFAULT_ENV}_${DATE}_env.yml
printf "\n%s env exported successfully to the environments folder!" $(echo $CONDA_DEFAULT_ENV)


conda activate infercnv
printf "\nExporting currently active env: %s" $(echo $CONDA_DEFAULT_ENV)
conda env export > ../environments/${CONDA_DEFAULT_ENV}_${DATE}_env.yml
printf "\n%s env exported successfully to the environments folder!" $(echo $CONDA_DEFAULT_ENV)
