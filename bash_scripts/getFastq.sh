#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=suharto.banerjee@mdc.berlin.de

# @param 1. outdir
# @param 2. SRR ID

fastq-dump --gzip --split-files --readids --outdir $1 $2
