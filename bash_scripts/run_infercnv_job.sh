#!/bin/bash
#
#SBATCH -o R_%10x_%j.out
#SBATCH -e R_%10x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=2G
#SBATCH --time=10:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=suharto.banerjee@mdc.berlin.de


source ~/.bashrc
conda activate infercnv
printf "\nConda env activated at: %s" $(echo $CONDA_DEFAULT_ENV)

Rscript ../R/run_infercnv.R \
    ../proc/pbmc3k_tall_ref_merged_counts.tsv.gz \
    "" \
    $SLURM_CPUS_PER_TASK \
    ../outputs/tall_scnova/infercnv_pbmc3k_tall

