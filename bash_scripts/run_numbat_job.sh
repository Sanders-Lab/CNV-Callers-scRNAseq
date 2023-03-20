#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=2G
#SBATCH --time=10:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=suharto.banerjee@mdc.berlin.de


source ~/.bashrc
conda activate numbat
printf "\nConda env activated at: %s" $(echo $CONDA_DEFAULT_ENV)


Rscript ../R/run_numbat.R \
    "" \
    "" \
    ../data/ega_data/tall_ega_data/10X_count/bamtofastq/outs/filtered_feature_bc_matrix \
    ../data/ega_data/tall_ega_data/pileup_phase_tall/tall_allele_counts.tsv.gz \
    $SLURM_CPUS_PER_TASK \
    ../outputs/tall_scnova/numbat_pbmc_tall/
