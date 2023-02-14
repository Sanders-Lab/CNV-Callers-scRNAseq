#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=10:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=suharto.banerjee@mdc.berlin.de

source ~/.bashrc
conda activate numbat
printf "\nME: Conda env activated: %s!\n" $(echo $CONDA_DEFAULT_ENV)


OUT_DIR=pileup_phase_output_tnbc1
mkdir $OUT_DIR
printf "\nME:Output dir created at ./%s" $(echo $OUT_DIR)


printf "\nME: Starting Numbat Pileup and Phasing...\n"

DATA_DIR=./tnbc1/outs
SAMPLE=tnbc1
#SAMPLE=$(ls $BAM_dir/*.bam | awk -F/ '{print $2}' | awk -F_ '{print $1}')

printf "\nME: Starting with SAMPLE: %s\n" "$SAMPLE"

Rscript /fast/work/users/sbanerj_m/miniconda3/envs/numbat/lib/R/library/numbat/bin/pileup_and_phase.R \
    --label $SAMPLE \
    --samples $SAMPLE \
    --bams $DATA_DIR/possorted_genome_bam.bam \
    --barcodes $DATA_DIR/filtered_feature_bc_matrix/barcodes.tsv.gz \
    --outdir $OUT_DIR/ \
    --gmap $g_tools/eagle/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf $g_data/1000g_snp_vcf/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
    --paneldir $g_data/1000g_snp_vcf/1000G_hg38_bcf/ \
    --ncores $SLURM_CPUS_PER_TASK

printf "\nME: Completed pileup and phasing for SAMPLE: %s\n" "\n$SAMPLE"


printf "\nME: Numbat Pileup and phase finished! Please check output folder and logs for errors!"
