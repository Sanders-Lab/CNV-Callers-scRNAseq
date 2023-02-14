#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=suharto.banerjee@mdc.berlin.de

FASTQ_DIR=tnbc_fastq

$g_tools/10x-BAM-to-Fastq/bamtofastq_linux --nthreads 32 BAM_TNBC1.bam $FASTQ_DIR 
