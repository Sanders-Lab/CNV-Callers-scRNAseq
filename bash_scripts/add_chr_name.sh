#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=suharto.banerjee@mdc.berlin.de

source ~/.bashrc
conda activate numbat
printf "Conda env activated: %s" $(echo $CONDA_DEFAULT_PREFIX)

# script to add chr before chromosome names in bam file
# and save the file to a temp.bam file
# first regex greps for any line starting with @SQ
# substitute mode is then switched on
# changes SN: to SN:chr
# SECOND regex similar
samtools view -h $1 | \
    sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | \
    samtools view -bS - > temp.bam

printf 'job complete!!'
