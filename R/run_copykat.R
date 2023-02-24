library(copykat)
library(tidyverse)


## ---- AIM ---- ##
# this script runs copykat cnv caller
# on the combined ref + exp scRNA-seq data
# combined because we want to use the ref baseline
# to call out cnvs
## ---- AIM ---- ##



## ---- sample data ---- ##
# ck_data <- exp.rawdata 
# ck_data %>% head()



## ---- own data ---- ##

# ref <- read.table('./proc/counts_atlas_pbmc_kiwi.tsv.gz')
# ref[1:5,1:5]


# ref_cells <- colnames(ref)
# ref_cells[1:5]
# length(ref_cells)

# count_mat <- read.table('./proc/infercnv_atc2_ref_merged_counts.tsv')
# count_matrix <- as.matrix(count_mat)
# annot_mat <- read.table('./infercnv_test/input_files/hg38_gencode_v27.txt')
# annot_mat %>% head()


## ---- tnbc1 fresh ---- ##

# reading in the tnbc1 pbmc combined mat
# making a matrix and supplying it as it is good practice
count_tnbc1 <- read.table('../proc/pbmc3k_epithelium_tnbc1_ref_merged_counts.tsv.gz')
count_tnbc1[1:5,1:5]
# count_tnbc1[1:5,6835:6842]
# count_tnbc1[1:5, 5740:5746]
count_mat <- as.matrix(count_tnbc1)


# get the cell ids of the normal cells
# 1. select the cols which have pbmc3k || SRR since it is cell prefix
# 2. get the colnames of this subsetted df
# SRR is for the cells from epithelium dataset
ref_mat_from_count <- count_tnbc1 %>% select(contains('pbmc3k') | contains('SRR'))
ref_mat_from_count %>% ncol()
ref_cells <- ref_mat_from_count %>% colnames()
ref_cells[1:5]
ref_cells %>% length()
# ref_cells[1:2] %in% colnames(count_tnbc1[,1:2])
# colnames(count_tnbc1[,1:2])


# dir where results will be saved
out_dir <- '../outputs/copykat_pbmc3k_epithelium_tnbc1/'


# setting the wd to the output dir
# so that the results are saved to the out_dir
setwd(out_dir)


# win.size is the window for smoothing the exp
# 25 is the recommended value in the vignette
# pasing the names of cells from the ref as the norm.cell.names
copykat_tnbc <- copykat(rawmat=count_mat,
                        id.type="S",
                        ngene.chr=5,
                        win.size=101,
                        KS.cut=0.1,
                        sam.name="tnbc1",
                        distance="euclidean",
                        norm.cell.names=ref_cells,
                        output.seg="FLASE",
                        plot.genes="TRUE",
                        genome="hg20",
                        n.cores=32)

