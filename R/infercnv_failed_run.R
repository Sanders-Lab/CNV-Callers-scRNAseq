library(infercnv)
library(tidyverse)
library(Seurat)
library(stringr)

load(".RData")

counts_atc2 <- read.table("./proc/infercnv_atc2_ref_merged_counts.tsv.gz")
sample_annot <- read.table("./proc/infercnv_atc2_sample_annot.tsv.gz")
counts_atc2[1:5, 1:5]


counts_uc <- read.table("./proc/infercnv_uc_kiwi_ref_merged_counts.tsv.gz")
counts_uc[1:5, 1:5]
counts_uc %>% tail()


## ---- using the fresh tnbc1 merged dataset ---- ##

counts_tnbc1 <- read.table("./proc/pbmc3k_tnbc1_ref_merged_counts.tsv", row.names = 1)
counts_tnbc1 <- read_tsv("./proc/pbmc3k_tnbc1_ref_merged_counts.tsv")


counts_tnbc1[1:5, 1:2]

counts_df <- as.data.frame(counts_tnbc1)
row.names(counts_df) <- counts_df$...1

counts_df %>% head()
counts_df$...1 <- NULL
counts_df %>% ncols()
counts_tnbc1[, 1] %>%
    is.nan() %>%
    table()




rownames(counts_tnbc1)
counts_tnbc1[, 1]
row.names(counts_tnbc1) <- counts_tnbc1[, 1]
counts_mat <- as.matrix(counts_tnbc1)
counts_mat[1:5, 1:5]
counts_mat %>% ncol()
counts_tnbc1[1:5, 1:5]
counts_tnbc1[1:5, 3790:3797]
counts_mat %>%
    colnames() %>%
    head()


colnames(counts_tnbc1) <- colnames(counts_tnbc1) %>% str_replace_all("[.]", "-")
counts_tnbc1[1:5, 1:5]

counts_mat <- as.matrix(counts_tnbc1)
counts_mat[1:5, 1:5]
nrow(counts_mat)
ncol(counts_mat)



## -- tnbc1 -- ##
filt_counts <- read.table("./proc/filtered_probGenesRemoved_pbmc3k_tnbc1_counts.tsv.gz")
filt_counts %>%
    rownames() %>%
    head()
filt_counts[1:5, 1:5]
colnames(filt_counts) <- colnames(filt_counts) %>% str_replace_all("[.]", "-")
counts_mat <- as.matrix(filt_counts)


sample_annot <- read.table("./proc/pbmc3k_tnbc1_sample_annot.tsv.gz")
sample_annot %>% head()
sample_annot %>% tail()
sample_annot %>% nrow()
sample_annot[250:260, ]


sample_annot$V1 <- sample_annot$V1 %>% str_replace_all("[.]", "_")
sample_annot %>% head()
sample_annot %>% tail()
row.names(sample_annot) <- sample_annot$V1
sample_annot %>% head()
all(rownames(sample_annot) %in% colnames(counts_mat)) %>% table()
colnames(counts_df)[1:5]

rownames(counts_mat)
"AC136612.1" %in% gene_order$V1


gene_order <- read.table("./infercnv_test/input_files/hg38_gencode_v27.txt")
gene_order <- read.table("./proc/hgnc_infercnv_gene_ranges.tsv.gz")
gene_order %>% head()




infercnvObj <- CreateInfercnvObject(
    raw_counts_matrix = counts_mat,
    annotations_file = sample_annot,
    gene_order_file = gene_order,
    ref_group_names = "normal"
)


# trying with null ref names
infercnvObj <- CreateInfercnvObject(
    raw_counts_matrix = counts_tnbc1,
    annotations_file = sample_annot,
    gene_order_file = gene_order,
    ref_group_names = NULL
)


out_dir <- "./outputs/infercnv/infercnv_cutoff_0-1"

# cutoff is 0.1 more suitable for 10X data
# 1 works well for smart-seq2
infercnv_obj <- infercnv::run(infercnvObj,
    cutoff = 0.1,
    out_dir = out_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE
)


out_dir_c1 <- "./outputs/infercnv/infercnv_cutoff_1"


# This did not work
infercnv_obj_c1 <- infercnv::run(infercnvObj,
    cutoff = 1,
    out_dir = out_dir_c1,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE
)
