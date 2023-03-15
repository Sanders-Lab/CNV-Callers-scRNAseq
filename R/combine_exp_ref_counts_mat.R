library(tidyverse)
library(Seurat)


## ---- Combine exp dataset with the ref dataet ---- ##
# RATIONALE: the combined counts mat will be used for
# cnv calling. The normal cells will be pointed out
# for ref.


# read in ref data
# counts_ref <- read.table('../proc/counts_atlas_pbmc_kiwi.tsv.gz')
# counts_ref %>% head()


## ---- PBMC 3K 10X ---- #
pbmc_3k <- Read10X('../data/reference/pbmc_3k')
pbmc_3k


# read in uc kiwi data
# uc <- Read10X('./data/P1835/sample_filtered_feature_bc_matrix')
# uc
# uc_mat <- as.matrix(uc)


# to change the hyphen to . to match the cell names in annot
# colnames(uc_mat) <- uc_mat %>% colnames() %>% str_replace(.,'-','.')
# uc_mat[1:5,1:5]


# reading in the epithelium dataset
epithelium <- Read10X('../data/reference/epithelium')
epithelium <- read.table('../proc/epithelium_hgnc_counts.tsv.gz',
                         row.names = 1)
epithelium[1:5,1:5]
# epithelium[,1] <- rownames(epithelium)
epithelium %>% rownames() %>% head()
epithelium %>% rownames() %>% is.na() %>% table()
epithelium %>% rownames() %>% str_detect('%') %>% table()
epi_mat <- as.matrix(epithelium)
epi_mat[1:5,1:5]
rownames(epi_mat)
rownames(epi_mat) %>% is.na() %>% table()
rownames(epi_mat) %>% str_detect('') %>% table()
epi_mat %>% nrow()
any(is.na(row.names(epi_mat)))


# One gene had a blank gene name which was 
# the culprit
# commands to sort out this shit
spurious_genes <- which(rownames(epithelium) == '')
spurious_genes
epithelium <- epithelium[-spurious_genes,]


## ---- tnbc1 fresh ---- ##
# SOURCE: dataset from copykat publication
# read in the saved file
tnbc1 <- read.table('../proc/tnbc1_gene_mat_raw_cleaned.tsv.gz')
tnbc1 %>% head(n = c(5,5))


# create seurat object for each of the counts
pbmc3k_so <- CreateSeuratObject(pbmc_3k, project = 'pbmc3k')
pbmc3k_so
epi_so <- CreateSeuratObject(epithelium, project = 'epi')
epi_so
tnbc1_so <- CreateSeuratObject(tnbc1, project = 'tnbc1')
tnbc1_so
# uc_so <- CreateSeuratObject(counts = uc_mat, project = 'uc')
# uc_so
epi_so@meta.data %>% rownames() %>% duplicated() %>% table()
colnames(epi_so)

all(colnames(pbmc3k_so) == colnames(epi_so))


# uc_so@assays$RNA@counts[1:5,1:5]


# merge the seurat objects
# cell ids are added as prefix to cell names
# so that it can be used to mark normal
and malignant cells easily when sample_annot is made
# had to merge in two steps because there was this error
# `please provide a cell identifier for each object provided to merge`
merged_so <- merge(pbmc3k_so,
                   tnbc1_so,
                   add.cell.ids = c('pbmc3k', 'tnbc1'),
                   project = 'combined')
merged_so
merged_so <- merge(epi_so,
                   merged_so)
merged_so
# merged_so <- merge(ref_so, uc_so, add.cell.ids = c('ref', 'uc'), project = 'combined')


# get the raw gene count matrix from the merged seurat object
merged_counts <- GetAssayData(merged_so, slot = 'counts')
merged_counts[1:5,1:5]
merged_counts[1:5,4420:4428]


# convert to matrix to save to file
merged_counts_mat <- as.matrix(merged_counts)
merged_counts_mat %>% head(n = c(5,5))
write.table(merged_counts_mat,
            gzfile('../proc/pbmc3k_epithelium_tnbc1_ref_merged_counts.tsv.gz'),
            sep = '\t')


# sanity file write check
# merged_mat <- read.table('./proc/infercnv_tnbc1_kiwi_ref_merged_counts.tsv.gz')
merged_mat %>% head()


# MAKE THE SAMPLE ANNOT FILE WITH RUN_INFERCNV
# make the sample annot file
# get the cell names
cells <- colnames(merged_so)
# cells <- colnames(merged_mat)
cells[1:5]
cells %>% tail()


# making a df
sample_annot <- data.frame(cells)
sample_annot %>% head()
sample_annot %>% tail()


# marking cells as normal if ref is prefix
# marking cells as malignant_tnbc1 if not
# sample_annot$group <- apply(sample_annot, 1, function(x) ifelse(str_split_i(x, '_', 1) == 'pbmc3k', 'normal', 'malignant_tnbc1'))
# sample_annot %>% head()
# sample_annot %>% tail()


# writing to file
# write.table(sample_annot, gzfile('./proc/pbmc3k_tnbc1_sample_annot.tsv.gz'), sep = '\t', row.names = F, col.names = F)
