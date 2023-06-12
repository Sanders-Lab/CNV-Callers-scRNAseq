library(tidyverse)
library(Seurat)

# atlas <- readRDS('/fast/groups/ag_sanders/work/data/scIBD/BloodAtlas/Granja/Seurat_Granja_scRNAseq.rds')
# atlas %>% head()

# atlas@meta.data$Group %>% unique()

# atlas <- subset(atlas, subset = Group == c('CD34_D2T1','CD34_D3T1','PBMC_D4T1','PBMC_D4T2'))

# counts <- GetAssayData(atlas, slot = 'counts')
# counts %>% head()

# write.table(counts, 'counts_atlas_pbmc_kiwi.tsv', sep = '\t')

# atlas@meta.data$cell_names <- rownames(atlas@meta.data)
# atlas@meta.data$cell_names %>% head()
# annot <- atlas@meta.data %>% select(cell_names, BioClassification)

# annot %>% head()

# annot <- annot %>% rename(cell = cell_names)
# annot <- annot %>% rename(group = BioClassification)

# annot %>% head()

# write.table(annot, 'cell_annot_kiwi.tsv', row.names = F, col.names = T, sep = '\t')




## ---- Read in matrix files and cell annotations ---- ##
pbmc3k <- Read10X("../data/reference/pbmc_3k/")
pbmc_cell_annot <- read.table("../proc/pbmc3k_cell_group_annot.tsv.gz",
    header = T
)
pbmc3k[1:5, 1:5]
row.names(pbmc3k) <- rownames(pbmc3k) %>% str_replace(., "_", "-")
pbmc3k[1:5, 1:5]
pbmc_cell_annot %>% head()
pbmc_cell_annot %>% nrow()
pbmc3k %>% ncol()
pbmc3k_mat <- as.matrix(pbmc3k)
pbmc3k_mat %>% nrow()
pbmc3k_mat %>% ncol()


# add -1 to the cell names in pbmc3k
# so that the cell names stay the same
pbmc_cell_annot$cell <- paste0("pbmc3k_", pbmc_cell_annot$cell, "-1")
pbmc_cell_annot$cell <- paste0(pbmc_cell_annot$cell, "-1")
pbmc_cell_annot %>% head()
pbmc_cell_annot %>% tail()
pbmc_cell_annot %>% nrow()


# read in epithelium data
epithelium <- Read10X("../data/reference/epithelium/")
epithelium[1:5, 1:5]




# making seurat objs
pbmc_so <- CreateSeuratObject(counts = pbmc3k)
pbmc_so
epi_so <- CreateSeuratObject(counts = epithelium)
epi_so <- CreateSeuratObject(counts = epithelium)
epi_so
epithelium[1:5, 1:5]
epithelium %>% nrow()
epithelium %>% rownames()


# merging seurat objs
merged_so <- merge(pbmc_so, epi_so, add.cell.ids = c("pbmc3k", "epi"))
merged_so <- merge(pbmc_so, epi_so)
merged_so



# get the merged counts
merged_counts <- GetAssayData(merged_so, slot = "counts")
merged_counts %>% head()
merged_counts_mat <- as.matrix(merged_counts)
merged_counts_mat[1:5, 1:5]
merged_counts_mat[3325:3331, 3325:3331]
merged_counts_mat %>% ncol()
merged_counts_mat %>%
    row.names() %>%
    head()



# saving to file
write.table(merged_counts_mat,
    gzfile("../proc/pbmc3k_epithelium_ref_counts.tsv.gz"),
    sep = "\t"
)




# getting the list of cell annot
merged_so@meta.data %>% head()
epi_cells <- epithelium %>% colnames()
epi_cells <- paste0("epi_", epi_cells)
epi_cell_annot <- as.data.frame(epi_cells)
epi_cell_annot$group <- "epithelium"
epi_cell_annot %>% head()
epi_cell_annot <- epi_cell_annot %>% rename(cell = epi_cells)




combined_cell_annot <- rbind(pbmc_cell_annot, epi_cell_annot)
combined_cell_annot %>% nrow()
combined_cell_annot %>% head()
combined_cell_annot %>% tail()
combined_cell_annot %>% ncol()


## ---- Numbat ref gen ---- ##

library(numbat)
library(Matrix)
# counts <- read.table('counts_atlas_pbmc_kiwi.tsv', sep = '\t')
# counts %>% head()



numbat_ref <- aggregate_counts(epithelium, epi_cell_annot)
numbat_ref <- aggregate_counts(pbmc3k, pbmc_cell_annot)
numbat_ref <- aggregate_counts(merged_counts_mat, combined_cell_annot)

numbat_ref %>% head()

write.table(numbat_ref,
    gzfile("../proc/numbat_ref_pbmc_epi.tsv.gz"),
    sep = "\t"
)
# write.table(numbat_ref, 'numbat_ref_granje_kiwi_pbmc.tsv', sep = '\t')

ref <- read.table("numbat_ref_granje_kiwi_pbmc.tsv")
ref %>% head()
ref %>% ncol()
