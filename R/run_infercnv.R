library(tidyverse)
library(infercnv)
library(stringr)


## ---- AIM ---- ##
# this script runs infercnv
# it also generates sample annot file
# required for the run here
# doing it here maintains consistent cell ids
## ---- AIM ---- ##



# read in the combined matrix.
# gene names need to be row.names for both count mat
# and the gene loc mat 
# convert the df to matrix? check if works without
# gene_mat <- read.table('../proc/pbmc3k_tnbc1_ref_merged_counts.tsv.gz',
#                        row.names = 1)


## ---- tnbc1 + pbmc3k + epi ---- ##
gene_mat <- read.table('../proc/pbmc3k_epithelium_tnbc1_ref_merged_counts.tsv.gz')
gene_order <- read.table('../infercnv_test/input_files/hg38_gencode_v27.txt',
                         row.names = 1)
# gene_order <- read.table('./proc/hgnc_infercnv_gene_ranges.tsv.gz',
#                          row.names = 1)
gene_order[1:5,1:3]
gene_mat[1:5,1:5]
# colnames(gene_mat) <- colnames(gene_mat) %>% str_replace('[.]', '_')
# colnames(gene_mat) %>% head
count_mat <- as.matrix(gene_mat)
count_mat[1:5,1:5]
count_mat %>% ncol()
count_mat[1:5,670:675]



# making the sample annot file
# check if the cell prefix has pbmc3k
# mark those rows as normal
# otherwise malignant_{sample}
cells <- colnames(gene_mat)
cells %>% length
sample_annot <- as.data.frame(cells)
sample_annot %>% head
row.names(sample_annot) <- sample_annot$cells
sample_annot$cells %>% str_detect('SRR') %>% table()

for (i in 1:nrow(sample_annot)){
    cell <- sample_annot$cell[i]

    if(str_detect(cell,'pbmc3k')){
        sample_annot$group[i] <- 'pbmc'
    }

    if(str_detect(cell, 'SRR')){
        sample_annot$group[i] <- 'epithelium'
    }

    else{
        sample_annot$group[i] <- 'malignant_tnbc1'
    }
}

sample_annot$group <- sample_annot %>% apply(., 
                                             1,
                                             function(x) ifelse(str_detect(x, 'pbmc3k'), 'pbmc', 'malignant_tnbc1'))
sample_annot$group <- sample_annot$cells %>% lapply(., 
                                             function(x) ifelse(str_detect(x, 'SRR'), 'epithelium', ()))


# sample_annot$cells  <-  sample_annot$cells %>% str_replace('[.]', '_')
# row.names(sample_annot) <- sample_annot$cells
sample_annot$cells = NULL
sample_annot %>% head
sample_annot[625:635,] 
sample_annot %>% tail
sample_annot %>% nrow
sample_annot$group %>% str_detect('pbmc') %>% table




# making the infercnv object
# pass the ref group names
# if we have identified cell groups, pass those as a chr vec
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = count_mat,
                                     annotations_file = sample_annot,
                                     gene_order_file = gene_order,
                                     ref_group_names = c('pbmc', 'epithelium'))


# set the output dir for the results
out_dir <- '../outputs/infercnv_pbmc3k_epithelium_tnbc1/'




# running infercnv
# FROM VIGNETTE:
# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
infercnv_obj <- infercnv::run(infercnv_obj,
                             cutoff = 0.1, 
                             out_dir = out_dir,
                             cluster_by_groups = TRUE,
                             denoise = TRUE,
                             HMM = TRUE,
                             num_threads = 32)


# check this obj properly
after_run <- readRDS('../outputs/infercnv_pbmc3k_epithelium_tnbc1/run.final.infercnv_obj')
after_run@reference_grouped_cell_indices
after_run@observation



