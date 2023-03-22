library(tidyverse)
library(infercnv)
library(Seurat)


## ---- AIM ---- ##
# this script runs infercnv
# it also generates sample annot file
# required for the run here
# doing it here maintains consistent cell ids
# 
# @param 1. count_mat_path: combined count mat
# @param 2. mat_type: the count mat type. 
#           for 10x, pass t (Default)
# @param 3: n_cores
# @param 4: out_dir
## ---- AIM ---- ##



# read in the combined matrix.
# gene names need to be row.names for both count mat
# and the gene loc mat 
# convert the df to matrix? check if works without
# gene_mat <- read.table("../proc/pbmc3k_tnbc1_ref_merged_counts.tsv.gz",
#                        row.names = 1)

ReadData <- function(count_mat_path,
                     gene_order_path = "../data/gene_chr_loc/hg38_gencode_v27.txt",
                     mat_type = "c") {

    if(mat_type == "t") {
    
        count_mat <- Read10X(count_mat_path)
    } 
    if(mat_type == "c") {

        count_mat <- read.table(count_mat_path)
    }

    gene_order <- read.table(gene_order_path,
                             row.names = 1)

    writeLines("Read in the Data!")
    return(list(count_mat, gene_order))
}




ProcAnnot <- function(count_mat) {

    # making the sample annot file
    # check if the cell prefix has pbmc3k
    # mark those rows as normal
    # otherwise malignant_{sample}
    sample_annot <- data.frame(cells = colnames(count_mat))
    row.names(sample_annot) <- sample_annot$cells

    sample_annot$groups <- sample_annot$cells %>%
        sapply(.,
               function(x) ifelse(str_detect(x, "pbmc"), "pbmc", "malignant_tall" ))
    sample_annot$cells <- NULL
    #     names(sample_annot) <- NULL

    writeLines("Sample annot processed successfully!")
    return(as.data.frame(sample_annot))
}


RunInfercnv <- function(count_mat,
                        sample_annot,
                        gene_order,
                        ref_group_names,
                        n_cores,
                        out_dir) {

    # making the infercnv object
    # pass the ref group names
    # if we have identified cell groups, pass those as a chr vec
    infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = count_mat,
                                         annotations_file = sample_annot,
                                         gene_order_file = gene_order,
                                         ref_group_names = ref_groups)




    # running infercnv
    # FROM VIGNETTE:
    # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    infercnv_obj <- infercnv::run(infercnv_obj,
                                  cutoff = 0.1, 
                                  out_dir = out_dir,
                                  cluster_by_groups = TRUE,
                                  denoise = TRUE,
                                  HMM = TRUE,
                                  num_threads = n_cores)

}




#### ---- MAIN SECTION ---- ####

cmd_args <- commandArgs(trailingOnly = T)

count_mat_path <- cmd_args[1]
mat_type <- cmd_args[2]
n_cores <- as.integer(cmd_args[3])
out_dir <- cmd_args[4]
ref_groups <- "pbmc"

writeLines(str_glue("The input args are:\n
                    count_mat_path: {count_mat_path}\n
                    mat_type: {mat_type}\n
                    n_cores: {n_cores}\n
                    out_dir: {out_dir}\n"))

# count_mat_path <- "../proc/pbmc3k_tall_ref_merged_counts.tsv.gz"
ret_list <- ReadData(count_mat_path)
count_mat <- data.frame(ret_list[[1]])
gene_order <- data.frame(ret_list[[2]])


sample_annot <- as.data.frame(ProcAnnot(count_mat))


RunInfercnv(count_mat,
            sample_annot,
            gene_order,
            ref_group_names,
            n_cores,
            out_dir)
