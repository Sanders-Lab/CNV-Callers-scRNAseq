library(tidyverse)
library(Seurat)
library(numbat)




#### ---- BRIEF ---- ####
# @param 1. ref_mat_path
# @param 2. ref_annot_path
# @param 3. count_mat_path
# @param 4. allele_count_path
# @param 5. n_cores
# @param 6. out_dir





#### ---- Read in matrix files and cell annotations ---- ####

ReadData <- function(ref_mat_path = "../data/reference/pbmc_3k/",
                     ref_annot_path = "../proc/pbmc3k_cell_group_annot.tsv.gz",
                     count_mat_path,
                     allele_count_path,
                     count_mat_type = "t") {
    pbmc3k <- Read10X(ref_mat_path)
    pbmc_cell_annot <- read.table(ref_annot_path,
        header = T
    )

    if (count_mat_type == "t") {
        count_mat <- Read10X(count_mat_path)
    }

    allele_count <- read.table(allele_count_path,
        sep = "\t",
        header = T
    )

    writeLines("\nRead in all the files required!\n")
    return(list(pbmc3k, pbmc_cell_annot, count_mat, allele_count))
}




#### ---- Make the reference for numbat ---- ####

MakeNumbatRef <- function(ref_mat, ref_mat_annot) {
    numbat_ref <- aggregate_counts(ref_mat, ref_mat_annot)

    writeLines("\nGenerated reference for numbat!\n")
    return(numbat_ref)
}


SaveNumbatRef <- function(numbat_ref, ref_save_path) {
    writeLines("\nSaving the numbat ref to file")
    write.table(numbat_ref,
        gzfile(save_path),
        sep = "\t"
    )
}




#### ---- run numbat ---- ####

RunNumbat <- function(count_mat,
                      ref_custom,
                      allele_count,
                      n_cores,
                      out_dir) {
    run_numbat(count_mat,
        ref_custom,
        allele_count,
        genome = "hg38",
        ncores = n_cores,
        plot = TRUE,
        out_dir = out_dir
    )
}




#### ---- MAIN SECTION ---- ####

cmd_args <- commandArgs(trailingOnly = T)

writeLines(str_glue("\nThe passed in args are:\n
                    ref_mat_path: {cmd_args[1]}\n
                    ref_annot_path: {cmd_args[2]}\n
                    count_mat_path: {cmd_args[3]}\n
                    allele_count_path: {cmd_args[4]}\n
                    n_cores: {cmd_args[5]}\n
                    out_dir: {cmd_args[6]}\n"))



# use the default pbmc loc if blank is provided
if (!cmd_args[1] == "") {
    read_list <- ReadData(
        ref_mat_path = cmd_args[1],
        ref_annot_path = cmd_args[2],
        count_mat_path = cmd_args[3],
        allele_count_path = cmd_args[4],
        count_mat_type = "t"
    )
} else {
    read_list <- ReadData(
        count_mat_path = cmd_args[3],
        allele_count_path = cmd_args[4],
        count_mat_type = "t"
    )
}

n_cores <- as.integer(cmd_args[5])
out_dir <- cmd_args[6]

ref_mat <- read_list[[1]]
ref_mat_annot <- read_list[[2]]
count_mat <- read_list[[3]]
allele_count <- read_list[[4]]


numbat_ref <- MakeNumbatRef(ref_mat, ref_mat_annot)
# SaveNumbatRef(numbat_ref, )


RunNumbat(
    count_mat,
    numbat_ref,
    allele_count,
    n_cores,
    out_dir
)

# QND TESTS
#   32 outputs/tall_scnova/numbat_pbmc_tall/
#
# read_list <- ReadData(count_mat_path = "../data/ega_data/tall_ega_data/10X_count/bamtofastq/outs/filtered_feature_bc_matrix",
#                       allele_count_path = "../data/ega_data/tall_ega_data/pileup_phase_tall/tall_allele_counts.tsv.gz",
#                       count_mat_type = "t")
# ref_mat <- read_list[[1]]
# ref_mat_annot <- read_list[[2]]
# count_mat <- read_list[[3]]
# allele_count <- read_list[[4]]
#
# ref_mat %>% head(n = c(5,5))
# ref_mat_annot %>% head()
#
# numbat_ref <- MakeNumbatRef(ref_mat, ref_mat_annot)
