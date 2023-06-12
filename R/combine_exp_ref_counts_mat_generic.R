library(tidyverse)
library(Seurat)


#### ---- Combine exp dataset with the ref dataet ---- ####
# RATIONALE: the combined counts mat will be used for
# cnv calling. The normal cells will be pointed out
# for ref.
###########################################################




#### ---- METHODS SECTION ---- ####

ReadData <- function(data_path, count_type) {
    if (count_type == "s") {
        return(raw_mat <- Read10X(data_path))
    }
    if (count_type == "r") {
        return(raw_mat <- read.table(data_path))
    }
}


CombineExpMat <- function(mat1, mat2, mat1_annot, mat2_annot) {
    mat1_so <- CreateSeuratObject(mat1,
        project = mat1_annot
    )

    mat2_so <- CreateSeuratObject(mat2,
        project = mat2_annot
    )

    merged_so <- merge(mat1_so,
        mat2_so,
        add.cell.ids = c(mat1_annot, mat2_annot),
        project = "combined"
    )


    # get the raw gene count matrix from the merged seurat object
    merged_counts <- GetAssayData(merged_so, slot = "counts")


    # convert to matrix to return
    return(merged_counts_mat <- as.matrix(merged_counts))
}


SaveCombinedMat <- function(merged_counts_mat, save_path) {
    write.table(merged_counts_mat,
        gzfile(save_path),
        sep = "\t"
    )
}

###################################




#### ---- MAIN SECTION ---- ####

pbmc_3k <- ReadData("../data/reference/pbmc_3k",
    count_type = "s"
)
pbmc_3k %>% head(n = c(5, 5))


tall <- ReadData("../data/ega_data/tall_ega_data/10X_count/bamtofastq/outs/filtered_feature_bc_matrix/",
    count_type = "s"
)
tall %>% head(n = c(5, 5))


merged_counts_mat <- CombineExpMat(
    mat1 = pbmc_3k,
    mat2 = tall,
    mat1_annot = "pbmc",
    mat2_annot = "tall"
)
merged_counts_mat %>% head(n = c(5, 5))
merged_counts_mat %>% tail(n = c(5, 5))



SaveCombinedMat(merged_counts_mat,
    save_path = "../proc/pbmc3k_tall_ref_merged_counts.tsv.gz"
)

################################
