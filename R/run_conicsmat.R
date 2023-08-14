library(CONICSmat)
library(tidyverse)
library(biomaRt)
library(data.table)




ReadData <- function(obs_ref_dir,
                     reg_dir = "../proc/chromosome_arm_positions_grch38.txt") {
    obs_ref_mat <- fread(obs_ref_dir)
    obs_ref_mat <- setnames(obs_ref_mat, "V1", "gene")


    regions <- fread(reg_dir,
        sep = "\t",
        header = T
    )


    return(list(obs_ref_mat, regions))
}




GetGeneLoc <- function(gene_names) {
    hmart <- useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl"
    )


    gene_pos <- biomaRt::getBM(
        attributes = c(
            "ensembl_gene_id",
            "hgnc_symbol",
            "chromosome_name",
            "start_position",
            "end_position"
        ),
        filters = "hgnc_symbol",
        values = gene_names,
        mart = hmart,
        useCache = F
    )


    return(gene_pos)
}




ProcConicsmatData <- function(obs_ref_mat,
                              gene_pos,
                              normal_prefix,
                              aneu_prefix) {

    obs_ref_mat <- as.data.frame(obs_ref_mat)
    rownames(obs_ref_mat) <- obs_ref_mat$gene
    obs_ref_mat$gene <- NULL

    obs_ref_mat <- filterMatrix(obs_ref_mat,
        gene_pos[, "hgnc_symbol"],
        minCells = 5
    )
    norm_factor <- calcNormFactors(obs_ref_mat)


    all_cells <- colnames(obs_ref_mat)
    normal_cells <- c(1:(colnames(obs_ref_mat) %>%
        str_detect(., normal_prefix) %>%
        sum(.)))
    names(normal_cells) <- colnames(obs_ref_mat) %>%
        grep(normal_prefix,
            .,
            value = T
        )

    aneu_cells <- c(1:(colnames(obs_ref_mat) %>%
        str_detect(., aneu_prefix) %>%
        sum(.)))
    names(aneu_cells) <- colnames(obs_ref_mat) %>%
        grep(aneu_prefix,
            .,
            value = T
        )


    return(list(
        obs_ref_mat,
        norm_factor,
        normal_cells,
        aneu_cells
    ))
}




RunConicsmat <- function(obs_ref_mat,
                         regions,
                         norm_factor,
                         gene_pos,
                         normal_cells,
                         aneu_cells,
                         save_dir) {
    cell_prefix <- c(rep("Ref", length(normal_cells)), rep("Observed", length(aneu_cells)))

    cs_all <- plotAll(obs_ref_mat,
        norm_factor,
        regions,
        gene_pos,
        str_glue("{save_dir}/supervised_postProb.pdf"),
        normal = normal_cells,
        tumor = aneu_cells
    )

    hist <- plotHistogram(cs_all,
        obs_ref_mat,
        clusters = 2,
        zscoreThreshold = 4,
        cell_prefix
    )

    bin_mat <- binarizeMatrix(cs_all,
        normal = normal_cells,
        tumor = tumour_cells,
        0.8
    )

    chr_heatmap <- plotChromosomeHeatmap(obs_ref_mat,
        normal = normal_cells,
        plotcells = all_cells,
        gene_pos = gene_pos,
        chr = T,
        windowsize = 101,
        expThresh = 0.2,
        thresh = 1
    )


    ggsave(
        str_glue("{save_dir}/all_cells_chr_heatmap.pdf"),
        chr_heatmap
    )

    ggsave(
        str_glue("{save_dir}/cna_hist.pdf"),
        hist
    )


    return(bin_mat)
}


## TESTS
# ret_list <- ReadData(obs_ref_dir = "../proc/pbmc3k_tall_ref_merged_counts.tsv.gz")
# obs_ref_mat <- ret_list[[1]]
# regions <- ret_list[[2]]
# obs_ref_mat[1:5, 1:5]
# obs_ref_mat %>% tail(n = c(5, 5))
# obs_ref_mat %>% nrow()
# all_cells %>% length()
# 
# 
# GetGeneLoc(obs_ref_mat$gene)[1:5, 1:5]
# colnames(obs_ref_mat[, -1]) %>%
#     str_detect(., "pbmc") %>%
#     sum(.)
# colnames(obs_ref_mat) %>%
#     grep("tall", ., value = T) %>%
#     length()
# colnames(obs_ref_mat) %>%
#         grep("tall",
#             .,
#             value = T
#         )
# 
#     obs_ref_mat <- as.data.frame(obs_ref_mat)
#     rownames(obs_ref_mat) <- obs_ref_mat$gene
#     obs_ref_mat$gene <- NULL
# gene_pos[1:5, "hgnc_symbol"]
# check <- filterMatrix(obs_ref_mat,
#         gene_pos[, "hgnc_symbol"],
#         minCells = 5
#     )
# 
# 1:(colnames(check) %>%
#     str_detect(., "pbmc") %>%
#     sum(.))
# 
# 
# 
# check <- as.data.frame(obs_ref_mat)[1:5,1:5]
# rownames(check) <- check$gene
# check$gene  <- NULL
# check
# gene_pos
# 
# ret_list <- ProcConicsmatData(obs_ref_mat,
#     gene_pos,
#     normal_prefix = "pbmc",
#     aneu_prefix = "tall"
# )
# obs_ref_mat <- ret_list[1]
# norm_factor <- ret_list[2]
# normal_cells <- ret_list[3]
# aneu_cells <- ret_list[4]
# 
# normal_cells <- as.character(normal_cells)
# aneu_cells <- as.vector(aneu_cells)
# typeof(normal_cells)
# normal_cells
# aneu_cells
# norm_factor
# gene_pos
# 
# RunConicsmat(obs_ref_mat,
#     regions,
#     norm_factor,
#     gene_pos,
#     normal_cells,
#     aneu_cells,
#     save_dir = "../outputs/tall_scnova/conicsmat_pbmc3k_tall"
# )
# 
# 
################################################################################
# MAIN SECTION

cmd_args <- commandArgs(trailingOnly = T)


if(length(cmd_args) < 5) {

    stop("\nPlease pass in the following inputs:\n
                        count_mat_path\n
                        region_dt (optional)\n
                        normal_prefix\n
                        aneu_prefix\n
                        out_dir\n")
} else {

    writeLines(str_glue("\nThe passed in args are:\n
                        count_mat_path: {cmd_args[1]}\n
                        region_dt: {cmd_args[2]}\n
                        normal_prefix: {cmd_args[3]}\n
                        aneu_prefix: {cmd_args[4]}\n
                        out_dir: {cmd_args[5]}\n"))
}



if (!dir.exists(cmd_args[5])) {

    dir.create(cmd_args[5])
    writeLines("output dir is created!")
}


# use the whole chr if the region dt is blank
if (cmd_args[2] == "") {
    read_list <- ReadData(obs_ref_dir = cmd_args[1])

} else {
    read_list <- ReadData(obs_ref_dir = cmd_args[1],
                          reg_dir = cmd_args[2])
}

obs_ref_mat <- read_list[[1]]
regions <- read_list[[2]]

gene_pos <- GetGeneLoc(gene_names = obs_ref_mat$gene)



ret_list <- ProcConicsmatData(obs_ref_mat,
    gene_pos,
    normal_prefix = cmd_args[3],
    aneu_prefix = cmd_args[4]
)

obs_ref_mat <- ret_list[1]
norm_factor <- ret_list[2]
normal_cells <- ret_list[3]
aneu_cells <- ret_list[4]




RunConicsmat(obs_ref_mat,
    regions,
    norm_factor,
    gene_pos,
    normal_cells,
    aneu_cells,
    save_dir = cmd_args[5]
)

