library(data.table)
library(magrittr)
library(tidyverse)
library(GenomicRanges)



setDTthreads(16)
nb_segs <- fread("../proc/numbat_tnbc_cna_segs.tsv.gz")
nb_segs


# one time tidy
# nb_segs$cell_name <- paste0("tnbc1_", str_sub(nb_segs$cell_name, 1, -3))
# nb_segs

# fwrite(nb_segs,
#        "../proc/numbat_tnbc_cna_segs.tsv.gz",
#        sep = "\t")

ck_segs <- fread("../proc/copykat_aneu_cnaSegs_vals_labels.tsv.gz")
ck_segs


# ck_segs$seqnames <- paste0("chr", ck_segs$seqnames)
# fwrite(ck_segs,
#        "../proc/copykat_aneu_cnaSegs_vals_labels.tsv.gz",
#        sep = "\t")



if_segs <- fread("../proc/infercnv_segs_cnv_label_allcells.tsv.gz")
if_segs


# create gro
nb_gro <- makeGRangesFromDataFrame(nb_segs,
                                   keep.extra.columns = T)
nb_gro


ck_gro <- makeGRangesFromDataFrame(ck_segs,
                                   keep.extra.columns = T)
ck_gro


if_gro <- makeGRangesFromDataFrame(if_segs,
                                   keep.extra.columns = T)
if_gro


gene_ranges <- fread("../data/gene_chr_loc/hg38_gencode_v27.txt",
                     fill = T,
                     sep = "\t")
gene_ranges[gene_ranges$V5 != NA,] 
gene_ranges$V5 %>% is.na() %>% length()
gene_ranges$V5 <- NULL
gene_ranges


gene_ranges <- gene_ranges %>% dplyr::rename(gene = "V1",
                                             seqnames = "V2",
                                             start = "V3",
                                             end = "V4")
gene_ranges$seqnames %>% unique()


# dropping chr Y and M
gene_ranges <- gene_ranges[seqnames != "chrY",]
gene_ranges <- gene_ranges[seqnames != "chrM",]
gene_ranges


gene_gro <- makeGRangesFromDataFrame(gene_ranges,
                                     keep.extra.columns = T)
gene_gro


gene_ranges <- as.data.table(gene_gro)
gene_ranges$fill <- "gene"
gene_ranges



tmp_plot <- "../plots/tmp/tmp_plot.pdf"
source("../../digital_karyotype/R/add_info_plt_layers.R")



col_gene <- "blue"
names(col_gene) <- "gene"
col_gene


gene_dk <- plotSkeletonH2(ideo_df)


gene_dk <- addPltLayers(gene_dk,
                        gene_ranges,
                        2,
                        "fill",
                        col_gene,
                        "Genes",
                        title = "Gene Locations")
ggsave(tmp_plot)


cna_glist <- GRangesList(nb_gro, if_gro, ck_gro) 
names(cna_glist) <- c("numbat_cna", "infercnv_cna", "copykat_cna")
cna_glist


# grouping by cells 
nb_segs_cell_list <- nb_segs %>% group_by(cell_name) %>% group_split(.) 
if_segs_cell_list <- if_segs %>% group_by(cell_name) %>% group_split(.)
ck_segs_cell_list <- ck_segs %>% group_by(cell_name) %>% group_split(.)


# getting the least common cells
common_cells <- Reduce(intersect, list(nb_segs$cell_name, 
                                       if_segs$cell_name,
                                       ck_segs$cell_name))
common_cells %>% length()
common_cells %>% head()


# subsetting the higher cell dfs by the common cells
nb_segs_sub <- nb_segs %>% dplyr::filter(cell_name %in% common_cells)
nb_segs_sub %>% head()
nb_segs_sub$cell_name %>% unique() %>% length()


if_segs_sub <- if_segs %>% dplyr::filter(cell_name %in% common_cells)
if_segs_sub$cell_name %>% unique() %>% length()


# QND TESTS WITH 1 CELL CNA
cell1 <- nb_segs_sub$cell_name[1]

nb1_cell <- nb_segs_sub[cell_name == cell1]
if1_cell <- if_segs_sub[cell_name == cell1]
ck1_cell <- ck_segs[cell_name == cell1]




# THIS FUNCTION WORKS!!!!
CalcOverlapProp <- function(tool1_cell, tool2_cell){

    tool1_cell_gro <- makeGRangesListFromDataFrame(tool1_cell,
                                                keep.extra.columns = T)
    tool2_cell_gro <- makeGRangesListFromDataFrame(tool2_cell,
                                                keep.extra.columns = T)



    # findOverlaps(if_1cell_gro, nb_1cell_gro)
    # I WAS HERE LAST!! AND IT WORKED!!!
    # now need to implement for all the cells.
    overlaps <- GenomicRanges::findOverlaps(tool1_cell_gro, tool2_cell_gro)

    tool1_hits <- tool1_cell[overlaps@from,]
    tool2_hits <- tool2_cell[overlaps@to,]


    # getting the unique ranges
    # from the common ranges
    tool1_hits_retain <- tool1_hits[cnv_state == tool2_hits$cnv_state] %>% unique()
    tool2_hits_retain <- tool2_hits[cnv_state == tool1_hits$cnv_state] %>% unique()


    # this gives the percent
    # NEED TO PARALLELIZE THIS!!!
    return(list(nrow(tool1_hits_retain) / nrow(tool1_cell),
                nrow(tool2_hits_retain) / nrow(tool2_cell)))

}

nb_no_loh <- nb1_cell[cnv_state != "loh"]
CalcOverlapProp(nb_no_loh, ck1_cell)[[1]]
ck1_cell


library(foreach)
library(parallel)


comp_cnv <- data.table(numbat_infercnv = double(length(common_cells)),
                       infercnv_numbat = double(length(common_cells)),
                       numbat_copykat = double(length(common_cells)),
                       copykat_numbat = double(length(common_cells)),
                       infercnv_copykat = double(length(common_cells)),
                       copykat_infercnv = double(length(common_cells)))



n_cores <- 15
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


nb_segs_sub <- nb_segs_sub[cnv_state != "loh"]
# very gross solution, took 2.5 mins
# but worked
start <- Sys.time()
foreach(c = 1:length(common_cells)) %do% {


    temp_list <- CalcOverlapProp(nb_segs_sub[cell_name == common_cells[c]],
                                 if_segs_sub[cell_name == common_cells[c]])

    comp_cnv$numbat_infercnv[c] <- temp_list[[1]]
    comp_cnv$infercnv_numbat[c] <- temp_list[[2]]

    temp_list <- CalcOverlapProp(nb_segs_sub[cell_name == common_cells[c]],
                                 ck_segs[cell_name == common_cells[c]])

    comp_cnv$numbat_copykat[c] <- temp_list[[1]]
    comp_cnv$copykat_numbat[c] <- temp_list[[2]]

    temp_list <- CalcOverlapProp(if_segs_sub[cell_name == common_cells[c]],
                                 ck_segs[cell_name == common_cells[c]])

    comp_cnv$infercnv_copykat[c] <- temp_list[[1]]
    comp_cnv$copykat_infercnv[c] <- temp_list[[2]]
}
Sys.time() - start
comp_cnv %>% head()

parallel::stopCluster(cluster)


fwrite(comp_cnv,
       "../proc/cnv_callers_comparison.tsv.gz",
       sep = "\t")



nb_cell_gro <- makeGRangesListFromDataFrame(nb_cell,
                                            keep.extra.columns = T)
if_cell_gro <- makeGRangesListFromDataFrame(if_cell,
                                            keep.extra.columns = T)
ck_cell_gro <- makeGRangesListFromDataFrame(ck_cell,
                                            keep.extra.columns = T)



# findOverlaps(if_1cell_gro, nb_1cell_gro)
# I WAS HERE LAST!! AND IT WORKED!!!
# now need to implement for all the cells.
nb_if_overlap <- GenomicRanges::findOverlaps(nb_cell_gro, if_cell_gro)
nb_ck_overlap <- GenomicRanges::findOverlaps(nb_cell_gro, ck_cell_gro)
if_ck_overlap <- GenomicRanges::findOverlaps(if_cell_gro, ck_cell_gro)

nb_hits <- nb_1cell[nb_if_overlap@from,]
if_hits <- if_1cell[nb_if_overlap@to,]


nb_hits %>% head()
if_hits %>% head()


nb_hits_common <- nb_hits[cnv_state == if_hits$cnv_state]
if_hits_common <- if_hits[cnv_state == nb_hits$cnv_state]
nb_hits_common %>% head()
if_hits_common %>% head()


# getting the unique ranges
nb_hits_retain <- nb_hits_common %>% unique()
if_hits_retain <- if_hits_common %>% unique()
nb_hits_retain
if_hits_retain



nb_no_loh <- nb_1cell[cnv_state != "loh"]


# this gives the percent
# NEED TO PARALLELIZE THIS!!!
nrow(nb_hits_retain) / nrow(nb_no_loh)
nrow(if_hits_retain) / nrow(if_1cell)


source("../../digital_karyotype/R/add_info_plt_layers.R")
tmp_plot <- "../plots/tmp/tmp_plot.pdf"

ncnv_col <- c("red", "blue", "green")
names(ncnv_col) <- c("amp", "del", "loh")

cnv_col <- c("orange", "violet")
names(cnv_col) <- c("amp", "del")


cnv_dk <- plotSkeletonH2(ideo_df)


cnv_dk <- addPltLayers(cnv_dk,
                        nb_hits_common,
                        2,
                        "cnv_state",
                        ncnv_col,
                        "Numbat CNV Calls")


cnv_dk <- addPltLayers(cnv_dk,
                        if_hits_common,
                        2,
                        "cnv_state",
                        cnv_col,
                        "InferCNV CNV Calls",
                        y_min = 1,
                        y_max = 1.5,
                        title = "CNV Callers Comparison 1 cell")

ggsave(tmp_plot)


# plot comparison between tools
comp_df <- fread("../proc/cnv_callers_comparison.tsv.gz")


avg_nb_if <- mean(comp_df$numbat_infercnv)
avg_nb_if


avg_if_nb <- mean(comp_df$infercnv_numbat)
avg_if_nb


avg_nb_ck <- mean(comp_df$numbat_copykat)
avg_nb_ck


avg_ck_nb <- mean(comp_df$copykat_numbat)
avg_ck_nb


df_melt <- melt(comp_df)
df_melt %>% head()

df_melt <- df_melt %>% rename(Tools_Compared = "comp")
df_melt



tmp_plot <- "../plots/tmp/tmp_plot.pdf"


ggplot(data = df_melt,
       aes(Tools_Compared[1], value, color = Tools_Compared)) +
geom_violin() + 
geom_point() +
geom_jitter() +
scale_color_manual(values = c("red",
                              "orange",
                              "salmon",
                              "blue",
                              "yellow",
                              "skyblue")) +
facet_wrap(~ Tools_Compared) +
theme(axis.ticks.x = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank()) +
ggtitle("Fraction of common calls between any two tools")





ggsave(tmp_plot)
