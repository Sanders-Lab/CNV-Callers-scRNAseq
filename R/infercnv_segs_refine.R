library(tidyverse)
library(infercnv)
library(parallel)
library(foreach)
library(data.table)
library(magrittr)


step_20 <- readRDS('../outputs/tnbc1/infercnv_pbmc3k_epithelium_tnbc1/20_HMM_pred.repr_intensitiesHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj')
expr_data_df <- as.data.frame(step_20@expr.data)
tnbc_cells <- expr_data_df %>% select(contains('tnbc1'))
tnbc_cells %>% ncol()
tnbc_cells <- cbind(step_20@gene_order, tnbc_cells)
tnbc_cells <- tnbc_cells %>% dplyr::rename(seqnames = 'chr',
                                           end = 'stop')
tnbc_cells %>% head(n = c(5,5))





GetSegDistBreaks <- function(pred_mat,
                             seg_interval = 5000000){

    segs_brks_dist <- c()
    for(i in 1:(nrow(pred_mat) - 1)){

        if((pred_mat$start[(i + 1)] - pred_mat$end[i]) > seg_interval &&
           pred_mat$seqnames[i] == pred_mat$seqnames[(i + 1)]){

            segs_brks_dist <- c(segs_brks_dist, pred_mat$idx[i])
        }

    }
    

    return(segs_brks_dist)
}


GetInfercnvSegs <- function(pred_mat){

    pred_mat <- as.data.frame(pred_mat)
    pred_mat <- pred_mat[pred_mat[,5] != "",]

    segs_brks_dist <- GetSegDistBreaks(pred_mat)


    if(nrow(pred_mat) > 0){

        start <- pred_mat$start[1]
        tmp_df <- data.frame(seqnames = character(1),
                             start = character(1),
                             end = character(1),
                             cell_name = character(1),
                             cnv_state = character(1))
        cnv_segs <- data.frame()


        for(row_num in 1:(nrow(pred_mat) - 1)){

            if (pred_mat[row_num,5] != pred_mat[(row_num + 1),5] |
                (pred_mat$idx[row_num] %in% segs_brks_dist) |
                pred_mat$seqnames[row_num] != pred_mat$seqnames[(row_num + 1)]) {

                if(pred_mat$seqnames[row_num] != 23){

                    tmp_df$seqnames <- pred_mat$seqnames[row_num]
                } else{

                    tmp_df$seqnames <- "chrX"
                }
                tmp_df$start <- start
                tmp_df$end <- pred_mat$end[row_num]
                tmp_df$cell_name <- colnames(pred_mat)[5]
                tmp_df$cnv_state <- pred_mat[row_num,5]

                cnv_segs <- rbind(cnv_segs, tmp_df)


                start <- pred_mat$start[(row_num + 1)]
            }
        }
    }

    return(cnv_segs)
}




segs_brks_dist <- integer()
segs_brks_dist <- GetSegDistBreaks(tnbc_cells, seg_interval = 5000000)
segs_brks_dist %>% length()
segs_brks_dist


# assigning labels
unique(tnbc_cells[,4])
call_mat <- tnbc_cells[,4:ncol(tnbc_cells)]
tnbc_cells_label <- call_mat %>% 
    apply(., 2, function(x) ifelse(x > 1, 'amp', ifelse(x < 1, 'del', '')))
tnbc_cells_label <- cbind(tnbc_cells[,1:3], tnbc_cells_label)
tnbc_cells_label %>% head(n = c(5,5))
ncol(tnbc_cells_label)
tnbc_cells_label[,1:4] %>% 
    filter(seqnames == "chr3")


# saving the gene level calls
fwrite(x = tnbc_cells_label,
       file = "../proc/infercnv_tnbc1_gene_calls.tsv.gz",
       sep = "\t")


tnbc_cells_label <- fread("../proc/infercnv_tnbc1_gene_calls.tsv.gz")
tnbc_cells_label[,1:5]


# appending a idx col
tnbc_cells_label$idx <- 1:nrow(tnbc_cells_label)
setcolorder(tnbc_cells_label,
            c(ncol(tnbc_cells_label), 1:(ncol(tnbc_cells_label) - 1)))
# tnbc_cells_label[,c(ncol(tnbc_cells_label), 1:(ncol(tnbc_cells_label) - 1))]
tnbc_cells_label[,1:5]



n_cores <- 7
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)



infercnv_cnv_segs <- data.frame()
system.time(infercnv_cnv_segs <- 
    foreach(cell_num = 5:ncol(tnbc_cells_label), .combine = "rbind") %dopar% {
    if(tnbc_cells_label[,cell_num] %>% grep ("amp", .) %>% length() > 0 |
       tnbc_cells_label[,cell_num] %>% grep ("del", .) %>% length() > 0){
        GetInfercnvSegs(tnbc_cells_label[, c(1:4, cell_num)], segs_brks_dist)

    }

})

infercnv_cnv_segs %>% 
    head()
infercnv_cnv_segs$cell_name %>% 
    unique() %>% 
    length()
parallel::stopCluster(cluster)


fwrite(x = infercnv_cnv_segs,
       file = "../proc/infercnv_tnbc1_segs_refined.tsv.gz",
       sep = "\t")


infercnv_cnv_segs <- as.data.table(infercnv_cnv_segs)
infercnv_cnv_segs[cell_name == infercnv_cnv_segs$cell_name[1] &
                  seqnames == "chr3",]


tnbc_cells[300:302,1:4]
tnbc_cells[,1:4] %>%
    filter(seqnames == "chr3") %>% 
    filter(.[,4] == 0.5)


load_cnv_segs <- fread( "../proc/infercnv_tnbc1_segs_refined.tsv.gz")
load_cnv_segs



source("../../digital_karyotype/R/utils.R")
cna_colors <- c('red', 'blue')
names(cna_colors) <- c('amp', 'del')


plot_digital_karyotype(plot_ideo_only = F,
                       layers_h2 = load_cnv_segs[cell_name == load_cnv_segs$cell_name[1],],
                       fill_arg = "cnv_state",
                       colors_arg_h2 = cna_colors,
                       legend_title_arg = "CNV Calls",
                       sub_title_arg = str_glue("InferCNV CNV calls - {load_cnv_segs$cell_name[1]}"),
                       plot_both_haplotypes = F,
                       save_digital_karyotype = T,
                       save_dir = "infercnv_digital_karyotypes/cnv_refine")



# gene level calls
gene_calls <- fread("../proc/infercnv_tnbc1_gene_calls.tsv.gz")
gene_calls[,1:5]
gene_calls$idx <- rep(1:nrow(gene_calls))
gene_calls$idx
setcolorder(gene_calls,
            c(ncol(gene_calls),1:(ncol(gene_calls) - 1)))
gene_calls[,1:5]
gene_calls[,c(1:4,7)]
gene_calls[, get("n_cores")] %>% grep ("del", .) %>% length() > 0
get("n_cores")
gene_calls_df <- as.data.frame(gene_calls)


n_cores <- 7
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)

# running segmentation using updated one arg
infercnv_cnv_segs <- data.frame()
system.time(infercnv_cnv_segs <- 
    foreach(cell_num = 5:ncol(gene_calls_df), .combine = "rbind") %dopar% {
    if(gene_calls_df[, cell_num] %>% grep ("amp", .) %>% length() > 0 |
       gene_calls_df[, cell_num] %>% grep ("del", .) %>% length() > 0){
        GetInfercnvSegs(gene_calls_df[, c(1:4, cell_num)])

    }

})

parallel::stopCluster(cluster)
infercnv_cnv_segs %>% 
    head()
infercnv_cnv_segs <- as.data.table(infercnv_cnv_segs)


fwrite(x = infercnv_cnv_segs,
       file = "../proc/infercnv_tnbc1_segs_refined.tsv.gz",
       sep = "\t")


source("../../digital_karyotype/R/utils.R")
cna_colors <- c('red', 'blue')
names(cna_colors) <- c('amp', 'del')


plot_digital_karyotype(plot_ideo_only = F,
                       layers_h2 = infercnv_cnv_segs[cell_name == infercnv_cnv_segs$cell_name[1],],
                       fill_arg = "cnv_state",
                       colors_arg_h2 = cna_colors,
                       legend_title_arg = "CNV Calls",
                       sub_title_arg = str_glue("InferCNV CNV calls - {infercnv_cnv_segs$cell_name[1]}"),
                       plot_both_haplotypes = F,
                       save_digital_karyotype = T,
                       save_dir = "infercnv_digital_karyotypes/cnv_refine")





gene_calls
one_cell_call <- gene_calls[,1:4]
one_cell_call %>% head()
one_cell_call$cell_name <- rep("tnbc1_AAACCTGCACCTTGTC", nrow(one_cell_call))
one_cell_call <- one_cell_call[,c(1:3,5,4)]
colnames(one_cell_call)[5] <- "cnv_state"
one_cell_call <- one_cell_call[cnv_state != "",]
one_cell_call %>% head()



dig_kar <- plot_digital_karyotype(plot_ideo_only = F,
                                  layers_h2 = one_cell_call,
                                  fill_arg = "cnv_state",
                                  colors_arg_h2 = cna_colors,
                                  legend_title_arg = "CNV Calls",
                                  sub_title_arg = str_glue("InferCNV CNV calls - tnbc1_AAACCTGCACCTTGTC"),
                                  plot_both_haplotypes = F,
                                  save_digital_karyotype = T,
                                  save_dir = "infercnv_digital_karyotypes/gene_calls")



one_cell_call
one_cell_call[seqnames == "chr1" &
              cnv_state == "del"] %>%
              print(n = 111)


one_cell_call
gene_calls[,1:3]
seg_brks <- GetSegDistBreaks(one_cell_call)
seg_brks <- GetSegDistBreaks(gene_calls)
seg_brks
seg_brks %>%
    length()


one_cell_call[seg_brks,]


gene_calls[70:73,1:5]



gene_calls$start[302] - gene_calls$end[301]
