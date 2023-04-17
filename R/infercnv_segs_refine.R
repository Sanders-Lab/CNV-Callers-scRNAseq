library(tidyverse)
library(infercnv)
library(parallel)
library(foreach)
library(data.table)


step_20 <- readRDS('../outputs/tnbc1/infercnv_pbmc3k_epithelium_tnbc1/20_HMM_pred.repr_intensitiesHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj')
expr_data_df <- as.data.frame(step_20@expr.data)
tnbc_cells <- expr_data_df %>% select(contains('tnbc1'))
tnbc_cells %>% ncol()
tnbc_cells <- cbind(step_20@gene_order, tnbc_cells)
tnbc_cells <- tnbc_cells %>% dplyr::rename(seqnames = 'chr',
                                           end = 'stop')
tnbc_cells %>% head(n = c(5,5))



seg_brk <- integer()
seg_interval <- 5000000
for(i in 2:nrow(tnbc_cells)){

    if(tnbc_cells$start[i] - tnbc_cells$end[(i - 1)] > seg_interval & tnbc_cells$seqnames[i] == tnbc_cells$seqnames[(i - 1)]){
        
        seg_brk <- c(seg_brk, (i - 1))
    }

}
seg_brk %>% length()
seg_brk



GetSegDistBreaks <- function(pred_mat,
                             seg_interval = 1000000){

    seg_brk <- integer()
    for(i in 2:nrow(pred_mat)){

        if(pred_mat$start[i] - pred_mat$end[(i - 1)] > seg_interval & pred_mat$seqnames[i] == pred_mat$seqnames[(i - 1)]){

            segs_brks_dist <- c(segs_brks_dist, (i - 1))
        }

    }
    

    return(segs_brks_dist)
}


GetInfercnvSegs <- function(per_chr_df, segs_brks_dist){

    per_chr_df <- as.data.frame(per_chr_df)
    per_chr_df <- per_chr_df[per_chr_df[,4] != "",]
    if(nrow(per_chr_df) > 0){

        start <- per_chr_df$start[1]
        tmp_df <- data.frame(seqnames = character(1),
                             start = character(1),
                             end = character(1),
                             cell_name = character(1),
                             cnv_state = character(1))
        per_chr_cnv_segs <- data.frame()


        for(row_num in 1:(nrow(per_chr_df) - 1)){

            if (per_chr_df[row_num,4] != per_chr_df[(row_num + 1),4] |
                (row_num %in% segs_brks_dist) |
                per_chr_df$seqnames[row_num] != per_chr_df$seqnames[(row_num + 1)]) {

                if(per_chr_df$seqnames[row_num] != 23){

                    tmp_df$seqnames <- per_chr_df$seqnames[row_num]
                } else{

                    tmp_df$seqnames <- "chrX"
                }
                tmp_df$start <- start
                tmp_df$end <- per_chr_df$end[row_num]
                tmp_df$cell_name <- colnames(per_chr_df)[4]
                tmp_df$cnv_state <- per_chr_df[row_num,4]

                per_chr_cnv_segs <- rbind(per_chr_cnv_segs, tmp_df)


                start <- per_chr_df$start[(row_num + 1)]
            }
        }
    }

    return(per_chr_cnv_segs)
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



n_cores <- 7
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)



infercnv_cnv_segs <- data.frame()
system.time(infercnv_cnv_segs <- foreach(cell_num = 4:ncol(tnbc_cells_label), .combine = "rbind") %dopar% {
    if(tnbc_cells_label[,cell_num] %>% grep ("amp", .) %>% length() > 0 |
       tnbc_cells_label[,cell_num] %>% grep ("del", .) %>% length() > 0){
        GetInfercnvSegs(tnbc_cells_label[, c(1:3, cell_num)], segs_brks_dist)

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



source("../../digital_karyotype/R/utils.R")

