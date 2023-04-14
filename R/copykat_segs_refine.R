library(tidyverse)
library(data.table)
library(foreach)
library(parallel)


raw_res <- fread('../outputs/tnbc1/copykat_pbmc3k_epithelium_tnbc1/tnbc1_copykat_CNA_raw_results_gene_by_cell.txt',
                 header = T)
idx <- c(1:nrow(raw_res))
raw_res <- cbind(idx, raw_res)
raw_res[1:5,1:10]
raw_res %>% nrow()
raw_res <- as.data.frame(raw_res)


cna_res <- fread('../outputs/tnbc1/copykat_pbmc3k_epithelium_tnbc1/tnbc1_copykat_CNA_results.txt',
                 header = T)
cna_res[1:5,1:10]
cna_res %>% nrow()
cna_res <- as.data.frame(cna_res)

copykat_preds <- fread('../outputs/tnbc1/copykat_pbmc3k_epithelium_tnbc1/tnbc1_copykat_prediction.txt',
                       header = T)
copykat_preds %>% head
copykat_preds <- as.data.frame(copykat_preds)
copykat_preds %>% grep("tnbc", .)

ck_preds <- copykat_preds[grep("tnbc", copykat_preds$cell.names),]

# binning is done by ck
# counting the number of instances when the 
# chrompos diff is more than 1MB
# will break segs here
count <- 0
seg_brk <- integer()
for(i in 2:nrow(cna_res)){

    if(cna_res$abspos[i] - cna_res$abspos[(i - 1)] > 1000000 & cna_res$chrom[i] == cna_res$chrom[(i - 1)]){
        
        count <- count + 1
        seg_brk <- c(seg_brk, i)
    }

}

# checking with raw res
count <- 0
seg_brk <- integer()
seg_interval <- 5000000
for(i in 2:nrow(raw_res)){

    if(raw_res$abspos[i] - raw_res$abspos[(i - 1)] > seg_interval & raw_res$chrom[i] == raw_res$chrom[(i - 1)]){
        
        count <- count + 1
        seg_brk <- c(seg_brk, (i - 1))
    }

}
seg_brk %>% length()


seg_brk <- integer()
seg_interval <- 5000000
for(i in 2:nrow(sub_raw_loc)){

    if(sub_raw_loc$start_position[i] - sub_raw_loc$end_position[(i - 1)] > seg_interval &&
       sub_raw_loc$chromosome_name[i] == sub_raw_loc$chromosome_name[(i - 1)]){
        
        seg_brk <- c(seg_brk, (i - 1))
    }
}
seg_brk %>% length()



raw_res$abspos[21] - raw_res$abspos[(21 - 1)]
cna_res$abspos[21] - cna_res$abspos[(21 - 1)]
count
raw_res[26:28, 1:5]
cna_res[idx,1:5]

cna_res[50:52,1:5]
cna_res[515:517,1:5]
cna_res$chrompos[5] - cna_res$chrompos[4]


table(ck_preds$copykat.pred == 'aneuploid')
aneu_cells <- ck_preds[ck_preds$copykat.pred == 'aneuploid',]
aneu_cells
raw_res %>% head(n = c(5,5))
mod_gene_mat <- raw_res[,aneu_cells$cell.names]
cna_gene_mat <- mod_gene_mat %>% apply(., 2, function(x) ifelse(x > 0.03, 'amp', ifelse(x < -0.03, 'del', '')))
mod_gene_mat[1:5,1:5]
cna_gene_mat[1:5,1:5]
cna_gene_mat %>% nrow
aneu_cells
mod_gene_mat

sub_raw_res <- raw_res[,8:=NULL]
sub_raw_res <- cbind(sub_raw_res[,1:7], cna_gene_mat)
sub_raw_res[1:5,1:8]
sub_raw_res[,1:8]
sub_raw_res %>% ncol() - 5
sub_raw_res %>% nrow()


sub_raw_loc <- sub_raw_res
sub_raw_loc[1:5,1:8]


fwrite(sub_raw_loc,
       "../proc/copykat_geneLevel_230331.tsv.gz",
       sep = "\t")



# READ THIS FILE
sub_raw_loc <- fread("../proc/copykat_geneLevel_230331.tsv.gz")
sub_raw_loc <- as.data.frame(sub_raw_loc)
sub_raw_loc %>% head(n = c(5,8))

# I WAS HERE
# BREAK SEGS AT LOCS WHERE IT IS MORE THAN 1MB
sub_raw_loc %>%
    group_by(chromosome_name) %>%
    group_split(.) %>%
    .[[20]] %>%
    .[, 1:10]

gene_calls_chr_split <- sub_raw_loc %>%
    group_by(chromosome_name) %>%
    group_split(.)
gene_calls_chr_split


# TESTS
chr1_df <- gene_calls_chr_split[[1]]
chr1_df <- chr1_df[,1:8]
chr1_df <- as.data.frame(chr1_df)
chr1_df %>% head
chr1_df[1:5,7:ncol(chr1_df)]
chr1_df %>% nrow()
seg_brk[1:20]




# only for tibble based
prev_call <- chr1_df[1,][[8]]
prev_call



# LOGIC START HERE
# TESTS
chr1_df <- chr1_df[chr1_df[,8] != "",]
chr1_df %>% nrow()
start <- chr1_df$start_position[1]
start
tmp_df <- data.frame(seqnames = character(1),
                     start = character(1),
                     end = character(1),
                     cnv_state = character(1))
cnv_segs <- data.frame()
for(row_num in 1:(nrow(chr1_df) - 1)){

    if (chr1_df[row_num,8] != chr1_df[(row_num + 1),8] || (chr1_df$idx[row_num] %in% seg_brk)) {


        tmp_df$seqnames <- paste0("chr", chr1_df$chromosome_name[row_num])
        tmp_df$start <- start
        tmp_df$end <- chr1_df$end_position[row_num]
        tmp_df$cnv_state <- chr1_df[row_num,8]

        cnv_segs <- rbind(cnv_segs, tmp_df)


        start <- chr1_df$start_position[(row_num + 1)]
    }
}
cnv_segs %>% nrow()
cnv_segs %>% head()
cnv_segs %>% tail()


source("../../digital_karyotype/R/add_info_plt_layers.R")
tmp_plot <- "../plots/tmp/tmp_plot.pdf"
cnv_colors <- c('red', 'blue')
names(cnv_colors) <- c('amp', 'del')



cnv_dig_kar <- plotSkeletonH2(ideo_df)
cnv_dig_kar <- addPltLayers(cnv_dig_kar,
                            cnv_segs,
                            2,
                            "cnv_state",
                            cnv_colors,
                            "CNV Calls",
                            title = "Copykat CNV Calls")

ggsave(tmp_plot)


GetSegDistBreaks <- function(pred_mat,
                             seg_interval = 1000000){

    segs_brks_dist <- integer()
    for(i in 2:nrow(pred_mat)){

        if(pred_mat$start_position[i] - pred_mat$end_position[(i - 1)] > seg_interval &&
           pred_mat$chromosome_name[i] == pred_mat$chromosome_name[(i - 1)]){

            segs_brks_dist <- c(segs_brks_dist, (i - 1))
        }
    }
    

    return(segs_brks_dist)
}


GetCopykatSegs <- function(per_chr_df, segs_brks_dist){

    per_chr_df <- as.data.frame(per_chr_df)
    per_chr_df <- per_chr_df[per_chr_df[,8] != "",]
    start <- per_chr_df$start_position[1]
    tmp_df <- data.frame(seqnames = character(1),
                         start = character(1),
                         end = character(1),
                         cell_name = character(1),
                         cnv_state = character(1))
    per_chr_cnv_segs <- data.frame()


    for(row_num in 1:(nrow(per_chr_df) - 1)){

        if (per_chr_df[row_num,8] != per_chr_df[(row_num + 1),8] |
            (per_chr_df$idx[row_num] %in% segs_brks_dist) |
            per_chr_df$chromosome_name[row_num] != per_chr_df$chromosome_name[(row_num + 1)]) {

            if(per_chr_df$chromosome_name[row_num] != 23){

                tmp_df$seqnames <- paste0("chr", per_chr_df$chromosome_name[row_num])
            } else{

                tmp_df$seqnames <- "chrX"
            }
            tmp_df$start <- start
            tmp_df$end <- per_chr_df$end_position[row_num]
            tmp_df$cell_name <- colnames(per_chr_df)[8]
            tmp_df$cnv_state <- per_chr_df[row_num,8]

            per_chr_cnv_segs <- rbind(per_chr_cnv_segs, tmp_df)


            start <- per_chr_df$start_position[(row_num + 1)]
        }
    }

    return(per_chr_cnv_segs)
}




pred_mat_grpd <- sub_raw_loc %>%
    group_by(chromosome_name) %>%
    group_split(.)
pred_mat_grpd %>% length()
pred_mat_grpd[[2]][1:8] %>% 
    head

n_cores <- 7
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)
segs_brks_dist


# test with 1 chr at a time
# only DO
tmp_df <- data.frame()
copykat_cnv_segs <- data.frame()
cnv_segs <- data.frame()
system.time(for(cell_num in 8:ncol(pred_mat_grpd[[2]])) {

    tmp_df <- GetCopykatSegs(as.data.frame(pred_mat_grpd[[2]][, c(1:7, 8)]), segs_brks_dist)

    cnv_segs  <- rbind(cnv_segs, tmp_df)

})
as.data.frame(pred_mat_grpd[[2]] %>% 
                   .[,c(1:7,8)])
pred_mat_grpd[[2]]$chromosome_name %>% unique()
tmp_df
cnv_segs
cnv_segs[cnv_segs$seqnames == "chr2",]



system.time(cnv_segs <- foreach(cell_num = 8:ncol(pred_mat_grpd[[12]]), .combine = "rbind") %dopar% {
    if(pred_mat_grpd[[12]][,cell_num] %>% grep ("amp", .) %>% length() == 0 &
       pred_mat_grpd[[12]][,cell_num] %>% grep ("del", .) %>% length() == 0){

    } else{

        GetCopykatSegs(pred_mat_grpd[[12]][, c(1:7, cell_num)], segs_brks_dist)
    }

    #     GetCopykatSegs(as.data.frame(pred_mat_grpd[[3]] %>% 
    #                    select(all_of(c(1:7, cell_num)))), segs_brks_dist)


})
cnv_segs %>% head()
cnv_segs %>% filter(seqnames == "chr12") %>% 
    head()


#     GetCopykatSegs(as.data.frame(pred_mat_grpd[[3]] %>% 
#                    select(all_of(c(1:7, cell_num)))), segs_brks_dist)


})

parallel::stopCluster(cluster)
copykat_cnv_segs[1:5,1:5] %>% as.data.frame()
copykat_cnv_segs %>% head()
copykat_cnv_segs %>% tail
copykat_cnv_segs$cell_name %>% unique() %>% length()



n_cores <- 7
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


# I WAS HERE.
# WORKS FOR 1 CHR!
# Something wrong with chr 3 df. Check
# UP CODE FOR 1 CHR
sub_raw_loc %>%
    group_by(chromosome_name) %>% 
    group_map(testFunc(.x, .y))



# THIS WORKED FINALLY FOR ALL CHR
segs_brks_dist <- GetSegDistBreaks(sub_raw_loc, seg_interval = 5000000)
segs_brks_dist %>% length()


n_cores <- 7
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)
tmp_df <- data.frame()
copykat_cnv_segs <- data.frame()
# start_time <- Sys.time()

# cell outside, chr inside
# system.time(copykat_cnv_segs <- foreach(cell_num = 8:ncol(pred_mat_grpd[[1]]), .combine = "rbind") %do% {
# 
#     foreach(chr = 1:length(pred_mat_grpd), .combine = "rbind") %do% {
# 
# 
#                 if(pred_mat_grpd[[chr]][,cell_num] %>% grep ("amp", .) %>% length() == 0 &
#                    pred_mat_grpd[[chr]][,cell_num] %>% grep ("del", .) %>% length() == 0){
# 
#                 } else{
# 
#                     GetCopykatSegs(pred_mat_grpd[[chr]][, c(1:7, cell_num)], segs_brks_dist)
#                 }
# }
# }
# )
# 
# 
pred_mat_grpd[[12]] %>% pull(8)
pred_mat_grpd[[12]][,1:8] %>% .[.[,8] != "",]
pred_mat_grpd[[12]] %>% pull(8) %>% grep ("amp", .) 

pred_mat_grpd[[12]][,8] %>% grep ("amp", .) %>% length() <= 0 &
    pred_mat_grpd[[12]][,8] %>% grep ("del", .) %>% length() <= 0


system.time(copykat_cnv_segs <- foreach(chr = 1:length(pred_mat_grpd), .combine = "rbind") %:% 

            foreach(cell_num = 8:ncol(pred_mat_grpd[[chr]]), .combine = "rbind") %dopar% {


                if(pred_mat_grpd[[chr]] %>% pull(cell_num) %>% grep ("amp", .) %>% length() <= 0 |
                   pred_mat_grpd[[chr]] %>% pull(cell_num) %>% grep ("del", .) %>% length() <= 0){
                } else{

                    GetCopykatSegs(pred_mat_grpd[[chr]] %>% as.data.frame(.) %>% .[, c(1:7, cell_num)], segs_brks_dist)
                }

})


# THIS WORKED FINALLY
# THIS
# THIS
# THIS
# THIS
system.time(copykat_cnv_segs <- foreach(cell_num = 8:ncol(sub_raw_loc), .combine = "rbind") %dopar% {
    if(sub_raw_loc[,cell_num] %>% grep ("amp", .) %>% length() > 0 |
       sub_raw_loc[,cell_num] %>% grep ("del", .) %>% length() > 0){
        GetCopykatSegs(sub_raw_loc[, c(1:7, cell_num)], segs_brks_dist)
    #         GetCopykatSegs(pred_mat_grpd[[12]][, c(1:7, cell_num)], segs_brks_dist)

    }

})

GetCopykatSegs(sub_raw_loc[, c(1:8)], segs_brks_dist)

copykat_cnv_segs %>% head()
copykat_cnv_segs %>% nrow()
copykat_cnv_segs$seqnames %>% unique()
copykat_cnv_segs[copykat_cnv_segs$seqnames == "chr2",] %>% head()
copykat_cnv_segs[copykat_cnv_segs$seqnames == "chr12",] %>% head()

one_cell_cnv <- copykat_cnv_segs[copykat_cnv_segs$cell_name == "tnbc1_AAACCTGCACCTTGTC",]
one_cell_cnv
one_cell_name <- one_cell_cnv$cell_name[2]
one_cell_name
one_cell_cnv[1:5,1:5]
one_cell_cnv$seqnames %>% 
    unique()

one_cell_cnv[one_cell_cnv$seqnames == "chr2",]
parallel::stopCluster(cluster)
sub_raw_loc[,1:8] %>% filter(chromosome_name == 12) %>% .[,8] %>% grep("amp",.)



source("../../digital_karyotype/R/utils.R")
source("../../digital_karyotype/R/digital_karyotype.R")
cnv_colors <- c('red', 'blue')
names(cnv_colors) <- c('amp', 'del')
tmp_plot <- "../plots/tmp/tmp_plot.pdf"

dig_kar_plt <- plot_digital_karyotype(type = "",
                                      layers_h2 = one_cell_cnv,
                                      fill_arg = "cnv_state",
                                      colors_arg_h2 = cnv_colors,
                                      legend_title_arg = "Copykat Calls",
                                      sub_title_arg = str_glue("Copykat CNV Calls: {one_cell_name}"),
                                      plot_both_haplotypes = FALSE)

ggsave(tmp_plot)


all_cells <- unique(copykat_cnv_segs$cell_name)
length(all_cells)
all_cells[1]
save_path <- "copykat_digital_karyotypes/cna_refine"
foreach(cell = 1:length(all_cells)) %dopar% {

    dig_kar_plt <- plot_digital_karyotype(type = "",
                                          layers_h2 = copykat_cnv_segs[copykat_cnv_segs$cell_name == all_cells[cell],],
                                          fill_arg = "cnv_state",
                                          colors_arg_h2 = cnv_colors,
                                          legend_title_arg = "Copykat Calls",
                                          sub_title_arg = str_glue("Copykat CNV Calls: {all_cells[cell]}"),
                                          plot_both_haplotypes = FALSE)

    saveDigitalKaryotype(dig_kar_plt, save_path, all_cells[cell])
}


