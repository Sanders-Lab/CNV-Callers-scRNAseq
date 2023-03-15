library(tidyverse)
library(foreach)
library(parallel)


segs_df <- read.table('../proc/infercnv_preds_mat_pbmc3k_tnbc1.tsv.gz',
                      header = T)
tumour_cells <- segs_df %>% select(gene,
                                   seqnames,
                                   start,
                                   end,
                                   contains('tnbc'))





source('../../digital_karyotype/R/add_info_plt_layers.R')
prj_id <- 'infercnv_digital_karyotypes'
tmp_plot <- '../plots/tmp/tmp_plot.pdf'
cna_colors <- c('red', 'blue', 'green')
names(cna_colors) <- c('amp', 'del', 'loh')


plotDigKar <- function(cell_mat, cell_name){
    
    cna_dk <- plotSkeletonH2(ideo_df)

    cna_dk <- addPltLayers(cna_dk,
                           cell_mat,
                           2,
                           cell_name,
                           cna_colors,
                           'InferCNV CNV Calls',
                           title = str_glue('InferCNV CNV Calls - {cell_name}'))

    saveDigitalKaryotype(cna_dk,
                         prj_id,
                         cell_name)
}


n_cores <- detectCores() - 4
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


start <- Sys.time()
foreach(c = 5:ncol(tumour_cells)) %dopar% {
    
    curr_cell_mat <- tumour_cells %>% select(seqnames,
                                             start,
                                             end,
                                             c)
    curr_cell_mat <- curr_cell_mat[!(curr_cell_mat[,4] == ''),]

    cell_name <- curr_cell_mat %>% colnames() %>% .[4]



    plotDigKar(curr_cell_mat, cell_name)

}
end <- Sys.time()
end - start
parallel::stopCluster(cluster)


