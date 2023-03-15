library(tidyverse)


joint_post <- read.table('../outputs/numbat_pbmc3k_epithelium_tnbc1/joint_post_2.tsv',
                         sep = '\t',
                         header = T)
joint_post %>% head()


segs_raw <- joint_post %>% select(cell,
                                  CHROM,
                                  cnv_state,
                                  seg_start,
                                  seg_end)

segs_raw <- segs_raw %>% dplyr::rename(seqnames = 'CHROM',
                                       start = 'seg_start',
                                       end = 'seg_end')
segs_raw$seqnames <- paste0('chr', segs_raw$seqnames)
segs_raw %>% head()

#### ---- PLOT DIGKARS ---- ####

segs_raw %>% head()

library(foreach)
library(parallel)

source('../../digital_karyotype/R/add_info_plt_layers.R')
prj_id <- 'numbat_digital_karyotypes'
tmp_plot <- '../plots/tmp/tmp_plot.pdf'
cna_colors <- c('red', 'blue', 'green')
names(cna_colors) <- c('amp', 'del', 'loh')


plotDigKar <- function(cell_mat, cell_name){
    
    cna_dk <- plotSkeletonH2(ideo_df)

    cna_dk <- addPltLayers(cna_dk,
                           cell_mat,
                           2,
                           'cnv_state',
                           cna_colors,
                           'Numbat CNV Calls',
                           title = str_glue('Numbat CNV Calls - {cell_name}'))

    saveDigitalKaryotype(cna_dk,
                         prj_id,
                         cell_name)
}


n_cores <- detectCores() - 4
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


u_cells <- unique(segs_raw$cell)
start <- Sys.time()
foreach(c = 1:length(u_cells)) %dopar% {
    
    cell_name <- u_cells[c]
    curr_cell_mat <- segs_raw %>% filter(cell == cell_name)

    plotDigKar(curr_cell_mat, cell_name)

}
end <- Sys.time()
end - start
parallel::stopCluster(cluster)
