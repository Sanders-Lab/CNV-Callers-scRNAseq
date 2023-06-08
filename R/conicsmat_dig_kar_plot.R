library(tidyverse)

label_mat <- read.table('../proc/conicsmat_tumour_cna_labels.tsv.gz',
                        header = T)




source('../../digital_karyotype/R/add_info_plt_layers.R')
prj_id <- 'conicsmat_pbmc3k_epithelium_tnbc1'
cnv_colors <- c('red', 'blue')
names(cnv_colors) <- c('amp', 'del')


plotCnvDigKar <- function(mat, cell_name, fill_param){
    cnv_dk <- plotSkeletonH2(ideo_df)

    cnv_dk <- addPltLayers(dig_kar_plt = cnv_dk,
                           info_df = mat,
                           2,
                           fill_param = fill_param,
                           colors_param = cnv_colors,
                           leg_titl = 'ConicSmat CNV Calls',
                           title = str_glue('ConicSmat CNV Calls - {cell_name}'))

    saveDigitalKaryotype(cnv_dk, prj_id, cell_name)
}


# loop through all the cells
library(foreach)
library(parallel)


ncores <- detectCores() - 4
ncores
cluster <- parallel::makeForkCluster(ncores)
doParallel::registerDoParallel(cluster)


start <- Sys.time()
foreach(c = 4:ncol(label_mat)) %dopar% {
    curr_cell <- label_mat %>% dplyr::select(1,
                                             2,
                                             3,
                                             c)

    curr_cell <- curr_cell[!curr_cell[,4] == '',]

    cell_name <- colnames(curr_cell)[4]

    plotCnvDigKar(curr_cell, cell_name, cell_name)
}
end <- Sys.time()
end - start
parallel::stopCluster(cluster)
