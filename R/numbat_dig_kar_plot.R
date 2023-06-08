library(tidyverse)
library(foreach)
library(parallel)




#### ---- BRIEF ---- ####

# @param 1. joint_post_dir: (str) the dir of joint post tsv
# @param 2. prj_id: (str) the dir name where the dig kars will 
#           be saved.
# @param 3. n_cores: (integer) the n_cores.






ReadData <- function(joint_post_dir){

    joint_post <- read.table(joint_post_dir,
                             sep = '\t',
                             header = T)

    writeLines("Joint Post data read in!")

    return(joint_post)

}


ProcData <- function(joint_post){

    cna_segs <- joint_post %>% select(cell,
                                      CHROM,
                                      cnv_state,
                                      seg_start,
                                      seg_end)

    cna_segs <- cna_segs %>% dplyr::rename(seqnames = 'CHROM',
                                           start = 'seg_start',
                                           end = 'seg_end')
    cna_segs$seqnames <- paste0('chr', cna_segs$seqnames)

    writeLines("Data Processed successfully!")
    return(cna_segs)

}


PlotDigKar <- function(cell_mat, cell_name){
    
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




#### ---- MAIN SECTION ---- #### 

source('../../digital_karyotype/R/add_info_plt_layers.R')
cmd_args <- commandArgs(trailingOnly = T)

joint_post_dir <- cmd_args[1]
prj_id <- cmd_args[2]
n_cores <- as.integer(cmd_args[3])

writeLines(str_glue("The input args are:\n
           joint_post_dir: {joint_post_dir}\n
           prj_id: {prj_id}\n
           n_cores: {n_cores}\n"))

joint_post <- ReadData(joint_post_dir)
cna_segs <- ProcData(joint_post)


cna_colors <- c('red', 'blue', 'green')
names(cna_colors) <- c('amp', 'del', 'loh')

cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


writeLines("Starting to plot Digital Karyotypes...")
u_cells <- unique(cna_segs$cell)
start <- Sys.time()
foreach(c = 1:length(u_cells)) %dopar% {

    cell_name <- u_cells[c]
    curr_cell_mat <- cna_segs %>% filter(cell == cell_name)

    PlotDigKar(curr_cell_mat, cell_name)

}
Sys.time() - start
parallel::stopCluster(cluster)

