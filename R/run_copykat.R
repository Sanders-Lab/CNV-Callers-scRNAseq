library(copykat)
library(tidyverse)


#### ---- AIM ---- ####
# this script runs copykat cnv caller
# on the combined ref + exp scRNA-seq data
# combined because we want to use the ref baseline
# to call out cnvs
# @param 1. count mat: (mat) merged count mat
# @param 2. n_cores: (int) the number of cores to use
# @param 3. out_dir: (str) the output dir
# @param 4. ref_cell_prefix: (str) the cell prefix for normal cells
#######################





preProcData <- function(cmd_args) {
    # reading in the count mat
    # making a matrix and supplying it as it is good practice
    count_raw <- read.table(cmd_args[1])
    count_mat <- as.matrix(count_raw)


    # getting the ncores
    n_cores <- cmd_args[2]


    # get the out dir
    out_dir <- cmd_args[3]


    # get the cell ids of the normal cells
    # 1. select the cols which have pbmc3k || SRR since it is cell prefix
    # 2. get the colnames of this subsetted df
    # SRR is for the cells from epithelium dataset
    count_raw <- as.data.frame(count_raw)
    ref_cells_names <- count_raw %>%
        select(contains(cmd_args[4])) %>%
        colnames()


    return(list(count_mat, n_cores, out_dir, ref_cells_names))
}






runCopykat <- function(count_mat, n_cores, out_dir, ref_cells_names) {
    # create the out_dir if not present
    # setting the wd to the output dir
    # so that the results are saved to the out_dir
    if (!dir.exists(out_dir)) {
        dir.create(out_dir)
    }
    setwd(out_dir)


    # win.size is the window for smoothing the exp
    # 25 is the recommended value in the vignette
    # pasing the names of cells from the ref as the norm.cell.names
    copykat_run <- copykat(
        rawmat = count_mat,
        id.type = "S",
        ngene.chr = 5,
        win.size = 101,
        KS.cut = 0.1,
        sam.name = "tnbc1",
        distance = "euclidean",
        norm.cell.names = ref_cells_names,
        output.seg = "FLASE",
        plot.genes = "TRUE",
        genome = "hg20",
        n.cores = n_cores
    )
}





#### ---- MAIN SECTION ---- ####

# read cmd line args
cmd_args <- commandArgs(trailingOnly = T)


ret_list <- preProcData(cmd_args)
count_mat <- ret_list[[1]]
n_cores <- as.integer(ret_list[[2]])
out_dir <- ret_list[[3]]
ref_cells_names <- ret_list[[4]]

print(str_glue("\n Running Copykat with {n_cores} and the output_dir is {out_dir}\n"))

runCopykat(count_mat, n_cores, out_dir, ref_cells_names)

################################
