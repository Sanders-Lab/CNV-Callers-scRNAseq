source("./run_copykat.R")

cmd_args[1] <- "../proc/ct_segs_ref_comb_filt_counts.tsv.gz"
cmd_args[2] <- as.integer(64)
cmd_args[3] <- "../outputs/tall_scnova/copykat_pbmc3k_regenotype_chr6"
cmd_args[4] <- "pbmc_"

cmd_args

ret_list <- preProcData(cmd_args)


count_mat <- ret_list[[1]]
n_cores <- as.integer(ret_list[[2]])
out_dir <- ret_list[[3]]
ref_cells_names <- ret_list[[4]]

colnames(count_mat) <- gsub(".1", "-1", colnames(count_mat))
ref_cells_names <- gsub(".1", "-1", ref_cells_names)
count_mat[1:5,1:5]
ref_cells_names[1:5]


getwd()
setwd("../../../R")
runCopykat(count_mat, n_cores, out_dir, ref_cells_names)



