source("./run_numbat.R")


read_list <- ReadData(count_mat_path = "../proc/ct_segs_ref_comb_filt_counts_numbat.tsv.gz",
                      count_mat_type = "m",
                      allele_count_path =  "../data/ega_data/tall_ega_data/pileup_phase_tall/tall_allele_counts.tsv.gz"
)

ref_mat <- read_list[[1]]
ref_mat_annot <- read_list[[2]]
count_mat <- read_list[[3]]
allele_count <- read_list[[4]]

# required to remove the prefix
# colnames(count_mat) <- gsub("pbmc_", "", colnames(count_mat))
colnames(count_mat) <- gsub("tall_", "", colnames(count_mat))
count_mat[1:5,1:5]



numbat_ref <- MakeNumbatRef(ref_mat, ref_mat_annot)

out_dir <- "../outputs/tall_scnova/numbat_chr6_regenotype_tall"
n_cores <- 64

RunNumbat(
    count_mat,
    numbat_ref,
    allele_count,
    n_cores,
    out_dir,
    t = 1,
    min_llr = -1
)
