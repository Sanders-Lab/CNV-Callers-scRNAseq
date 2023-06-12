library(tidyverse)
library(foreach)
library(parallel)


segs_df <- read.table("../proc/infercnv_segs_cnv_label_allcells.tsv.gz",
    header = T
)
segs_df %>% head()






source("../../digital_karyotype/R/add_info_plt_layers.R")
prj_id <- "infercnv_digital_karyotypes"
tmp_plot <- "../plots/tmp/tmp_plot.pdf"
cna_colors <- c("red", "blue", "green")
names(cna_colors) <- c("amp", "del", "loh")


PlotDigKar <- function(cell_mat) {
    cna_dk <- plotSkeletonH2(ideo_df)

    cna_dk <- addPltLayers(cna_dk,
        cell_mat,
        2,
        "cnv_state",
        cna_colors,
        "InferCNV CNV Calls",
        title = str_glue("InferCNV CNV Calls - {cell_mat$cell_name[1]}")
    )

    saveDigitalKaryotype(
        cna_dk,
        prj_id,
        cell_mat$cell_name[1]
    )
}

segs_df <- segs_df[!segs_df$cnv_state == "", ]
segs_df %>% head()


segs_grpd_df <- segs_df %>%
    group_by(cell_name) %>%
    group_split(.)
segs_grpd_df[[1]]


n_cores <- 63
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


start <- Sys.time()
foreach(cell_num = 1:length(segs_grpd_df)) %dopar% {
    PlotDigKar(segs_grpd_df[[cell_num]])
}
Sys.time() - start
parallel::stopCluster(cluster)
