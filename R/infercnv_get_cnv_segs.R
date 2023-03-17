library(tidyverse)
library(infercnv)
library(parallel)
library(foreach)


step_20 <- readRDS('../outputs/infercnv_pbmc3k_epithelium_tnbc1/20_HMM_pred.repr_intensitiesHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj')
expr_data_df <- as.data.frame(step_20@expr.data)
tnbc_cells <- expr_data_df %>% select(contains('tnbc1'))
tnbc_cells %>% ncol()
tnbc_cells <- cbind(step_20@gene_order, tnbc_cells)
tnbc_cells <- tnbc_cells %>% dplyr::rename(seqnames = 'chr',
                                           end = 'stop')
tnbc_cells %>% head(n = c(5,5))
tnbc_cells[,4] %>% unique()



# QnD Tests
cell_tib <- tnbc_cells[,1:4] %>% group_by(seqnames) %>% group_split(.)
cell_tib
segs_df <- data.frame(seqnames = character(),
                      start = integer(),
                      end = integer(),
                      cell_name = character(),
                      cnv_state = character()) 
segs_df <- do.call(rbind, lapply(cell_tib, function(x) SegsChromCell(x)))
segs_df


cell_tib[[1]] %>% pull(4) %>% rle()
r <- tnbc_cells[,1:4] %>% group_by(seqnames) %>% group_map(~ rle(.$'tnbc1_AAACCTGCACCTTGTC'))




#### ---- PROD CODE ---- #### 
# proud of this! :D

# get the segs for each chr and the curr cell
SegsChromCell <- function(cell_chr_df){

    rle_chr_cell <- cell_chr_df %>% pull(4) %>% rle()


    tmp_df <- data.frame(seqnames = character(length(rle_chr_cell$lengths)),
                         start = integer(length(rle_chr_cell$lengths)),
                         end = integer(length(rle_chr_cell$lengths)),
                         cell_name = character(length(rle_chr_cell$lengths)),
                         cnv_state = character(length(rle_chr_cell$lengths)))


    tmp_df$start <- cumsum(c(1, rle_chr_cell$lengths[-length(rle_chr_cell$lengths)])) %>% 
        cell_chr_df$start[.]
    tmp_df$end <- cumsum(rle_chr_cell$lengths) %>% 
        cell_chr_df$end[.]
    tmp_df$seqnames <- rep(cell_chr_df$seqnames[1], length(tmp_df$start))
    tmp_df$cell_name <- rep(colnames(cell_chr_df)[4], length(tmp_df$start))
    tmp_df$cnv_state <- rle_chr_cell$values


    return(tmp_df)
}


segs_df <- data.frame(seqnames = character(),
                      start = integer(),
                      end = integer(),
                      cell_name = character(),
                      cnv_state = character())

ncores <- 63
cluster <- parallel::makeForkCluster(ncores)
doParallel::registerDoParallel(cluster)


start <- Sys.time()
segs_df <- foreach(cell_col = 4:ncol(tnbc_cells), .combine = 'rbind') %dopar% {
    
    rle_cell_chr <- tnbc_cells %>% select(1,
                                          2,
                                          3,
                                          cell_col) %>%
    group_by(seqnames) %>%
    group_split(.)

    segs_df <- do.call(rbind, lapply(rle_cell_chr, function(x) SegsChromCell(x)))


    return(segs_df)
}
Sys.time() - start
segs_df %>% head()
segs_df$cell_name %>% unique() %>% length()
parallel::stopCluster(cluster)


segs_label_df <- segs_df
segs_label_df$cnv_state <- segs_df$cnv_state %>% sapply(., function(x) ifelse(x > 1, 'amp', ifelse(x < 1, 'del', '')))
segs_label_df %>% head()

write.table(segs_label_df,
            gzfile("../proc/infercnv_segs_cnv_label_allcells.tsv.gz"),
            sep = "\t",
            col.names = T)
