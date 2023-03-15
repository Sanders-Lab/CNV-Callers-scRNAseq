library(tidyverse)
library(stringr)
library(data.table)


raw_res <- read.table('../outputs/copykat_pbmc3k_epithelium_tnbc1/tnbc1_copykat_CNA_raw_results_gene_by_cell.txt',
                      header = T)
raw_res[1:5,1:10]


cna_res <- read.table('../outputs/copykat_pbmc3k_epithelium_tnbc1/tnbc1_copykat_CNA_results.txt',
                      header = T)
cna_res[1:5,1:10]



copykat_preds <- read.table('../outputs/copykat_pbmc3k_epithelium_tnbc1/tnbc1_copykat_prediction.txt',
                            header = T)



# change vals to amp or del
# if val is > 0.05 or < -0.05

table(ck_preds$copykat.pred == 'aneuploid')
aneu_cells <- ck_preds[ck_preds$copykat.pred == 'aneuploid',]
mod_gene_mat <- raw_res[,aneu_cells$cell.names]
cna_gene_mat <- mod_gene_mat %>% apply(., 2, function(x) ifelse(x > 0.05, 'amp', ifelse(x < -0.05, 'del', '')))
mod_gene_mat[1:5,2900:2909]
cna_gene_mat[1:5,1:5]


gene_loc <- raw_res %>% select(6,
                               2,
                               3,
                               4)
gene_loc %>% head()
cna_gloc <- cbind(gene_loc, cna_gene_mat)
cna_gloc[1:5,1:10]


write.table(cna_gloc,
            gzfile('../proc/copykat_aneu_labels.tsv.gz'),
            sep = '\t',
            row.names = T)



## ---- Plotting dig kar ---- ##
source('../../digital_karyotype/R/add_info_plt_layers.R')
prj_id <- 'copykat_digital_karyotypes'
tmp_plot <- '../plots/tmp/tmp_plot.pdf'
ck_colors <- c('red', 'blue')
names(ck_colors) <- c('amp', 'del')


pred_ck_genes_loc <- read.table('../proc/copykat_aneu_labels.tsv.gz',
                                header = T)
pred_ck_genes_loc[1:5,1:10]



# one time tidy 
# changing the col names for easy plotting on digKar
# adding chr before chromosome_name
# changing chr23 to chrX
pred_ck_genes_loc <- pred_ck_genes_loc %>% dplyr::rename(gene = 'hgnc_symbol',
                                                         seqnames = 'chromosome_name',
                                                         start = 'start_position',
                                                         end = 'end_position')

pred_ck_genes_loc$seqnames <- paste0('chr', pred_ck_genes_loc$seqnames)
tmp <- 1
pred_ck_genes_loc$seqnames[pred_ck_genes_loc$seqnames == 'chr23'] <- 'chrX'
unique(pred_ck_genes_loc$seqnames)


# saving the table with the updated col names
write.table(pred_ck_genes_loc,
            gzfile('../proc/copykat_aneu_labels.tsv.gz'),
            sep = '\t',
            row.names = T)


# making a grangeslist of length = ncells
# - 4 to remove the gro ncols
ck_preds_cell_grl <- vector(mode = 'list',
                            length = (ncol(pred_ck_genes_loc) - 4))
cell_num <- 1
for(cell_idx in 5:ncol(pred_ck_genes_loc)){



    # getting each cell
    cell_preds <- pred_ck_genes_loc %>% select(gene,
                                               seqnames,
                                               start,
                                               end,
                                               cell_idx)

    # get the cell name from the 5th col
    cell_name <- colnames(pred_ck_genes_loc[,cell_idx, drop = F]) 


    # subset to only CNA cols
    cna_idx <- which(cell_preds[,5] == 'amp' | cell_preds[, 5] == 'del')
    cell_preds_gloc <- cell_preds[cna_idx,]



    # making the digital karyotype
    ck_cnv_dk <- plotSkeletonH2(ideo_df)

    # add layers only if cna are there
    if(nrow(cell_preds_gloc) > 0){

        # Add the kar info per cell
        ck_cnv_dk <- addPltLayers(dig_kar_plt = ck_cnv_dk,
                                  info_df = cell_preds_gloc,
                                  2,
                                  fill_param = cell_name,
                                  colors_param = ck_colors,
                                  leg_titl = 'Copykat CNV Calls',
                                  title = str_glue('Digital Karyotype - {cell_name}'))

        # making gro only if cna present
        cell_gro <- makeGRangesFromDataFrame(cell_preds_gloc,
                                             keep.extra.columns = T)
        ck_preds_cell_grl[[cell_num]] <- cell_gro
    }

    # saving the plot
    saveDigitalKaryotype(ck_cnv_dk, prj_id, cell_name)

    cell_num <- cell_num + 1 
}

# save the granges list obj
saveRDS(ck_preds_cell_grl,
        file = '../proc/ck_cna_cells_granges_list.rds')



## ---- tidying up cna preds ---- ##

cna_res <- read.table('../outputs/copykat_pbmc3k_epithelium_tnbc1/tnbc1_copykat_CNA_results.txt',
                      header = T)

copykat_preds <- read.table('../outputs/copykat_pbmc3k_epithelium_tnbc1/tnbc1_copykat_prediction.txt',
                            header = T)

cna_res[1:5,1:10]
cna_res[1:5,2905:2912]

# checking if chrompos is the same as abspos
cna_res[cna_res$chrompos != cna_res$abspos,]  %>% head(5,5)
cna_res[cna_res$chrompos != cna_res$abspos,][1:5,1:5]
# chrom pos is the pos relative to each chrom
# abs pos is for the copykat panel plot


# checking the spread of the max and min vals
# for aneuploid regions
ck_preds <- copykat_preds[copykat_preds$copykat.pred != 'not.defined',]
alt_hyp <- vector(mode = 'integer', length = 796)
max_aneu <- vector(mode = 'integer', length = 796)
min_aneu <- vector(mode = 'integer', length = 796)


# checking the number of pred values of gene exp
# to set a cutoff for amp / del 
cell_num <- 1
for (i in 1:nrow(ck_preds)){

    if(ck_preds$copykat.pred[i] == 'aneuploid'){
        min_aneu[cell_num] <- min(cna_res[, ck_preds$cell.names[i]])
        max_aneu[cell_num] <- max(cna_res[, ck_preds$cell.names[i]])
        
        if(min(cna_res[, ck_preds$cell.names[i]]) > -0.05){
            alt_hyp[cell_num] <- i
        }
        if(max(cna_res[, ck_preds$cell.names[i]]) < 0.05){
            alt_hyp[cell_num] <- i
        }
        cell_num <- cell_num + 1
    }
}
which(alt_hyp != NULL)
alt_hyp[1:50]
unique(alt_hyp)
min_aneu[1:5]
min(max_aneu)
max(min_aneu)




# ALT HYP STANDS!!!


# change vals to amp or del
# if val is > 0.05 or < -0.05
table(ck_preds$copykat.pred == 'aneuploid')
aneu_cells <- ck_preds[ck_preds$copykat.pred == 'aneuploid',]
mod_gene_mat <- cna_res[,aneu_cells$cell.names]
cna_gene_mat <- mod_gene_mat %>% apply(., 2, function(x) ifelse(x > 0.05, 'amp', ifelse(x < -0.05, 'del', '')))
mod_gene_mat[1:5,1:10]
cna_gene_mat[1:5,1:10]


# taking the chrompos
cna_res[1:5,1:5]
gene_loc <- cna_res %>% select(1,
                               2)



id <- tnbc_cna_segs[,3] %>% consecutive_id()
r <- rle(id)

cumsum(r$lengths)



## ---- inspired by chat GPT ---- ##
# rle
# THIS WORKED!!!
tnbc_cna_segs <- cbind(cna_res[,c(1,2)], mod_gene_mat)
tnbc_cna_segs[1:4,1:10]


cna_segs_df <- data.frame(seqnames = character(),
                          start = integer(),
                          end = integer(),
                          mod_expr_value = double(),
                          cell_name = character())

ck_segs_grl <- vector(mode = 'list', length = nrow(tnbc_cna_segs))
cell_num <- 1

for(cell_col in 3:ncol(tnbc_cna_segs)){

    cna_cell <- tnbc_cna_segs %>% select('chrom',
                                         'chrompos',
                                         cell_col)
    cell_name <- colnames(cna_cell)[3]

    # rle of chr
    rle_chr <- rle(cna_cell[,1])
    chr_end <- cumsum(rle_chr$lengths)
    
    cell_segs_df <- data.frame(seqnames = character(),
                               start = integer(),
                               end = integer(),
                               mod_expr_value = double(),
                               cell_name = character())


    for (chr_idx in 1:23){

        if(chr_idx == 1){
            rle_cna <- rle(cna_cell[1:chr_end[chr_idx], 3])
        } else{
            rle_cna <- rle(cna_cell[(chr_end[chr_idx - 1] + 1):chr_end[chr_idx], 3])
        }

        # find the start and end of each
        start_idx <- cumsum(c(1, rle_cna$lengths[-length(rle_cna$lengths)]))
        end_idx <- cumsum(rle_cna$lengths)

        # alloc the tmp_df
        tmp_df <- data.frame(seqnames = character(length(start_idx)),
                             start = integer(length(start_idx)),
                             end = integer(length(start_idx)),
                             mod_expr_value = double(length(start_idx)),
                             cell_name = character(length(start_idx)))


        for(i in 1:length(start_idx)){
            if(chr_idx == 1){

                tmp_df$start[i] <- cna_cell$chrompos[start_idx[i]]
            } else{

                tmp_df$start[i] <- cna_cell$chrompos[(chr_end[(chr_idx - 1)] + start_idx[i])]
            }
        }

        for(j in 1:length(end_idx)){
            if(chr_idx == 1){
    
                tmp_df$end[j] <- cna_cell$chrompos[end_idx[j]]
            } else{
    
                tmp_df$end[j] <- cna_cell$chrompos[(chr_end[(chr_idx - 1)] + end_idx[j])]
            }
        }

        tmp_df$seqnames <- rep(chr_idx, length(start_idx))
        tmp_df$mod_expr_value <- rle_cna$values
        tmp_df$cell_name <- rep(cell_name, length(start_idx))

        if(nrow(tmp_df) > 0){

            cell_segs_df <- rbind(cell_segs_df, tmp_df)
        }
    }

    # making the gro if the cell_segs_df is not empty
    # making the gro and appending to the grl
    # adding to cna_segs_df
    cna_segs_df <- rbind(cna_segs_df, cell_segs_df)
    cell_gro <- makeGRangesFromDataFrame(cell_segs_df,
                                         keep.extra.columns = T)
    ck_segs_grl[[cell_num]] <- cell_gro
    cell_num <- cell_num + 1

}

dim(cna_segs_df)
cell_segs_df[1:5,1:5]
cna_segs_df[1:5,1:5]
ck_segs_grl[[1]]
length(ck_segs_grl)
is.null(ck_segs_grl)
tmp_df[1:5,1:4]
cell_gro
rle_cna
unique(cna_segs_df[,5])











cna_segs_df[cna_segs_df$mod_expr_value]
cna_segs_label <- cna_segs_df$mod_expr_value %>% sapply(., function(x) ifelse(x > 0.05, 'amp', ifelse(x < -0.05, 'del', '')))
is.numeric(cna_segs_label)
cna_segs_label


cna_segs_df <- cbind(cna_segs_df, cna_segs_label)
cna_segs_df[1:5,1:6]


write.table(cna_segs_df,
            gzfile('../proc/copykat_aneu_cnaSegs_vals_labels.tsv.gz'),
            sep = '\t',
            row.names = T)


cna_segs <- cna_segs_df[cna_segs_df$cna_segs_label != "",]
cna_segs$seqnames <- paste0('chr', cna_segs$seqnames)
cna_segs %>% head()
dim(cna_segs)
cna_segs %>% filter(seqnames == 'chr23') %>% select(seqnames)
cna_segs$seqnames[cna_segs$seqnames == 'chr23'] <- 'chrX'
cna_segs %>% tail()



write.table(cna_segs,
            gzfile('../proc/copykat_aneu_cnaSegs_labels_OnlyCNA.tsv.gz'),
            sep = '\t',
            row.names = T)



## ---- Plotting the digKars ---- ##

source('../../digital_karyotype/R/add_info_plt_layers.R')
prj_id <- 'copykat_digital_karyotypes'
tmp_plot <- '../plots/tmp/tmp_plot.pdf'
cna_colors <- c('red', 'blue')
names(cna_colors) <- c('amp', 'del')




plotCkCNA <- function(x, prj_id, cell_name){

    cna_segs_label <- x$cna_segs_label

    # making the digital karyotype
    ck_cnv_dk <- plotSkeletonH2(ideo_df)

    # add layers only if cna are there
    if(nrow(x) > 0){

        # Add the kar info per cell
        ck_cnv_dk <- addPltLayers(dig_kar_plt = ck_cnv_dk,
                                  info_df = x,
                                  2,
                                  fill_param = 'cna_segs_label',
                                  colors_param = cna_colors,
                                  leg_titl = 'Copykat CNV Calls',
                                  title = str_glue('Digital Karyotype - {cell_name}'))

        saveDigitalKaryotype(ck_cnv_dk, prj_id, cell_name)
    }
}


cell_names <- cna_segs$cell_name %>% unique()
cell_names
cna_segs  %>% head()




# parallel computing for making digKar
library(foreach)
library(parallel)


# setting up parallel backend
n_cores <- detectCores() - 4
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


start <- Sys.time()
foreach(i = 1:length(cell_names)) %dopar% {
    cell_n <- cell_names[i]

    cell_cna_segs <- cna_segs %>% filter(cell_name == cell_n)
    cell_cna_segs %>% head()

    plotCkCNA(cell_cna_segs, prj_id, cell_name) 
}
end <- Sys.time()
end-start

parallel::stopCluster(cluster)
