# Dependencies:
# - beanplot
# - mixtools
# - pheatmap
# - zoo
# - squash
# - biomaRt
# for issues with the stdc++ lib not found error: 
# install r-rsqlite from conda-forge
# install grimbough/biomaRt from biocManager

# reproduced from the tutorial:
# SmartSeq2 scRNA seq of Oligodendroglioma
# Paper in ../docs/
# The github repo is here -> https://github.com/diazlab/CONICS 




# .libPaths("/fast/work/users/sbanerj_m/miniconda/envs/baseR/lib/R/library/")

# BiocManager::install('grimbough/biomaRt')

# BiocManager::install("biomaRt")

# BiocManager::install("scran")

devtools::install_github("diazlab/CONICS/CONICSmat", dep = FALSE)



# ---- Begin analysis ----


library(CONICSmat)
library(tidyverse)
library(biomaRt)


# tnbc1_counts = as.matrix(read.table(paste("data/GSE70630_OG_processed_data_v2.txt.gz", sep = ""), 
#                                  sep="\t",
#                                  header=T,
#                                  row.names=1,
#                                  check.names=F))
# tnbc1_counts [which(is.na(tnbc1_counts ))]=0


# dim(tnbc1_counts)


# tnbc1_counts[1:5,1:5]


# patients=unlist(strsplit(colnames(tnbc1_counts),
#                          "_",
#                          fixed=TRUE))[seq(1,(3*ncol(tnbc1_counts))-1,3)]
# unique(patients)


# patients[which(patients=="93")]="MGH93"
# patients[which(patients=="97")]="MGH97"




## ---- TNBC1 ---- ##
tnbc1_counts <- read.table('../proc/pbmc3k_epithelium_tnbc1_ref_merged_counts.tsv.gz',
                           header = T,
                           row.names = 1)
tnbc1_counts[1:5,1:5]
tnbc1_counts %>% head(n = c(5,5))
tnbc1_counts %>% nrow
tnbc1_counts %>% ncol


patients <- tnbc1_counts %>% colnames() %>% grep('tnbc', .)
patients <- tnbc1_counts %>% colnames()
patients %>% length()


# log transforming the counts mat
# might not be required!
# norm_mat <- log2((tnbc1_counts / 10) + 1)
# norm_mat[1:5,]


# This matrix contains the chromosomal coordinates of all chromosome arms 
# on autosomes. If no information on chromosomal alterations based on DNA 
# sequencing (e.g. exome-seq) is available, this file can be used to test for
# chromosome arm-scale copy number alterations. If exome-seq is available,
# this file should be replaced with a file holding large-scale CNVs.

regions=read.table("../conicsmat_test/data/chromosome_arm_positions_grch38.txt",
                   sep="\t",
                   row.names = 1, 
                   header = T)
head(regions,n=5)



# getting the gene positions from ensembl
# gene_pos <- getGenePositions(rownames(tnbc1_counts))
# gene_pos %>% head(n = 5)


# trying with the local
# trying manually
hmart <- useMart(biomart =  "ENSEMBL_MART_ENSEMBL",
                 dataset = "hsapiens_gene_ensembl")


gene_pos <- biomaRt::getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),
                           filters ='hgnc_symbol',
                           values = rownames(tnbc1_counts),
                           mart = hmart)
gene_pos %>% nrow


# gene_pos=getGenePositions(gene_names = rownames(tnbc1_counts),
#                           ensembl_version = "apr2019.archive.ensembl.org")
# gene_pos=getGenePositions(gene_names = rownames(tnbc1_counts))
# gene_pos <- read.table('../infercnv_test/input_files/hg38_gencode_v27.txt')
# gene_pos[1:5,1:4]
# gene_pos <- gene_pos %>% rename(hgnc_symbol = V1,
#                                 chr = V2,
#                                 start = V3,
#                                 end = V4)


# Subsequently, we can filter uniformative genes. 
# These are genes which are expressed in only 
# very few cells (default >5 cells).
# tnbc1_counts=filterMatrix(tnbc1_counts,gene_pos[,"hgnc_symbol"],minCells=5)
tnbc1_counts <- filterMatrix(tnbc1_counts,gene_pos[,"hgnc_symbol"],minCells=5)
tnbc1_counts %>% nrow


# Now, we calculate a normalization factor for each cell. 
# Because the average gene expression in each cell depends on
# the number of expressed genes 
# (the more genes expressed in one cell, 
# the less reads are "available" per gene), the normalization factor 
# centers the gene expression in each cell around the mean.
normFactor=calcNormFactors(tnbc1_counts)


# defining normal an tumour cells
# all_cells <- colnames(tnbc1_counts)
# normal_cells <- colnames(tnbc1_counts[str_detect(all_cells, 'pbmc3k')])
# normal_cells %>% length
# normal_cells %>% head
# tumour_cells <- colnames(tnbc1_counts[!str_detect(all_cells, 'pbmc3k')])
# tumour_cells %>% length
# tumour_cells %>% head

all_cells <- colnames(tnbc1_counts)
normal_cells <- c(1:3331)
names(normal_cells) <- colnames(tnbc1_counts[str_detect(all_cells, 'pbmc3k') | str_detect(all_cells, 'SRR')])
normal_cells %>% head()
normal_cells %>% tail()
normal_cells %>% length()


tumour_cells <- c(3332:4428)
names(tumour_cells) <- colnames(tnbc1_counts[str_detect(all_cells, 'tnbc')])
tumour_cells %>% head()
tumour_cells %>% length()
tumour_cells %>% tail()


# all_cells <- c(normal_cells, tumour_cells)
all_cells[1:5]


# CHECK the numbers before doing
cell_prefix <- c(rep('pbmc3k',2700), rep('tnbc1', 1097))
cell_prefix %>% head()
cell_prefix %>% length()
ncol(tnbc1_counts)


# The next step is to determine if the average gene expression 
# any of the regions show a bimodal distribution across cells.
# First, the average expression in each cell is centered using
# the previously calculated normalization factor. 
# Then, the z-score of the centered gene expression across all
# cells is calculated. Based on these z-scores, a Gaussian mixture model
# is calculated with the mixtools package. 
# Important: We only calculate results for regions harboring more 
# than 100 expressed genes (as defined by the initial filtering step) 
# to make sure the predictions are not influenced by a few deferentially 
# expressed genes in a small region.
# l=plotAll(tnbc1_counts,normFactor,regions,gene_pos,"../outputs/conicsmat_pbmc3k_tnbc1/tnbc1_cnv") w

# plotting cnv with the knowledge of normal and tumour cells
# l <- plotAll(tnbc1_counts,
#                  normFactor,
#                  regions,
#                  gene_pos,
#                  "../outputs/conicsmat_pbmc3k_tnbc1/tnbc1_cnv")




# hi=plotHistogram(l,tnbc1_counts,clusters=2,zscoreThreshold=4)

pdf("../outputs/conicsmat_pbmc3k_tnbc1/tnbc1_hist.pdf")
hi <- plotHistogram(l,
                    tnbc1_counts,
                    clusters=2,
                    zscoreThreshold=4,
                    cell_prefix)
dev.off()



# To obtain a final assignmant as malignant or non-malignant cells, 
# we first filter uninformative, noisy regions based on the
# results of the likelihood ratio test and the BIC for each region.
# lrbic=read.table("../outputs/conicsmat_pbmc3k_tnbc1/tnbc1_cnv_BIC_LR.txt",
#                  sep="\t",
#                  header=T,
#                  row.names=1,
#                  check.names=F)
# colnames(lrbic)
# lrbic %>% head()
# candRegions=rownames(lrbic)[which(lrbic[,"BIC difference"]>200 & 
#                                     lrbic[,"LRT adj. p-val"]<0.01)]
# candRegions %>% head()
# l[1:5,candRegions]


# For these regions we again generate a heatmap of posterior probabilities.
# hi=plotHistogram(l[,candRegions],
#                  tnbc1_counts,
#                  clusters=4,
#                  zscoreThreshold=4,
#                  patients)
# 
# pdf("plots/posterior_prob_final_assignment.pdf")
# plotHistogram(l[,candRegions],
#                  tnbc1_counts,
#                  clusters=4,
#                  zscoreThreshold=4,
#                  patients)
# dev.off()




# Since these results support our previous assignments, 
# we can now assign a label as malignant or non-malignant to 
# each cell based on the clusters returned by the heatmap function. 
# The cluster ID for each cell is returned by the plotHistogram () function
# and the leftmost cluster has ID 1
# (as given by the ID information on the bottom right).
# normal= which(hi==1)
# tumor=which(hi!=1)
# normal %>% length()




# Now, we can plot the posterior probabilities again, 
#but with statistics for normal and tumor cells:
# redu=plotAll(tnbc1_counts,
#              normFactor,
#              regions[candRegions,],
#              gene_pos,"plots/SUVA_CNVs_with_info.pdf",
#              normal=normal,
#              tumor=tumor)

# Now, we can plot the posterior probabilities again, 
#but with statistics for normal and tumor cells:
redu=plotAll(tnbc1_counts,
             normFactor,
             regions,
             gene_pos,
             "../outputs/conicsmat_pbmc3k_tnbc1/supervised_postProb",
             normal=normal_cells,
             tumor=tumour_cells)
redu %>% head()
regions %>% head()



# By thresholding on the posterior probabilities we can
# next generate a binary matrix, where 1 indicates the 
# presence of a CNV and 0 the absence. Based on the average
# expression in the normal cells, we know if an alteration is
# either a copy number gain or a loss.
bin_mat=binarizeMatrix(redu,
                       normal = normal_cells,
                       tumor = tumour_cells,
                       0.8)
bin_mat %>% head()
bin_mat %>% tail()
bin_mat %>% dim()
bin_mat %>% is.na() %>% table()
bin_mat %>% is.infinite() %>% table()
bin_mat %>% is.nan() %>% table()
bin_mat[is.na(bin_mat)] <- 0



## ---- DIGKAR PLOTTING ---- ##

# preprocessing
bin_mat_tr <- bin_mat %>% t()
bin_mat_tr %>% dim()
bin_mat_tr %>% head(n = c(5,5))


# subsetting only for tumour cells
bin_mat_tr <- as.data.frame(bin_mat_tr)
tumour_bin_mat_tr <- bin_mat_tr %>% dplyr::select(contains('tnbc1'))


# changing 1 to either del or amp acc to rownames
# changing 0 to blank
# also adding glocs
tumour_label_mat_tr <- data.frame()
for(r in 1:nrow(tumour_bin_mat_tr)){
    
    r_name_split <- rownames(tumour_bin_mat_tr)[r] %>% str_split(., '_')
    cna_state_row <-  r_name_split[[1]][1]
    chr_arm <- r_name_split[[1]][2]


    curr_row <- tumour_bin_mat_tr[r,]
    curr_row[curr_row == 1] <- cna_state_row
    curr_row[curr_row == 0] <- ''
    curr_row <- as.data.frame(curr_row)
    curr_row <- cbind(as.data.frame(regions[chr_arm,]), curr_row)

    tumour_label_mat_tr <- rbind(tumour_label_mat_tr, curr_row)
}
tumour_label_mat_tr$Length <- NULL
tumour_label_mat_tr %>% head(n = c(5,5))


# save to file
write.table(tumour_label_mat_tr,
            gzfile('../proc/conicsmat_tumour_cna_labels.tsv.gz'),
            sep = '\t',
            row.names = T)

write.table(label_mat,
            gzfile('../proc/conicsmat_tumour_cna_labels.tsv.gz'),
            sep = '\t',
            row.names = T)

## digkar plotting begins here
label_mat <- read.table('../proc/conicsmat_tumour_cna_labels.tsv.gz',
                        header = T)
# label_mat <- label_mat %>% dplyr::rename(seqnames = 'Chrom',
#                                   start = 'Start',
#                                   end = 'End')
label_mat %>% head(n = c(5,5))
label_mat$seqnames <- paste0('chr', label_mat$seqnames)



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


## ---- DIGKAR PLOTTING ---- ##

# replacing





pdf('../outputs/conicsmat_pbmc3k_tnbc1/bin_mat.pdf')
bin_mat_plot <- plotBinaryMat(bin_mat,
                              normal = normal_cells,
                              tumor = tumour_cells,
                              patients = patients)
dev.off()
ggsave('../outputs/conicsmat_pbmc3k_tnbc1/bin_mat.pdf',
       device = 'pdf',
       plot = bin_mat_plot)




# To detect individual chr. Not mandatory.
# Alternative is to use the plotAllChromosomes function
# detectBreakPoints (tnbc1_counts,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=1,patients=patients,patient="MGH36",breakpoints=regions)
# detectBreakPoints (tnbc1_counts,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=4,patients=patients,patient="MGH36",breakpoints=regions)
# detectBreakPoints (tnbc1_counts,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=7,patients=patients,patient="MGH36",breakpoints=regions)
# detectBreakPoints (tnbc1_counts,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=8,patients=patients,patient="MGH36",breakpoints=regions)
# detectBreakPoints (tnbc1_counts,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=10,patients=patients,patient="MGH36",breakpoints=regions)
# detectBreakPoints (tnbc1_counts,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=12,patients=patients,patient="MGH36",breakpoints=regions)




plotAllChromosomes (mat = tnbc1_counts,
                    normal = normal_cells,
                    tumor = tumour_cells,
                    windowsize = 101,
                    gene_pos = gene_pos,
                    patients = cell_clusters,
                    breakpoints = regions)


tnbc1_counts %>% nrow()
tnbc1_counts %>% ncol()
# tnbc one here
plotAllChromosomes(mat = tnbc1_counts,
                    normal = normal_cells,
                    tumor = tumour_cells,
                    windowsize = 101,
                    gene_pos = gene_pos,
                    fname = "tnbc1",
                    patients = cell_prefix,
                    patient = 'tnbc1',
                    breakpoints = regions)



# To obtain a matrix of corrected p-values for each cell and each CNV, we identify the component representing the cell population without CNVs for each region. Subsequently, we utilize the estimated parameters (mean and standard deviation) of this component to calculate a p-value for each of the tumor cells based on its z-scored expression with the pnorm() function.

# The matrix r holds all adjusted p-values (Benjamini-Hochberg) for each cell and each CNV candidate region. The matrix binr is a binarized version of r, where the presence of a CNV is thresholded on an adjusted p-val<0.1.
r=generatePvalMat(tnbc1_counts,
                  regions[candRegions,],
                  normFactor,
                  normal,
                  tumor,
                  gene_pos,
                  threshold=0.8)
binr=ifelse(r>0.1,0,1)

pdf("plots/adj_pValues_cnv_presence.pdf")
boxplot(r)
dev.off()


# Last, we can generate a visualization of chromosomal alterations in each 
# single cell across the genome for each patient. In the underlying function we compute the average expression in normal cells in a smoothing window (again, genes are sorted by genomic position). To reduce noise we filter for genes that are robustly measured in normal and tumor cells. In each cell, we center the expression ratio to this normal control by the mean, assuming that less than 50% of loci will be altered by copy number.

# tnbc here
# this can be directly done after 
pdf('../outputs/conicsmat_pbmc3k_tnbc1/all_chr_tnbc1.pdf')
plotChromosomeHeatmap(tnbc1_counts,
                      normal = normal_cells,
                      plotcells = tumour_cells,
                      gene_pos = gene_pos,
                      chr = T,
                      windowsize = 101,
                      expThresh=0.2,
                      thresh = 1)
dev.off()


# plotting the ref cells too
pdf('../outputs/conicsmat_pbmc3k_tnbc1/all_chr_all_cells.pdf')
plotChromosomeHeatmap(tnbc1_counts,
                      normal = normal_cells,
                      plotcells = all_cells,
                      gene_pos = gene_pos,
                      chr = T,
                      windowsize = 101,
                      expThresh=0.2,
                      thresh = 1)
dev.off()
