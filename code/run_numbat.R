library(tidyverse)
library(numbat)

counts <- read.table('data//GSM4476492_combined_UMIcount_CellTypes_ATC2.txt.gz')
counts %>% head()
counts %>% tail()

# first two rows had preds. Removed them
raw_counts <- counts[3:nrow(counts),]
raw_counts %>% head()
raw_counts <- as.data.frame(raw_counts)
raw_counts <- apply(raw_counts, 2, as.numeric)
raw_counts <- as.matrix(raw_counts, sparse = T)

write.table(raw_counts, 'data/atc2_raw_counts_cleaned.tsv', sep = '\t')

raw_counts <- read.table('data/atc2_raw_counts_cleaned.tsv')
counts_mat <- as.matrix(raw_counts)

count_mat_ATC2 = readRDS(url('http://pklab.med.harvard.edu/teng/data/count_mat_ATC2.rds'))
allele_df <- readRDS(url('http://pklab.med.harvard.edu/teng/data/df_allele_ATC2.rds'))
allele_df %>% head()

ref_custom <- read.table('numbat_ref_granje_kiwi_pbmc.tsv')
ref_custom %>% head()
ref_custom <- as.matrix(ref)




## ---- tnbc1 ---- ##
# convert to matrix for each of the matrices
tnbc1_counts <- read.table('../proc/tnbc1_gene_mat_raw_cleaned.tsv.gz',
                           sep = '\t')
tnbc1_counts[1:5,1:5]

tnbc1_count_mat <- as.matrix(tnbc1_counts)
tnbc1_count_mat %>% dim()
tnbc1_count_mat %>% head()
colnames(tnbc1_count_mat) <- paste0(colnames(tnbc1_count_mat), '-1')


# as matrix is key here
ref_custom <- read.table('../proc/numbat_ref_pbmc_epi.tsv.gz')
ref_custom %>% dim()
ref_custom %>% head
ref_custom %>% is.na() %>% table
ref_custom_mat <- as.matrix(ref_custom)


# there were 3 NA in the 5th col. 
# substituted them for 0.0
tnbc1_allele_df <- read.table('../data/tnbc1_ck/tnbc_fastq/ART13preTX_0_1_HFYMLDMXX/pileup_phase_output_tnbc1/tnbc1_allele_counts.tsv.gz',
                        sep = '\t',
                        header = T)
tnbc1_allele_df %>% head()
tnbc1_allele_df %>% is.na() %>% which(arr.ind = T)
tnbc1_allele_df[828054:828056,5] <- 0.0
tnbc1_allele_df[828054:828056,5] 


out <- run_numbat(tnbc1_count_mat,
                  ref_custom_mat,
                  tnbc1_allele_df,
                  genome = "hg38",
                  ncores = 16,
                  plot = TRUE,
                  out_dir = '../outputs/numbat_tnbc1_pbmcEpiRef/')


## ---- p1794 kiwi ---- ##

library(numbat)
library(Seurat)
library(stringr)

sample <- 'P1794_GEX_HD1'

feat_mat_dir <- str_glue("/fast/groups/ag_sanders/work/projects/suharto/cnv_callers/numbat_test/p1794_kiwi/count_output/{sample}/outs/filtered_feature_bc_matrix")
allele_count_dir <- str_glue("/fast/groups/ag_sanders/work/projects/suharto/cnv_callers/numbat_test/p1794_kiwi/pileup_phase/{sample}/{sample}_allele_counts.tsv.gz")
nout_dir <- str_glue("outputs/numbat_{sample}")


feat_mat <- Read10X(feat_mat_dir)
allele_count <- read.table(allele_count_dir, sep = '\t', header=T)
feat_mat %>% head()
allele_count %>% head()



out <- run_numbat(feat_mat,
                  ref_custom,
                  allele_count,
                  genome = "hg38",
                  ncores = 16,
                  plot = TRUE,
                  out_dir = nout_dir)
