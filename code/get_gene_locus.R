library(biomaRt)
library(tidyverse)


raw_gene <- read.table('./proc/pbmc3k_tnbc1_ref_merged_counts.tsv.gz')
raw_gene[1:5,1:2]

hgnc_gene_names <- read_tsv('./hgnc_gene_names.tsv')
hgnc_gene_names %>% nrow()
hgnc_gene_names %>% head()
hgnc_gene_names['HGNC ID'] %>% is.na() %>% table()
hgnc_id  <- hgnc_gene_names[['HGNC ID']]
hgnc_id %>% length



gene_list <- raw_gene %>% rownames()
gene_list %>% length()
gene_list %>% head()


# loading the mart for humans
h_mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")



gene_order <- getBM(attributes = c('hgnc_symbol','chromosome_name', 'start_position', 'end_position'),
                    filter = 'hgnc_symbol',
                    values = gene_list,
                    mart = h_mart,
                    useCache = F,
                    uniqueRows = T)
gene_order %>% nrow

gene_order$hgnc_symbol  %in% gene_list  %>% which(arr.ind = T)
filt_count <- raw_gene[(rownames(raw_gene) %in% gene_order$hgnc_symbol),]
rownames(raw_gene) %in% gene_order$hgnc_symbol %>% table
filt_count %>% head
filt_count %>% nrow

write.table(filt_count,
            gzfile('pbmc3k_tnbc1_filtered_counts.tsv.gz'))
gene_order %>% head
write.table(gene_order,
            gzfile('./proc/hgnc_sample_coords.tsv.gz'))


# get the details from the hgnc id
hgnc_order <- getBM(attributes = c('hgnc_id', 'hgnc_symbol','chromosome_name', 'start_position', 'end_position'),
                    filter = 'hgnc_id',
                    values = hgnc_id,
                    mart = h_mart,
                    useCache = F)
hgnc_order %>% nrow()
hgnc_order %>% head()

# removing entries with weird chr names like contigs and stuff
bad_chr_name_idx <- hgnc_order$chromosome_name %>% sapply(., function(x) str_detect(x,'_')) %>% which(arr.ind = T)
filt_hgnc_order <- hgnc_order[-bad_chr_name_idx,]
filt_hgnc_order %>%  nrow()
filt_hgnc_order %>% head()
filt_hgnc_order$hgnc_symbol %>% is.na() %>%  table

# writing to file with hgnc id for future use
write.table(filt_hgnc_order, 
          gzfile('./proc/hgnc_withID_gene_ranges.tsv.gz'),
          sep = '\t',
          row.names = F,
          col.names = T)


infercnv_gene_range <- filt_hgnc_order
infercnv_gene_range$hgnc_id = NULL
infercnv_gene_range$chromosome_name %>%  is.na() %>% table
infercnv_gene_range$start_position %>%  is.na() %>% table
infercnv_gene_range$end_position %>%  is.na() %>% table
# writing to file for infercnv
write.table(infercnv_gene_range, 
          gzfile('./proc/hgnc_infercnv_gene_ranges.tsv.gz'),
          sep = '\t',
          row.names = F,
          col.names = T)





# get the attributes required
gene_order <- getBM(attributes = c('hgnc_id', 'chromosome_name', 'start_position', 'end_position'),
                    filter = 'hgnc_id',
                    values = gene_list,
                    mart = h_mart,
                    useCache = F)


# get the ensembl gene ids
# this gave 24198 hits only
gene_names <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id'),
                    filter = 'hgnc_symbol',
                    values = gene_list,
                    mart = h_mart,
                    useCache = F)
gene_names %>% nrow()


# get the ensembl exon ids
# this gave all the genes
# but that leads to multiple exons per gene
exon_names <- getBM(attributes = c('hgnc_symbol','ensembl_exon_id'),
                    filter = 'hgnc_symbol',
                    values = gene_list,
                    mart = h_mart,
                    useCache = F)
exon_names %>% nrow()
exon_names %>% head()
gene_list %in% exon_names$hgnc_symbol %>% length()



# getting the exonchrom start and end and chromosome name
# using the exon id obtained from the previous step
exon_loci <- getBM(attributes = c('ensembl_exon_id', 'chromosome_name', 'exon_chrom_start', 'exon_chrom_end'),
                    filter = 'ensembl_exon_id',
                    values = exon_names$ensembl_exon_id,
                    mart = h_mart,
                    useCache = F)
exon_loci %>% nrow()
exon_loci %>% head()
exon_loci$ensembl_exon_id %in% exon_names$ensembl_exon_id %>% length()
exon_loci$ensembl_exon_id %>% duplicated() %>% which(arr.ind = T)
exon_loci[682430:682442,]



# making the exon_details df
# removing duplicate exon name col
exon_ranges <- rbind(exon_names, exon_loci)
exon_ranges[,1]
exon_ranges[,1] = NULL





# this worked for some
gene_order <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
                    filter = 'hgnc_symbol',
                    values = gene_list,
                    mart = h_mart,
                    useCache = F)


gene_order <- getBM(attributes = c('uniprot_gn_symbol', 'chromosome_name', 'start_position', 'end_position'),
                    filter = 'uniprot_gn_symbol',
                    values = gene_list,
                    mart = h_mart,
                    useCache = F)

gene_order %>% nrow()
gene_order %>%  head()
gene_order$hgnc_symbol %in% gene_list %>% table()
gene_order$hgnc_symbol[24198]
gene_list[24198]



gene_order$uniprot_gn_symbol %in% gene_list %>% table()






## ---- getting hgnc gene id for genes in the mat ---- ##
count_data <- read.table('./proc/pbmc3k_tnbc1_ref_merged_counts.tsv.gz')
gene_list <- rownames(count_data)
gene_list %>% head()
gene_list %>% length



prob_transcripts <- read.table('./proc/problematic_gene_names.txt')
prob_transcripts %>% nrow()
prob_transcripts %>%  head()


transcript_order <- getBM(attributes = c('external_transcript_name', 'chromosome_name', 'start_position', 'end_position'),
                    filter = 'external_transcript_name',
                    values = prob_transcripts$V1,
                    mart = h_mart,
                    useCache = F)
transcript_order %>% nrow()

counts_mat <- read.table('./proc/pbmc3k_tnbc1_ref_merged_counts.tsv.gz')
prob_genes <- read.table('./proc/problematic_gene_names.txt')
prob_genes <- prob_genes$V1
prob_genes %>% length
rownames(counts_mat)


(rownames(counts_mat) %in% prob_genes) %>% which(arr.ind = T)
filt_mat <- counts_mat[!(rownames(counts_mat) %in% prob_genes),]
filt_mat %>% head

write.table(filt_mat,
            gzfile('filtered_probGenesRemoved_pbmc3k_tnbc1_counts.tsv.gz'),
            sep = '\t')
