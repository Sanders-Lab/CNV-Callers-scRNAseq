library(biomaRt)
library(tidyverse)
library(Seurat)


epithelium <- Read10X('../data/reference/epithelium')
epithelium[1:5,1:5]
epithelium %>% row.names() %>% length()


gene_list <- row.names(epithelium)

# loading the mart for humans
h_mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")



gene_order <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol') ,
                    filter = 'ensembl_gene_id',
                    values = gene_list,
                    mart = h_mart,
                    useCache = F,
                    uniqueRows = T)


gene_order %>% nrow()
gene_order %>% head()
gene_order <- gene_order[!(gene_order$hgnc_symbol %>% duplicated()),]


rownames(epithelium) %in% gene_order$ensembl_gene_id %>% table
epi_subset <- epithelium[rownames(epithelium) %in% gene_order$ensembl_gene_id,]
epi_subset %>% nrow()
row.names(epi_subset) <-  gene_order$hgnc_symbol
epi_subset %>% head()
epi_subset %>% colnames()
epi_subset %>% rownames() %>% duplicated() %>% table() 



epi_subset_mat <- as.matrix(epi_subset)

row.names(epithelium)
write.table(epi_subset_mat,
            gzfile('../proc/epithelium_hgnc_counts.tsv.gz'),
            row.names = T,
            col.names = T)
