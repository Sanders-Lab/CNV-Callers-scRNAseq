library(data.table)
library(GenomicRanges)


#########
# making dt tidy and long form
#########

ck_calls <- fread("../proc/copykat_geneLevel_230331.tsv.gz")
ck_calls <- ck_calls[,!c("abspos",
             "ensembl_gene_id",
             "hgnc_symbol",
             "abspos")]
setnames(ck_calls,
         c("chromosome_name", "start_position", "end_position"),
         c("chr", "start_loc", "end_loc")
)
ck_calls[1:5,1:5]

# making into long form
ck_tidy <- melt(ck_calls, id = c("idx", "chr", "start_loc", "end_loc"))
setnames(ck_tidy,
         c("chr", "variable", "value"),
         c("chr_name", "cell_name", "cnv_state")
)
ck_tidy[, chr_name := paste0("chr", chr_name)]
ck_tidy

fwrite(ck_tidy,
       "../proc/copykat_geneLevel_230331.tsv.gz",
       sep = "\t"
)




inf_calls <- fread("../proc/infercnv_tnbc1_gene_calls.tsv.gz")
setnames(inf_calls,
         c("seqnames", "start", "end"),
         c("chr", "start_loc", "end_loc")
)
inf_calls[1:5,1:5]
fwrite(inf_calls,
       "../proc/infercnv_tnbc1_gene_calls.tsv.gz",
       sep = "\t"
)
inf_tidy <- melt(inf_calls, id = c("idx", "chr", "start_loc", "end_loc"))
setnames(inf_tidy,
         c("chr", "variable", "value"),
         c("chr_name", "cell_name", "cnv_state")
)
inf_tidy

fwrite(inf_tidy,
       "../proc/infercnv_tnbc1_gene_calls.tsv.gz",
       sep = "\t"
)



nb_calls <- fread("../proc/numbat_tnbc_cna_segs.tsv.gz")
nb_calls
setnames(nb_calls,
         c("seqnames", "start", "end"),
         c("chr", "start_loc", "end_loc")
)

setnames(nb_calls,
         "chr",
         "chr_name"
)

fwrite(nb_calls,
       "../proc/numbat_tnbc_cna_segs.tsv.gz",
       sep = "\t"
)

#####################




ck_calls <- fread("../proc/copykat_geneLevel_230331.tsv.gz")
inf_calls <- fread("../proc/infercnv_tnbc1_gene_calls.tsv.gz")
nb_calls <- fread("../proc/numbat_tnbc_cna_segs.tsv.gz")
