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
ck_calls[1:5,1:6]

# making into long form
ck_tidy <- melt(ck_calls, id = c("idx", "chr_name", "start_loc", "end_loc"))
setnames(ck_tidy,
         c("chr", "variable", "value"),
         c("chr_name", "cell_name", "cnv_state"),
         skip_absent = T
)
ck_tidy[, chr_name := paste0("chr", chr_name)]
ck_tidy[chr_name == "chr23", chr_name := "chrX"]
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


inf_tidy <- inf_calls[cell_name %in% ck_calls[, unique(cell_name)]]
inf_tidy[, length(unique(cell_name))]

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


nb_calls <- nb_calls[cell_name %in% ck_calls[, unique(cell_name)]]
nb_calls[, length(unique(cell_name))]

fwrite(nb_calls,
       "../proc/numbat_tnbc_cna_segs.tsv.gz",
       sep = "\t"
)

#####################




ck_calls <- fread("../proc/copykat_geneLevel_230331.tsv.gz")
inf_calls <- fread("../proc/infercnv_tnbc1_gene_calls.tsv.gz")
nb_calls <- fread("../proc/numbat_tnbc_cna_segs.tsv.gz")


ck_calls
inf_calls
nb_calls



inf_calls[, length(unique(cell_name))]
nb_calls[, length(unique(cell_name))]

ck_calls[chr_name == "chr23", chr_name := "chrX"]
fwrite(ck_calls,
       "../proc/copykat_geneLevel_230331.tsv.gz",
       sep = "\t"
)


ck_gro <- makeGRangesFromDataFrame(ck_calls,
                                   start.field = "start_loc",
                                   end.field = "end_loc",
                                   seqnames.field = "chr_name",
                                   keep.extra.columns = T
)
inf_gro <- makeGRangesFromDataFrame(inf_calls,
                                   start.field = "start_loc",
                                   end.field = "end_loc",
                                   seqnames.field = "chr_name",
                                   keep.extra.columns = T
)
nb_gro <- makeGRangesFromDataFrame(nb_calls,
                                   start.field = "start_loc",
                                   end.field = "end_loc",
                                   seqnames.field = "chr_name",
                                   keep.extra.columns = T
)


ck_gro
inf_gro


ck_dis <- disjoin(ck_gro)
inf_dis <- disjoin(inf_gro)

ck_dis
inf_dis


union_ranges <- union(ck_dis, inf_dis)
union_ranges$bin_id <- 1:length(union_ranges)
union_ranges
union_dt <- as.data.table(union_ranges)
union_dt



hits <- findOverlaps(union_ranges, nb_gro)
nb_seg <- cbind(nb_calls[subjectHits(hits)],
      union_dt[queryHits(hits), "bin_id"])
nb_seg



hits <- findOverlaps(union_ranges, ck_gro)
ck_seg <- cbind(ck_calls[subjectHits(hits)],
      union_dt[queryHits(hits), "bin_id"])
ck_seg



hits <- findOverlaps(union_ranges, inf_gro)
inf_seg <- cbind(inf_calls[subjectHits(hits)],
      union_dt[queryHits(hits), "bin_id"])
inf_seg



all_bins <- union_dt[, bin_id]


comp <- merge(inf_seg[, cnv_state, by = .(cell_name, bin_id)],
              ck_seg[, cnv_state, by = .(cell_name, bin_id)],
              by = c("cell_name", "bin_id"))

comp <- merge(comp,
              nb_seg[, cnv_state, by = .(cell_name, bin_id)],
              by = c("cell_name", "bin_id"))
)
comp
setnames(comp,

)



for(bin in all_bins) {
    nb_seg[bin_id == bin,
           cnv_state]
}




