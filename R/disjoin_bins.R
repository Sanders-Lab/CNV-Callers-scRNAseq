library(data.table)
library(GenomicRanges)
library(UpSetR)
library(ggplot2)


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


# ck_calls
# inf_calls
# nb_calls



# inf_calls[, length(unique(cell_name))]
# nb_calls[, length(unique(cell_name))]
# 
# ck_calls[chr_name == "chr23", chr_name := "chrX"]
# fwrite(ck_calls,
#        "../proc/copykat_geneLevel_230331.tsv.gz",
#        sep = "\t"
# )


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


# ck_gro
# inf_gro


ck_dis <- disjoin(ck_gro)
inf_dis <- disjoin(inf_gro)

# ck_dis
# inf_dis


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

setnames(comp,
         c("cnv_state.x", "cnv_state.y", "cnv_state"),
         c("infercnv", "copykat", "numbat")
)
comp
comp[duplicated(comp)]
unique(comp)


# comp[bin_id == 2 &
#      cell_name == comp[1, cell_name]]


# Adding neu tag
comp[infercnv != "amp" &
     infercnv != "del", infercnv := "neu"]
comp[copykat != "amp" &
     copykat != "del", copykat := "neu"]
comp[numbat != "amp" &
     numbat != "del", numbat := "neu"]
comp



# upset_list <- list(infercnv = comp[, infercnv],
#                    copykat = comp[, copykat],
#                    numbat = comp[, numbat]
# )
# str(upset_list)
upset_dt <- comp[, .(infercnv,
                        copykat,
                        numbat)]
# upset_dt[upset_dt$"infercnv" == "amp", "infercnv"]


upset_dt[infercnv == "del", infercnv := 1]
upset_dt[infercnv == "neu", infercnv := 2]
upset_dt[infercnv == "amp", infercnv := 3]

upset_dt[copykat == "del", copykat := 1]
upset_dt[copykat == "neu", copykat := 2]
upset_dt[copykat == "amp", copykat := 3]

upset_dt[numbat == "del", numbat := 1]
upset_dt[numbat == "neu", numbat := 2]
upset_dt[numbat == "amp", numbat := 3]
upset_dt

# binary
upset_dt[infercnv == "del", infercnv := 1]
upset_dt[infercnv == "neu", infercnv := 0]
upset_dt[infercnv == "amp", infercnv := 1]

upset_dt[copykat == "del", copykat := 1]
upset_dt[copykat == "neu", copykat := 0]
upset_dt[copykat == "amp", copykat := 1]

upset_dt[numbat == "del", numbat := 1]
upset_dt[numbat == "neu", numbat := 0]
upset_dt[numbat == "amp", numbat := 1]

upset_df <- data.frame(infercnv = as.numeric(as.vector(upset_dt[, infercnv])),
                       copykat = as.numeric(as.vector(upset_dt[, copykat])),
                       numbat = as.numeric(as.vector(upset_dt[, numbat]))
)
str(upset_df)



upset_dt[, infercnv := as.integer(infercnv)]
upset_dt[, copykat := as.integer(copykat)]
upset_dt[, numbat := as.integer(numbat)]
str(upset_dt)

upset_df <- as.data.frame(upset_dt)
upset_list <- as.list(upset_dt)

attr(upset_df, "index") <- NULL
attr(upset_list, "index") <- NULL


row.names(upset_df) <- NULL
str(upset_df)
str(upset_list)
upset_df[1:5,1:3]
nrow(upset_dt[infercnv != 2])

# Example data
data <- data.frame(
  Set1 = c(1, 0, 1, 1, 0, 0),
  Set2 = c(1, 1, 1, 0, 0, 1),
  Set3 = c(1, 0, 1, 0, 1, 0)
)
str(data)

# Create UpSet plot
pdf("../plots/upset_plots/upset_plot_chatgpt.pdf")
upset(data)
dev.off()



pdf("../plots/upset_plots/upset_plot.pdf")
upset(upset_df,
      order.by = "freq"
)
dev.off()

str(obj)



pdf("../plots/upset_plots/upset_plot.pdf")
upset(fromList(movies))
dev.off()

str(up_plt)

set_dt <- data.table()
set_dt$all_overlap <- upset_dt[infercnv == copykat &
         infercnv == numbat &
         numbat == copykat, .N]
set_dt$infercnv_copykat <- upset_dt[infercnv == copykat
         , .N]
set_dt$numbat_copykat <- upset_dt[numbat == copykat
         , .N]
set_dt$infercnv_numbat <- upset_dt[infercnv == numbat
         , .N]
set_dt

plot_set_dt <- melt(set_dt)
plot_set_dt


movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
    header = T, sep = ";")
str(movies)

pdf("../plots/upset_plots/manual_upset.pdf")
ggplot(data = plot_set_dt,
       aes(x = variable,
           y = value
       )) +
      geom_col()
dev.off()



comp
comp[cell_name == "tnbc1_AAACCTGCACCTTGTC" &
     bin_id == 4212]
comp <- unique(comp)
comp




melt_comp <- melt(comp,
     id.vars = c("cell_name", "bin_id"),
     measure.vars = c("infercnv", "copykat", "numbat")
)
melt_comp[value == "neu"]
melt_comp[value == "del"]
melt_comp[value == "amp"]


# converting cnv_state to values
# melt_comp[, value := gsub("amp", 1, value)]
# melt_comp[, value := gsub("neu", 0, value)]
# melt_comp[, value := gsub("del", -1, value)]
# 
# 
# 
# melt_comp[value == 1]
# melt_comp[value == -1]
# melt_comp[value == 0]
# melt_comp


str(melt_comp)
melt_comp



melt_comp <- unique(melt_comp)

melt_comp[value == "neu", n_overlaps := .N, by = .(cell_name,
                                               bin_id)]
melt_comp[value == "amp", n_overlaps := .N, by = .(cell_name,
                                               bin_id)]
melt_comp[value == "del", n_overlaps := .N, by = .(cell_name,
                                               bin_id)]
melt_comp[, unique(n_overlaps)]
melt_comp

melt_comp[n_overlaps == 2,
          variable[1],
          by = .(cell_name, bin_id)]
melt_comp[n_overlaps == 2,
          variable[2],
          by = .(cell_name, bin_id)]
melt_comp[n_overlaps == 2,
          overlap_pairs := paste(variable[1], variable[2], sep = "-"),
          by = .(cell_name, bin_id)]
melt_comp[is.na(overlap_pairs), overlap_pairs := variable]
melt_comp[overlap_pairs == "infercnv-infercnv", .SD, by = .(cell_name, bin_id)]



plot_dt <- as.data.table(melt_comp[n_overlaps == 3, round(table(value) / .N, 2)])
plot_dt <- rbind(plot_dt,
                 as.data.table(melt_comp[n_overlaps == 2, round(table(value) / .N, 2),
                                            by = .(overlap_pairs)]),
                 fill = T)
plot_dt
melt_comp[n_overlaps == 2, table(value) / .N]
melt_comp[n_overlaps == 2, round(table(value) / .N, 2), by = .(overlap_pairs)]


as.data.table(melt_comp[, table(n_overlaps)])
