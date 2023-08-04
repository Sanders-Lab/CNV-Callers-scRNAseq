library(data.table)
library(GenomicRanges)
library(UpSetR)
library(ggplot2)
library(ComplexUpset)
source("../../digital_karyotype/R/utils.R", chdir = T)


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

################################################################################




ck_calls <- fread("../proc/copykat_geneLevel_230331.tsv.gz")
inf_calls <- fread("../proc/infercnv_tnbc1_gene_calls.tsv.gz")
nb_calls <- fread("../proc/numbat_tnbc_cna_segs.tsv.gz")


# removing n_genes col
nb_calls[, n_genes := NULL]
ck_calls[, idx := NULL]
inf_calls[, idx := NULL]

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

nb_dt <- as.data.table(nb_gro)
ck_dt <- as.data.table(ck_gro)
inf_dt <- as.data.table(inf_gro)


nb_dt


# union of the three ranges
union_ranges <- makeGRangesFromDataFrame(inf_calls,
                                         start.field = "start_loc",
                                         end.field = "end_loc",
                                         seqnames.field = "chr_name"
)
union_ranges <- inf_gro
union_ranges <- c(ck_gro, inf_gro, nb_gro)
union_ranges <- GenomicRanges::union(ck_gro, inf_gro)
GenomicRanges::merge(ck_gro, 
                     inf_gro
)
isDisjoint(union_ranges)
union_ranges  <- disjoin(union_ranges)
union_ranges$sv_state  <- "bg"
union_ranges$cell_name  <- "inf"
union_ranges
union_dt <- as.data.table(union_ranges)
setnames(union_dt,
         c("seqnames", "start", "end"),
         c("chrom", "start_loc", "end_loc")
)
union_dt


# visualising the union ranges on digKar
str(Plot_Digital_Karyotype)
Plot_Digital_Karyotype(cell_sv_dt = union_dt,
                       plot_both_haplotypes = F,
                       save_digital_karyotype = T,
                       plot_dir = "union_ranges"

)






union_ranges$bin_id <- 1:length(union_ranges)
union_dt <- as.data.table(union_ranges)
union_dt
union_ranges
isDisjoint(union_ranges)


# ck_gro
# inf_gro


# ck_dis <- disjoin(ck_gro)
# inf_dis <- disjoin(inf_gro)
# nb_dis <- disjoin(nb_gro)
# 
# # ck_dis
# # inf_dis
# nb_dis
# 
# 
# union_ranges <- union(ck_dis, inf_dis, nb_dis, ignore.strand = T)
# union_ranges$bin_id <- 1:length(union_ranges)
# union_ranges
# union_dt <- as.data.table(union_ranges)
# union_dt
# union_dt[bin_id == 1342]



# function to get the overlaps
Get_Overlaps <- function(union_ranges,
                         call_gro, 
                         caller) {

    .union_dt <- as.data.table(union_ranges)
    .call_dt <- as.data.table(call_gro)

    hits <- findOverlaps(call_gro, union_ranges)
    call_seg <- cbind(.call_dt[queryHits(hits)],
                    .union_dt[subjectHits(hits), "bin_id"])
    call_seg <- as.data.table(call_seg)
    call_seg[, caller := caller]


    return(call_seg)
}

# hits <- findOverlaps(ck_gro, union_ranges)
# hits
# cbind(ck_dt[queryHits(hits)],
#       union_dt[subjectHits(hits), "bin_id"])[bin_id == 30]
#                   by = c("seqnames", "start", "end"))
# subsetByOverlaps(ck_gro, union_ranges)



nb_seg <- Get_Overlaps(union_ranges, nb_gro, "numbat")
inf_seg <- Get_Overlaps(union_ranges, inf_gro, "infercnv")
ck_seg <- Get_Overlaps(union_ranges, ck_gro, "copykat")



comp <- c()
comp <- rbind(nb_seg,
              ck_seg,
              inf_seg
)
comp[cnv_state == "", cnv_state := "neu"]
comp


comp_dt <- comp[, .(cell_name,
                    cnv_state,
                    bin_id,
                    caller)]
comp_dt 




all_bins <- union_dt[, bin_id]


# comp <- merge(inf_seg[, cnv_state, by = .(cell_name, bin_id)],
#               ck_seg[, cnv_state, by = .(cell_name, bin_id)],
#               by = c("cell_name", "bin_id"))
# 
# comp <- merge(comp,
#               nb_seg[, cnv_state, by = .(cell_name, bin_id)],
#               by = c("cell_name", "bin_id"))
# 
# setnames(comp,
#          c("cnv_state.x", "cnv_state.y", "cnv_state"),
#          c("infercnv", "copykat", "numbat")
# )


# Adding neu tag
# comp_dt[cnv_state == "", cnv_state := "neu"]
# comp_dt
# comp_dt[cell_name == "tnbc1_AAACCTGCACCTTGTC" &
#           bin_id == 2299]


# this dt is required
comp_wide <- dcast(comp_dt,
      cell_name + bin_id + cnv_state ~ caller,
      fun = length
#       value.var = "cnv_state"
)
comp_wide

str(comp_wide)
comp_wide[, unique(copykat)]
comp_wide[, unique(infercnv)]
comp_wide[, unique(numbat)]


# getting more than 1
curr_caller <- "numbat"
dups_dt <- comp_wide[get(curr_caller) > 1, .(cell_name, bin_id)]
dups_dt

test_dt <- c()
for(nr in 1:nrow(dups_dt)) {

    if(comp[cell_name == dups_dt[nr, cell_name] &
       bin_id == dups_dt[nr, bin_id] &
       caller == curr_caller][, length(unique(cnv_state))] > 1) {

        tmp_dt <- data.table(cell_name = dups_dt[nr, cell_name],
                             bin_id = dups_dt[nr, bin_id],
                             cnv_state = T
        )
        test_dt <- rbind(test_dt, tmp_dt)
#         tmp_dt <- data.table(cell_name = dups_dt[nr, cell_name],
#                              bin_id = dups_dt[nr, bin_id],
#                              cnv_state = F
#         )
#         print(paste0("cell_name : ",
#                      dups_dt[nr, cell_name],
#                      "\nbin_id : ",
#                      dups_dt[nr, bin_id]))
    }


}
test_dt
inf_test <- test_dt
inf_test


comp[cell_name == "tnbc1_AAACCTGCACCTTGTC" &
     bin_id == 30 &
     caller == "copykat"]


comp_wide[copykat == 2 &
          cell_name == "tnbc1_AAACCTGCACCTTGTC"]



all_instances <- comp_wide[, length(unique(cell_name))] *
    comp_wide[, length(unique(bin_id))]
comp_wide[, length(unique(bin_id))]

all_instances

comp_wide[cell_name == "tnbc1_AAACCTGCACCTTGTC" &
          bin_id == 2299]
comp_wide_df <- data.frame(comp_wide)
head(comp_wide_df)
str(comp_wide_df)
typeof(comp_wide_df)


comp[cell_name == "tnbc1_AAACCTGCACCTTGTC" &
          bin_id == 29 &
          caller == "copykat"]

comp[cell_name == "tnbc1_AAACCTGCACCTTGTC" &
          bin_id == 1316 &
          caller == "copykat"]

comp[cell_name == "tnbc1_AAACCTGCACCTTGTC" &
     start == 2391775 &
     seqnames == "chr1" &
     caller == "infercnv"]

union_dt[bin_id == 29]
union_dt[bin_id == 30]
union_dt[bin_id == 1316]

comp[cell_name == "tnbc1_AAACGGGTCCAGAGGA" &
     bin_id == 3869]


################################################################################

str(upset)
pdf("../plots/upset_plots/upset_plot_overlaps.pdf")
upset(data = comp_wide_df,
      sets = c("copykat", "infercnv", "numbat"),
      order.by = "freq"
)
dev.off()



cnv_colors <- c("neu" = "gray",
                "del" = "dodgerblue2",
                "amp" = "tomato2",
                "loh" = "springgreen3"

)

callers <- c("copykat", "infercnv", "numbat")


pdf("../plots/upset_plots/upset_plot_stacked.pdf")
ComplexUpset::upset(data = comp_wide_df,
      callers,
      base_annotations=list(
                            'Intersection size'=intersection_size(
                                                                  counts=FALSE,
                                                                  mapping=aes(fill=cnv_state)
                            ) +
            scale_fill_manual(values = cnv_colors)
      )
)
dev.off()

head(comp_wide_df, n = 15)




pdf("../plots/upset_plots/upset_plot_overlaps_prop.pdf")
ComplexUpset::upset(data = comp_wide_df,
      callers,
      base_annotations = list(
        'Intersection size'=(
            intersection_size(
                text_mapping=aes(
                    label=paste0(round(!!get_size_mode('exclusive_intersection') / nrow(comp_wide_df) * 100, 2), '')
            ),
                mapping=aes(fill=cnv_state)) +
            ylab('Intersection %') +
            scale_y_continuous(labels=scales::percent_format(scale=100 / nrow(comp_wide_df))) +
            scale_fill_manual(values = cnv_colors)
        )
    ),
    set_sizes=(FALSE)
)
dev.off()


pdf("../plots/upset_plots/upset_plot_overlaps_fill.pdf")
ComplexUpset::upset(data = comp_wide_df,
      callers,
      base_annotations = list('Proportion'=list(
        aes=aes(x=intersection, fill=cnv_state),
        geom=list(
            geom_bar(stat='count', position='fill'),
            scale_fill_manual(values = cnv_colors)
        )
    )

    )
)
dev.off()

comp_wide[copykat == 0 &
          infercnv == 0 &
          numbat == 1] 

comp_wide[copykat == 1 &
          infercnv == 1 &
          numbat == 1
      ][, unique(cnv_state)]





################################################################################
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



################################################################################
# melting begins
melt_comp <- melt(comp,
     id.vars = c("cell_name", "bin_id"),
     measure.vars = c("infercnv", "copykat", "numbat")
)
melt_comp
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


neu_calls <- melt_comp[value == "neu" ]
neu_calls[,unique(value)]
neu_calls_wide <- dcast(neu_calls,
                        cell_name + bin_id ~ variable,
                        value.var = "value"
)
neu_calls_wide


# stackoverflow solution to replace all NA with 0
f_dowle2 = function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)), (i):=as.numeric(0)]
}

binarize_calls = function(DT, call) {
  for (i in names(DT))
    DT[get(i) == call, (i):=as.numeric(1)]
}

make_numeric = function(DT) {
  for (i in names(DT))
    DT[, get(i):= as.numeric(i)]
}
binarize_calls(neu_calls_wide, call = "neu")
f_dowle2(neu_calls_wide)

upset_dt <- neu_calls_wide[, .(infercnv,
                   copykat,
                   numbat)]
upset_dt

upset_dt[, infercnv := as.numeric(infercnv)]
upset_dt[, copykat := as.numeric(copykat)]
upset_dt[, numbat := as.numeric(numbat)]

upset_dt[, infercnv := as.numeric(infercnv)]
upset_dt[, copykat := as.numeric(copykat)]
upset_dt[, numbat := as.numeric(numbat)]

# make_numeric(upset_dt)
neu_calls_wide
str(upset_dt)

upset_df <- as.data.frame(upset_dt)
upset_df
str(upset_df)


(head(upset_df))

pdf("../plots/upset_plots/upset_neu.pdf")
upset(upset_df,
      order.by = "freq"
)
dev.off()






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
