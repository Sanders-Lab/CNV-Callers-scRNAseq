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

# nb_dt <- as.data.table(nb_gro)
# ck_dt <- as.data.table(ck_gro)
# inf_dt <- as.data.table(inf_gro)




# union of the three ranges
# union_ranges <- makeGRangesFromDataFrame(inf_calls,
#                                          start.field = "start_loc",
#                                          end.field = "end_loc",
#                                          seqnames.field = "chr_name"
# )


# removing the cell_name and cnv_state as I want the 
# unique of only the ranges
ck_gro$cell_name <- NULL
ck_gro$cnv_state <- NULL

inf_gro$cell_name <- NULL
inf_gro$cnv_state <- NULL

union_ranges <- c(ck_gro, inf_gro)
union_ranges <- unique(union_ranges)
union_ranges

# union_ranges <- GenomicRanges::union(ck_gro, inf_gro)
# GenomicRanges::merge(ck_gro, 
#                      inf_gro
# )
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
                       plot_dir = "upset_plots"

)






# assigning bin ids
union_ranges$bin_id <- 1:length(union_ranges)
union_dt <- as.data.table(union_ranges)
union_dt
union_ranges


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



# ATTENTION: Please rerun the gro creating steps above
# since cell_name and cnv_state was removed earlier
# to make the union ranges.
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
      fun = function(x) {
          length(unique(x))
      },
      value.var = "cnv_state"
)
comp_wide

str(comp_wide)
comp_wide[, unique(copykat)]
comp_wide[, unique(infercnv)]
comp_wide[, unique(numbat)]


# getting more than 1
curr_caller <- "infercnv"
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
     caller == "infercnv"]


comp_wide[copykat == 2 &
          cell_name == "tnbc1_AAACCTGCACCTTGTC"]



all_instances <- comp_wide[, length(unique(cell_name))] *
    comp_wide[, length(unique(bin_id))]
comp_wide[, length(unique(bin_id))]

all_instances

comp_wide[cell_name == "tnbc1_AAACCTGCACCTTGTC" &
          bin_id == 2299]
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
# plotting upset plots.
comp_wide_df <- data.frame(comp_wide)


str(upset)
pdf("../plots/upset_plots/upset_plot_overlaps.pdf")
UpSetR::upset(data = comp_wide_df,
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


# checking how many neu_calls go to other calls
neu_calls <- comp_wide[cnv_state == "neu"]
neu_calls_df <- data.frame(neu_calls)


pdf("../plots/upset_plots/neu_calls_overlap.pdf")
UpSetR::upset(data = neu_calls_df,
      sets = c("copykat", "infercnv", "numbat"),
      order.by = "freq"
)
dev.off()


# prop of neu calls overlap
pdf("../plots/upset_plots/neu_calls_prop.pdf")
ComplexUpset::upset(data = neu_calls_df,
      callers,
      base_annotations = list(
        'Intersection size'=(
            intersection_size(
                text_mapping=aes(
                    label=paste0(round(!!get_size_mode('exclusive_intersection') / nrow(neu_calls_df) * 100, 2), '')
            )) +
            ylab('Intersection %') +
            scale_y_continuous(labels=scales::percent_format(scale=100 / nrow(neu_calls_df)))
        )
    ),
    set_sizes=(FALSE)
)
dev.off()


# checking what calls are called in place of neu
neu_calls

# removing the rows which have all overlap
par_overlaps <- neu_calls[rowSums(neu_calls[, .(copykat, infercnv, numbat)]) < 3]
par_overlaps


# way to get those cell_name and bin_id of par overlaps in comp_dt
new_df <- merge(par_overlaps, comp_dt, all = F, by = c("cell_name", "bin_id"))
new_df

# sanity check
new_df[rowSums(new_df[, .(copykat, infercnv, numbat)]) < 3]
new_df[cnv_state.y != "neu", .N, by = .(cell_name, bin_id)
       ][,unique(N)]


alt_calls <- unique(new_df[cnv_state.y != "neu", .SD, by = .(cell_name, bin_id)])
alt_calls[, cnv_state.x := NULL]
alt_calls[, copykat := NULL]
alt_calls[, infercnv := NULL]
alt_calls[, numbat := NULL]
alt_calls


alt_wide <- dcast(alt_calls,
      cell_name + bin_id + cnv_state.y ~ caller,
      fun = function(x) {
          length(unique(x))
      },
      value.var = "cnv_state.y"
)
alt_wide

alt_wide[, unique(copykat)]
alt_wide[, unique(infercnv)]
alt_wide[, unique(numbat)]

################################################################################
# plotting alt calls (calls where any one caller calls neu)
alt_wide_df <- data.frame(alt_wide)


pdf("../plots/upset_plots/upset_plot_neu_alt_overlaps.pdf")
UpSetR::upset(data = alt_wide_df,
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


pdf("../plots/upset_plots/upset_plot_neu_alt_stacked.pdf")
ComplexUpset::upset(data = alt_wide_df,
      callers,
      base_annotations=list(
                            'Intersection size'=intersection_size(
                                                                  counts=FALSE,
                                                                  mapping=aes(fill=cnv_state.y)
                            ) +
            scale_fill_manual(values = cnv_colors)
      )
)
dev.off()





pdf("../plots/upset_plots/upset_plot_overlaps_neu_alt_prop.pdf")
ComplexUpset::upset(data = alt_wide_df,
      callers,
      base_annotations = list(
        'Intersection size'=(
            intersection_size(
                text_mapping=aes(
                    label=paste0(round(!!get_size_mode('exclusive_intersection') / nrow(alt_wide_df) * 100, 2), '')
            ),
                mapping=aes(fill=cnv_state.y)) +
            ylab('Intersection %') +
            scale_y_continuous(labels=scales::percent_format(scale=100 / nrow(alt_wide_df))) +
            scale_fill_manual(values = cnv_colors)
        )
    ),
    set_sizes=(FALSE)
)
dev.off()


pdf("../plots/upset_plots/upset_plot_overlaps_neu_alt_fill.pdf")
ComplexUpset::upset(data = alt_wide_df,
      callers,
      base_annotations = list('Proportion'=list(
        aes=aes(x=intersection, fill=cnv_state.y),
        geom=list(
            geom_bar(stat='count', position='fill'),
            scale_fill_manual(values = cnv_colors)
        )
    )

    )
)
dev.off()



new_df[cell_name == "tnbc1_TTTGTCATCTTGTATC" &
       bin_id == 7640]



comp_dt[cell_name %in% par_overlaps[, cell_name] &
        bin_id %in% par_overlaps[, bin_id]]


comp_dt[cell_name == "tnbc1_AAACCTGCACCTTGTC" &
        bin_id == 5
    ][cnv_state != "neu"]


################################################################################

# plotting the loc of all overlaps
all_overlaps <- comp_wide[rowSums(comp_wide[, .(copykat, infercnv, numbat)]) == 3]
all_overlaps


two_overlaps <- comp_wide[rowSums(comp_wide[, .(copykat, infercnv, numbat)]) == 2]
two_overlaps


all_overlap_bins <- all_overlaps[, unique(bin_id)]
two_overlap_bins <- two_overlaps[, unique(bin_id)]


all_overlaps_loc <- union_dt[bin_id %in% all_overlap_bins
     ][, .(seqnames, start, end)]
all_overlaps_loc


two_overlaps_loc <- union_dt[bin_id %in% two_overlap_bins
     ][, .(seqnames, start, end)]
two_overlaps_loc


non_overlaps <- union_dt[!bin_id %in% c(all_overlap_bins, two_overlap_bins)
     ][, .(seqnames, start, end)]
non_overlaps

non_overlaps <- union_dt[!bin_id %in% all_overlap_bins
                         ][, .(seqnames, start, end)]

all_overlaps_loc[, cell_name := "all_overlaps"]
all_overlaps_loc[, sv_state := "all_agreement"]
all_overlaps_loc


two_overlaps_loc[, cell_name := "all_overlaps"]
two_overlaps_loc[, sv_state := "two_agreement"]
two_overlaps_loc



non_overlaps[, cell_name := "all_overlaps"]
non_overlaps[, sv_state := "non_agreement"]
non_overlaps

plot_all_overlaps <- rbind(all_overlaps_loc,
#                            two_overlaps_loc,
                           non_overlaps)

plot_all_overlaps <- rbind(all_overlaps_loc,
                           two_overlaps_loc,
                           non_overlaps)

setnames(plot_all_overlaps,
         c("seqnames", "start", "end"),
         c("chrom", "start_loc", "end_loc")
)
plot_all_overlaps


agree_col <- c("springgreen3", "darkorchid", "tomato2")
names(agree_col) <- c("all_agreement", "two_agreement", "non_agreement")


Plot_Digital_Karyotype(cell_sv_dt = plot_all_overlaps,
                       plot_both_haplotypes = F,
                       save_digital_karyotype = T,
                       color_pal = agree_col,
                       plot_dir = "upset_plots",
                       kar_label = "Agreement loci between three callers"

)


# changing cell name to save to a diff file
plot_all_overlaps[, cell_name := "all_combo_overlaps"]
Plot_Digital_Karyotype(cell_sv_dt = plot_all_overlaps,
                       plot_both_haplotypes = F,
                       save_digital_karyotype = T,
                       color_pal = agree_col,
                       plot_dir = "upset_plots",
                       kar_label = "Agreement loci between three callers"

)

