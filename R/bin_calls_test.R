library(data.table)
library(ggplot2)
library(GenomicRanges)
library(foreach)
library(parallel)
library(primatR)



# reading in the segs
copykat_segs <- fread("../proc/copykat_tnb1_segs_refined.tsv.gz")
copykat_segs


infercnv_segs <- fread("../proc/infercnv_tnbc1_segs_refined.tsv.gz")
infercnv_segs

infercnv_segs <- infercnv_segs[cell_name %in% copykat_segs$cell_name]


numbat_segs <- fread("../proc/numbat_tnbc_cna_segs.tsv.gz")
numbat_segs <- numbat_segs[cell_name %in% copykat_segs$cell_name]
numbat_segs <- numbat_segs[cnv_state != ""]
numbat_segs


ideo_dt <- fread("../../digital_karyotype/proc/ideogram_scaffold.tsv")
ideo_width <- ideo_dt[, width]
names(ideo_width) <- ideo_dt[, chrom]
ideo_width


bin_size <- 1e6
bin_genome <- tileGenome(ideo_width,
    tilewidth = bin_size,
    cut.last.tile.in.chrom = T
)
length(bin_genome)
bin_genome$bin_id <- 1:length(bin_genome)
bin_genome




# making gro for each
ck_gro <- makeGRangesFromDataFrame(copykat_segs,
    keep.extra.columns = T
)
inf_gro <- makeGRangesFromDataFrame(infercnv_segs,
    keep.extra.columns = T
)
nb_gro <- makeGRangesFromDataFrame(numbat_segs,
    keep.extra.columns = T
)

# nb_temp_gro <- makeGRangesFromDataFrame(numbat_segs[cell_name == numbat_segs[1,cell_name]],
#                                         keep.extra.columns = T
# )
# nb_temp_gro


# findOverlaps(bin_genome,
#                             nb_temp_gro
# )


nb_gro
# nb_perc_overlaps <- as.data.table(getSegDupOverlaps(nb_bins_gro,
#                                                     nb_gro
# )
# )
# nb_perc_overlaps


nb_overlaps <- findOverlaps(
    bin_genome,
    nb_gro
)
nb_overlaps
# overlap <- IRanges::intersect(bin_genome,
#                    nb_gro
# )
# overlap
# percOverlap <- (sum(width(overlap))/width(nb_gro))*100
# percOverlap

nb_gro
nb_gro[nb_overlaps@to]
nb_overlaps
bin_genome[nb_overlaps@from]$bin_id

nb_bins <- as.data.table(nb_gro[nb_overlaps@to])
nb_bins[, bin_id := bin_genome[nb_overlaps@from]$bin_id]
nb_bins






# width(nb_gro[5])
# width(bin_genome[(nb_overlaps@from) == 5])
#
# width(nb_gro[5]) /
#     width(bin_genome[(nb_overlaps@from) == 5])
#
# queryHits(nb_overlaps) |>
#     head()
# width(bin_genome[(nb_overlaps@from) == 5])
# width(bin_genome[nb_overlaps@from][5]) |>
#     sum()



# nb_overlaps |>
#     str()
# width(nb_bins_gro)
#
# binned_genome_length <- sum(width(binned_granges))
# overlap_length <-
#
# sum(width(bin_genome[queryHits(nb_overlaps)])) /sum(width(bin_genome))
#
#
#
# subsetByOverlaps(bin_genome,
#                  nb_gro
# )
#
#
# width(bin_genome[nb_overlaps@from]) /
#     width(nb_gro[nb_overlaps@to]) *100
#
# nb_gro[nb_overlaps@to]




inf_overlaps <- findOverlaps(
    bin_genome,
    inf_gro
)
inf_overlaps
inf_bins <- as.data.table(inf_gro[inf_overlaps@to])
inf_bins[, bin_id := bin_genome[inf_overlaps@from]$bin_id]
inf_bins


ck_overlaps <- findOverlaps(
    bin_genome,
    ck_gro
)
ck_gro
ck_overlaps@to
ck_gro[1]
ck_gro[65]
ck_bins <- as.data.table(ck_gro[ck_overlaps@to])
ck_bins

bin_genome[ck_overlaps@from]$bin_id

ck_bins[, !which(duplicated(.SD)), by = cell_name]
ck_bins[, bin_id := bin_genome[ck_overlaps@from]$bin_id]
ck_bins


# mode func copied from SO
# Mode <- function(x) {
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }
#
#
#
# nb_bins_mode <- nb_bins[, Mode(cnv_state), by = .(bin_id, cell_name)]
# nb_bins_mode
# setnames(nb_bins_mode,
#          "V1",
#          "mode_call"
# )
# nb_bins_mode
#
# inf_bins_mode <- inf_bins[, Mode(cnv_state), by = bin_id]
# setnames(inf_bins_mode,
#          "V1",
#          "mode_call"
# )
# inf_bins_mode
#
# ck_bins_mode <- ck_bins[, Mode(cnv_state), by = bin_id]
# setnames(ck_bins_mode,
#          "V1",
#          "mode_call"
# )
# ck_bins_mode





inf_bins[, .(.N, duplicated(bin_id)), by = cell_name]
ck_bins[, .SD, by = cell_name]


temp <- union(
    nb_bins[
        cell_name == nb_bins[1, cell_name],
        bin_id
    ],
    inf_bins[
        cell_name == nb_bins[1, cell_name],
        bin_id
    ]
)
temp <- union(
    temp,
    ck_bins[
        cell_name == nb_bins[1, cell_name],
        bin_id
    ]
)
temp




nb_bins[bin_id == 2999]
all_cells <- unique(copykat_segs$cell_name)
all_cells
bin_gen_dt <- as.data.table(bin_genome)
all_bins <- bin_gen_dt$bin_id
length(all_bins)




nb_bins[
    cell_name == all_cells[1],
    which(duplicated(bin_id))
]
cn

bin_gen_dt[bin_id == 120]
nb_bins[cell_name == all_cells[1]][121, bin_id]

bin_gen_dt[bin_id == cn[1], start]
nb_bins[cell_name == all_cells[1]][cn[1]][, start] -
    bin_gen_dt[bin_id == cn[1], start] > 0
(bin_gen_dt[bin_id == 121, end] -
    nb_bins[cell_name == all_cells[1]][121][, end]) /
    1e6 * 100

nb_bins[cell_name == all_cells[1]][121][, start] -
    bin_gen_dt[bin_id == 120, start]


(bin_gen_dt[bin_id == cn[1], start] - nb_bins[cell_name == all_cells[1]][cn[1] + 1][, start]) / 1e6
bin_gen_dt[bin_id == bins_nb[1], end] - nb_bins[cell_name == all_cells[1]][bins_nb[1] + 1][, end]






# c <- nb_bins[cell_name == all_cells[1],
#         which(duplicated(bin_id))]
# c
# nb_bins[cell_name == all_cells[1]][!c]
#         !c]
#

bin_gen_dt[bin_id == 120, start]
up_part <-
    (nb_bins[cell_name == all_cells[1]][120][, end] - bin_gen_dt[bin_id == 120, start]) / 1e6 * 100
nb_bins[cell_name == all_cells[1]][120][, start] - bin_gen_dt[bin_id == 120, start]
abs(nb_bins[cell_name == all_cells[1]][120][, end] - bin_gen_dt[bin_id == 120, end])


nb_bins[cell_name == all_cells[1]][120]



ck_bins[cnv_state == "confused"]
inf_bins[cnv_state == "confused"]
# comp_bins[nb_calls != ""]





# 2023-06-09 15:46
# PercentOverlap working!
PercentOverlap <- function(calls_dt) {
    for (entry in 1:nrow(calls_dt)) {
        curr_bin <- calls_dt[entry, bin_id]

        # for 2nd event bigger than 1st
        if (calls_dt[entry, start] -
            bin_gen_dt[bin_id == curr_bin, start] < 0 &
            calls_dt[entry, end] -
                bin_gen_dt[bin_id == curr_bin, end] > 0) {
            calls_dt[entry, perc_overlap := 100]
            next
        }

        # for 2nd event right shifted
        if (calls_dt[entry, start] -
            bin_gen_dt[bin_id == curr_bin, start] > 0) {
            calls_dt[entry, perc_overlap := ((calls_dt[entry, start] -
                bin_gen_dt[bin_id == curr_bin, start]) / 1e6 * 100)]
        }

        # for 2nd event left shifted
        if (calls_dt[entry, start] -
            bin_gen_dt[bin_id == curr_bin, start] < 0) {
            calls_dt[entry, perc_overlap := ((bin_gen_dt[bin_id == curr_bin, end] -
                calls_dt[entry, end]) / 1e6 * 100)]
        }
    }

    return(calls_dt)
}


cellname <- nb_bins[1, cell_name]
nb_bins
nb_bins[cell_name == cellname]
inf_bins[cell_name == cellname]
ck_bins[cell_name == cellname]
check <- PercentOverlap(nb_bins[cell_name == cellname])
check


n_cores <- 63
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


# for nb calls
system.time(
    nb_perc <- foreach(celln = seq_along(all_cells), .combine = "rbind") %dopar% {
        PercentOverlap(nb_bins[cell_name == all_cells[celln]])
    }
)
nb_perc
nb_perc[cell_name == cellname][120:121]

parallel::stopCluster(cluster)





n_cores <- 63
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)

# for inf calls
system.time(
    inf_perc <- foreach(celln = seq_along(all_cells), .combine = "rbind") %dopar% {
        PercentOverlap(inf_bins[cell_name == all_cells[celln]])
    }
)
inf_perc
parallel::stopCluster(cluster)



n_cores <- 63
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)

# for inf calls
system.time(
    ck_perc <- foreach(celln = seq_along(all_cells), .combine = "rbind") %dopar% {
        PercentOverlap(ck_bins[cell_name == all_cells[celln]])
    }
)
ck_perc
parallel::stopCluster(cluster)


fwrite(nb_perc,
    "../proc/numbat_binned_calls_perc_overlap.tsv.gz",
    sep = "\t"
)

fwrite(inf_perc,
    "../proc/infercnv_binned_calls_perc_overlap.tsv.gz",
    sep = "\t"
)

fwrite(ck_perc,
    "../proc/copykat_binned_calls_perc_overlap.tsv.gz",
    sep = "\t"
)


nb_perc[, unique(perc_overlap)]
nb_perc
nb_perc[perc_overlap > 50][, which(duplicated(bin_id)), by = cell_name]

inf_perc[perc_overlap > 50][, which(duplicated(bin_id)), by = cell_name]

ck_perc[perc_overlap > 50][, which(duplicated(bin_id)), by = cell_name]

ck_perc[
    cell_name == cellname,
    which(duplicated(bin_id))
]

dups <- integer()
temp <- c(ck_perc[
    cell_name == cellname,
    which(duplicated(bin_id))
])
temp
dups[cellname] <- c(temp)
c(ck_perc[
    cell_name == cellname,
    which(duplicated(bin_id))
])
dups

ck_perc[cell_name == cellname][190:191]
ck_perc[cell_name == cellname][758:759]

ck_perc[, idx := 1:nrow(ck_perc)]
ck_perc



ck_perc[which(duplicated(bin_id)), .SD, by = cell_name]

remove_dups <- function(sub_dt) {
    for (dups in 1:nrow(sub_dt)) {
        curr_idx <- sub_dt[dups, idx]

        if (ck_perc[curr_idx, perc_overlap] >
            ck_perc[(curr_idx - 1), perc_overlap]) {
            ck_perc <- ck_perc[!(curr_idx - 1)]
            next
        }
        if (ck_perc[curr_idx, perc_overlap] <
            ck_perc[(curr_idx - 1), perc_overlap]) {
            ck_perc <- ck_perc[!curr_idx]
        }
        if (ck_perc[curr_idx, perc_overlap] ==
            ck_perc[(curr_idx - 1), perc_overlap]) {
            ck_perc <- ck_perc[!(curr_idx - 1)]
        }
    }

    return(1)
}


copy <- ck_perc
copy[which(duplicated(bin_id)), remove_dups(.SD), by = cell_name]


ck_perc
ck_dedup
sub_dt <- ck_perc[cell_name %like% "AAACCTGCACCTTGTC" &
    idx %in% which(duplicated(bin_id))]
ck_perc[cell_name %like% "AAACCTGCACCTTGTC" &
    idx %in% which(duplicated(bin_id))]



ck_perc[which(duplicated(bin_id)), .SD, by = cell_name]
ck_perc[189:191]
ck_perc[1336:1338]
ck_perc[perc_overlap > 50][which(duplicated(bin_id)), .SD, by = cell_name]

ck_perc[1338:1340, bin_id]



ck_dedup[
    cell_name == "tnbc1_AAACGGGTCCAGAGGA",
    which(duplicated(bin_id))
]
ck_dedup[cell_name == "tnbc1_AAACGGGTCCAGAGGA", ][100:102]
ck_dedup[cell_name == "tnbc1_AAACGGGTCCAGAGGA", ][313:315]


ck_dedup[1664:1666]


fwrite(ck_dedup,
    "../proc/copykat_binned_calls_perc_overlap.tsv.gz",
    sep = "\t"
)


all_bins <- bin_gen_dt$bin_id
all_bins
comp_dt <- data.table()

# for each cell
compare_bins <- function(call_dt, caller) {
    foreach(bin = all_bins) %do% {
        temp <- data.table(
            bin_id = integer(1),
            cell_name = character(1),
            caller = character(1),
            cnv_state = character(1)
        )

        if (bin %in% call_dt[, bin_id]) {
            temp[1, bin_id := bin]
            temp[1, cell_name := (call_dt[1, cell_name])]
            temp[1, caller := caller]
            temp[1, cnv_state := (call_dt[bin_id == bin, cnv_state])]
        }

        comp_dt <- rbind(comp_dt, temp)
    }


    return(comp_dt)
}

nb_perc <- fread("../proc/numbat_binned_calls_perc_overlap.tsv.gz")

nb_perc[, .SD, by = cell_name][, cell_name]
comp_dt <- nb_perc[, compare_bins(.SD, "numbat"), by = cell_name]
comp_dt

nb_perc <-
    nb_perc[, which(duplicated(nb_perc[, .(seqnames, start, end)])), by = cell_name]
nb_perc


nb_perc[10041:10043]
nb_perc[cell_name == "tnbc1_ACGCAGCCAACACCCG"]


nb_perc <- nb_perc[perc_overlap > 50]
nb_perc


inf_perc


comp_dt <- nb_perc[, .(
    cell_name,
    bin_id,
    cnv_state,
    perc_overlap
)][, caller := "Numbat"]
comp_dt

comp_dt <- rbind(
    comp_dt,
    inf_perc[, .(
        cell_name,
        bin_id,
        cnv_state,
        perc_overlap
    )][, caller := "InferCNV"]
)



comp_dt <- rbind(
    comp_dt,
    ck_perc[, .(
        cell_name,
        bin_id,
        cnv_state,
        perc_overlap
    )][, caller := "Copykat"]
)

comp_dt

comp_dt[, cnv_numeric := ifelse(cnv_state == "amp", 1, -1)]
comp_dt[cnv_state == "amp"]


tmp_plot <- "../plots/tmp/tmp_plot.pdf"
pdf("../plots/cnv_callers_stats/cell1_binned_calls_callers.pdf")
ggplot(data = comp_dt[cell_name == cellname]) +
    geom_col(aes(
        x = bin_id,
        y = cnv_numeric,
        fill = caller
    ))
dev.off()

ggsave(tmp_plot)



bin_gen_dt[, median := start + median(c(start, end))]
bin_gen_dt


sapply(
    bin_gen_dt,
    median()
)


comp_dt <- merge(comp_dt, bin_gen_dt, by = "bin_id")



ggplot(data = comp_dt[cell_name == cellname &
    seqnames == "chr1"]) +
    geom_col(aes(
        x = median,
        y = cnv_numeric,
        fill = caller
    )) +
    coord_fixed(ratio = 1000000) +
    xlab("") +
    ylab("")
ggsave(tmp_plot)



####
# STATISTICS
####


nb_perc <- fread("../proc/numbat_binned_calls_perc_overlap.tsv.gz")
inf_perc <- fread("../proc/infercnv_binned_calls_perc_overlap.tsv.gz")
ck_perc <- fread("../proc/copykat_binned_calls_perc_overlap.tsv.gz")

nb_perc


comp_dt <- nb_perc[, .(
    cell_name,
    bin_id,
    cnv_state,
    perc_overlap
)][, caller := "Numbat"]

comp_dt <- rbind(
    comp_dt,
    inf_perc[, .(
        cell_name,
        bin_id,
        cnv_state,
        perc_overlap
    )][, caller := "InferCNV"]
)


comp_dt <- rbind(
    comp_dt,
    ck_perc[, .(
        cell_name,
        bin_id,
        cnv_state,
        perc_overlap
    )][, caller := "Copykat"]
)

comp_dt <- comp_dt[perc_overlap > 50]



cellname <- comp_dt[1, cell_name]
comp_dt
comp_dt[, .N, by = .(cell_name, bin_id)]
comp_dt[, .N, by = .(cell_name, bin_id, cnv_state)]
comp_dt[, .N, by = .(cell_name, bin_id, cnv_state)][N == 3]



# VERY IMPORTANT STEP
# grouping by cells, binid and cnv_state
caller_coverage <- comp_dt[, .SD, by = .(cell_name, bin_id, cnv_state)]
caller_coverage



caller_coverage[hits == 3, length(unique(cell_name))]
a <- paste(caller_coverage[hits == 2, .SD, by = .(cell_name, bin_id, caller)][cell_name == cellname & bin_id == 21][, caller], collapse = ",")
length(a)
a




called_bins <- caller_coverage[, length(unique(bin_id))]
caller_coverage[, hits := length(unique(caller)), by = .(bin_id, cell_name, cnv_state)]
caller_coverage[, .SD, by = .(hits, cell_name, cnv_state)]
caller_coverage[hits == 3]
caller_coverage



caller_coverage[, length(unique(bin_id)), by = .(hits,  cnv_state)]
caller_coverage[, round(length(unique(bin_id)) / called_bins * 100, 2), by = .(hits, cell_name, cnv_state)]
n_hits <- caller_coverage[, round(length(unique(bin_id)) / called_bins * 100, 2), by = .(hits, cnv_state)]
n_hits


# tmp_plot <- "../plots/tmp/tmp_plot.pdf"

pdf("../plots/cnv_callers_stats/agreement_bin_prop_between_callers.pdf")

ggplot(data = n_hits,
       aes(
           x = hits,
           y = V1,
           fill = cnv_state,
           label = V1
       ) 
       ) +
                 geom_col(position = "dodge") +
                 geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
        scale_fill_manual(values = cnv_colors) +
        xlab("Number of callers") +
        ylab("Proportion of called bins (%)") +
        ggtitle("Extent of agreement of cnv calls between callers per cell per bin")

dev.off()


# prop of cells hit 
caller_coverage[, n_cells := length(unique(cell_name)), by = .(hits, bin_id, cnv_state)]
caller_coverage[, unique(n_cells), by = .(hits, cnv_state)]

caller_coverage[, length(unique(cell_name)), by = .(hits, cnv_state)]

caller_coverage[, length(unique(cell_name)), by = .(hits, cnv_state)]
n_cells <- caller_coverage[, round(length(unique(cell_name)) / 796 * 100, 2), by = .(hits, cnv_state)]
n_cells

pdf("../plots/cnv_callers_stats/agreement_cells_between_callers.pdf")

ggplot(data = n_cells,
       aes(
           x = hits,
           y = V1,
           fill = cnv_state,
           label = V1
       ) 
       ) +
                 geom_col(position = "dodge") +
                 geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
                 scale_fill_manual(values = cnv_colors) +
                 xlab("Number of callers with the same call") +
                 ylab("Proportion of cells (%)") +
                 ggtitle("Extent of agreement of cnv calls between callers per cell per bin")

dev.off()



# prop of bins in 2 caller coverage
# getting the pairs of callers
caller_pairs <- caller_coverage[hits == 2,
                                paste(caller, sep = "-", collapse = "-"),
                                by = .(cell_name, bin_id, cnv_state)]
setnames(caller_pairs,
         "V1",
         "caller_pairs"
)
caller_pairs
caller_pairs[, unique(caller_pairs)]
caller_pairs_stats <- caller_pairs[,round(length(unique(bin_id)) / called_bins * 100, 2), 
                                   by = .(cnv_state, caller_pairs)]
setnames(caller_pairs_stats,
         "V1",
         "prop_bins"
)
caller_pairs_stats




cnv_colors <- c('tomato2', 'dodgerblue2', "springgreen3")
names(cnv_colors) <- c('amp', 'del', "loh")
cnv_colors

# plotting the caller pairs 
pdf("../plots/cnv_callers_stats/agreement_caller_pairs.pdf")

ggplot(data = caller_pairs_stats,
       aes(x = caller_pairs,
           y = prop_bins,
           fill = cnv_state,
           label = prop_bins
       )
       )+
                 geom_col(position = "dodge") +
                 geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
         scale_fill_manual(values = cnv_colors) +
         xlab("Caller pairs") +
         ylab("Proportion of bins (%)") +
         ggtitle("Caller pairs in agreement per cell per bin per cnv call")

dev.off()



bin_gen_dt <- as.data.table(bin_genome)
all_bins <- bin_gen_dt$bin_id
caller_coverage

caller_coverage <- merge.data.table(caller_coverage, 
                                    bin_gen_dt[, .(seqnames,
                                                   start,
                                                   end,
                                                   bin_id)],
                                    by = "bin_id"
)
caller_coverage


nbins_chr <- caller_coverage[, length(unique(bin_id)), by = seqnames]
setnames(nbins_chr,
         "V1",
         "nbins_chr"
)
nbins_chr
caller_coverage[, length(unique(bin_id))]
nbins_chr[,sum(nbins_chr)]

caller_coverage <- merge(caller_coverage,
                         nbins_chr,
                         by = "seqnames"
)
caller_coverage


caller_coverage[, length(unique(bin_id)),
                             by = .(seqnames,
                                    cnv_state,
                                    caller)]

# prop of chr bins covered
calls_chr <- caller_coverage[, unique(round(length(unique(bin_id)) / nbins_chr * 100, 2)),
                             by = .(seqnames,
                                    cnv_state,
                                    caller)]
setnames(calls_chr,
         "V1",
         "prop_bins"
)
# mean prop of chr coverage by caller
calls_chr[, mean_chr_coverage := round(mean(prop_bins), 2), by = .(seqnames,
                                                                   caller)]
calls_chr[, unique(mean_chr_coverage), by = .(seqnames, caller)] 
calls_chr




# plotting calls across chr
pdf("../plots/cnv_callers_stats/calls_across_chr.pdf",
    width = 29.7,
    height = 21
)

ggplot(data = calls_chr,
       aes(x = seqnames,
           y = prop_bins,
           fill = cnv_state,
           label = prop_bins,
           width = 0.6
       )
       ) +
                 geom_col(position = "dodge") +
                 geom_text(position = position_dodge(width = 1), vjust = -0.5, size = 5) +
    geom_line(aes(x = seqnames,
                  y = mean_chr_coverage,
                  group = 1),
              linetype = "dashed") +
    facet_grid(rows = vars(caller)) +
    scale_fill_manual(values = cnv_colors) +
    labs(fill = "CNV State",
         x = "chromosomes",
         y = "Proportion of chr bins covered (%)") +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 0.5),
          axis.text = element_text(size = 22),
          axis.title = element_text(size = 22),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(size = 26),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22,
                                      face = "bold")
    ) 

dev.off()



caller_coverage[, med_loc := start + median(c(start, end))]
caller_coverage
caller_coverage[, .SD, by = .(cell_name, seqnames)]


bin_gen_dt$bin_id


caller_coverage[, cnv_numeric := ifelse(cnv_state == "amp", 1, -1)]



ggplot(data = caller_coverage[cell_name == cellname &
    seqnames == "chr1"]) +
    geom_col(aes(
        x = med_loc,
        y = cnv_numeric,
        fill = caller
    )) +
    coord_fixed(ratio = 1000000) +
    xlab("") +
    ylab("")



# number of neutral bins per cell
all_bins <- bin_gen_dt$bin_id
all_bins
caller_coverage[, bin_id, by = cell_name]
caller_coverage[cell_name == cellname, bin_id]
caller_coverage[cell_name == cellname, length(unique(bin_id))]

all_bins[all_bins %in% caller_coverage[cell_name == cellname, bin_id]] |>
    length()
all_bins[!all_bins %in% caller_coverage[cell_name == cellname, bin_id]] |>
    length()
length(all_bins)
all_bins[!all_bins %in% caller_coverage[cell_name == cellname, bin_id]]
caller_coverage[cell_name == cellname & bin_id == 248]


# the neutral bins
neu_calls <- caller_coverage[, all_bins[!all_bins %in% bin_id],
                             by = .(cell_name,
                                    caller)]
setnames(neu_calls,
         "V1",
         "neu_bins"
)
neu_calls
neu_calls[, .SD, by = .(neu_bins, cell_name)]
neu_calls[, hits := length(unique(caller)), by = .(neu_bins, cell_name)]
neu_calls[, unique(hits)]
neu_calls
neu_calls[neu_bins == 249]
neu_calls[neu_bins == 1]
neu_calls[neu_bins == 249 &
          cell_name == neu_calls[1, cell_name]]


total_bins <- max(bin_gen_dt$bin_id)
total_bins
# prop of bins called neu by number of callers
neu_calls[, length(unique())]
neu_calls[, unique(neu_bins)]
neu_calls[, .SD, by = .(hits)]
neu_overlap <- neu_calls[, round(length(unique(neu_bins)) / total_bins * 100, 2), by = .(hits)]


pdf("../plots/cnv_callers_stats/neutral_calls_overlap_callers.pdf")

ggplot(data = neu_overlap,
       aes(x = hits,
           y = V1,
           label = V1
       )
) +
    geom_col() +
    labs(x = "Number of callers",
         y = "Proportion of all bins (%)",
         title = "Extent of overlap of neutral called bins per cell"
    ) +
    geom_text(vjust = -0.5)

dev.off()


# for two caller overlap of neu
neu_calls[hits == 2, .SD, by = .(cell_name, neu_bins)][order(caller)]
neu_calls[neu_bins == 79 &
          like(cell_name, "AGAGGA")]


neu_caller_pairs <- neu_calls[order(caller)
                              ][hits == 2, 
                              paste(caller, sep = "-", collapse = "-"),
                              by = .(cell_name, neu_bins)]
neu_caller_pairs[, unique(V1)]
neu_caller_pairs
setnames(neu_caller_pairs,
         "V1",
         "caller_pairs"
)


neu_caller_pairs
neu_caller_stats <- neu_caller_pairs[, round(length(unique(neu_bins)) / total_bins * 100, 2), by = caller_pairs]
neu_caller_stats



pdf("../plots/cnv_callers_stats/neutral_calls_caller_pairs.pdf")

ggplot(data = neu_caller_stats,
       aes(x = caller_pairs,
           y = V1,
           label = V1
       )
) +
    geom_col() +
    geom_text(vjust = -0.5) +
    labs(x = "Caller Pairs",
         y = "Proportion of all bins (%)",
         title = "Caller pairs in agreement: Neutral Bins"
    )

dev.off()


loh_calls <- caller_coverage[cnv_state == "loh"]
loh_calls


caller_coverage[cell_name == "tnbc1_CCCAATCTCCGTACAA" &
                bin_id == 1 &
                cnv_state != "loh"]


loh_comp <- data.table()
for(entry in 1:nrow(loh_calls)) {

    curr_entry <- caller_coverage[cell_name == loh_calls[entry, cell_name] &
                                  bin_id == loh_calls[entry, bin_id] &
                                  cnv_state != "loh"]
    if(nrow(curr_entry) > 0) {
        loh_comp <- rbind(loh_comp,
                          curr_entry
        )
    }
}
entry
loh_comp

loh_bins <- loh_calls[, length(unique(bin_id))]


loh_calls[bin_id == 594]
loh_comp[, nrow(.SD) > 1, by = .(bin_id)]
loh_comp[bin_id == 594]


copy <- caller_coverage
caller_coverage

loh_cells <- loh_calls[, .SD, by = .(bin_id, cell_name)][, .(bin_id, cell_name)]
loh_cells[, loh_call := 1]
loh_cells


copy <- merge(copy,
              loh_cells,
              by = c("bin_id", "cell_name")
)
copy[cnv_state != "loh", .SD, by = .(bin_id, cell_name, caller)]
copy[cnv_state != "loh", hits := length(unique(caller)), by = .(bin_id, cell_name, cnv_state)]


copy[hits == 1 &
     caller == "InferCNV"]


 # calc bin prop of loh in other calls
loh_miscalls_stats <- copy[cnv_state != "loh", round(length(unique(bin_id)) / loh_bins * 100, 2), by = .(hits, cnv_state, caller)]

loh_miscalls_stats[hits == 2, caller := "Copykat-InferCNV"]
loh_miscalls_stats <- unique(loh_miscalls_stats)
loh_miscalls_stats



pdf("../plots/cnv_callers_stats/loh_miscalled_bins.pdf")

ggplot(data = loh_miscalls_stats,
       aes(
           x = caller,
           y = V1,
           label = V1,
           fill = cnv_state
       )
       ) +
                 geom_col(position = "dodge") +
                 geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
              labs(x = "Caller",
                   y = "Proportion of loh called bins (%)",
                   title = "Proportion of bins with a cnv call for LOH calls"
              ) +
              scale_fill_manual(values = cnv_colors)

dev.off()

