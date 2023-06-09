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


nb_overlaps <- findOverlaps(bin_genome,
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




inf_overlaps <- findOverlaps(bin_genome,
                            inf_gro
)
inf_overlaps
inf_bins <- as.data.table(inf_gro[inf_overlaps@to])
inf_bins[, bin_id := bin_genome[inf_overlaps@from]$bin_id]
inf_bins


ck_overlaps <- findOverlaps(bin_genome,
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


temp <- union(nb_bins[cell_name == nb_bins[1, cell_name],
        bin_id],
inf_bins[cell_name == nb_bins[1, cell_name],
        bin_id])
temp <- union(temp,
              ck_bins[cell_name == nb_bins[1, cell_name],
                      bin_id] 
)
temp

|>
    sort()



nb_bins[bin_id == 2999]
all_cells <- unique(copykat_segs$cell_name)
all_cells
bin_gen_dt <- as.data.table(bin_genome)
all_bins <- bin_gen_dt$bin_id
length(all_bins)




nb_bins[cell_name == all_cells[1],
                 which(duplicated(bin_id))]
nb_bins[cell_name == ]
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


    for(entry in 1:nrow(calls_dt)) {

        curr_bin <- calls_dt[entry, bin_id]

        # for 2nd event bigger than 1st
        if(calls_dt[entry, start] - 
           bin_gen_dt[bin_id == curr_bin, start] < 0 &
           calls_dt[entry, end] - 
           bin_gen_dt[bin_id == curr_bin, end] > 0) {

            calls_dt[entry, perc_overlap := 100]
            next
        }

        # for 2nd event right shifted
        if(calls_dt[entry, start] - 
           bin_gen_dt[bin_id == curr_bin, start] > 0) {


            calls_dt[entry, perc_overlap := ((calls_dt[entry, start] - 
                                             bin_gen_dt[bin_id == curr_bin, start]) / 1e6 * 100)]

        } 

        # for 2nd event left shifted
        if (calls_dt[entry, start] - 
           bin_gen_dt[bin_id == curr_bin, start] < 0){

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
            })
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
            })
inf_perc
parallel::stopCluster(cluster)



n_cores <- 63
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)

# for inf calls
system.time(
            ck_perc <- foreach(celln = seq_along(all_cells), .combine = "rbind") %dopar% {
                PercentOverlap(ck_bins[cell_name == all_cells[celln]])
            })
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
nb_perc[perc_overlap > 50
        ][, which(duplicated(bin_id)), by = cell_name]

inf_perc[perc_overlap > 50
        ][, which(duplicated(bin_id)), by = cell_name]

ck_perc[perc_overlap > 50
        ][, which(duplicated(bin_id)), by = cell_name]

ck_perc[cell_name == cellname,
        which(duplicated(bin_id))] 

dups <- integer()
temp <- c(ck_perc[cell_name == cellname,
                          which(duplicated(bin_id))]) 
temp
dups[cellname] <- c(temp)
    c(ck_perc[cell_name == cellname,
                          which(duplicated(bin_id))]) 
dups

ck_perc[cell_name == cellname][190:191]
ck_perc[cell_name == cellname][758:759]

ck_perc[, idx := 1:nrow(ck_perc)]
ck_perc



ck_perc[which(duplicated(bin_id)), .SD, by = cell_name]

remove_dups <- function(sub_dt) {

    for(dups in 1:nrow(sub_dt)) {
        curr_idx <- sub_dt[dups, idx]

        if(ck_perc[curr_idx, perc_overlap] > 
           ck_perc[(curr_idx - 1), perc_overlap]) {
            ck_perc <- ck_perc[!(curr_idx - 1)]
            next
        } 
        if (ck_perc[curr_idx, perc_overlap] < 
           ck_perc[(curr_idx - 1), perc_overlap]){
            ck_perc <- ck_perc[!curr_idx]
        }
        if(ck_perc[curr_idx, perc_overlap] == 
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
ck_perc[perc_overlap > 50
        ][which(duplicated(bin_id)), .SD, by = cell_name]

ck_perc[1338:1340, bin_id]



ck_dedup[cell_name == "tnbc1_AAACGGGTCCAGAGGA",
         which(duplicated(bin_id))]
ck_dedup[cell_name == "tnbc1_AAACGGGTCCAGAGGA",
         ][100:102]
ck_dedup[cell_name == "tnbc1_AAACGGGTCCAGAGGA",
         ][313:315]


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

        temp <- data.table(bin_id = integer(1),
                           cell_name = character(1),
                           caller = character(1),
                           cnv_state = character(1)
        )

        if(bin %in% call_dt[, bin_id]) {

            temp[1, bin_id := bin]
            temp[1, cell_name := (call_dt[1, cell_name])]
            temp[1, caller :=  caller]
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


comp_dt <- nb_perc[, .(cell_name, 
                       bin_id,
                       cnv_state,
                       perc_overlap
)
][, caller := "Numbat"]
comp_dt

comp_dt <- rbind(comp_dt,
                 inf_perc[, .(cell_name, 
                             bin_id,
                             cnv_state,
                             perc_overlap
                 )
                 ][, caller := "InferCNV"]
)



comp_dt <- rbind(comp_dt,
                 ck_perc[, .(cell_name, 
                             bin_id,
                             cnv_state,
                             perc_overlap
                 )
                 ][, caller := "Copykat"]
)

comp_dt

comp_dt[, cnv_numeric := ifelse(cnv_state == "amp", 1, -1)]
comp_dt[cnv_state == "amp"]


tmp_plot <- "../plots/tmp/tmp_plot.pdf"
pdf("../plots/cnv_callers_stats/cell1_binned_calls_callers.pdf")
ggplot(data = comp_dt[cell_name == cellname]) +
    geom_col(aes(x = bin_id,
                 y = cnv_numeric,
                 fill = caller
                 ))
dev.off()

ggsave(tmp_plot)



bin_gen_dt[,median := start + median(c(start,end))]
bin_gen_dt


sapply(bin_gen_dt,
       median()
)


comp_dt <- merge(comp_dt, bin_gen_dt, by = "bin_id")



ggplot(data = comp_dt[cell_name == cellname &
       seqnames == "chr1"]) +
    geom_col(aes(x = median,
                 y = cnv_numeric,
                 fill = caller
                 )) +
       coord_fixed(ratio = 1000000) +
       xlab("") +
       ylab("")
ggsave(tmp_plot)




