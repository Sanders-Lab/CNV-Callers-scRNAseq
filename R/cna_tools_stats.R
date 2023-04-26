library(data.table)
# library(tidyverse)
library(ggplot2)
# library(magrittr)


# reading in the segs
copykat_segs <- fread("../proc/copykat_tnb1_segs_refined.tsv.gz")
copykat_segs


infercnv_segs <- fread("../proc/infercnv_tnbc1_segs_refined.tsv.gz")
infercnv_segs

infercnv_segs <- infercnv_segs[cell_name %in% copykat_segs$cell_name]


numbat_segs <- fread("../proc/numbat_tnbc_cna_segs.tsv.gz")
numbat_segs <- numbat_segs[cell_name %in% copykat_segs$cell_name]
setnames(numbat_segs,
         "cnv_state_map",
         "cnv_state")
numbat_segs <- numbat_segs[cnv_state != ""]
numbat_segs






# number of amps and dels per cell
# num_calls <- data.table(cell_name = character(1),
#                         n_amps = integer(1),
#                         n_dels = integer(1),
#                         tool = character(1),
#                         cnv_call = character(1))
# 



infercnv_segs_df <- as.data.frame(infercnv_segs)

numbat_segs[numbat_segs$cell_name == numbat_segs$cell_name[200],] %>% 
    .$cnv_state %>% 
    grep("del", .) %>%
    length(.)


infercnv_segs[, .(count = length(grep("del", cnv_state)),
                  call = "del"), by = cell_name]
infercnv_segs



# calls from here
# INFERCNV
count_dt <- infercnv_segs[, .(count = length(grep("del", cnv_state)),
                              cnv_call = "del",
                              caller = "infercnv"), by = cell_name]
count_dt <- rbind(count_dt,
                  infercnv_segs[, .(count = length(grep("amp", cnv_state)),
                                    cnv_call = "amp",
                                    caller = "infercnv"), by = cell_name])
count_dt



# COPYKAT
count_dt <- rbind(count_dt,
                  copykat_segs[, .(count = length(grep("del", cnv_state)),
                                   cnv_call = "del",
                                   caller = "copykat"), by = cell_name])
count_dt <- rbind(count_dt,
                  copykat_segs[, .(count = length(grep("amp", cnv_state)),
                                    cnv_call = "amp",
                                    caller = "copykat"), by = cell_name])
count_dt


# NUMBAT
count_dt <- rbind(count_dt,
                  numbat_segs[, .(count = length(grep("del", cnv_state)),
                                  cnv_call = "del",
                                  caller = "numbat"), by = cell_name])
count_dt <- rbind(count_dt,
                  numbat_segs[, .(count = length(grep("amp", cnv_state)),
                                  cnv_call = "amp",
                                  caller = "numbat"), by = cell_name])
count_dt <- rbind(count_dt,
                  numbat_segs[, .(count = length(grep("loh", cnv_state)),
                                  cnv_call = "loh",
                                  caller = "numbat"), by = cell_name])
count_dt


# Plotting
tmp_plot <- "../plots/tmp/tmp_plot.pdf"

cnv_colors <- c('tomato2', 'dodgerblue2', "springgreen3")
names(cnv_colors) <- c('amp', 'del', "loh")
cnv_colors


pdf("../plots/cnv_callers_stats/cnv_counts_violin.pdf")
ggplot(data = count_dt,
       mapping = aes(x = cnv_call,
                     y = count,
                     colour = cnv_call)) +
geom_violin() +
geom_point() +
geom_jitter() +
stat_summary(fun = "mean",
             geom = "crossbar",
             width = 1,
             colour = "black") +
facet_grid(~ caller) +
scale_color_manual(values = cnv_colors) +
ggtitle("CNV Segment Counts")
dev.off()

# ggsave(tmp_plot)




# mean count of segs
mean_count <- count_dt[, .(mean_count = mean(count)), by = .(caller, cnv_call)]
mean_count


# plotting
tmp_plot <- "../plots/tmp/tmp_plot.pdf"

cnv_colors <- c('tomato2', 'dodgerblue2', "springgreen3")
names(cnv_colors) <- c('amp', 'del', "loh")
cnv_colors


pdf("../plots/cnv_callers_stats/cnv_mean_counts_bar.pdf")
ggplot(data = mean_count,
       mapping = aes(x = cnv_call,
                     y = mean_count,
                     fill = cnv_call)) + 
geom_col() +
facet_grid(~ caller) +
scale_fill_manual(values = cnv_colors) +
    ggtitle("Mean Counts of CNV Segments")
dev.off()
# ggsave(tmp_plot)


# mean length
numbat_segs
numbat_segs[, width := (end - start)]
infercnv_segs[, width := (end - start)]
copykat_segs[, width := (end - start)]

len <- numbat_segs[, .(width, cnv_call = cnv_state)
                      ][, caller := "numbat"]
len <- rbind(len, infercnv_segs[, .(width, cnv_call = cnv_state)
             ][, caller := "infercnv"])
len <- rbind(len, copykat_segs[, .(width, cnv_call = cnv_state)
             ][, caller := "copykat"])
len

mean_len <- len[, .(mean_length = mean(width)), by = .(caller, cnv_call)]
mean_len


# plotting
pdf("../plots/cnv_callers_stats/cnv_mean_length.pdf")
ggplot(data = mean_len,
       mapping = aes(x = cnv_call,
                     y = mean_length,
                     fill = cnv_call)) + 
geom_col() +
facet_grid(~ caller) +
scale_fill_manual(values = cnv_colors) +
ggtitle("CNV Segments Mean Lengths")

dev.off()

# ggsave(tmp_plot)


# number of overlaps
nb_overlaps <- fread("../proc/numbat_other_tools_overlap.tsv.gz")
nb_overlaps

inf_overlaps <- fread("../proc/infercnv_other_tools_overlap.tsv.gz")
inf_overlaps

ck_overlaps <- fread("../proc/copykat_other_tools_overlap.tsv.gz")
ck_overlaps


# NUMBAT
nb_overlaps
nb_sig_overlaps <- nb_overlaps[SDTrackPercOverlap > 50]
nb_sig_overlaps


nb_overlap_count <- nb_sig_overlaps[cnv_state == "amp",
                .(count = .N),
                by = .(cell_name, caller)
                ][, cnv_call := "amp"] 
nb_overlap_count

nb_overlap_count <- rbind(nb_overlap_count, 
                          nb_sig_overlaps[cnv_state == "del",
                                          .(count = .N),
                                          by = .(cell_name, caller)
                                          ][, cnv_call := "del"]) 
nb_overlap_count


loh_overlaps <- nb_sig_overlaps[cnv_state == "loh",
                                .(.N, matched_call),
                                by = .(cell_name, caller)
                                ][, cnv_call := "loh"]
setnames(loh_overlaps, "N", "count")
print(loh_overlaps)


# plotting
pdf("../plots/cnv_callers_stats/cnv_overlap_numbat_count.pdf")
ggplot(data = nb_overlap_count,
       mapping = aes(x = cnv_call,
                     y = count,
                     colour = cnv_call)) +
geom_violin() +
geom_point() +
geom_jitter() +
stat_summary(fun = "mean",
             geom = "crossbar",
             width = 1,
             colour = "black") +
facet_grid(~ caller) +
scale_colour_manual(values = cnv_colors) +
ggtitle("Numbat - Count of Overlapping (50% overlap) CNV Segments")
dev.off()

# ggsave(tmp_plot)



# plotting the loh calls
pdf("../plots/cnv_callers_stats/cnv_overlap_loh_count.pdf")
ggplot(data = loh_overlaps,
       mapping = aes(x = matched_call,
                     y = count,
                     color = matched_call)) +
geom_violin() +
geom_point() +
geom_jitter() +
stat_summary(fun = "mean",
             geom = "crossbar",
             width = 1,
             colour = "black") +
facet_grid(~ caller) +
scale_colour_manual(values = cnv_colors) +
ggtitle("Numbat - Count of Overlapping (50% overlap) loh Segments")
dev.off()

ggsave(tmp_plot)


# INFERCNV
inf_overlaps
inf_sig_overlaps <- inf_overlaps[SDTrackPercOverlap > 50]
inf_sig_overlaps

inf_overlap_count <- inf_sig_overlaps[cnv_state == "amp",
                .(count = .N),
                by = .(cell_name, caller)
                ][, cnv_call := "amp"] 
inf_overlap_count

inf_overlap_count <- rbind(inf_overlap_count, 
                          inf_sig_overlaps[cnv_state == "del",
                                          .(count = .N),
                                          by = .(cell_name, caller)
                                          ][, cnv_call := "del"]) 
inf_overlap_count


# plotting infercnv overlaps
pdf("../plots/cnv_callers_stats/cnv_overlap_infercnv_count.pdf")
ggplot(data = inf_overlap_count,
       mapping = aes(x = cnv_call,
                     y = count,
                     colour = cnv_call)) +
geom_violin() +
geom_point() +
geom_jitter() +
stat_summary(fun = "mean",
             geom = "crossbar",
             width = 1,
             colour = "black") +
facet_grid(~ caller) +
scale_colour_manual(values = cnv_colors) +
ggtitle("InferCNV - Count of Overlapping (50% overlap) CNV Segments")
dev.off()


# COPYKAT
ck_overlaps
ck_sig_overlaps <- ck_overlaps[SDTrackPercOverlap > 50]
ck_sig_overlaps

ck_overlap_count <- ck_sig_overlaps[cnv_state == "amp",
                .(count = .N),
                by = .(cell_name, caller)
                ][, cnv_call := "amp"] 
ck_overlap_count

ck_overlap_count <- rbind(ck_overlap_count, 
                          ck_sig_overlaps[cnv_state == "del",
                                          .(count = .N),
                                          by = .(cell_name, caller)
                                          ][, cnv_call := "del"]) 
ck_overlap_count


# plotting infercnv overlaps
pdf("../plots/cnv_callers_stats/cnv_overlap_copykat_count.pdf")
ggplot(data = ck_overlap_count,
       mapping = aes(x = cnv_call,
                     y = count,
                     colour = cnv_call)) +
geom_violin() +
geom_point() +
geom_jitter() +
stat_summary(fun = "mean",
             geom = "crossbar",
             width = 1,
             colour = "black") +
facet_grid(~ caller) +
scale_colour_manual(values = cnv_colors) +
ggtitle("Copykat - Count of Overlapping (50% overlap) CNV Segments")
dev.off()


# proportion of calls dropped
prop_drop <- data.table(caller = character(7),
                        prop = numeric(7))
prop_drop[1, prop := ((nrow(ck_overlaps[cnv_state == "amp"]) - nrow(ck_sig_overlaps[cnv_state == "amp"]))
                       / nrow(ck_overlaps[cnv_state == "amp"]))
          ][1, caller := "copykat"
          ][1, cnv_call := "amp"]
prop_drop[2, prop := ((nrow(ck_overlaps[cnv_state == "del"]) - nrow(ck_sig_overlaps[cnv_state == "del"]))
                       / nrow(ck_overlaps[cnv_state == "del"]))
          ][2, caller := "copykat"
          ][2, cnv_call := "del"]
prop_drop[3, prop := ((nrow(inf_overlaps[cnv_state == "amp"]) - nrow(inf_sig_overlaps[cnv_state == "amp"]))
                       / nrow(inf_overlaps[cnv_state == "amp"]))
          ][3, caller := "infercnv"
          ][3, cnv_call := "amp"]
prop_drop[4, prop := ((nrow(ck_overlaps[cnv_state == "del"]) - nrow(ck_sig_overlaps[cnv_state == "del"]))
                       / nrow(ck_overlaps[cnv_state == "del"]))
          ][4, caller := "infercnv"
          ][4, cnv_call := "del"]
prop_drop[5, prop := ((nrow(nb_overlaps[cnv_state == "amp"]) - nrow(nb_sig_overlaps[cnv_state == "amp"]))
                       / nrow(ck_overlaps[cnv_state == "amp"]))
          ][5, caller := "numbat"
          ][5, cnv_call := "amp"]
prop_drop[6, prop := ((nrow(nb_overlaps[cnv_state == "del"]) - nrow(nb_sig_overlaps[cnv_state == "del"]))
                       / nrow(nb_overlaps[cnv_state == "del"]))
          ][6, caller := "numbat"
          ][6, cnv_call := "del"]
prop_drop[7, prop := ((nrow(nb_overlaps[cnv_state == "del"]) - nrow(nb_sig_overlaps[cnv_state == "del"]))
                       / nrow(nb_overlaps[cnv_state == "del"]))
          ][7, caller := "numbat"
          ][7, cnv_call := "del"]


prop_drop


# plotting prop_drop
pdf("../plots/cnv_callers_stats/cnv_prop_drop_calls.pdf")
ggplot(data = prop_drop,
       mapping = aes(x = cnv_call,
                     y = prop,
                     fill = cnv_call)) + 
geom_col() +
facet_grid(~ caller) +
scale_fill_manual(values = cnv_colors) +
ggtitle("Proportions of calls retained after 50% overlap")

dev.off()


