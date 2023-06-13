library(data.table)
# library(tidyverse)
library(ggplot2)
# library(magrittr)


# reading in the segs
copykat_segs <- fread("../proc/copykat_tnb1_segs_refined.tsv.gz")
copykat_segs
nsegs_ck <- nrow(copykat_segs)


infercnv_segs <- fread("../proc/infercnv_tnbc1_segs_refined.tsv.gz")
infercnv_segs

infercnv_segs <- infercnv_segs[cell_name %in% copykat_segs$cell_name]
infercnv_segs
nsegs_inf <- nrow(infercnv_segs)


numbat_segs <- fread("../proc/numbat_tnbc_cna_segs.tsv.gz")
numbat_segs <- numbat_segs[cell_name %in% copykat_segs$cell_name]
numbat_segs <- numbat_segs[cnv_state != ""]
numbat_segs
nsegs_nb <- nrow(numbat_segs)





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

nb_sig_overlaps[, N := .N, by = cnv_state]
nb_sig_overlaps[, .SD, by = .(seqnames, start, end)]
nb_sig_overlaps[, mean_width := mean(width), by = cnv_state]
nb_sig_overlaps


nb_overlap_stats <- nb_sig_overlaps[cnv_state != "loh", .SD, by = .(cell_name, caller, cnv_state)
                                    ][, unique(round(.N / (N + log(mean_width)) * 100, 2)), by = .(cnv_state, caller)]
nb_overlap_stats

cnv_colors <- c('tomato2', 'dodgerblue2', "springgreen3")
names(cnv_colors) <- c('amp', 'del', "loh")
cnv_colors


# plotting
pdf("../plots/cnv_callers_stats/cnv_overlap_numbat_count.pdf")

ggplot(data = nb_overlap_stats,
       mapping = aes(x = caller,
                     y = V1,
                     fill = cnv_state,
                     label = V1)) +
geom_col(position = "dodge") +
geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
scale_fill_manual(values = cnv_colors) +
labs(y = "Proportion of Segments (%)",
     title = "Numbat - Count of Overlapping (> 50% overlap) CNV Segments"
)

dev.off()

# ggsave(tmp_plot)



# plotting the loh calls
nb_sig_overlaps[cnv_state == "loh"]
nb_overlap_stats <- nb_sig_overlaps[cnv_state == "loh", .SD, by = .(cell_name, caller, cnv_state)
                                    ][, unique(round(.N / (N + log(mean_width)) * 100, 2)), by = .(matched_call, caller)]
nb_overlap_stats


pdf("../plots/cnv_callers_stats/cnv_overlap_loh_count.pdf")

ggplot(data = nb_overlap_stats,
       mapping = aes(x = caller,
                     y = V1,
                     fill = matched_call,
                     label = V1)) +
geom_col(position = "dodge") +
geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
scale_fill_manual(values = cnv_colors) +
labs(y = "Proportion of Segments (%)",
     title = "Numbat - Count of Overlapping loh (> 50% overlap) CNV Segments"
)

dev.off()

ggsave(tmp_plot)


# INFERCNV
inf_overlaps
inf_overlaps[, width := (end - start)]
inf_sig_overlaps <- inf_overlaps[SDTrackPercOverlap > 50]
inf_sig_overlaps

inf_sig_overlaps[, N := .N, by = cnv_state]
inf_sig_overlaps[, mean_width := mean(width), by = cnv_state]
inf_sig_overlaps


inf_overlap_stats <- inf_sig_overlaps[, .SD, by = .(cell_name, caller, cnv_state)
                                    ][, unique(round(.N / (N + log(mean_width)) * 100, 2)), by = .(cnv_state, caller)]
inf_overlap_stats




# plotting infercnv overlaps
pdf("../plots/cnv_callers_stats/cnv_overlap_infercnv_count.pdf")

ggplot(data = inf_overlap_stats,
       mapping = aes(x = caller,
                     y = V1,
                     fill = cnv_state,
                     label = V1)) +
geom_col(position = "dodge") +
geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
scale_fill_manual(values = cnv_colors) +
labs(y = "Proportion of Segments (%)",
     title = "InferCNV - Count of Overlapping loh (> 50% overlap) CNV Segments"
)

dev.off()


# COPYKAT
ck_overlaps
ck_overlaps[, width := (end - start)]
ck_sig_overlaps <- ck_overlaps[SDTrackPercOverlap > 50]
ck_sig_overlaps

ck_sig_overlaps[, N := .N, by = cnv_state]
ck_sig_overlaps[, mean_width := mean(width), by = cnv_state]
ck_sig_overlaps


ck_overlap_stats <- ck_sig_overlaps[, .SD, by = .(cell_name, caller, cnv_state)
                                    ][, unique(round(.N / (N + log(mean_width)) * 100, 2)), by = .(cnv_state, caller)]
ck_overlap_stats


# plotting infercnv overlaps
pdf("../plots/cnv_callers_stats/cnv_overlap_copykat_count.pdf")

ggplot(data = ck_overlap_stats,
       mapping = aes(x = caller,
                     y = V1,
                     fill = cnv_state,
                     label = V1)) +
geom_col(position = "dodge") +
geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
scale_fill_manual(values = cnv_colors) +
labs(y = "Proportion of Segments (%)",
     title = "Copykat - Count of Overlapping loh (> 50% overlap) CNV Segments"
)

dev.off()


# proportion of calls dropped
prop_drop <- data.table(caller = character(7),
                        prop = numeric(7))



prop_drop[1, prop := round((nrow(ck_overlaps[cnv_state == "amp"]) - nrow(ck_sig_overlaps[cnv_state == "amp"]))
                       / nrow(ck_overlaps[cnv_state == "amp"]) * 100, 2)
          ][1, caller := "copykat"
          ][1, cnv_call := "amp"]
prop_drop[2, prop := round((nrow(ck_overlaps[cnv_state == "del"]) - nrow(ck_sig_overlaps[cnv_state == "del"]))
                       / nrow(ck_overlaps[cnv_state == "del"]) * 100, 2)
          ][2, caller := "copykat"
          ][2, cnv_call := "del"]
prop_drop[3, prop := round((nrow(inf_overlaps[cnv_state == "amp"]) - nrow(inf_sig_overlaps[cnv_state == "amp"]))
                       / nrow(inf_overlaps[cnv_state == "amp"]) * 100, 2)
          ][3, caller := "infercnv"
          ][3, cnv_call := "amp"]
prop_drop[4, prop := round((nrow(ck_overlaps[cnv_state == "del"]) - nrow(ck_sig_overlaps[cnv_state == "del"]))
                       / nrow(ck_overlaps[cnv_state == "del"]) * 100, 2)
          ][4, caller := "infercnv"
          ][4, cnv_call := "del"]
prop_drop[5, prop := round((nrow(nb_overlaps[cnv_state == "amp"]) - nrow(nb_sig_overlaps[cnv_state == "amp"]))
                       / nrow(ck_overlaps[cnv_state == "amp"]) * 100, 2)
          ][5, caller := "numbat"
          ][5, cnv_call := "amp"]
prop_drop[6, prop := round((nrow(nb_overlaps[cnv_state == "del"]) - nrow(nb_sig_overlaps[cnv_state == "del"]))
                       / nrow(nb_overlaps[cnv_state == "del"]) * 100, 2)
          ][6, caller := "numbat"
          ][6, cnv_call := "del"]
prop_drop[7, prop := round((nrow(nb_overlaps[cnv_state == "del"]) - nrow(nb_sig_overlaps[cnv_state == "del"]))
                       / nrow(nb_overlaps[cnv_state == "del"]) * 100, 2)
          ][7, caller := "numbat"
          ][7, cnv_call := "del"]


prop_drop


# plotting prop_drop
pdf("../plots/cnv_callers_stats/cnv_prop_drop_calls.pdf")

ggplot(data = prop_drop,
       mapping = aes(x = caller,
                     y = prop,
                     fill = cnv_call,
                     label = prop)) + 
geom_col(position = "dodge") +
geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
scale_fill_manual(values = cnv_colors) +
labs(y = "Proportion of calls (%)",
     title = "Proportions of calls retained after (> 50%) overlap"
)

dev.off()



# n_genes vs n_segs
numbat_segs
genes_segs <- data.table()
genes_segs <- rbind(genes_segs,
                    infercnv_segs[, .(n_genes = sum(n_genes)), by = .(cell_name, cnv_state)
                                  ][, caller := "infercnv"],
                    copykat_segs[, .(n_genes = sum(n_genes)), by = .(cell_name, cnv_state)
                                 ][, caller := "copykat"])
genes_segs

numbat_genes_segs <- numbat_segs[, .(n_genes = sum(n_genes)), by = .(cell_name, cnv_state)
                                 ][, caller := "numbat"]
nb_segs_tmp <- numbat_segs[, .(n_segs = .N), by = .(cell_name, cnv_state)
                                ][, caller := "infercnv"]
numbat_genes_segs
numbat_genes_segs <- cbind(numbat_genes_segs, nb_segs_tmp$n_segs)
setnames(numbat_genes_segs, "V2", "n_segs")

segs_tmp <- data.table()
segs_tmp <- rbind(segs_tmp,
                  infercnv_segs[, .(n_segs = .N), by = .(cell_name, cnv_state)
                                ][, caller := "infercnv"],
                  copykat_segs[, .(n_segs = .N), by = .(cell_name, cnv_state)
                               ][, caller := "copykat"])
segs_tmp

genes_segs <- cbind(genes_segs, segs_tmp$n_segs)
setnames(genes_segs, "V2", "n_segs")
genes_segs

ck_genes_segs <- genes_segs[caller == "copykat"]
inf_genes_segs <- genes_segs[caller == "infercnv"]

#     infercnv_segs[, .(n_genes = sum(n_genes)), by = .(cell_name, cnv_state)
                                ][, caller := "infercnv"],


# plotting
pdf("../plots/cnv_callers_stats/cnv_copykat_corr_genes_segs.pdf")
ggplot(data = ck_genes_segs,
       mapping = aes(x = n_segs,
                     y = n_genes,
                     colour = cnv_state)) +
geom_point() +
geom_jitter() +
geom_smooth(method = "lm",
            show.legend = T) +
scale_color_manual(values = cnv_colors) +
ggtitle("Copykat - Correlation between n_genes and n_segs")
dev.off()

ggsave(tmp_plot)


pdf("../plots/cnv_callers_stats/cnv_infercnv_corr_genes_segs.pdf")
ggplot(data = inf_genes_segs,
       mapping = aes(x = n_segs,
                     y = n_genes,
                     colour = cnv_state)) +
geom_point() +
geom_jitter(width = 10, height = 10) +
geom_smooth(method = "lm",
            show.legend = T) +
scale_color_manual(values = cnv_colors) +
ggtitle("Infercnv - Correlation between n_genes and n_segs")
dev.off()


pdf("../plots/cnv_callers_stats/cnv_numbat_corr_genes_segs.pdf")
ggplot(data = numbat_genes_segs,
       mapping = aes(x = n_segs,
                     y = n_genes,
                     colour = cnv_state)) +
geom_point() +
geom_jitter() +
geom_smooth(method = "lm",
            show.legend = T) +
scale_color_manual(values = cnv_colors) +
ggtitle("Numbat - Correlation between n_genes and n_segs")
dev.off()
