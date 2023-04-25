library(data.table)
library(tidyverse)
library(ggplot2)
library(magrittr)


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
num_calls <- data.table(cell_name = character(1),
                        n_amps = integer(1),
                        n_dels = integer(1),
                        tool = character(1),
                        cnv_call = character(1))




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

cna_colors <- c('tomato2', 'dodgerblue2', "springgreen3")
names(cna_colors) <- c('amp', 'del', "loh")
cna_colors


pdf("../plots/cnv_callers_stats/cnv_counts_violin.pdf")
ggplot(data = count_dt,
       mapping = aes(x = cnv_call,
                     y = count,
                     colour = cnv_call)) +
       geom_violin() +
       geom_point() +
       geom_jitter() +
       facet_grid(~ caller) +
       scale_color_manual(values = cna_colors)
dev.off()
ggsave(tmp_plot)





