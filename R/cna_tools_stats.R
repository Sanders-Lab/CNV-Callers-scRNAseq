library(data.table)
library(tidyverse)
library(ggplot2)
library(magrittr)


# reading in the segs
infercnv_segs <- fread("../proc/infercnv_tnbc1_segs_refined.tsv.gz")
infercnv_segs

infercnv_segs <- infercnv_segs[cell_name %in% copykat_segs$cell_name]


copykat_segs <- fread("../proc/copykat_tnb1_segs_refined.tsv.gz")
copykat_segs


numbat_segs <- fread("../proc/numbat_tnbc_cna_segs.tsv.gz")
numbat_segs

numbat_segs <- numbat_segs[cell_name %in% copykat_segs$cell_name]

# checking if cell order is the same
numbat_segs[,cell_name] %>% 
    unique() == unique(copykat_segs[,cell_name]) %>% 
    




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


infercnv_segs[, .(length(grep("del", cnv_state)),
                              length(grep("amp", cnv_state))), by = cell_name] %>% 
.[, c("V1", "V2")]




# calls from here
inf_calls <- infercnv_segs[, .(length(grep("del", cnv_state)),
                               length(grep("amp", cnv_state))), by = cell_name] %>% 
.[, c("V1", "V2")]
setnames(inf_calls, c("V1", "V2"), c("n_dels", "n_amps"))
inf_calls



ck_calls <- copykat_segs[, .(length(grep("del", cnv_state)),
                               length(grep("amp", cnv_state))), by = cell_name] %>% 
.[, c("V1", "V2")]
setnames(ck_calls, c("V1", "V2"), c("n_dels", "n_amps"))
ck_calls




nb_calls <- numbat_segs[, .(length(grep("del", cnv_state)),
                            length(grep("amp", cnv_state)),
                            length(grep("loh", cnv_state))), by = cell_name] %>% 
.[, c("V1", "V2", "V3")]
setnames(nb_calls, c("V1", "V2", "V3"), c("n_dels", "n_amps", "n_loh"))
nb_calls



# 2023-04-19
# rbinding calls
# next make a col and rep it 3 times
# the nrow of this rbind dt, with the names
# of the tools.
# then melt the three col into a count col
# 
all_calls <- rbind(inf_calls,
                   ck_calls,
                   nb_calls[,1:2])
all_calls
all_calls %>% 
    nrow() / 3


all_calls$tool <- "tool"
all_calls


all_calls[1:796, "tool"] <- "InferCNV"
all_calls[797:1592, "tool"] <- "Copykat"
all_calls[1593:2388, "tool"] <- "Numbat"


all_calls


plot_table <- data.table()
plot_table$tool <- all_calls$tool %>% rep(., 2)
plot_table


plot_table$count <- c(all_calls$n_amps, all_calls$n_dels)
plot_table


plot_table$cnv_call <- c(rep("n_amps", (nrow(plot_table) / 2)), rep("n_dels", (nrow(plot_table) / 2)))
plot_table


# adding numbat loh calls
loh_calls <- data.table(tool = rep("Numbat", nrow(nb_calls)),
                     count = nb_calls$n_loh,
                     cnv_call = rep("n_loh", nrow(nb_calls)))
loh_calls

plot_table <- rbind(plot_table, loh_calls)
plot_table


tmp_plot <- "../plots/tmp/tmp_plot.pdf"

plot_table_df <- as.data.frame(plot_table)
cna_colors <- c('red', 'blue', "green")
names(cna_colors) <- c('n_amps', 'n_dels', "n_loh")
cna_colors


pdf("../plots/cnv_callers_stats/cnv_counts_violin.pdf")
ggplot(data = plot_table_df,
       mapping = aes(x = cnv_call,
                     y = count,
                     colour = cnv_call)) +
       geom_violin() +
       geom_point() +
       facet_grid(~ tool) +
       scale_color_manual(values = c('red', 'blue', "green"))
dev.off()
ggsave(tmp_plot)






