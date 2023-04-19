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
                   nb_calls,
                   fill = T)
all_calls
all_calls %>% 
    nrow() / 3



setnames(all_calls,
         "V1",
         "cell_name")
all_calls

melt(all_calls[,2:ncol(all_calls)],
     id.vars = c(colnames(all_calls)[2:7]))



plot_table <- data.table()
plot_table$cell_name <- all_calls$cell_name %>% 
    rep(.,7)
plot_table


plot_table





