library(devtools)


# installing primatR
install_github("daewoooo/primatR", force=TRUE)




library(primatR)
library(data.table)
library(magrittr)
library(parallel)
library(foreach)



ck_calls <- fread("../proc/copykat_tnb1_segs_refined.tsv.gz")
ck_calls
ck_calls$cell_name %>% 
    unique() %>% 
    length()


inf_calls <- fread("../proc/infercnv_tnbc1_segs_refined.tsv.gz")
inf_calls

inf_calls <- inf_calls[cell_name %in% ck_calls$cell_name]


inf_calls$cell_name %>% 
    unique() %>% 
    length()


nb_calls <- fread("../proc/numbat_tnbc_cna_segs.tsv.gz")
nb_calls


nb_calls <- nb_calls[cnv_state_map != ""]


nb_calls <- nb_calls[cell_name %in% ck_calls$cell_name]
nb_calls$cell_name %>% 
    unique() %>% 
    length()


inf_gro <- makeGRangesFromDataFrame(inf_calls,
                                    keep.extra.columns = T)
inf_gro


ck_gro <- makeGRangesFromDataFrame(ck_calls,
                                   keep.extra.columns = T)
ck_gro


nb_gro <- makeGRangesFromDataFrame(nb_calls,
                                   keep.extra.columns = T)
nb_gro



# test with one cell
one_cell_inf <- inf_calls[cell_name == inf_calls$cell_name[1]]
one_cell_inf

one_cell_inf_gro <- makeGRangesFromDataFrame(one_cell_inf,
                                             keep.extra.columns = T)
one_cell_inf_gro


one_cell_nb <- nb_calls[cell_name == inf_calls$cell_name[1]]
one_cell_nb


one_cell_nb_gro <- makeGRangesFromDataFrame(one_cell_nb,
                                            keep.extra.columns = T)
one_cell_nb_gro


# system.time(inf_nb_perc_overlap <- getRangesOverlaps(one_cell_inf_gro,
#                                                      one_cell_nb_gro,
#                                                      "inf"))
system.time(nb_inf_perc_overlap <- getRangesOverlaps(one_cell_nb_gro,
                                                     one_cell_inf_gro,
                                                     index = "idx"))
as.data.table(nb_inf_perc_overlap)
t <- as.data.table(getSegDupOverlaps(one_cell_nb_gro,
                                one_cell_inf_gro))
t$tool <- "inf"
t



nb_ck_perc_overlap <- getSegDupOverlaps(one_cell_nb_gro,
nb_ck_perc_overlap <- getSegDupOverlaps(one_cell_nb_gro,

inf_nb_perc_overlap
nb_inf_perc_overlap
nb_inf_perc_overlap[SDTrackPercOverlap >= 50]
overlap_dt <- as.data.table(nb_inf_perc_overlap)


# to get the calls above 50% overlap
overlap_dt
overlap_dt[, unique(SDTrackPercOverlap)]
overlap_dt[SDTrackPercOverlap >= 50, .N]


# LOGIC BEGINS HERE
# 2023-04-24 15:58
# using the lcd of the cells
all_cells <- ck_calls[, unique(cell_name)]
all_cells
tools
tools[1]



# REF = numbat
n_cores <- 63
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


# to get the overlaps with numbat
# as the ref
# tmp_dt <- data.table()
tools <- c("nb_calls",
           "inf_calls",
           "ck_calls")
nb_overlaps <- data.table()
system.time(
            nb_overlaps <- foreach(curr_tool = 2:length(tools),
                                   .combine = "rbind") %:%

            foreach(cell_n = 1:length(all_cells),
                    .combine = "rbind") %dopar%
                {

                    as.data.table(getSegDupOverlaps(makeGRangesFromDataFrame(get(tools[1])[cell_name == all_cells[cell_n]],
                                                                             keep.extra.columns = T),
                                                    makeGRangesFromDataFrame(get(tools[curr_tool])[cell_name == all_cells[cell_n]],
                                                                             keep.extra.columns = T)
                    )
                    ) %>%
                    .[, tool := tools[curr_tool]]
                }
)
parallel::stopCluster(cluster)
nb_overlaps
nb_overlaps
nb_overlaps[, unique(tool)]
nb_overlaps[, length(unique(cell_name))]


nb_overlaps[, unique(.N), by = cell_name]


# checking the num of segs with >= 50 overlap
nb_overlaps[SDTrackPercOverlap >= 50, .(overlaps = .N), by = cell_name]
nb_overlaps[SDTrackPercOverlap >= 50, .N, by = cell_name
              ][, mean(N)]
nb_overlaps[SDTrackPercOverlap >= 50, .N, by = cell_name
              ][, unique(N)]


# removing pointless cols
nb_overlaps[, ":="(width = NULL,
                   strand = NULL)]
nb_overlaps


fwrite(x = nb_overlaps,
       file = "../proc/numbat_other_tools_overlap.tsv.gz",
       sep = "\t")




# ref = INF
n_cores <- 63
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


tools <- c("inf_calls",
           "nb_calls",
           "ck_calls")
inf_overlaps <- data.table()
system.time(
            inf_overlaps <- foreach(curr_tool = 1:length(tools),
                                    .combine = "rbind") %:%

            foreach(cell_n = 1:length(all_cells),
                    .combine = "rbind") %dopar%
                {

                    as.data.table(getSegDupOverlaps(makeGRangesFromDataFrame(get(tools[1])[cell_name == all_cells[cell_n]],
                                                                             keep.extra.columns = T),
                                                    makeGRangesFromDataFrame(get(tools[curr_tool])[cell_name == all_cells[cell_n]],
                                                                             keep.extra.columns = T)
                    )
                    ) %>%
                    .[, tool := tools[curr_tool]]

                #                     tmp_dt$tool <- tools[curr_tool]
                }
)
parallel::stopCluster(cluster)
inf_overlaps


# removing pointless cols
inf_overlaps[, ":="(width = NULL,
                    strand = NULL)]


fwrite(x = inf_overlaps,
       file = "../proc/infercnv_other_tools_overlap.tsv.gz",
       sep = "\t")






# ref = CK
n_cores <- 63
cluster <- makeForkCluster(n_cores)
doParallel::registerDoParallel(cluster)


tools <- c("ck_calls",
           "nb_calls",
           "inf_calls")
ck_overlaps <- data.table()
system.time(
            ck_overlaps <- foreach(curr_tool = 1:length(tools),
                                   .combine = "rbind") %:%

            foreach(cell_n = 1:length(all_cells),
                    .combine = "rbind") %dopar%
                {

                    as.data.table(getSegDupOverlaps(makeGRangesFromDataFrame(get(tools[1])[cell_name == all_cells[cell_n]],
                                                                             keep.extra.columns = T),
                                                    makeGRangesFromDataFrame(get(tools[curr_tool])[cell_name == all_cells[cell_n]],
                                                                             keep.extra.columns = T)
                    )
                    ) %>%
                    .[, tool := tools[curr_tool]]

                #                     tmp_dt$tool <- tools[curr_tool]
                }
)
parallel::stopCluster(cluster)
ck_overlaps


# removing pointless cols
ck_overlaps[, ":="(width = NULL,
                   strand = NULL)]


fwrite(x = ck_overlaps,
       file = "../proc/copykat_other_tools_overlap.tsv.gz",
       sep = "\t")

