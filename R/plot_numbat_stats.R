library(tidyverse)
library(patchwork)

joint_post <- read.table('../outputs/numbat_pbmc3k_epithelium_tnbc1/joint_post_2.tsv',
                         header = T,
                         sep = '\t')
joint_post %>% head()


tnbc_cells <- unique(joint_post$cell)


# cell1 <- joint_post %>% filter(cell == tnbc_cells[1])
# cell1 %>% head()
# cell1 %>% nrow()


# counting CNVs: logic here
# doing for one cell
amp_count <- 0
del_count <- 0
loh_count <- 0
neu_count <- 0
# for (i in 1:nrow(cell1)){
#
#   if(cell1$cnv_state[i] == 'amp'){
#     amp_count <- amp_count + 1
#   }
#
#   if(cell1$cnv_state[i] == 'del'){
#     del_count <- del_count + 1
#   }
#
#   if(cell1$cnv_state[i] == 'loh'){
#     loh_count <- loh_count + 1
#   }
#
#   if(cell1$cnv_state[i] == 'neu'){
#     neu_count <- neu_count + 1
#   }
# }


# doing for all the cells

cnv_count_df <- data.frame()
for (cell_num in 1:length(tnbc_cells)){

  # resetting the counts
  amp_count <- 0
  del_count <- 0
  loh_count <- 0
  neu_count <- 0

  # getting the curr_cell
  curr_cell <- joint_post %>% filter(cell == tnbc_cells[cell_num])

  # looping from the 2nd cell
  for (i in 1:nrow(curr_cell)){

    amp_count <- str_count(curr_cell$cnv_state_map, 'amp') %>% sum()

    del_count <- str_count(curr_cell$cnv_state_map, 'del') %>% sum()

    loh_count <- str_count(curr_cell$cnv_state_map, 'loh') %>% sum()

    neu_count <- str_count(curr_cell$cnv_state_map, 'neu') %>% sum()
  }

  # saving a temp df with the same col format
  # will use this to rbind to the master df
  temp_df <- data.frame(amp_count,
                        del_count,
                        neu_count,
                        loh_count)
  cnv_count_df <- rbind(cnv_count_df, temp_df)
}

# ggplot() +
#   geom_bar(
#     aes(x = amp_count),
#     data = cnv_count_df
#   )
#
#
# ggplot(aes(x = 'amp', y = amp_count)) +
#   geom_violin(data = cnv_count_df) +
#   geom_point()


# # these plots worked
# (amp_plot <- ggplot(cnv_count_df, aes('',amp_count)) +
#   geom_violin() +
#   geom_point() +
#   ylab('Count') +
#   xlab('amp') +
#   geom_jitter())
#
#
# (del_plot <- ggplot(cnv_count_df, aes('',del_count)) +
#     geom_violin() +
#     geom_point() +
#     ylab('') +
#     xlab('del') +
#     geom_jitter())
#
#
# (loh_plot <- ggplot(cnv_count_df, aes('',loh_count)) +
#     geom_violin() +
#     geom_point() +
#     ylab('') +
#     xlab('loh') +
#     geom_jitter())
#
#
# (neu_plot <- ggplot(cnv_count_df, aes('',neu_count)) +
#     geom_violin() +
#     geom_point() +
#     ylab('') +
#     xlab('neu') +
#     geom_jitter())

# # par(mfrow = c(2,2))
# (numbat_cnv_analysis_plt <-
#     amp_plot + del_plot + neu_plot + loh_plot +
#     plot_annotation('Numbat CNV Calls Per Cell') +
#     plot_layout(nrow = 1,
#                 ncol = 4))

cnv_count_df <- as.data.frame(cnv_count_df)

df_melt <- melt(cnv_count_df)

# THIS
ggplot(df_melt, aes(variable[1],
                    value,
                    color = variable)) +
  geom_violin() +
  geom_point() +
  scale_colour_manual(values = c('red', 'blue', 'gray', 'green')) +
  ylab('Segment Counts') +
  xlab('') +
  geom_jitter() +
  facet_grid(~ variable)



