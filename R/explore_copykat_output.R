library(tidyverse)
# library(infercnv)
library(reshape2)

copykat_preds <- read.table('../outputs/copykat_pbmc3k_epithelium_tnbc1/tnbc1_copykat_prediction.txt',
                            header = T)

aneu_cells <- copykat_preds[grep('tnbc1',copykat_preds$cell.names), ]

copykat_segs <- read.table('../outputs/copykat_pbmc3k_epithelium_tnbc1/tnbc1_copykat_CNA_results.txt',
                           header = T)
copykat_segs[1:5,1:5]
tnbc_cells <- copykat_segs %>% select(chrom,
                                      abspos,
                                      contains('tnbc1'))
tnbc_cells[1:5,1:5]
tnbc_cells %>% ncol()

curr_cell <- tnbc_cells[3]


# logic to get counts per cell

cnv_count_df <- data.frame()

# as cells start from 3rd col
for (i in 3:ncol(tnbc_cells)){

  amp_count <- 0
  del_count <- 0
  neu_count <- 0
  curr_cell <- tnbc_cells[,i]

  repeat_patterns <- curr_cell %>% rle()

      # store the current value of score
      pos_scores <- which(repeat_patterns$values > 0)
      neg_scores <- which(repeat_patterns$values < 0)
      neu_scores <- which(repeat_patterns$values == 0)

      # prev_score <- pos_scores[1] + 1
      for (k in 2:length(pos_scores)){
        # curr_score <- pos_scores[k]

        if(pos_scores[k] - 1 != pos_scores[k-1]){

          if(k == 2 | k == length(pos_scores)){
            amp_count <- amp_count + 1
          }
          amp_count <- amp_count + 1
        }
      }


      # prev_score <- neg_scores[1] + 1
      for (h in 2:length(neg_scores)){
        # curr_score <- neg_scores[k]

        if(neg_scores[h] - 1 != neg_scores[h-1]){
          if(h == 2 | h == length(neg_scores)){
            del_count <- del_count + 1
          }
          del_count <- del_count + 1
        }
      }



      if(length(neu_scores) > 0){
        neu_count <- neu_count + 1

        for (k in 2:length(neu_scores)){
          # curr_score <- neg_scores[k]

          if(neu_scores[k] - 1 != neu_scores[k-1]){
            if(k == 2 | k == length(neu_scores)){
              neu_count <- neu_count + 1
            }
            neu_count <- neu_count + 1
          }
        }
      }

    #   if (j == 1){
    #     prev_score <- cur
    #   }
    #
    # if (curr_score > 0){
    #   # checking for negative values: DEL
    #   # should give a positive num when mul with -1
    #   if(prev_score < 0){
    #     del_count <- del_count + 1
    #   }
    #   if(prev_score == 0)
    #   {
    #     neu_count <- neu_count + 1
    #   }
    # }
    #
    # if (curr_score < 0){
    #   # checking for positive values: AMP
    #   # should give a negative num when mul with -1
    #   if (prev_score > 0){
    #     amp_count <- amp_count + 1
    #   }
    #   if(prev_score == 0)
    #   {
    #     neu_count <- neu_count + 1
    #   }
    # }
    #
    #   if(curr_score == 0)
    #   {
    #     neu_count <- neu_count + 1
    #   }

      # updating the prev_score
      # prev_score <- curr_score


  temp_df <- data.frame(amp_count,
                        del_count,
                        neu_count)
  cnv_count_df <- rbind(cnv_count_df, temp_df)
}

# adding 1 for the lapse in logic
neu_count <- cnv_count_df[,3]
cnv_count_df <- cnv_count_df[,1:2] %>% apply(2, function(x) x+1)
cnv_count_df <- cbind(cnv_count_df, neu_count)
cnv_count_df <- as.data.frame(cnv_count_df)

# Plotting
(amp_plot <- ggplot(cnv_count_df, aes('',amp_count)) +
    geom_violin() +
    geom_point() +
    ylab('Count') +
    xlab('amp') +
    geom_jitter())


(del_plot <- ggplot(cnv_count_df, aes('',del_count)) +
    geom_violin() +
    geom_point() +
    ylab('') +
    xlab('del') +
    geom_jitter())


(neu_plot <- ggplot(cnv_count_df, aes('',neu_count)) +
    geom_violin() +
    geom_point() +
    ylab('') +
    xlab('neu') +
    geom_jitter())

# facets <- colnames(cnv_count_df)



# this!!!
df_melt <- melt(cnv_count_df)
# df_melt <- df_melt %>% rename(neu_count = neu_counts_df)

# (copykat_cnv_analysis_plt <-
#     amp_plot + del_plot + neu_plot +
#     plot_annotation('Copykat CNV Calls Per Cell') +
#     plot_layout(nrow = 1,
#                 ncol = 3) +
#     facet_grid(~ variable))



# TIHIS
ggplot(df_melt, aes(variable[1],
                    value,
                    color = variable)) +
  geom_violin() +
  geom_point() +
  scale_colour_manual(values = c('red', 'blue', 'gray')) +
  ylab('Segment Counts') +
  xlab('') +
  geom_jitter() +
  facet_grid(~ variable)

