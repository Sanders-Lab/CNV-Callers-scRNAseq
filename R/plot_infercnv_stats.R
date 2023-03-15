library(tidyverse)
library(infercnv)
library(patchwork)
library(reshape2)


# after_run <- readRDS('../outputs/infercnv_pbmc3k_epithelium_tnbc1/run.final.infercnv_obj')


# step_22 <- readRDS('../outputs/infercnv_pbmc3k_epithelium_tnbc1/22_denoiseHMMi6.leiden.NF_NA.SD_1.5.NL_FALSE.infercnv_obj')

step_20 <- readRDS('../outputs/infercnv_pbmc3k_epithelium_tnbc1/20_HMM_pred.repr_intensitiesHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj')

expr_data_df <- as.data.frame(step_20@expr.data)
tnbc_cells <- expr_data_df %>% select(contains('tnbc1'))
tnbc_cells %>% ncol()


cnv_count_df <- data.frame()
for (i in 1:ncol(tnbc_cells)){
   curr_cell <- tnbc_cells[,i]
   amp_count <- 0
   del_count <- 0
   neu_count <- 0

    repeat_patterns <- curr_cell %>% rle()
    
          # store the current value of score
      pos_scores <- which(repeat_patterns$values > 1)
            neg_scores <- which(repeat_patterns$values < 1)
            neu_scores <- which(repeat_patterns$values == 1)
                  
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
      
   
     # saving a temp df with the same col format
     # will use this to rbind to the master df
     temp_df <- data.frame(amp_count,
                           del_count,
                           neu_count)
    cnv_count_df <- rbind(cnv_count_df, temp_df)
}

# checking if the min if 0
cnv_count_df$neu_count %>% min()
cnv_count_df <- cnv_count_df %>% apply(2, function(x) x+1)
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


pdf('../Infercnv_CNV_Calls_Per_Cell.pdf')
(infercnv_cnv_analysis_plt <- 
        amp_plot + del_plot + neu_plot +
            plot_annotation('Infercnv CNV Calls Per Cell') +
                plot_layout(nrow = 1,
                                            ncol = 3))
dev.off()


# THIS WORKS
df_melt <- melt(cnv_count_df)



pdf('../Infercnv_CNV_Calls_Per_Cell.pdf')
# TIHIS
(infercnv_cnv_analysis_plt <- ggplot(df_melt, aes(variable[1],
                                                                      value,
                                                                                          color = variable)) +
  geom_violin() +
    geom_point() + 
      scale_colour_manual(values = c('red', 'blue', 'gray')) +
        ylab('Segment Counts') +
          xlab('') +
            geom_jitter() +
              facet_grid(~ variable))

dev.off()


# ggsave(infercnv_cnv_analysis_plt,
#        '../Infercnv_CNV_Calls_Per_Cell.pdf',
#        device = 'pdf')
