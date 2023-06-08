library(tidyverse)



check <- read.table('/fast/groups/ag_sanders/work/projects/suharto/cnv_callers/proc/infercnv_atc2_ref_merged_counts.tsv')

check %>% head()

t=seq(0,10,0.1)
y=sin(t)
plot(t,y,type="l", xlab="time", ylab="Sine wave")
