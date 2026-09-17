library(tidyverse)
source('utilities.R')
HDs  <-  readRDS ('data/sensitivity.outcomes2.RDS') %>% data.frame
class(HDs)
HDs$y_axis  <-  1:nrow(HDs)

colnames(HDs)  <-  c('outcome', 'estimate', 'low_ci', 'high_ci', 'y_axis')
rownames(HDs)
#add suffix _bool to rownames of HDs
rownames(HDs)  <-  paste0(rownames(HDs), '_bool')
g  <-  make.HD.plot(HDs, label_list2)
ggsave(g, width=8, height=5, filename = sprintf('figs/sensitivity.analysis2.pdf'))

