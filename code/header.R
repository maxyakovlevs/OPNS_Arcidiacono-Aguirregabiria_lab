my.library <- 'C:/Users/yakov/Documents/R/win-library/3.6'
.libPaths(my.library)

library('tidyverse')
c('reshape2', 'stringr', 'magrittr', 'doParallel', 'chebpol', 'purrr', 'plyr', 'EnvStats', 'dplyr', 'purrrlyr', 'np') %>%
  walk(~library(., character.only=TRUE))

dir('modules') %>% 
  walk(~source(paste('./modules/', ., sep="")))

registerDoParallel(cores=28)
