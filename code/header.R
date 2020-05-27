my.library <- 'C:/Users/yakov/Documents/R/win-library/3.6'
.libPaths(my.library)

library('tidyverse')
c('reshape2', 'stringr', 'magrittr', 'doParallel', 'fExtremes','dplyr') %>%
  walk(~library(., character.only=TRUE))

dir('modules') %>% 
  walk(~source(paste('./modules/', ., sep="")))

var_save <- '../variables/'

registerDoParallel(cores=28)

