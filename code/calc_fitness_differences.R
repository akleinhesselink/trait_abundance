rm(list = ls())
library(stringr)
library(tidyverse)

source('code/helper_functions.R')

demo_pars <- read_csv('output/demo_pars.csv')

A <- as.matrix(demo_pars %>% select(demo_pars$prior_code)) 
A[is.na(A)] <- 0 

g <- as.numeric(demo_pars$g)

lambda <- as.numeric(demo_pars$lambda)
s <- as.numeric(demo_pars$s)

all_pairs <- as.matrix(expand.grid( X1= 1:length(lambda), X2 = 1:length(lambda)))
all_pairs <- all_pairs[ all_pairs[,1] != all_pairs[,2],  ]

fdiff <- rep(NA, nrow(all_pairs))

for( i in 1:nrow(all_pairs)){ 
  myset <- all_pairs[i,]
  fdiff[i] <- fitness_difference(A[myset, myset], lambda[myset], g[myset], s[myset])
}


fit_diffs <- 
  data.frame( all_pairs, fdiff ) %>% 
  mutate( X1 = factor(X1, labels = colnames(A))) %>% 
  mutate( X2 = factor(X2, labels = colnames(A))) %>% 
  mutate( X2_wins = fdiff > 1) 


X1_diffs <- 
  fit_diffs %>% 
  arrange( X1, X2) %>% 
  filter( is.finite(fdiff) ) %>% 
  group_by( X1 ) %>% 
  summarise( mean_fd = mean(fdiff, na.rm = TRUE) )

fit_diffs <- 
  fit_diffs %>% 
  group_by(X1) %>% 
  summarise( losses = sum(X2_wins, na.rm = T), n_nas = sum(is.na(X2_wins))) %>% 
  arrange(losses)

fit_diffs %>% 
  left_join(X1_diffs, by = c('X1')) %>% 
  rename('prior_code' = X1) %>% 
  dplyr::select(prior_code, losses, mean_fd) %>% 
  rename('pred_rank_fd' = losses)  %>% 
  write_csv('output/fitness_diffs.csv')



