rm(list = ls())
library(stringr)
library(sedgwickspecies)
library(tidyverse)

source('code/helper_functions.R')

alphas <- read_csv('output/alpha_estimates_row_is_target.csv')
demo_pars <- read_csv('data/old_data/species_rates.csv')

### Fix mispelled species 
demo_pars <- 
  demo_pars %>% 
  mutate( species = str_replace(species, 'heretophylla', 'heterophylla'))

demo_pars <- 
  demo_pars %>% 
  separate( lambda, c('lambda', 'lambda_sd'), sep = '\xb1') %>% 
  separate( g, c('g', 'g_sd'), sep = '\xb1') %>% 
  separate( s, c('s', 's_sd'), sep = '\xb1') 

demo_pars <- 
  demo_pars %>% 
  gather(stat, value, lambda:s_sd ) %>% 
  mutate(value = as.numeric(value)) %>% 
  spread(stat, value)

demo_pars <- 
  demo_pars %>% 
  left_join(sedgwick_plants, by = c('species' = 'prior_name')) %>% 
  dplyr::select(prior_code, lambda, g, s) %>% 
  left_join( ( alphas %>% 
                 rename('prior_code' = X1) ) )

demo_pars <- 
  demo_pars %>% 
  arrange(prior_code) 

demo_pars <- 
  demo_pars %>% 
  dplyr::select(prior_code, lambda, g, s, demo_pars$prior_code) 

demo_pars %>%
  gather( competitor, alpha, AGHE:SACA) %>% 
  left_join(sedgwick_plants, by = 'prior_code') %>% 
  select( USDA_symbol, lambda, g, s, competitor, alpha) %>% 
  left_join(sedgwick_plants, by = c('competitor' = 'prior_code')) %>% 
  select( USDA_symbol.x, lambda, g, s, USDA_symbol.y, alpha) %>% 
  rename( 'focal' = USDA_symbol.x, competitor = 'USDA_symbol.y') %>% 
  arrange( focal, competitor) %>% 
  distinct() %>%
  spread( competitor, alpha) %>% 
  write_csv('output/demo_pars.csv')

