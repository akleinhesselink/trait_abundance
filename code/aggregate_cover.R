rm(list = ls())

library(tidyverse)
library(sedgwickenv)
library(sedgwickspecies)
library(sedgwickcover)

species_cover <- 
  site_cover %>%
  left_join(sedgwick_plants, by = c('species' = 'calflora_binomial')) %>% 
  dplyr::select(USDA_symbol, site, cover) %>% 
  filter( !is.na(USDA_symbol)) %>% 
  distinct() %>%
  left_join(sedgwickenv, by = 'site') %>% 
  dplyr::select( site, USDA_symbol, cover, site_name, type, microsite) %>% 
  spread( USDA_symbol, cover, fill = 0 ) %>% 
  gather( USDA_symbol, cover, -c(site:microsite)) 

global_cover <- 
  species_cover %>% 
  group_by( USDA_symbol) %>% 
  summarise( abu = mean(cover))

type_cover <- 
  species_cover %>% 
  group_by( USDA_symbol, type) %>% 
  summarise( abu = mean(cover))
  
microsite_cover <- 
  species_cover %>% 
  group_by( USDA_symbol, type, microsite) %>% 
  summarise( abu = mean(cover))

save(global_cover, type_cover, microsite_cover, file = 'output/avg_cover.rda')

