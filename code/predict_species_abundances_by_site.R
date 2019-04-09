rm(list = ls())

library(tidyverse)
library(sedgwickenv)
library(sedgwickspecies)
library(sedgwickcover)
library(sedgwicktraits)

seed_mass <- 
  sedgwicktraits %>% 
  left_join(sedgwick_plants, by = c('species' = 'calflora_binomial')) %>% 
  ungroup() %>% 
  dplyr::select(USDA_symbol, seed_mass) %>% 
  group_by( USDA_symbol ) %>% 
  summarise(seed_mass = mean(seed_mass)) %>% 
  distinct()

fa <- readRDS('output/site_feasible_abundances.rds')

fa <- 
  fa %>% 
  left_join(sedgwick_plants, by = 'prior_code') %>% 
  dplyr::select( site, USDA_symbol, stat, value) %>% 
  left_join(seed_mass, by = 'USDA_symbol')

fa <- 
  fa %>% 
  spread( stat, value) %>% 
  mutate( pred_abu2 = pred_abu*seed_mass, pred_abu3 = (pred_abu/g), pred_abu4 = pred_abu3*seed_mass) %>% 
  gather( stat, value, g:pred_abu4) %>% 
  dplyr::select(-seed_mass)

species_cover <- 
  site_cover %>%
  left_join(sedgwick_plants, by = c('species' = 'calflora_binomial')) %>% 
  dplyr::select(USDA_symbol, site, cover) %>% 
  filter( !is.na(USDA_symbol)) %>% 
  distinct() %>%
  left_join(sedgwickenv, by = 'site') %>% 
  dplyr::select( site, USDA_symbol, cover, site_name, type, microsite) %>% 
  spread( USDA_symbol, cover, fill = 0 ) %>% 
  gather( USDA_symbol, cover, -c(site:microsite)) %>% 
  group_by( site ) %>% 
  mutate( pcov = cover/sum(cover), 
          rank = rank(cover)) 

global <- 
  species_cover %>% 
  group_by( USDA_symbol ) %>% 
  summarise( obs_cover = mean(cover), obs_pcov = mean(pcov))  

type <- 
  species_cover %>% 
  group_by(USDA_symbol, type) %>% 
  summarise( obs_cover = mean(cover), obs_pcov = mean(pcov)) %>% 
  arrange( USDA_symbol, type)


fa_wide <- 
  fa %>% 
  left_join(sedgwickenv %>% dplyr::select(site, type, microsite), by = 'site') %>% 
  unite( microsite, c(type, microsite), sep = '_')  %>% 
  group_by(microsite,stat, USDA_symbol) %>% 
  summarise( value = mean(value))  %>% 
  spread(stat, value )

microsite <- 
  species_cover %>% 
  group_by( USDA_symbol, type, microsite) %>% 
  summarise( obs_cover = mean(cover), obs_pcov = mean(pcov), obs_rank = mean(rank))%>% 
  arrange( USDA_symbol, type, microsite) %>% 
  unite( microsite, c(type, microsite), sep = '_') %>% 
  left_join(fa_wide, by = c('microsite', 'USDA_symbol')) %>% 
  gather(type_y, y, obs_cover:obs_rank) %>% 
  gather(type_x, x, c(pred_abu, pred_abu2, pred_abu3, pred_abu4, lambda, g, s)) %>% 
  dplyr::select(USDA_symbol, microsite, type_x, x, type_y, y) %>% 
  filter( !is.na(y), !is.na(x)) %>% 
  group_by(microsite, type_x, type_y) %>% 
  mutate( plot_label = ifelse(y == max(y), USDA_symbol, '')) 

cor_df <- 
  microsite %>% 
  #filter( USDA_symbol != 'HECO7') %>% 
  group_by( microsite, type_x, type_y ) %>% 
  summarise(r = cor.test(x = x, y = y)$estimate, 
         p.val = cor.test(x = x, y = y)$p.value) %>%   
  mutate( r_label = paste0('r=', round(r, 2)), p.val_label = paste0('p=',round(p.val, 3))) %>% 
  mutate( p.val_label = ifelse(p.val > 0.05, 'n.s.', p.val_label)) %>% 
  left_join(microsite) %>% 
  filter( !is.na(y), !is.na(x)) %>% 
  group_by( microsite, type_x, type_y, r, p.val, r_label, p.val_label) %>% 
  summarise(x = max(x, na.rm = T), y = max(y, na.rm = T))

microsite %>% 
  #filter( USDA_symbol != 'HECO7') %>% 
  filter(type_y == 'obs_cover') %>% 
  ggplot(aes( x = x, y = y, color = microsite)) +
  geom_point() +
  geom_smooth(se = F, method = 'lm') + 
  geom_text( data = cor_df %>% 
               filter(type_y == 'obs_cover'), 
             aes( x = x*0.1, y = y*0.9, label = r_label), hjust=0)  + 
  geom_text( data = cor_df %>% 
               filter(type_y == 'obs_cover'), 
             aes( x = x*0.1, y = y*0.7, label = p.val_label), hjust=0) + 
  geom_text(aes(label = plot_label)) + 
  facet_grid(microsite ~ type_x, scales = 'free') + 
  guides(color = F) + 
  ylab('observed cover (%)')

microsite %>% 
  #filter( USDA_symbol != 'HECO7') %>% 
  filter(type_y == 'obs_pcov') %>% 
  ggplot(aes( x = x, y = y, color = microsite)) +
  geom_point() +
  geom_smooth(se = F, method = 'lm') + 
  geom_text( data = cor_df %>% 
               filter(type_y == 'obs_pcov'), 
             aes( x = x*0.1, y = y*0.9, label = r_label), hjust=0)  + 
  geom_text( data = cor_df %>% 
               filter(type_y == 'obs_pcov'), 
             aes( x = x*0.1, y = y*0.7, label = p.val_label), hjust=0) + 
  geom_text(aes(label = plot_label)) + 
  facet_grid(microsite ~ type_x, scales = 'free') + 
  guides(color = F) + 
  ylab('observed proportion cover (%)')

saveRDS(microsite, file = 'output/avg_microsite_cover.rds')


