rm(list = ls())

library(tidyverse)
library(sedgwickspecies)
library(sedgwicktraits)

fds <- read_csv('output/fitness_diffs.csv')
demo_pars <- read_csv('output/demo_pars.csv')
ndd <- readRDS('output/ndd.rds')
fa <- readRDS('output/species_feasible_abundances.rds')

load('output/avg_cover.rda')

predictor_traits <- 
  sedgwicktraits %>% 
  left_join(sedgwick_plants, by = c('species' = 'calflora_binomial')) %>% 
  group_by( prior_code) %>%
  summarise( max_height = mean(max_height), seed_mass = mean(seed_mass), CN_ratio = mean(CN_ratio), leaf_size = mean(leaf_size)) %>% 
  distinct() %>% 
  filter( !is.na(max_height), !is.na(seed_mass), !is.na(CN_ratio), !is.na(leaf_size))
  
alpha_ii <- diag(as.matrix( demo_pars[, demo_pars$prior_code]))
demo_pars$alpha_ii <- alpha_ii

predictors <- demo_pars %>% 
  select( -c(AGHE:SACA)) %>% 
  left_join(fds, by = 'prior_code') %>% 
  left_join(ndd, by = 'prior_code') %>% 
  left_join(fa, by = 'prior_code') %>% 
  left_join(predictor_traits, by = 'prior_code')

reserve_corr <- 
  global_cover %>% 
  left_join(sedgwick_plants) %>% 
  distinct(USDA_symbol, prior_code, global_cover) %>% 
  filter( ! global_cover == 0 ) %>% 
  left_join(predictors) %>% 
  filter( !is.na(prior_code), !is.na(lambda)) %>% 
  gather( stat, value, lambda:leaf_size, -n_comm) %>% 
  mutate( stat = factor(stat, levels = c('lambda', 'g', 's', 'alpha_ii', 'pred_abu', 'pred_abu2', 'ndd_mean', 'pred_rank_fd', 'mean_fd', 'max_height', 'seed_mass', 'CN_ratio', 'leaf_size'))) %>% 
  group_by( stat)  


type_corr <- 
  type_cover %>% 
  left_join(sedgwick_plants) %>% 
  distinct(USDA_symbol, prior_code, type, type_cover) %>% 
  filter( ! type_cover == 0 ) %>% 
  left_join(predictors) %>% 
  filter( !is.na(prior_code), !is.na(lambda)) %>% 
  group_by( type ) %>% 
  gather( stat, value, lambda:leaf_size, -n_comm) %>% 
  mutate( stat = factor(stat, levels = c('lambda', 'g', 's', 'alpha_ii', 'pred_abu', 'pred_abu2', 'ndd_mean', 'pred_rank_fd', 'mean_fd', 'max_height', 'seed_mass', 'CN_ratio', 'leaf_size'))) %>% 
  group_by( stat, type )


microsite_corr <- 
  microsite_cover %>% 
  left_join(sedgwick_plants) %>% 
  distinct(USDA_symbol, prior_code, microsite, microsite_cover) %>% 
  filter( ! microsite_cover == 0 ) %>% 
  left_join(predictors) %>% 
  filter( !is.na(prior_code), !is.na(lambda)) %>% 
  group_by( microsite) %>% 
  gather( stat, value, lambda:leaf_size, -n_comm) %>% 
  mutate( stat = factor(stat, levels = c('lambda', 'g', 's', 'alpha_ii', 'pred_abu', 'pred_abu2', 'ndd_mean', 'pred_rank_fd', 'mean_fd', 'max_height', 'seed_mass', 'CN_ratio', 'leaf_size'))) %>% 
  group_by(stat, type, microsite ) 


reserve_cortest <- 
  reserve_corr %>% 
  summarise(r = cor.test(global_cover, value)$estimate, 
            p = cor.test(global_cover, value)$p.value, 
            label_x = 0.3*(max(value)), 
            label_y = 0.3*(max(global_cover))) %>% 
  mutate( r_label = paste0('r=', round(r, 2)), p.val_label = paste0('p=',round(p, 3))) %>% 
  mutate( p.val_label = ifelse(p > 0.05, 'n.s.', p.val_label)) 

type_cortest <- 
  type_corr %>% 
  summarise(r = cor.test(type_cover, value)$estimate, 
            p = cor.test(type_cover, value)$p.value, 
            label_x = 0.3*(max(value)), 
            label_y = 0.3*(max(type_cover))) %>% 
  mutate( r_label = paste0('r=', round(r, 2)), p.val_label = paste0('p=',round(p, 3))) %>% 
  mutate( p.val_label = ifelse(p > 0.05, 'n.s.', p.val_label)) 

microsite_cortest <- 
  microsite_corr %>% 
  summarise(r = cor.test(microsite_cover, value)$estimate, 
            p = cor.test(microsite_cover, value)$p.value, 
            label_x = 0.3*(max(value)), 
            label_y = 0.3*(max(microsite_cover))) %>% 
  mutate( r_label = paste0('r=', round(r, 2)), p.val_label = paste0('p=',round(p, 3))) %>% 
  mutate( p.val_label = ifelse(p > 0.05, 'n.s.', p.val_label)) 


gg_reserve <- 
  reserve_corr %>% 
  ungroup() %>% 
  filter( stat %in% c('lambda', 'alpha_ii', 'pred_abu', 'pred_rank_fd', 'mean_fd')) %>% 
  ggplot( aes( x = value, y = global_cover)) + 
  facet_wrap(~stat, scales = 'free') + 
  geom_point() +
  geom_smooth(method = 'lm', se = F)  + 
  geom_text( data = reserve_cortest %>% filter(stat %in% c('lambda', 'alpha_ii', 'pred_abu', 'pred_rank_fd', 'mean_fd') ), 
             aes( x = label_x, y = label_y, label = paste(r_label, p.val_label)), hjust=0, color = 'darkblue') + 
  ylab('Average Cover %') + 
  ggtitle('All Sites') 



gg_upper <- 
  type_corr %>% 
  filter( type == 'upper') %>%
  ungroup() %>% 
  filter( stat %in% c('lambda', 'alpha_ii', 'pred_abu', 'pred_rank_fd', 'mean_fd')) %>% 
  ggplot( aes( x = value, y = type_cover)) + 
  facet_wrap(~stat, scales = 'free') + 
  geom_point() +
  geom_smooth(method = 'lm', se = F)  + 
  geom_text( data = type_cortest %>% filter(type == 'upper', stat %in% c('lambda', 'alpha_ii', 'pred_abu', 'pred_rank_fd', 'mean_fd')), 
             aes( x = label_x, y = label_y, label = paste(r_label, p.val_label)), hjust=0, color = 'darkblue') +
  ylab('Average Cover %') + 
  ggtitle('Upper Sites Only') 


gg_lower <- 
  type_corr %>% 
  filter( type == 'lower') %>%
  ungroup() %>% 
  filter( stat %in% c('lambda', 'alpha_ii', 'pred_abu', 'pred_rank_fd', 'mean_fd')) %>% 
  ggplot( aes( x = value, y = type_cover)) + 
  facet_wrap(~stat, scales = 'free') + 
  geom_point() +
  geom_smooth(method = 'lm', se = F)  + 
  geom_text( data = type_cortest %>% filter(type == 'lower', stat %in% c('lambda', 'alpha_ii', 'pred_abu', 'pred_rank_fd', 'mean_fd')), 
             aes( x = label_x, y = label_y, label = paste(r_label, p.val_label)), hjust=0, color = 'darkblue')  + 
  ylab('Average Cover %') + 
  ggtitle('Lower Sites Only') 


gg_hummock <- 
  microsite_corr %>% 
  filter( type == 'upper', microsite == 'hummock') %>%
  ungroup() %>% 
  filter( stat %in% c('lambda', 'alpha_ii', 'pred_abu', 'pred_rank_fd', 'mean_fd')) %>% 
  ggplot( aes( x = value, y = microsite_cover)) + 
  facet_wrap(~stat, scales = 'free') + 
  geom_point() +
  geom_smooth(method = 'lm', se = F)  + 
  geom_text( data = microsite_cortest %>% filter( type == 'upper', microsite == 'hummock', stat %in% c('lambda', 'alpha_ii', 'pred_abu', 'pred_rank_fd', 'mean_fd')), 
             aes( x = label_x, y = label_y, label = paste(r_label, p.val_label)), hjust=0, color = 'darkblue') + 
  ylab('Average Cover %') + 
  ggtitle('Upper Hummock Sites Only') 


gg_upper_nonhummock <- 
  microsite_corr %>% 
  filter( type == 'upper', microsite == 'not_hummock') %>%
  ungroup() %>% 
  filter( stat %in% c('lambda', 'alpha_ii', 'pred_abu', 'pred_rank_fd', 'mean_fd')) %>% 
  ggplot( aes( x = value, y = microsite_cover)) + 
  facet_wrap(~stat, scales = 'free') + 
  geom_point() +
  geom_smooth(method = 'lm', se = F)  + 
  geom_text( data = microsite_cortest %>% filter( type == 'upper', microsite == 'not_hummock', stat %in% c('lambda', 'alpha_ii', 'pred_abu', 'pred_rank_fd', 'mean_fd')), 
             aes( x = label_x, y = label_y, label = paste(r_label, p.val_label)), hjust=0, color = 'darkblue') + 
  ylab('Average Cover %') + 
  ggtitle('Upper Non-Hummock Sites Only') 

ggsave(gg_reserve, filename = 'figures/cor_raw_vals_reserve.png', width = 8, height = 5, units = 'in', dpi = 'print')

library(gridExtra)
ggsave( grid.arrange(gg_upper, gg_lower, nrow = 2), filename = 'figures/cor_raw_vals_upper_lower.png', width = 8, height = 5, units = 'in', dpi = 'print')

ggsave( grid.arrange(gg_hummock, gg_upper_nonhummock, nrow = 2), filename = 'figures/cor_raw_vals_hummock.png', width = 8, height = 5, units = 'in', dpi = 'print')

save(gg_reserve, gg_lower, gg_upper, gg_hummock, gg_upper_nonhummock, file = 'output/correlation_rawvals_plots.rda')
