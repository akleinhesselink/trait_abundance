rm(list = ls())
library(tidyverse)
library(sedgwicktraits)
library(sedgwickspecies)
library(glmnet)

fa <- readRDS('output/species_feasible_abundances.rds')

fa <- 
  fa %>% 
  left_join(sedgwick_plants, by = 'prior_code') %>% 
  left_join(sedgwicktraits, by = c('calflora_binomial' = 'species')) %>% 
  filter( dataset == 'TAPIOCA') %>%
  select( -sp  )

fa %>% 
  select( -c(calflora_binomial:USDA_symbol), -c(notes:dataset) ) %>% 
  gather( trait, trait_val, leaf_size:leaf_pH) %>% 
  ggplot( aes( x = trait_val, y = value) ) + 
  geom_point() + 
  geom_smooth(se = F, method = 'lm') + 
  facet_grid(stat ~ trait, scales = 'free')


fa <- 
  fa %>% 
  spread(stat, value) %>% 
  select( USDA_symbol, aboveground_abundance, avg_abundance, lambda, leaf_size:leaf_pH)

m1 <- lm(data = fa, aboveground_abundance ~ leaf_size + SLA + foliar_N + CN_ratio + LDMC + phenology + seed_mass + max_height)
summary(m1)

fit1 <- glmnet(x = scale(m1$model)[,-1], y = scale(m1$model[,1]))
print(fit1)

coef(fit1,s=0.01) # extract coefficients at a single value of lambda

saveRDS(fit1, file = 'output/lasso_fit.rds')
