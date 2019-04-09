rm(list = ls())

library(tidyverse)
source('code/helper_functions.R')

demo_pars <- read_csv('output/demo_pars.csv')

A <- as.matrix(demo_pars %>% select(demo_pars$prior_code)) 
A[is.na(A)] <- 0 
g <- as.numeric(demo_pars$g)
lambda <- as.numeric(demo_pars$lambda)
s <- as.numeric(demo_pars$s)
a_intra <- diag(A)

eta <- (g*lambda)/(1 - s*(1 -g))

mono_eq <- (eta - 1)/a_intra

out <- matrix(NA, length(lambda), length(lambda))

for ( i in 1:length(lambda)){ 
  for(j in 1:length(lambda)){ 
    
    out[i,j] <- (1-g[i])*s[i] + g[i]*lambda[i]/(1 + A[i, j]*mono_eq[j] )
    
  }
}

ndd_df <- 
  data.frame(out) %>% 
  mutate( focal_species = demo_pars$prior_code) %>% 
  gather( competitor, r_hat, starts_with('X')) %>% 
  mutate( competitor = factor(competitor, labels = demo_pars$prior_code)) %>% 
  filter( focal_species != competitor ) %>% 
  left_join(data.frame(prior_code = demo_pars$prior_code, mono_eq = mono_eq), by = c('focal_species' = 'prior_code'))

ndd_mean <- apply( out, 1, mean)/mono_eq 
demo_pars$ndd_mean <- ndd_mean
demo_pars <- demo_pars[order(demo_pars$ndd_mean), ]
demo_pars$prior_code <- factor(demo_pars$prior_code, levels = demo_pars$prior_code, ordered = T)
ndd <- demo_pars %>% select(prior_code, ndd_mean)

max(out)

ndd_labels <- 
  ndd %>% 
  mutate( x = 0.5*max(mono_eq)) %>% 
  mutate( y = 0.5*max(out)) %>% 
  mutate( label = round(ndd_mean, 4)) %>% 
  mutate( focal_species = prior_code)

ndd_labels

gg_ndd <- 
  ndd_df %>% 
  gather( type, value, r_hat:mono_eq) %>% 
  mutate( x = 1e-3, y = 1e-5) %>% 
  mutate( y = ifelse(type == 'r_hat', value, y)) %>% 
  mutate( x = ifelse(type == 'mono_eq', value, x )) %>% 
  ggplot( aes( x = x, y = y, group = interaction(focal_species, competitor))) + 
    geom_point() + 
    geom_line() + 
    stat_summary(aes( group = focal_species), fun.y = 'mean', geom = 'line', color = 'red') + 
    geom_text(data = ndd_labels, aes( x = x, y = y, label = label, group = prior_code), col = 'red') + 
    facet_wrap(~focal_species, 3, 6)  + 
    scale_x_continuous(breaks = c(0, 1000, 2000, 3000)) + 
    xlab( 'Abundance') + 
    ylab( 'Per Capita Invasion Growth Rate')

ggsave(gg_ndd, filename = 'figures/gg_ndd.png', width = 7, height = 5)


saveRDS(ndd, 'output/ndd.rds')
