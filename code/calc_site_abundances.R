rm(list = ls())
library(stringr)
library(sedgwickspecies)
library(tidyverse)

source('code/helper_functions.R')

alphas <- read_csv('output/alpha_estimates_row_is_target.csv')
lambdas <- read_csv('output/lambda_estimates.csv')
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

A <- as.matrix(demo_pars %>% select(AGHE:SACA)) 
A[is.na(A)] <- 0 
g <- as.numeric(demo_pars$g)
lambda <- as.numeric(demo_pars$lambda)
s <- as.numeric(demo_pars$s)

site_pools <- 
  sedgwickcover::site_cover %>% 
  left_join(sedgwick_plants, by = c('species' = 'calflora_binomial')) %>% 
  filter( prior_code %in% demo_pars$prior_code) %>% 
  distinct(site, species, prior_code, cover) %>% 
  group_by( site ) %>% 
  summarise(my_list = list( unique(prior_code[!is.na(prior_code)])))


solution <- as.list( rep(NA, length(site_pools$my_list)) )

for(p in 1:length(site_pools$my_list)){
  
  temp_site <- site_pools$site[p]
  temp_pool <- site_pools$my_list[[p]]
  
  temp_pool <- as.numeric(factor(demo_pars$prior_code)[demo_pars$prior_code %in% temp_pool])
  
  mysets <- list(NA)
  for( i in 1:length(temp_pool)){ 
    mysets[[i]] <- combn(temp_pool, i)
  }

  site_solution <- lapply( mysets, function(x){ x[] <- NA; x} )

  for( i in 1:length(mysets)){ 
    
    temp_set <- mysets[[i]]
    
    for( j in 1:ncol(temp_set)){
      
      my_spp <- temp_set[,j]
      g_temp <- g[my_spp]
      lambda_temp <- lambda[my_spp]
      s_temp <- s[my_spp]
      A_temp <- A[my_spp, my_spp]
      
      eta <- (g_temp*lambda_temp)/(1 - s_temp*(1 -g_temp))
      site_solution[[i]][,j] <- solve(A_temp) %*% ( eta - 1)
      
      #x <- solution[[i]][,j]
      #Jac[[i]][[j]] <- get_Jacobian(x, A_temp, lambda_temp, g_temp, s_temp)
    }
    
  }
  
  quick_gather <- function(x, set = 'set', sp = 'sp') { data.frame(x) %>% gather( !!set, !!sp ) }

  all_sets <- do.call(bind_rows, 
                      c(lapply( mysets, 
                                FUN = quick_gather, 
                                set = 'set', 
                                sp = 'sp'), 
                        .id  = "nsp" ))

  all_eq <- do.call(bind_rows, 
                    c(lapply(site_solution, 
                             FUN = quick_gather, 
                             set = 'set', 
                             sp = 'abundance'), 
                      .id = 'nsp'))

  all_sets <- 
    all_sets %>% 
    bind_cols(all_eq) %>% 
    group_by( nsp, set ) %>% 
    mutate( feasible = all(abundance > 0))

  name_df <- data.frame( sp = 1:length(demo_pars$prior_code), prior_code = demo_pars$prior_code)

  sp_summary <- 
    all_sets %>% 
    dplyr::select( nsp, set, sp, abundance, feasible) %>% 
    left_join(name_df, by = 'sp') %>% 
    ungroup() %>% 
    unite(comm, c(nsp, set), sep = '-') %>%   
    filter( feasible) %>% 
    dplyr::select( comm, prior_code, abundance) %>% 
    group_by( comm ) %>% 
    group_by( prior_code ) %>% 
    summarise( pred_abu = mean(abundance))

  sp_summary <- 
    sp_summary %>% 
    left_join((demo_pars %>% select(prior_code:s)))

  sp_summary %>% 
    select( -prior_code ) %>% 
    pairs()

  sp_summary <- 
    sp_summary %>% 
    gather( stat, value, c(pred_abu, lambda, g, s)) %>% 
    mutate( site = !!temp_site)
  
  solution[[p]] <- sp_summary
}

site_fa <- do.call( rbind, solution ) 

saveRDS(site_fa, 'output/site_feasible_abundances.rds')
