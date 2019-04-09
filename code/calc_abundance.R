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
eta <- (g*lambda)/(1 - s*(1 -g))

mysets <- list(NA)

for( i in 1:18){ 
  mysets[[i]] <- combn(1:18, i)
}

solution <- lapply( mysets, function(x){ x[] <- NA; x} )
#Jac <- lapply( mysets, function(x){ y <- rep(list(NA), ncol(x)); y} )

for( i in 1:length(mysets)){ 
  
  temp_set <- mysets[[i]]
  
  for( j in 1:ncol(temp_set)){
    
    my_spp <- temp_set[,j]
    A_temp <- A[my_spp, my_spp]
    eta_temp <- eta[my_spp]
    
    solution[[i]][,j] <- solve(A_temp) %*% ( eta_temp - 1)
    
    #x <- solution[[i]][,j]
    #Jac[[i]][[j]] <- get_Jacobian(x, A_temp, lambda_temp, g_temp, s_temp)
  }

}

saveRDS(mysets, file =  'output/species_combos.rda')
saveRDS(solution, file = 'output/solutions.rda')

quick_gather <- function(x, set = 'set', sp = 'sp') { data.frame(x) %>% gather( !!set, !!sp ) }

all_sets <- do.call(bind_rows, 
        c(lapply( mysets, 
                  FUN = quick_gather, 
                  set = 'set', 
                  sp = 'sp'), 
          .id  = "nsp" ))

all_eq <- do.call(bind_rows, 
                  c(lapply(solution, 
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

two_sp_comm <- all_sets %>% filter(nsp == 2)

two_sp_comm %>% 
  left_join(name_df, by = 'sp') %>% 
  mutate( abundance = ifelse(abundance < 1, 0, abundance)) %>% 
  group_by(set) %>% 
  mutate( rel_abu = abundance/sum(abundance))  %>% 
  group_by( prior_code ) %>% 
  summarise( rel_abu = mean(rel_abu))  %>% 
  write_csv('output/pred_rel_abu.csv')


sp_summary <- 
  all_sets %>% 
  dplyr::select( nsp, set, sp, abundance, feasible) %>% 
  left_join(name_df, by = 'sp') %>% 
  ungroup() %>% 
  unite(comm, c(nsp, set), sep = '-', remove = F) %>%   
  filter( feasible) %>% 
  dplyr::select(nsp, comm, prior_code, abundance) %>% 
  group_by( comm ) %>% 
  mutate( rank = rank(abundance)) %>% 
  group_by( prior_code ) %>% 
  summarise( pred_abu = mean(abundance), 
             pred_abu2= abundance[nsp == 1],
             n_comm = n_distinct(comm[abundance > 0]))



saveRDS(sp_summary, 'output/species_feasible_abundances.rds')
