library(tidyverse)
source('code/helper_functions.R')

demo <- read_csv('output/demo_pars.csv')

load('output/avg_cover.rda')

sp_pool <- global_cover$USDA_symbol

demo <- 
  demo %>% 
  filter( focal %in% sp_pool) %>% 
  select( focal:s, .$focal)

all_pairs <- 
  expand.grid( focal = demo$focal, comp = demo$focal) %>% 
  filter( focal != comp ) %>% 
  arrange( focal, comp) %>% 
  mutate( focal = as.character(focal), comp = as.character(comp))


all_pairs$fd <- NA
all_pairs$rho <- NA


for( i in 1:nrow(all_pairs)) { 
  temp_pair <- as.character(all_pairs[i, c('focal', 'comp') ])
  
  demo_temp <- 
    demo %>% 
    select( focal:s , temp_pair) %>% 
    filter( focal %in% temp_pair) 

  A <- as.matrix(demo_temp[, temp_pair])
  
  g <- demo_temp$g
  s <- demo_temp$s
  lambda <- demo_temp$lambda
  
  A[ is.na(A) ] <- 0 # set NA to zero 
  
  all_pairs$fd[i] <- fitness_difference(A, lambda, g, s )
  all_pairs$rho[i] <- rho(A) 
}

all_pairs$rho


all_pairs %>% 
  filter( fd > 1 ) %>% 
  mutate( rho = ifelse( rho < 0, 0, rho)) %>% 
  mutate( rho = ifelse( rho > 1, 1, rho)) %>% 
  mutate( coexist = ifelse(rho < 1/fd, 'coexist', 'exclude' )) %>% 
  ggplot( aes(x = 1 - rho, y = fd, color = coexist)) + 
  scale_y_log10() + 
  geom_point()  + 
  geom_line(data = test, aes( x = x, y= y, color = NA))

test <- data.frame( x = seq(0, 1, 0.001)) %>% 
  mutate( y = 1/(1-x))





