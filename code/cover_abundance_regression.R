rm(list = ls())

library(tidyverse)

abu_list <- read_lines('data/cover_abundance.csv')

groups <- seq_along(abu_list) %% 3 == 1
split_at <- seq_along(abu_list)[groups] - 1

out <- list()
for ( i in 1:length(split_at)){ 
  out[[i]] <- read_csv( 'data/cover_abundance.csv', skip = split_at[i], n_max = 2)
  out[[i]]$type <- c('cover', 'abu')
}

out <- do.call( bind_rows, lapply(out, function(x) {names(x)[1:2] <- c('site', 'plot'); x} ) )

out <- 
  out %>% 
  data.frame() %>% 
  select(-starts_with('X')) %>% 
  gather( species, value, -c(site, plot, type)) 

out <- 
  out %>% 
  spread(type, value) 

out <- 
  out %>% 
  mutate(species = ifelse( species == 'trifolium', 'TRIF', species)) %>% 
  group_by( site, plot, species) %>% 
  summarise( abu = sum(abu, na.rm = T), cover = sum(cover, na.rm = T)) %>% 
  ungroup() %>% 
  mutate( abu = ifelse( cover > 0 & abu == 0 , NA, abu))

grass_cover <- 
  out %>% 
  filter( species == 'grass' ) %>% 
  select( - abu )

out <- 
  out %>% 
  filter( species != 'grass')

temp_list <- 
  out %>% 
  split(. , .$species)

library(sedgwickspecies)

View( sedgwick_plants %>% 
  distinct(USDA_symbol, prior_code) )


names( temp_list ) <- c('LOWR2', 'EUSP', 'HECO7', 'LACA7', 'MICA', 'MIDO', 'NAAT', 'PLER3', 'TRAL5')

cover_regressions <- 
  lapply( temp_list, function(x) lm(abu ~ -1 + cover, data = x))

out %>% 
  filter( species != 'grass') %>%
  ggplot( aes( x = cover, y = abu)) + 
  geom_point() + 
  geom_smooth(se = F, method = 'lm', formula = 'y ~ -1 + x') + 
  facet_wrap(~species, scales = 'free')


plants_per_cover <- data.frame( lapply(cover_regressions, coef) )

plants_per_cover <- 
  plants_per_cover %>% 
  gather( spp, n_per_cover)

save(plants_per_cover, cover_regressions, file = 'output/plants_per_cover.rda')







