rm(list = ls())

library(tidyverse)
library(sedgwickspecies)

load('output/avg_cover.rda')

pairwise_nd <- read_csv('data/old_data/for_will.csv')

pairwise_nd <- 
  pairwise_nd %>% 
  separate(species_pair, c('focal', 'competitor')) %>% 
  rename( 'nd' = `niche difference`, 'fd' = `fitness difference`) 

sedgwick_codes <- 
  sedgwick_plants %>% 
  ungroup() %>% 
  select(prior_code, USDA_symbol) %>% 
  filter( !is.na(prior_code)) %>% 
  distinct()

fd <- 
  pairwise_nd %>% 
  left_join(sedgwick_codes, by = c('focal' = 'prior_code')) %>% 
  mutate( focal = USDA_symbol) %>% 
  select( - USDA_symbol ) %>% 
  left_join(sedgwick_codes, by = c('competitor' = 'prior_code')) %>% 
  mutate( competitor = USDA_symbol) %>% 
  select( - USDA_symbol) %>% 
  unite( 'pair', c(focal, competitor)) %>% 
  select( pair, fd )

comp_grid <- 
  expand.grid(focal = sedgwick_codes$USDA_symbol, competitor = sedgwick_codes$USDA_symbol) %>% 
  filter( competitor != focal)


join_cover_fd <- function( abu, pair_grid, fd_pairs)  { 

  abu_ratio <- 
    pair_grid %>% 
    left_join(abu, by = c('focal' = 'USDA_symbol')) %>% 
    rename( 'focal_abu' = abu) %>% 
    left_join(abu, by = c('competitor' = 'USDA_symbol')) %>% 
    rename( 'competitor_abu' = abu) %>% 
    mutate( abu_ratio = competitor_abu/focal_abu  ) %>% 
    unite( 'pair', c(focal, competitor)) %>% 
    select(pair, abu_ratio)

  all_pairs <- 
    fd_pairs %>% 
    left_join(abu_ratio)  %>% 
    filter( ! is.na(abu_ratio))
  
  return( all_pairs) 
} 

x <- rnorm(20)
y <- rnorm(20)
m <- lm(x ~ y)
sm <- summary(m)

sm$coefficients[1, 4]

plot_correlation <- function( x  ){ 
  
  x_pos <- 0.5*max( log(x$abu_ratio))
  y_pos <- 0.9*max( log(x$fd))
  x_pos2 <- 1.1*max(log(x$abu_ratio))
  y_pos2 <- min(log(x$abu_ratio)) + 0.2*(max(log(x$abu_ratio)))
  
  x_cor <- cor.test( log( x$fd), log(x$abu_ratio))
  mod <- summary( lm(log(x$abu_ratio) ~ log(x$fd)))$coefficients 
  
  x %>% 
    ggplot(aes( x = log(fd), y = log(abu_ratio)) ) + 
    geom_point() + 
    geom_smooth(se = F, method = 'lm') + 
    annotate( x = x_pos, y = y_pos, geom = 'text', 
              label = paste('r =', round(x_cor$estimate, 2), 
                            '; p =', round(x_cor$p.value, 2)), color = 'darkblue') + 
    annotate( x = x_pos2, y = y_pos2, geom = 'text', 
              label = paste( 'Intcpt =', round(mod[1, 1], 2), 
                             '; p-val:', round( mod[1, 4], 3), '; \n', 
                             'slope = ', round(mod[2, 2],2), 
                             '; p-val:', round(mod[2, 4], 3)), hjust = 1)
}



global_pairs <- join_cover_fd(global_cover, comp_grid, fd)

upper_pairs <- join_cover_fd(type_cover %>% filter( type == 'upper', abu > 0), comp_grid, fd)

lower_pairs <- join_cover_fd(type_cover %>% filter( type == 'lower', abu > 0), comp_grid, fd)

hummock_pairs <- join_cover_fd(microsite_cover %>% filter( type == 'upper', microsite == 'hummock', abu > 0), comp_grid, fd)

non_hummock_pairs <- join_cover_fd(microsite_cover %>% filter( type == 'upper', microsite == 'not_hummock', abu > 0), comp_grid, fd)

## global comparison ----------------------------------------------- # 

plot_correlation(global_pairs)

## upper sites comparison ----------------------------------------------- # 

plot_correlation(upper_pairs)

## lower sites comparison ----------------------------------------------- # 

plot_correlation(lower_pairs)

## hummock sites comparison ----------------------------------------------- # 

plot_correlation(hummock_pairs)

## non-hummock upper sites comparison ----------------------------------------------- # 

plot_correlation(non_hummock_pairs)

##  Divide by seed size to rescale cover as individual abundance 

library(sedgwicktraits)

size <- 
  sedgwicktraits %>% 
  ungroup() %>% 
  select(USDA_symbol, seed_mass, projected_area) 

global_abu <- 
  global_cover %>% 
  left_join(size) %>% 
  mutate( abu = abu/seed_mass) %>% 
  select( USDA_symbol, abu)

type_abu <- 
  type_cover %>% 
  ungroup %>% 
  left_join(size) %>% 
  mutate( abu = abu/seed_mass) %>% 
  select(USDA_symbol, type, abu)

microsite_abu <- 
  microsite_cover %>% 
  ungroup %>% 
  left_join(size) %>% 
  mutate( abu = abu/seed_mass) %>% 
  select(USDA_symbol, type, microsite, abu)

global_pairs2 <- join_cover_fd(global_abu, comp_grid, fd)
upper_pairs2 <- join_cover_fd(type_abu %>% filter( type == 'upper', abu > 0 ) ,comp_grid, fd)
lower_pairs2 <- join_cover_fd(type_abu %>% filter( type == 'lower', abu > 0 ), comp_grid, fd)
hummock_pairs2 <- join_cover_fd(microsite_abu %>% filter( type == 'upper', microsite == 'hummock', abu > 0 ), comp_grid, fd)
non_hummock_pairs2 <- join_cover_fd(microsite_abu %>% filter( type == 'upper', microsite == 'not_hummock', abu > 0 ), comp_grid, fd)

plot_correlation(global_pairs2)
plot_correlation(upper_pairs2)
plot_correlation(lower_pairs2)
plot_correlation(hummock_pairs2)
plot_correlation(non_hummock_pairs2)

## Divide by canopy area to rescale cover as number of plants 

global_abu <- 
  global_cover %>% 
  left_join(size) %>% 
  mutate( abu = abu/projected_area) %>% 
  select( USDA_symbol, abu)

type_abu <- 
  type_cover %>% 
  ungroup %>% 
  left_join(size) %>% 
  mutate( abu = abu/projected_area) %>% 
  select(USDA_symbol, type, abu)

microsite_abu <- 
  microsite_cover %>% 
  ungroup %>% 
  left_join(size) %>% 
  mutate( abu = abu/projected_area) %>% 
  select(USDA_symbol, type, microsite, abu)


global_pairs3 <- join_cover_fd(global_abu, comp_grid, fd)
upper_pairs3 <- join_cover_fd(type_abu %>% filter( type == 'upper', abu > 0 ) ,comp_grid, fd)
lower_pairs3 <- join_cover_fd(type_abu %>% filter( type == 'lower', abu > 0 ), comp_grid, fd)
hummock_pairs3 <- join_cover_fd(microsite_abu %>% filter( type == 'upper', microsite == 'hummock', abu > 0 ), comp_grid, fd)
non_hummock_pairs3 <- join_cover_fd(microsite_abu %>% filter( type == 'upper', microsite == 'not_hummock', abu > 0 ), comp_grid, fd)

plot_correlation(global_pairs3)
plot_correlation(upper_pairs3)
plot_correlation(lower_pairs3)
plot_correlation(hummock_pairs3)
plot_correlation(non_hummock_pairs3)

# use emperical relationship to rescale cover to abundance 

load('output/plants_per_cover.rda')

global_abu <- 
  global_cover %>% 
  left_join(plants_per_cover, by = c('USDA_symbol' = 'spp')) %>%
  mutate( abu = abu*n_per_cover) %>% 
  select( USDA_symbol, abu)

type_abu <- 
  type_cover %>% 
  ungroup %>% 
  left_join(plants_per_cover, by = c('USDA_symbol' = 'spp')) %>%
  mutate( abu = abu*n_per_cover) %>% 
  select(USDA_symbol, type, abu)

microsite_abu <- 
  microsite_cover %>% 
  ungroup %>% 
  left_join(plants_per_cover, by = c('USDA_symbol' = 'spp')) %>%
  mutate( abu = abu*n_per_cover) %>% 
  select(USDA_symbol, type, microsite, abu)


global_pairs4 <- join_cover_fd(global_abu, comp_grid, fd)
upper_pairs4 <- join_cover_fd(type_abu %>% filter( type == 'upper', abu > 0 ) ,comp_grid, fd)
lower_pairs4 <- join_cover_fd(type_abu %>% filter( type == 'lower', abu > 0 ), comp_grid, fd)
hummock_pairs4 <- join_cover_fd(microsite_abu %>% filter( type == 'upper', microsite == 'hummock', abu > 0 ), comp_grid, fd)
non_hummock_pairs4 <- join_cover_fd(microsite_abu %>% filter( type == 'upper', microsite == 'not_hummock', abu > 0 ), comp_grid, fd)


plot_correlation(global_pairs4)
plot_correlation(upper_pairs4)
plot_correlation(lower_pairs4)
plot_correlation(hummock_pairs3)
plot_correlation(non_hummock_pairs4)

plot_correlation( upper_pairs[ upper_pairs$pair %in% upper_pairs4$pair, ] )
plot_correlation( upper_pairs2[ upper_pairs2$pair %in% upper_pairs4$pair, ] )
plot_correlation( upper_pairs3[ upper_pairs3$pair %in% upper_pairs4$pair, ] )
plot_correlation(upper_pairs4)


# print all plots to pdf 

pdf( file = 'figures/pairwise_abundance_v_fitness_diff.pdf',  width = 7, height = 5 )

print( 
plot_correlation(global_pairs) + 
  ggtitle('All Plots') + 
  xlab('log ratio of abundances (cover)') + 
  ylab('log ratio of fitness differences') 

) 
print( 
plot_correlation(upper_pairs) + 
  ggtitle('Upper Plots') + 
  xlab('log ratio of abundances (cover)') + 
  ylab('log ratio of fitness differences')
)
print( 
  
plot_correlation(lower_pairs) + 
  ggtitle('Lower Plots') + 
  xlab('log ratio of abundances (cover)') + 
  ylab('log ratio of fitness differences')

)
print( 
  
plot_correlation(hummock_pairs) + 
  ggtitle('Hummock Plots') + 
  xlab('log ratio of abundances (cover)') + 
  ylab('log ratio of fitness differences')
)
print( 
  
plot_correlation(non_hummock_pairs) + 
  ggtitle('Non-Hummock Plots') + 
  xlab('log ratio of abundances (cover)') + 
  ylab('log ratio of fitness differences')
)
print( 
  

plot_correlation(global_pairs2) + 
  ggtitle('All Plots') + 
  xlab('log ratio of abundances (cover/seed_mass)') + 
  ylab('log ratio of fitness differences')
)
print( 
  
plot_correlation(upper_pairs2) + 
  ggtitle('Upper Plots') + 
  xlab('log ratio of abundances (cover/seed_mass)') + 
  ylab('log ratio of fitness differences')
)
print( 
  
plot_correlation(lower_pairs2) + 
  ggtitle('Lower Plots') + 
  xlab('log ratio of abundances (cover/seed_mass)') + 
  ylab('log ratio of fitness differences')
)
print( 
  
plot_correlation(hummock_pairs2) + 
  ggtitle('Hummock Plots') + 
  xlab('log ratio of abundances (cover/seed_mass)') + 
  ylab('log ratio of fitness differences')
)
print( 
  
plot_correlation(non_hummock_pairs2) + 
  ggtitle('Non-Hummock Plots') + 
  xlab('log ratio of abundances (cover/seed_mass)') + 
  ylab('log ratio of fitness differences')
)
print( 
  

plot_correlation(global_pairs3) + 
  ggtitle('All Plots') + 
  xlab('log ratio of abundances (cover/canopy_area)') + 
  ylab('log ratio of fitness differences')
)
print( 
  
plot_correlation(upper_pairs3) + 
  ggtitle('Upper Plots') + 
  xlab('log ratio of abundances (cover/canopy_area)') + 
  ylab('log ratio of fitness differences')
)
print( 
  
plot_correlation(lower_pairs3) + 
  ggtitle('Lower Plots') + 
  xlab('log ratio of abundances (cover/canopy_area)') + 
  ylab('log ratio of fitness differences')
)
print( 
  
plot_correlation(hummock_pairs3) + 
  ggtitle('Hummock Plots') + 
  xlab('log ratio of abundances (cover/canopy_area)') + 
  ylab('log ratio of fitness differences')
)
print( 
  
plot_correlation(non_hummock_pairs3) + 
  ggtitle('Non-Hummock Plots') + 
  xlab('log ratio of abundances (cover/canopy_area)') + 
  ylab('log ratio of fitness differences') 
)


dev.off()
