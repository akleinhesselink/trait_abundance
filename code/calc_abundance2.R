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

A_mean <- list(NA)
eta_i <- list(NA)

for( i in 1:ncol(A)){ 
  a_ij <- mean(A[i, -i])
  a_ji <- mean(A[-i, i])
  a_ii <- diag(A)[i]
  a_jj <- mean(diag(A)[-i])
  A_mean[[i]] <- matrix( c(a_ii, a_ij, a_ji, a_jj), 2, 2, byrow = T)
  g_i <- c(g[i], mean(g[-i]))
  lambda_i <- c(lambda[i], mean(lambda[-i]))
  s_i <- c(s[i], mean(s[-i]))
  eta_i[[i]] <- (g_i*lambda_i)/(1 - s_i*(1 - g_i))
}

solution <- list(NA)

for( i in 1:length(A_mean)){ 
  solution[[i]] <- solve(A_mean[[i]]) %*% (eta_i[[i]] - 1)
}

solution

solutions_df <- data.frame( do.call(rbind, lapply( solution, function(x) t(x) ) ) )

solutions_df %>% 
  mutate( feasible = X1 > 0 & X2 > 0 ) %>% 
  mutate( focal_wins = X1 > X2, 
          comp_wins = X2 > X1) %>%
  mutate( prior_code = demo_pars$prior_code)

#Jac <- lapply( mysets, function(x){ y <- rep(list(NA), ncol(x)); y} )

mag_fac <- exp(seq(0,10, by = 0.5))

all_sp <- all_sp2 <- list(NA)

plot(log(eta))

scale_eta <- scale(log(eta))

attributes(scale_eta)[1]

scale_reduction <- 1/exp(seq(0,3, by = 0.1))

for( i in 1:length(scale_reduction)){ 
  
  log_eta <- scale_eta*attributes(scale_eta)[[3]]*(red_factor[i]) + attributes(scale_eta)[[2]]
  eta_temp <- exp(log_eta)
  A_temp <- A[1:18, 1:18]
  all_sp[[i]] <- solve(A_temp) %*% ( eta_temp - 1)
  
}

matplot( t(do.call( cbind, all_sp)) , type = 'l')

for( i in 1:length(mag_fac)){ 
  A_temp <- A[1:18, 1:18]
  eta_temp <- eta[1:18]
  diag(A_temp) <- diag(A_temp)*mag_fac[i]
  all_sp2[[i]] <- solve(A_temp) %*% ( eta_temp - 1)
}

matplot( t(do.call( cbind, all_sp2)) , type = 'l')

