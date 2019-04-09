rm(list = ls())
library(stringr)
library(sedgwickspecies)
library(tidyverse)

source('helper_functions.R')

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
  left_join(sedgwick_plants, by = c('species' = 'prior_name')) %>% 
  select(prior_code, lambda, g, s) %>% 
  left_join( ( alphas %>% 
                 rename('prior_code' = X1) ) )


N <- rep(1, 18)
nrow(demo_pars)

A <- as.matrix(demo_pars %>% select(AGHE:SIGA)) 
A[is.na(A)] <- 0 
g <- as.numeric(demo_pars$g)
lambda <- as.numeric(demo_pars$lambda)
s <- as.numeric(demo_pars$s)

myset <- c(2,3)
nsp <- length(myset)
init_abu <- c(1, 5000)
time <- 200
A_temp <- A[myset, myset]
lambda_temp <- lambda[myset]
g_temp <- g[myset]
s_temp <- s[myset]
N1 <- matrix(NA, time, length(myset))
N1[1, ] <- init_abu
N2 <- N1

for( i in 2:time){ 
  N1[i,] <- N1[i-1,]*(s_temp*(1 - g_temp) + g_temp*(lambda_temp/(1 + A_temp%*%(g_temp*N1[i-1, ]) ) ))
  for( j in 1:length(myset))  { 
    N2[i,j] <- N2[i-1,j]*(s_temp[j]*(1 - g_temp[j]) + g_temp[j]*(lambda_temp[j]/(1 + sum(A_temp[j,]*g_temp*N2[i-1,]))))
  }
}

matplot(N1, type = 'l')
matplot(N2, type = 'l')

eta <- (g_temp*lambda_temp)/(1 - s_temp*(1 -g_temp))
x <- solve(A_temp) %*% ( eta - 1)
x
tail(N1)[1,]*g_temp
N <- x*(g_temp)^-1
tail(N1)[1,]
x

J11 <- (1-g_temp[1])*s_temp[1] - (-lambda_temp[1]*A_temp[1,1]*g_temp[1]^2*N[1])/(1 + sum(A_temp[1,]*g_temp*N))^2 + (lambda_temp[1]*g_temp[1])/(1 + sum(A_temp[1,]*g_temp*N))
J12 <- (-lambda_temp[1]*A_temp[1,2]*g_temp[2]*g_temp[1]*N[1])/(1 + sum(A_temp[1,]*g_temp*N))^2
J21 <- (-lambda_temp[2]*A_temp[2,1]*g_temp[1]*g_temp[2]*N[2])/(1 + sum(A_temp[2,]*g_temp*N))^2
J22 <- (1-g_temp[2])*s_temp[2] - (-lambda_temp[2]*A_temp[2,2]*g_temp[2]^2*N[2])/(1 + sum(A_temp[2,]*g_temp*N))^2 + (lambda_temp[2]*g_temp[2])/(1 + sum(A_temp[2,]*g_temp*N))

Jmat <- matrix(c(J11, J12, J21, J22), 2,2, byrow = T)

J <- get_Jacobian(x = tail(N1)[1,], A_temp, lambda_temp, g_temp, s_temp)
eigen(Jmat)
eigen(J)




tail(N2)

myset <- 1:2
nsp <- length(myset)
init_abu <- c(1,1)
time <- 10000
alpha <- matrix( c(1e-2, 5.2e-3, 
                   4.75e-3, 1.1e-2), 2, 2, byrow = T)

lambda <- rep(2,nsp)
N1 <- matrix(NA, time, length(myset))
N1[1, ] <- init_abu
N2 <- N1

for( i in 2:time){ 
  N1[i, ] <- N1[i-1,]*lambda/(1 + alpha%*%N1[i-1, ])
  
  for( j in 1:length(myset))  { 
    N2[i,j] <- N2[i-1,j]*lambda[j]/(1 + sum(alpha[j, ]*N2[i-1, ]))
  }
}

matplot(N1, type = 'l')
matplot(N2, type = 'l')
tail(N1)

x <- solve(alpha) %*% (1 - 1/lambda)
lambda*x

lv_jacobian <- function(x, A, lambda){ 
  
  J <- matrix(NA, length(x), length(x))
  d <- diag(J)  
  for( i in 1:length(x)){
    for( j in 1:length(x)){ 
      J[i,j] <- (-lambda[i]*A[i,j]*x[i])/((1 + sum(A[i,]*x))^2)
    }
  }
  
  for( i in 1:length(x)){ 
    d[i] <-  (-lambda[i]*A[i,i]*x[i])/((1 + sum(A[i,]*x))^2) + (lambda[i])/(1 + sum(A[i,]*x))
  }
  
  diag(J) <- d
  return(J)
}  
x[,1]

J <- lv_jacobian(x[,1]*lambda, alpha, lambda )
eigen(J)

