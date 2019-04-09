
get_Jacobian <- function(x, A, lambda, g, s){
  
  # x is vector of equilibrium abundances
  # A is interaction matrix for set
  # lambda is per capita seed production vector
  # g is germination vector
  # s is survival vector
  
  J <- matrix(NA, length(x), length(x))
  d <- diag(J)  
  for( i in 1:length(x)){
    for( j in 1:length(x)){ 
      J[i,j] <- (-lambda[i]*A[i,j]*g[j]*g[i]*x[i])/(1 + sum(A[i,]*g*x))^2
    }
  }
  
  for( i in 1:length(x)){ 
    d[i] <-  (1-g[i])*s[i] - (-lambda[i]*A[i,i]*g[i]^2*x[i])/(1 + sum(A[i,]*g*x))^2 + (lambda[i]*g[i])/(1 + sum(A[i,]*g*x))
  }
  
  diag(J) <- d
  
  return(J) 
}




fitness_difference <- function(alpha, lambda, g, s){ 
  
  eta <- (g*lambda)/(1 - (1 - g)*s)
  
  as.numeric(((eta[2] - 1)/(eta[1] - 1))*sqrt((alpha[1,2]*alpha[1,1])/(alpha[2,2]*alpha[2,1])))
}


rho <- function( alpha ) { 
  
  sqrt( (alpha[1,2]*alpha[2,1] )/(alpha[1,1]*alpha[2,2]) )  
  
}
