library(tidyverse)
library(mvtnorm)
library(mgcv)
#source("feas_dom.R")  #load the toolbox from Guideline

#function that generates a uniform initial conditons and solve the LV model
#inputs: alpha = interaction matrix; S = species numbers; r = vector of intrinsic growth rates
#output: out = indexes of survived species
##>>LV is revised to include specified initial conditions
LV <- function(alpha,N,r) {
  survival <- 0
    delta_t <- 0.01 #time step
    time_step <- seq(0,300,by=delta_t)
    threshold <- 0.0000001
    
    N0 <- N/ sum(N)
    
    d <- sqrt(sum(r^2))
    r <- r / d
    
    parms <- list(r=r,alpha = alpha) ##ODE
    model <- function(t,N,parms){
    dN <- N * (parms$r + parms$alpha %*% N);
    list(dN)
    }
    
    e <- try(sol <- ode(N0,time_step,model,parms), silent = F) ## T is short for TRUE?
    index <- 0
    
    if(class(e)!="try-error"){
      for(z in 1:nrow(alpha)){   #replace S with nrow(alpha)
          if (sol[nrow(sol),z+1] > threshold){
            #survived <- survived + 1
            index <- c(index,z)
          }
      }
    }
    
  out <- index[-1]
  return(out)
}

# function that computes the normalized feasibility from an interaction matrix
# inputs: alpha = interaction matrix
# output: out = the normalized feasibility
Omega <- function(alpha) {
  S <- nrow(alpha)
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    out <- d[1]^(1 / S)
    return(out)
  }
  #   if (length(which(diag(alpha) == 0)) == 0) {
  #     Sigma <- chol2inv(alpha, size = NCOL(alpha), LINPACK = FALSE)
  #     return(omega(S, Sigma))
  #   }
  #   else {
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  if (f(alpha) == FALSE) {
    return(0)
  }
  else {
    Sigma <- solve(t(alpha) %*% alpha)
    return(omega(S, Sigma))
  }
}

#function that generates a random interaction matrix
#inputs: num = community size; stren = scale of interaction strength; conne = connectance
#output: Inte = random interaction matrix
generate_Interaction_matrix <- function(num, stren, conne){
  ## initialize interaction matrix
  Inte <- (rnorm(num*num, mean = 0, sd = stren))
  
  ## sample from several ones and several zeros, generates random permutation
  zeroes <- sample(c(rep.int(1,floor(num*num*conne)), rep.int(0,(num*num-floor(num*num*conne))))) 
  
  ## assign zeros
  Inte[which(zeroes==0)] <- 0
  
  Inte <- matrix(Inte, ncol = num, nrow = num)
  diag(Inte) <- -1
  return(Inte)
}

#function that checks whether the chosen parameterization satisfies feasibility
#inputs: A = interaction matrix; r = vector of intrinsic growth rates
#output: unfeasible = 0, feasible = 1
check_feasibility <- function(A, r){
  if(sum(solve(-A, r)<0)==0) return(1)
  else return(0)
}

#function that checks whether the chosen parameterization satisfies local stability
#inputs: A = interaction matrix; r = vector of intrinsic growth rates
#output: unstable = 0; stable= 1
check_stability <- function(A, r){
  N <- solve(A, -r)
  jacobian <- A %*% diag(N) 
  auxi <- sum(Re(eigen(jacobian)$values) < -1e-8)
  if(auxi==nrow(A)) return(1) #To comveniently check sta for subcommnuity, change ==num to ==nrow(A)
  else return(0)
}

#function that generates the spanning vectors of the feasible cone
#inputs: A = interaction matrix; num = number of species
#output: generates the matrix where each column is a spanning vector
spanned_vectors <- function(A, num){
  G <- matrix(0, ncol=ncol(A), nrow=nrow(A))
  for(k in 1:ncol(A)) G[,k] <- -A[,k]/sqrt(sum(A[,k]^2))
  return(G) ## modified: from G to return(G)
}

#function that generates all the elements of the vector of intrinsic growth rate with the same value
#inputs: num = number of species
#output: vector of intrinsic growth rates with 'fixed' parameterization 
parameterization_fixed <- function(num){
  rep.int(1, num)
}

#function that generates the vector of intrinsic growth rates randomly with no other constraints
#inputs: num = number ofspecies
#output: vector of intrinsic growth rates with 'random' parameterization 
parameterization_random <- function(num){
  runif(num, -1, 1)
}

#function that generates the vector of intrinsic growth rates randomly inside the feasibility domain
#inputs: A = interaction matrix; num = number of species
#output: vector of intrinsic growth rates with 'feasible' parameterization 
parameterization_feasible <- function(A, num){
  G <- spanned_vectors(A, num)
  lambda <- runif(num, min = 0, max = 1)
  lambda <- lambda/sum(lambda)
  growth <- G %*% matrix(lambda, ncol=1) %>%
    as.vector()
}

#function that generates the vector of intrinsic growth rates located at the geometric centroid of the feasibility domain
#inputs: A = interaction matrix; num = number of species
#output: vector of intrinsic growth rate with 'centroid' parameterization 
parameterization_center <- function(A, num){
  G <- spanned_vectors(A, num)
  lambda <- rep.int(1, num)
  lambda <- lambda/sum(lambda)
  growth <- G %*% matrix(lambda, ncol=1) %>%
    as.vector()
}

convert2names <- function(vec){
  str = ""
  i = 1
  while (!is.na(vec[i])){
    str = paste0(str,as.character(vec[i]))
    i = i+1
  }
  str
}