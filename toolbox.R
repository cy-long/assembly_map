library(tidyverse)
library(mvtnorm)
library(mgcv)
library(binaryLogic)
##------------ Original funcs ------------

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

# function that computes the raw Omega value of feasibility from an interaction matrix
# inputs: alpha = interaction matrix
# output: out = the raw Omega value of feasibility 
Omega <- function(alpha) {
  if(is.null(nrow(alpha))) { return(1/2) }
  else{
    S <- nrow(alpha)
    omega <- function(S, Sigma) {
      m <- matrix(0, S, 1)
      a <- matrix(0, S, 1)
      b <- matrix(Inf, S, 1)
      d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
      #out <- d[1]^(1 / S)
      out <- d[1]
      return(out)
    }
    f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
    if (f(alpha) == FALSE) {
      return(0)
    } else {
      Sigma <- solve(t(alpha) %*% alpha) #Inverse
      return(omega(S, Sigma))
    }
  }
}

# function that normalizes a vector in the L2 norm
# inputs: a = the orignal vector
# output: the normalized vector
normalization <- function(a) {
  a / sqrt(sum(a^2))
}

# function that normalizes the spanning vectors of the feasibility domain in the L2 norm
# inputs: alpha = interaction matrix
# output: Span = the normalized spanning vectors
span_vectors <- function(alpha) {
  Span <- matrix(0, ncol = ncol(alpha), nrow = nrow(alpha))
  for (k in 1:ncol(alpha)) {
    Span[, k] <- -alpha[, k] / sqrt(sum(alpha[, k]^2))
  }
  Span
}

# function that computes all the extreme points that belong to original vertexes
# inputs: A = one interaction matrix, B = another interaction matrix
# output: inside_vertex = all the extreme points that belong to original vertexes
inside_vertex_detection <- function(A, B) {
  SpanA <- span_vectors(A)
  SpanB <- span_vectors(B)
  # to determine whether a vertex of one cone is inside another cone or not.
  inside_detection <- function(Span, vector) {
    lambda <- solve(Span, vector)
    if (sum(lambda >= -1e-10) == length(lambda)) {
      return(1)
    } else {
      return(0)
    }
  }
  inside_vertex <- list()
  l <- 1
  for (i in 1:ncol(B)) {
    auxi <- inside_detection(SpanA, SpanB[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanB[, i]
      l <- l + 1
    }
  }
  for (i in 1:ncol(A)) {
    auxi <- inside_detection(SpanB, SpanA[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanA[, i]
      l <- l + 1
    }
  }
  return(inside_vertex)
}

# function that computes all the extreme points generated that are generated by the intersections of the cones
# inputs: S = one interaction matrix, M = another interaction matrix
# output: intersection_vertex = all the extreme points that are generated by the intersections of the cones
intersection_vertex_detection <- function(S, M) {
  num <- ncol(S)
  combination_S <- combn(1:ncol(S), 2)
  combination_M <- combn(1:ncol(S), (num - 1))
  Span_S <- span_vectors(S)
  Span_M <- span_vectors(M)
  
  border_M <- list()
  extreme_point_M <- list()
  
  for (i in 1:ncol(M)) {
    coeff_matrix <- matrix(1, ncol = num, nrow = num)
    for (j in 1:(num - 1))
      coeff_matrix[j, ] <- Span_M[, combination_M[j, i]]
    coeff_vector <- c(rep(0, num - 1), 1)
    border_M[[i]] <- solve(coeff_matrix, coeff_vector)
    extreme_point_M[[i]] <- t(coeff_matrix)[1:(num - 1), 1:(num - 1)]
  }
  
  inside_face_detection <- function(extreme_point, test_vector) {
    lambda <- solve(extreme_point, test_vector)
    if (sum(lambda >= -1e-10) == length(lambda)) {
      return(1)
    } else {
      return(0)
    }
  }
  
  l <- 1
  intersection_vertex <- list()
  side <- c()
  for (i in 1:ncol(combination_S)) {
    vertex_1 <- Span_S[, combination_S[1, i]]
    vertex_2 <- Span_S[, combination_S[2, i]]
    
    for (j in 1:length(border_M)) {
      n1 <- sum(vertex_1 * border_M[[j]])
      n2 <- sum(vertex_2 * border_M[[j]])
      
      auxi <- n1 * n2
      if (auxi < -1e-10) {
        lambda <- n2 / (n2 - n1)
        possible <- lambda * vertex_1 + (1 - lambda) * vertex_2
        #>det() only accept matrix, for 2*2, calculation must be in a degenerated form
        #>just change here, is it legal?
        emj <- extreme_point_M[[j]]
        if ((!is.null(ncol(emj)) && det(emj) != 0) || (is.null(ncol(emj)) && emj != 0)) {
          auxi2 <- inside_face_detection(extreme_point_M[[j]], possible[1:(num - 1)])
          if (auxi2 == 1) {
            intersection_vertex[[l]] <- possible
            side[l] <- j
            l <- l + 1
          }
        }
      }
    }
  }
  
  if (length(intersection_vertex) > 0) {
    for (i in 1:length(intersection_vertex)) {
      intersection_vertex[[i]] <- normalization(intersection_vertex[[i]])
    }
  }
  
  return(intersection_vertex)
}

# function that computes all the extreme points
# inputs: A = one interaction matrix, B = another interaction matrix
# output: out = all the extreme points that generate the intersection region
vertex_detection <- function(A, B) {
  num <- ncol(A)
  inside_vertex <- inside_vertex_detection(A, B)
  intersection_vertex <- intersection_vertex_detection(A, B)
  
  # combine the two vertex lists
  if (length(inside_vertex) > 0) {
    vertex <- matrix(unlist(inside_vertex), nrow = num, byrow = FALSE)
  } else {
    vertex <- matrix(0, nrow = num, ncol = 2)
  }
  if (length(intersection_vertex) > 0) {
    vertex <- cbind(vertex, matrix(unlist(intersection_vertex), nrow = num, byrow = FALSE))
  }
  
  # delete the points that are nonzero due to numerical error
  delete_zeroes <- c()
  for (i in 1:ncol(vertex)) {
    if (near(sum(vertex[, i]^2), 0)) {
      delete_zeroes <- c(delete_zeroes, i)
    }
  }
  if (length(delete_zeroes) > 0) vertex <- vertex[, -delete_zeroes]
  
  
  # delete the same ones
  if (length(vertex) > num) {
    for (test in 1:ncol(vertex)) {
      vertex[, test] <- normalization(vertex[, test])
    }
    delete_duplicates <- c()
    for (i in 1:(ncol(vertex) - 1)) {
      for (j in (i + 1):ncol(vertex)) {
        if (sum(near(vertex[, i], vertex[, j])) == nrow(vertex)) {
          delete_duplicates <- c(delete_duplicates, j)
        }
      }
    }
    if (length(delete_duplicates) > 0) vertex <- vertex[, -unique(delete_duplicates)]
  }
  return(vertex)
}

# function that computes all the extreme points
# inputs: p = the list of all extreme points
# output: out = all the extreme points that generate the intersection region
partitionize <- function(p) {
  triangulize <- function(p) {
    set.seed(100)
    inside_detection <- function(Span, vector) {
      lambda <- solve(Span, vector)
      if (sum(lambda >= -1e-10) == length(lambda)) {
        return(1)
      } else {
        return(0)
      }
    }
    border <- function(position) {
      coeff_matrix <- t(p[, position])
      coeff_matrix <- rbind(coeff_matrix, 1)
      coeff_vector <- c(rep(0, NR - 1), 1)
      if (abs(det(coeff_matrix)) > 1e-10) {
        border_M <- solve(coeff_matrix, coeff_vector)
      } else {
        border_M <- rep(0, NR)
      }
      return(border_M)
    }
    side_determine <- function(all_del) {
      # choose a candidate
      candidate <- sample(setdiff(c(1:NC), all_del), NR - 1)
      candidate_face <- border(candidate)
      
      # check which side other points are
      side <- rep(0, ncol(p))
      for (j in 1:length(side)) {
        side[j] <- sum(p[, j] * candidate_face)
      }
      side <- ifelse(abs(side) < 1e-10, 0, side)
      left_side <- which(side < 0)
      right_side <- which(side > 0)
      left_side <- setdiff(left_side, all_del)
      right_side <- setdiff(right_side, all_del)
      
      abundant_side <- which(side == 0)
      abundant_side <- setdiff(abundant_side, all_del)
      abundant_side <- setdiff(abundant_side, candidate)
      
      if (length(left_side) == 1) {
        all_del <- c(all_del, left_side)
        all_del <- unique(c(all_del, abundant_side))
        return(list(all_del, c(candidate, left_side)))
      }
      else if (length(right_side) == 1) {
        all_del <- c(all_del, right_side)
        all_del <- unique(c(all_del, abundant_side))
        return(list(all_del, c(candidate, right_side)))
      }
      else {
        return(c())
      }
    }
    isolation_finder <- function(all_del) {
      try <- side_determine(all_del)
      while (length(try) == 0) {
        try <- side_determine(all_del)
      }
      return(try)
    }
    NR <- nrow(p)
    NC <- ncol(p)
    partition <- matrix(0, nrow = 1, ncol = NR)
    all_del <- c()
    while (length(setdiff(c(1:NC), all_del)) >= NR) {
      result <- isolation_finder(all_del)
      all_del <- result[[1]]
      if (length(result) > 1) partition <- rbind(partition, result[[2]])
    }
    # partition
    partition <- partition[-1, ]
    partition <- matrix(partition, ncol = NR)
    # partition <- partition[!duplicated(partition[,]),]
    return(partition)
  }
  # if the number of vertexes is the same as the dimension
  if (nrow(p) == ncol(p)) {
    matrix(1:(nrow(p)), nrow = 1)
  }
  # else it is a triangulization problem
  else {
    return(triangulize(p))
  }
}

# function that computes the full volume
# inputs: partition = triangulation of the intersection region, vertex = all the extreme points
# output: the normalize feasibility of the intersection region
total_volume <- function(partition, vertex) {
  vol <- c()
  num <- ncol(partition)
  if (length(partition) == num) {
    auxi <- vertex[, partition]
    if (near(det(auxi), 0)) {
      vol <- 0
    } else {
      vol <- Omega(auxi)
    }
  } else {
    for (i in 1:nrow(partition)) {
      auxi <- vertex[, partition[i, ]]
      if (near(det(auxi), 0)) {
        vol[i] <- 0
      } else {
        vol[i] <- Omega(auxi)
      }
    }
  }
  sum(vol^num)^(1 / num)
}

sampling_overlap <- function(A,B) {
  mat <- solve(B) %*% A
  Nsample <- 10^5
  abundance_all <- rmvnorm(n = Nsample, mean = rep(0, nrow(A))) %>% 
    {abs(./sqrt(rowSums(.^2)))}
  get_feasibility <- function(N_A){
    N_B <- mat %*% matrix(N_A,ncol=1) %>% c()
    if_else(sum(N_B >= -1e-10) == length(N_B), 1, 0) 
  }
  percent <- 1:Nsample %>% 
    map_dbl(~get_feasibility(abundance_all[.x,])) %>% 
    mean()
  
  #Omega(A) * percent^(1/nrow(A))
  #Raw overlap, not the normailzed one
  return(Omega(A) * percent) 
}

# function that computes the overlap of two feasibility domains
# inputs: A = one interaction matrix, B = another interaction matrix
# output: volume_overlap = the raw feasibility of the intersection region
Omega_overlap <- function(A, B) {
  num <- nrow(A)
  
  overlap_vertex <- vertex_detection(A, B)
  
  if (qr(overlap_vertex)$rank < num) {
    volume_overlap <- 0
  } else {
    partition <- partitionize(overlap_vertex)
    volume_overlap <- total_volume(partition, overlap_vertex)
  }
  
  volume_overlap <- sampling_overlap(A,B)
  volume_overlap
}

##------------ Added funcs ------------

# function that convert numeric vector to strings
# input: vec = numeric vector, such as c(1,2), oputput: string, such as "12"
convert2names <- function(vec) {
  str = ""
  i = 1
  while (!is.na(vec[i])){
    str = paste0(str,as.character(vec[i]))
    i = i+1
  }
  str
}

convert2sets <- function(str) {
  as.numeric(strsplit (str,"")[[1]] )
}

# function that computes the overlap of two feasibility domains with different dimensions
# inputs: A = one interaction matrix, B = another interaction matrix, 
# comp_A = composition of community A, comp_B = composition of community B
# output: volume_overlap = the normalized feasibility of the intersection region
Omega_overlap_ext <- function(A, B, comp_A, comp_B) {
  if (all(is.element(comp_A, comp_B))) {
    
    #create a location vector to store indexes where insertion happens
    loc <- c(rep(0, length((comp_B))))
    for (i in 1:length(comp_B)) {
      if(!is.element(comp_B[i], comp_A)) {
        loc[i] <- 1
      }
    }
    
    #create a matrix of all possible diagnal combinations of +-
    generate_diag_string <- function(N){
      record <- matrix(0, 2^N, N)
      for (i in 1:(2^N-1)){
        record[i+1, (N+1-length(as.binary(i))):N] <- as.binary(i)
      }
      record = apply(record,c(1,2), function(x) {if(x <= 0){return(-1)} else return(1)})
      return(record)
    }
    diag_mat <- generate_diag_string(sum(loc))
    
    #compute Omega_overlap
    #Omega_ext <- 0
    Omega_ext <- c(rep(0, (2 ^ sum(loc))))
    for (i in 1:(2 ^ sum(loc))){
      A_ext <- matrix_scatter(A, loc, diag_mat[i, ])
      #Omega_ext = Omega_ext + Omega_overlap(B,A_ext)
      #Omega_ext[i] = Omega_overlap(A_ext, B)
      Omega_ext[i] = Omega_overlap(B, A_ext)
    }
    return(sum(Omega_ext))
    
  } else if (all(is.element(comp_B, comp_A))) {
    Omega_overlap_ext(B, A, comp_B, comp_A) #change the variations and do recursion
  } else {
    warning("connecting two sub-communities that are not subset of each other")
  }
}

# function that rearrange matrix elements based on a control vector
# input: mat = one interaction matrix, usually at a lower dimension;
# input: loc = location vector indicating the index of the extended matrix
#where zeros and diagonal elements are inserted, coded in sequence of 0/1
# input: apnd_diag = diagonal elements series
# output: the extended matrix
matrix_scatter <- function(mat, loc, apnd_diag) {
  if(is.null(ncol(mat))) {
    n_mat <- 1
  } else {
    n_mat <- ncol(mat)
  }
  mov <- c(rep(0, n_mat))
  j <- 0
  #sum all movement steps needed to insert the new species
  for (i in 1:length(loc)) {
    if(loc[i] == 0) {
      j <- j + 1
      mov[j] <- sum(loc[1:i])
    }
  }
  if(j != n_mat) {
    warning("matrix_scatter: discrepancy between ncol(mat) and length(loc)")
  }
  
  mat_ext <- matrix(0, ncol = length(loc), nrow = length(loc))
  # scatter the original elements
  if(n_mat > 1){
    for (i in 1:n_mat) {
      for (j in 1:n_mat) {
        mat_ext[i+mov[i], j+mov[j]] <- mat[i,j]
      }
    }
  } else {mat_ext[1+mov[1], 1+mov[1]] <- mat}
  # set new diagonal element
  ##>for other modification policy, rewrite this part
  d <- 1
  for (k in 1:length(loc)) {
    if(loc[k] == 1){
      mat_ext[k,k] <- apnd_diag[d]
      d <- d + 1
    }
  }
  return(mat_ext)
}

# function that adds one species to a sub-commnutiy
# inputs: comm = one vector that represents community components, num = total number of species, 
#augm = adding number of species
# output: extended communities in a matrix, where each col has one more species
extend_communities <- function(comm, num, augm) {
  if(isTRUE(augm > num - length(comm))) warning("augmentation exceeds total number of species")
  else{
    s <- length(comm)
    candid_sp <- setdiff(c(combn(num, 1)), comm)
    
    if(length(candid_sp) == 1) comm_ext <- matrix(sort(c(comm, candid_sp)), ncol = 1)
    else{
      comm_ext <- matrix(0, nrow = length(comm)+augm, ncol = choose(length(candid_sp), augm))
      candid_com <- combn(candid_sp, augm)
      for (i in 1:ncol(candid_com)) {
        comm_ext[, i] <- sort(c(comm, candid_com[,i]))
      }
    }
    return(comm_ext)
  }
}

# function that check if a vector is from one column of one matrix
# input: vec = the vector, mat = the matrix
# output: TRUE or FALSE
is.vec_in_mat <- function(vec,mat){
  out <- apply(mat, 2, function(x, y){isTRUE(all.equal(x, y))}, vec)
  any(out)
}

# function that normailze each row of a matrix
# input: mat = the original matrix
# output: the normalized one, with each row sum equals to one
norm_row_sum <- function(mat){
  t(apply(mat,1,function(x) x/sum(x)))
}

# function that returns the cartesian product from sets (filter zero elements)
# input: mat = matrix whose rows are original sets, using zeros as placeholders
# output: out = matrix whose rows are possible combinations
cartesian_prod <- function(mat){
  if(is.null(nrow(mat))) {
    out <- mat[!mat %in% 0]
  }
  else {
    mat_rows <- list()
    for (s in 1:nrow(mat)){
      mat_rows[[s]] <- mat[s, !mat[s,] %in% 0]
    }
    out <- expand.grid(mat_rows)
    colnames(out) <- NULL
  }
  return(as.matrix(out))
}

# function that computes pathwise probability of specified path
# input: path = nodes of specified path, vector type; mat_sto = stochastic matrix
# output: pathwise probability (the product of conditional probabilities) 
prob_path <- function(path, mat_sto){
  path <- c(1, path, 2 ^ num)
  ind <- matrix(0,nrow = (length(path) - 1), ncol = 2)
  for(i in 1:(length(path) - 1)){
    ind[i,] <- c(path[i],path[i+1])
  }
  return(prod(mat_sto[ind]))
}
