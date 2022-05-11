# nolint start
library(tidyverse)
library(mvtnorm)
library(mgcv)
library(binaryLogic)

library(feasoverlap)
# library(dplyr)
# library(geometry)
# library(uniformly)
# library(pracma)
# source("overlap.R")


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
Omega_overlap_ext <- function(A, B, comp_A, comp_B, raw = TRUE, order = 1) {
  if (all(is.element(comp_A, comp_B))) {
    # create a location vector to store indexes where insertion happens
    loc <- c(rep(0, length((comp_B))))
    for (i in 1:length(comp_B)) {
      if(!is.element(comp_B[i], comp_A)) {
        loc[i] <- 1
      }
    }
    
    # create a matrix of all possible diagonal combinations of +-
    generate_diag_string <- function(N){
      record <- matrix(0, 2^N, N)
      for (i in 1:(2^N-1)){
        record[i+1, (N+1-length(as.binary(i))):N] <- as.binary(i)
      }
      record = apply(record,c(1,2), function(x) {if(x <= 0){return(-1)} else return(1)})
      return(record)
    }
    diag_mat <- generate_diag_string(sum(loc))
    
    # compute Omega_overlap
    Omega_overlap_ts <- c(rep(0, (2 ^ sum(loc))))
    if (length(loc) == 1) return(1/2)
    for (i in 1:(2 ^ sum(loc))){
      A_ext <- matrix_scatter(A, loc, diag_mat[i, ])
      if(order == 1) {
        Omega_overlap_ts[i] = calculate_omega_overlap(A_ext, B, raw = TRUE)
      } else {
        Omega_overlap_ts[i] = calculate_omega_overlap(B, A_ext, raw = TRUE)
      }
    }
    
    if (raw) {
      return(sum(Omega_overlap_ts))
    } else {
      return(sum(Omega_overlap_ts)^(1/length(comp_B)))
    } 
  } else if (all(is.element(comp_B, comp_A))) {
    Omega_overlap_ext(B, A, comp_B, comp_A) #change the variations and do recursion
  } else {
    warning("overlapping two sub-communities that are not subset of each other")
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

  #sum all movement steps needed to insert the new species
  mov <- c(rep(0, n_mat))
  j <- 0
  for (i in 1:length(loc)) {
    if(loc[i] == 0) {
      j <- j + 1
      mov[j] <- sum(loc[1:i])
    }
  }
  if(j != n_mat) {
    warning("matrix_scatter: discrepancy between ncol(mat) and length(loc)")
  }
  
  # scatter the original elements
  mat_ext <- matrix(0, ncol = length(loc), nrow = length(loc))
  if(n_mat > 1){
    for (i in 1:n_mat) {
      for (j in 1:n_mat) {
        mat_ext[i+mov[i], j+mov[j]] <- mat[i,j]
      }
    }
  } else if (n_mat == 1) {
    mat_ext[1+mov[1], 1+mov[1]] <- mat
  }

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
# augm = adding number of species
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
  t(apply(mat,1,function(x) x/sum(x, na.rm = TRUE)))
}

# function that returns the cartesian product from sets (filter zero elements)
# input: mat = matrix whose rows are original sets, using zeros as placeholders
# output: out = matrix whose rows are possible combination()
cartesian_prod <- function(mat){
  if(is.null(nrow(mat))) {
    out <- mat[!mat %in% NA]
  }
  else {
    mat_rows <- list()
    for (s in 1:nrow(mat)){
      mat_rows[[s]] <- mat[s, !mat[s,] %in% NA]
    }
    out <- expand.grid(mat_rows)
    colnames(out) <- NULL
  }
  return(as.matrix(out))
}

# function
# input
# output
path_filter <- function(paths_raw,ti,tf){
  comp_of_t <- function(t){
    if(t == 1){
      return(NULL)
    } else {
      return(sub_coms[[l_ind[t,1]]][,l_ind[t,2]])
    }
  }

  paths_fil <- matrix(NA, nrow = 1, ncol = ncol(paths_raw))
  for (pa in 1:nrow(paths_raw)){
    testpath <- c(ti, paths_raw[pa,], tf)
    pathresult <- c()
    for (no in 1:(length(testpath)-1)){
      s1 <- comp_of_t(testpath[no])
      s2 <- comp_of_t(testpath[no+1])
      pathresult <- append(pathresult, all(s1 %in% s2))
    }
    if(all(pathresult)){
      paths_fil <- rbind(paths_fil, paths_raw[pa,])
    }
  }
    return(as.matrix(paths_fil[-1,]))
}

# function that generate potential paths of assembly
# input: ti = t_ind of starting node; tf = t_ind of ending node
# output: list of potential paths represented by t_inds
find_path <- function(ti, tf){
  s <- l_ind[ti,1]; i <- l_ind[ti,2]
  p <- l_ind[tf,1]; j <- l_ind[tf,2]
  
  if (s == 0){
    set_i <- c()
  } else {
    set_i <- sub_coms[[s]][ ,i]
  }
  set_f <- sub_coms[[p]][ ,j]

  if(!(all(set_i %in% set_f))){
    stop("Current version only deals with sub-communities")
  }

  paths <- list()
  distance <- abs(p-s-1)
  if(distance == 0){
    stop("There is no path within one layer")
  }
  for (k in 1:distance){
    if(k == 1){
      paths_k <- as.matrix((ti+1):(tf-1))
      paths[[k]] <- path_filter(paths_k,ti,tf)
    } else {
      paths_k <- matrix(nrow = 1, ncol = k)
      for (r in 1:choose(distance,k)){
        t_rows <- combn((s+1):(p-1),k)[,r] #select layers
        paths_k <- rbind(paths_k, cartesian_prod(t_ind[t_rows, ]))
      }
      paths_k <- paths_k[-1,]
      paths[[k]] <- path_filter(paths_k, ti, tf)
    }
  }
  return(paths)
}

# function that computes pathwise probability of specified path
# input: path = t_ind vec of specified path; nodes = raw omega value of nodes;
# input: overlaps = raw omega overlaps between two nodes; method = "r"(random), "e"(environ), "s"(species)
# output: pathwise probability
prob_path <- function(path, nodes, overlaps, method){
  ti <- path[1]; tf <- path[length(path)]
  ts <- ncol(Overlap) # May be specified in params
  chain <- matrix(0, nrow = (length(path) - 1), ncol = 2)

  if(method == "r"){
    prob <- prod(nodes[path])
  }
  else if (method == "e"){
    prob <- prod(overlaps[path,ts])/((nodes[ts])^(length(path)))
  }
  else if (method == "s"){
    for(i in 1:(length(path) - 1)){
      chain[i,] <- c(path[i],path[i+1])
    }
    pr_chains <- prod(overlaps[chain])/(prod(nodes[chain[,1]]))
    pr_init <- overlaps[ti,ts]/nodes[ts]
    prob <- pr_chains * pr_init
  }
  return(prob)
}

# function that
# input
# output
interaction_matrix_ill <- function(num, stren, conne, epsilon, threshold = 0) {
  inte <- interaction_matrix_random(num, stren, conne)
  new_col <- floor(num/2) + 1
  #inte[new_col,1] <- (max(inte[new_col,1], threshold))
  factor <- (-1)/inte[new_col,1]
  inte[,new_col] <- factor*inte[,1] + epsilon * rnorm(num)
  inte[new_col,new_col] <- -1
  return(inte)
}

# ----- Compute raw/norm/overlap omega values ------
# Compute raw/norm omega for each node
evaluate_nodes <- function(A, raw){
  A <- as.matrix(A); num <- ncol(A)
  omega_single<- c(0.5, rep(0, 2^num-1))

  for (s in 1:num){
    for (i in 1:choose(num, s)){
      ori <- sub_coms[[s]][, i]
      omega_single[t_ind[s,i]] <- calculate_omega(A[ori, ori], raw)
    }
  }
  return(omega_single)
}


evaluate_overlaps <- function(A, raw){
  # Compute matrices and omega_overlaps
  Overlap <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
  for (s in 0:(num - 1)) {
    for (i in 1:choose(num, s)) {
      if (s == 0){
        ori <- NULL
      } else {
        ori_layer <- sub_coms[[s]]
        ori <- ori_layer[, i]
      }
      ori_mat <- A[ori, ori]

      for (p in (s + 1):num){
        tar_layer <- sub_coms[[p]]
        for (j in 1:ncol(tar_layer)) {
          tar <- tar_layer[, j]
          tar_mat <- A[tar, tar]

          # filter those that completely have ori as their subset
          if(is.vec_in_mat(tar, extend_communities(ori, num, p-s))) {
            Overlap_value <- Omega_overlap_ext(ori_mat, tar_mat, ori, tar, raw)
            ti <- t_ind[s, i]
            tf <- t_ind[p, j]
            if(s == 0){
              Overlap[1, tf] <- Overlap_value
            } else {
              Overlap[ti, tf] <- Overlap_value
            }
          }
        }
      }
    }
  }
  diag(Overlap) <- evaluate_nodes(A,raw)
  return(Overlap)
}


#nolint end