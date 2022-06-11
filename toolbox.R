# nolint start
library(tidyverse)
library(mvtnorm)
library(mgcv)
library(binaryLogic)
library(igraph)

# ----- feasoverlap pkg or codefile -----
# library(feasoverlap)

library(dplyr)
library(geometry)
library(uniformly)
library(pracma)
source("overlap.R")

# ----- indexing functions ------

# functions that generate indices for the subcommunities
generate_subcoms <- function(num) {
  # Original list: location to composition
  sub_coms = list()
  for (s in 1:num) {
    sub_coms[[s]] <- combn(num, s)
  }
  sub_coms
}
generate_t_index <- function(num) {
  # t_index: location to order
  t_ind <- matrix(NA, nrow = num, ncol = max(choose(num, 1:num)))
  for (s in 1:num) {
    for (i in 1:choose(num, s)) {
      t_ind[s,i] <- sum(choose(num,1:s-1)) + i
    }
  }
  t_ind
}
generate_l_index <- function(num) {
  # l_index: order to location
  l_ind <- matrix(0, 2^num, 2)
  for (t in 2:nrow(l_ind)){
    l_ind[t,] <- which(t_ind == t, arr.ind = TRUE)
  }
  l_ind
}

# ----- computational functions -----

# function that computes raw/norm omega for single communities
evaluate_solo <- function(A, raw){
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

# function that computes overlap/extended omega for pairs of communities
evaluate_dual <- function(A, raw, state){
  Overlap <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
  Exten <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
  
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
            Overlap_value <- Omega_overlap_ext(ori_mat, tar_mat, ori, tar, raw, 2)
            Exten_value <- Omega_overlap_ext(ori_mat, tar_mat, ori, tar, raw, 1)
            ti <- t_ind[s, i]
            tf <- t_ind[p, j]
            if(s == 0){
              Overlap[1, tf] <- Overlap_value
              Exten[1, tf] <- Exten_value
            } else {
              Overlap[ti, tf] <- Overlap_value
              Exten[ti, tf] <- Exten_value
            }
          }
        }
      }
    }
  }
  # diag(Overlap) <- evaluate_solo(A,raw)
  if(state == 2) {return(Overlap)}
  else if (state == 1) {return(Exten)}
  else {stop("no output")}
}


# function that computes the overlap of two feasibility domains with different dimensions
# inputs: A = one interaction matrix, B = another interaction matrix, 
# comp_A = composition of community A, comp_B = composition of community B
# output: volume_overlap = the normalized feasibility of the intersection region
Omega_overlap_ext <- function(A, B, comp_A, comp_B, raw = TRUE, state = 2) {
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
    
    # compute Overlap values or extented omega values
    Omega_comb <- c(rep(0, (2 ^ sum(loc))))
    if (length(loc) == 1) return(1/2)
    for (i in 1:(2 ^ sum(loc))){
      A_ext <- matrix_scatter(A, loc, diag_mat[i, ])
      if(state == 2) {
        Omega_comb[i] = calculate_omega_overlap(A_ext, B, raw = TRUE)
      } else {
        Omega_comb[i] = calculate_omega(A_ext, raw = TRUE)
      }
    }
    if (raw) {
      return(sum(Omega_comb))
    } else {
      return(sum(Omega_comb)^(1/length(comp_B))) #? true for state == 1?
    }

  } else if (all(is.element(comp_B, comp_A))) {
    Omega_overlap_ext(B, A, comp_B, comp_A) #change the variations and do recursion
  } else {
    warning("overlapping two sub-communities that are not subset of each other")
  }
}

# function that rearranges matrix elements based on a control vector
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
  #? for other modification methods, rewrite this part
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
# input: comm = one vector that represents community components, num = total number of species, 
# input: augm = adding number of species
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

# ----- pathwise functions -----

# function that computes pathwise probabilities
# input: paths = list of paths generate between initial and final points
# output: a tibble with characters and probabilities for paths
quantify_paths <- function(paths, r_omega, s_omega, Overlap){
  entire_pr <- tibble(
    n_step = c(0), n_path = c(0),
    random_pr = c(0), environ_pr = c(0), species_pr = c(0)
  )

  for (k in 1:length(paths)){
    paths_k <- cbind(ti, paths[[k]], tf)
    for (p in 1:nrow(paths_k)){
      rp <- prob_path(paths_k[p,], r_omega, s_omega, Overlap, "r")
      ep <- prob_path(paths_k[p,], r_omega, s_omega, Overlap, "e")
      sp <- prob_path(paths_k[p,], r_omega, s_omega, Overlap, "s")
      #? Do some filtering
      entire_pr <- entire_pr %>% 
        add_row(n_step=k, n_path=p, random_pr=rp, environ_pr=ep, species_pr=sp)
    }
  }
  return(entire_pr)
}


# function that computes pathwise probability of specified path
# input: path = t_ind vec of specified path; nodes = raw omega value of nodes;
# input: overlaps = raw omega overlaps between two nodes; 
# input: method = "r"(random), "e"(environ), "s"(species)
# output: pathwise probability
prob_path <- function(path, nodes, extens, overlaps, method){
  ti <- path[1]; tf <- path[length(path)]
  ts <- ncol(Overlap) #? May be specified directly with params
  chain <- matrix(0, nrow = (length(path) - 1), ncol = 2)

  if(method == "r"){
    prob <- prod(2*nodes[path]) 
  } #? Do we need to rescale it within [0,1]

  else if (method == "e") {
    prob <- prod(pmin(overlaps[path,ts],nodes[ts]))/((nodes[ts])^(length(path)))
  } #? use pmin to avoid ratio over 1

  else if (method == "s") {
    for(i in 1:(length(path) - 1)){
      chain[i,] <- c(path[i],path[i+1])
    }
    pr_chains <- prod(overlaps[chain])/prod(extens[chain])
    pr_init <- overlaps[ti,ts]/nodes[ts]
    prob <-  pr_init * pr_chains
  }

  return(prob)
}


# function that generates potential paths of assembly
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


# function that filters only paths forming a chain of subsets
# input: paths_raw = matrix of candidate paths
# input: ti = initial t_ind; tf = final t_ind
# output: matrix of paths that are chains of subsets
path_filter <- function(paths_raw, ti, tf){
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

# ----- matrix operation functions -----

# function that normalize matrix accroding to markov process
# input: mat = original mat, weights = normalizing factor
# output: row normalized matrix
Markov_norm <- function(mat, weights){
  mat[is.na(mat)] <- 0
  norm_row_sum <- function(mat){
    t(apply(mat,1,function(x) x/sum(x, na.rm = TRUE)))
  }
  mat_markov <- norm_row_sum(mat)
  mat <- (1-weights) * mat
  diag(mat) <- weights
  return(mat_markov)
}

# function that computes stationary distribution of Markov process
# input: mat_trans = transitional matrix of Markov process
# output: probability distribution 
Stat_dist <- function(mat_trans){
  mat_trans <- t(mat_trans)
  vec_w <- eigen(mat_trans)$vectors[,1]
  vec_w <- vec_w/sum(vec_w)
  return(vec_w)
}

# function that computes entropy of a Markov transitional matrix
# input: mat_trans = transitional matrix of Markov process
# output: entropy value accroding to (Hill et al. 2004)
Entropy <- function(mat_trans){
  mat_trans <- t(mat_trans)
  logsum <- function(x){
    sum <- 0
    for (i in 1:length(x)){
      if(x[i] > 0) sum <- sum + x[i]*log(x[i], base = exp(1))
    }
    return(sum)
  }
  vec_p <- apply(mat_trans, 2, logsum)
  vec_w <- eigen(mat_trans)$vectors[,1] #the dominated vector
  vec_w <- vec_w/(sum(vec_w))
  return(-sum(vec_p * vec_w))
}

# function that computes species level probabilities
# input: num = total number of species, pd = probability distri of states;
# input: names = names of subcommunities
# output: species level probabilities
Species_pr <- function(num, pd, names){
  Pr <- c(rep(0, num))
  for (i in 1:num){
    for (j in 1:2^num){
      if(is.element(i, convert2sets(names[j]))){
        Pr[i] <- Pr[i] + pd[j]
      }
    }
  }
  return(Pr)
}

# ----- tool sets -----

# function that checks if a vector is from one column of one matrix
# input: vec = the vector, mat = the matrix
# output: TRUE or FALSE
is.vec_in_mat <- function(vec,mat){
  out <- apply(mat, 2, function(x, y){isTRUE(all.equal(x, y))}, vec)
  any(out)
}


# function that converts numeric vector to strings or vice versa
# input: vec = numeric vector, such as c(1,2)
# output: string, such as "12"
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


# nolint end