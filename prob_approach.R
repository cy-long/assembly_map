# R-code that generate a web of nodes and paths
#to illustrate possible path for the assembly of a given community
#in a probablistic approach

rm(list = ls())
source("toolbox.R")
library(deSolve)
library(igraph)

# ---- Initialize ----
num <- 4
stren <- 1
conne <- 0.6
A <- generate_Interaction_matrix(num, stren, conne)

# ---- Create the list of sub-communities ----
##>all subcmmunities of the same species number (defined as a group) is stored in one combn matrix
##>all combn matrix is stored in one list
sub_coms = list()
for (s in 1:num) {
  sub_coms[[s]] <- combn(num, s)
}

# ---- Compute transition possibilities via feasibility----
# Create the transition matrix
T <- matrix(0, ncol = 2 ^ num, nrow = 2 ^ num)

##>For every node, compute all the possibilities (omega_overlap for different dim),
##>normalize again for each row?

# index transformation
t <- matrix(0, nrow = num, ncol = max(choose(num, 1:num)))
for (s in 1:num) {
  for (i in 1:choose(num, s)) {
    t[s,i] <- sum(choose(num,1:s-1)) + i
  }
}


# discuss s=0, i.e. the initial state of nothing
##>just randomly set to 1/num, for we have no a priori knowledge on them?
T[1, 2:(num+1)] <- 1/num

# discuss s>0
for (s in 1:(num-1)) {
  for (i in 1:choose(num, s)) {
    ori <- sub_coms[[s]][, i]
    ori_mat <- A[ori, ori]
    
    targ1 <- sub_coms[[s]]
    targ2 <- sub_coms[[s+1]]  
    targ3 <- extend_communities(ori, num)
    
    if (s >= 2){
      for (j in (1:ncol(targ1))[-i]) {
        targ_mat <- A[targ1[, j], targ1[, j]]
        T[t[s, i], t[s, j]] <- Omega_overlap(ori_mat, targ_mat) / Omega(ori_mat)
      }
      for (j in 1:ncol(targ2)) {
        # filter those that completely have ori as their subset
        if(is.vec_in_mat(targ2[, j], targ3)) {
          targ_mat <- A[targ2[, j], targ2[, j]]
          T[t[s, i] , t[s + 1, j]] <- Omega_overlap_ext(ori_mat, targ_mat, ori, targ2[ , j]) / Omega(ori_mat)
        }
      }
    } else if (s == 1) {
      for (j in 1:ncol(targ2)) {
        if(is.vec_in_mat(targ2[, j], targ3)) {
          targ_mat <- A[targ2[, j], targ2[, j]]
          T[t[s, i] , t[s + 1, j]] <- Omega_overlap_ext(ori_mat, targ_mat, ori, targ2[ , j]) / (1/2)
        }
      }
    }
    
  }
}

# ---- Visualize ----
