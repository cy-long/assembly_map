# nolint start
# Test if feasoveralp library works well in the assembly process
rm(list = ls())

source("toolbox.R")
#library(feaoverlap)
source("overlap.R")

# ---- Initialize ----
num <- 5; stren <- 1; conne <- 0.5

A <- interaction_matrix_random(num, stren, conne)

sub_coms = list()
for (s in 1:num) {
  sub_coms[[s]] <- combn(num, s)
}

# index transformation
t <- matrix(0, nrow = num, ncol = max(choose(num, 1:num)))
for (s in 1:num) {
  for (i in 1:choose(num, s)) {
    t[s,i] <- sum(choose(num,1:s-1)) + i
  }
}

# matL & matR computes ADO in two different orders
# matO computes min of original Omgea values
matL <- matrix(0, ncol = 2 ^ num, nrow = 2 ^ num)
matR <- matrix(0, ncol = 2 ^ num, nrow = 2 ^ num)
matO <- matrix(0, ncol = 2 ^ num, nrow = 2 ^ num)

# Compute transitional probabilities
for (s in 1:(num - 1)) {
  targ1 <- sub_coms[[s]]
  
  for (i in 1:choose(num, s)) {
    ori <- sub_coms[[s]][, i]
    ori_mat <- A[ori, ori]
    
    for (k in (s + 1):num){
      targ2 <- sub_coms[[k]]
      targ3 <- extend_communities(ori, num, k-s)
      
      # currently works only for >= 2 species 
      if (s >= 2){
        for (j in (1:ncol(targ1))[-i]) {
          targ_mat <- A[targ1[, j], targ1[, j]]
          matL[t[s, i], t[s, j]] <- calculate_omega_overlap(ori_mat, targ_mat)
          matR[t[s, i], t[s, j]] <- calculate_omega_overlap(targ_mat, ori_mat)
          matO[t[s, i], t[s, j]] <- min(calculate_omega(ori_mat), calculate_omega(targ_mat))
        }
      }
      
      for (j in 1:ncol(targ2)) {
        #filter those that completely have ori as their subset
        if(is.vec_in_mat(targ2[, j], targ3)) {
          targ_mat <- A[targ2[, j], targ2[, j]]
          print(c(s,i,k,j))
          matL[t[s, i], t[k, j]] <- Omega_overlap_ext(ori_mat, targ_mat, ori, targ2[, j], 1)
          matR[t[s, i], t[k, j]] <- Omega_overlap_ext(ori_mat, targ_mat, ori, targ2[, j], 2)
          matO[t[s, i], t[k, j]] <- min(calculate_omega(ori_mat), calculate_omega(targ_mat))
        }
      }
    }
  }
}

# Visualize in matrix_heatmap.R

# nolint end
