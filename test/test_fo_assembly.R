# nolint start
source("toolbox.R") # feasoverlap pkg/source included
source("wrapper.R") # change defaults to raw omega if using feaoverlap

num <- 4;
sub_coms = list()
for (s in 1:num) {
  sub_coms[[s]] <- combn(num, s)
}
t_ind <- matrix(0, nrow = num, ncol = max(choose(num, 1:num)))
for (s in 1:num) {
  for (i in 1:choose(num, s)) {
    t_ind[s,i] <- sum(choose(num,1:s-1)) + i
  }
}

g1 <- ggplot(); g2 <- ggplot()
stren <- 1; conne <- 0.5; order <- 1
i <- 0
  set.seed(331+i^2-10*i)
  A <- interaction_matrix_random(num, stren, conne)
   A <- interaction_matrix_ill(num, stren, conne, 10^(-5))
  
  # matL & matR computes ADO in two different orders
  # matO computes min of original Omgea values
  matL <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
  matR <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
  matO <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
  matS <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
  
  # Compute transitional probabilities
  for (s in 1:(num - 1)) {
    targ1 <- sub_coms[[s]]
    
    for (i in 1:choose(num, s)) {
      ori <- sub_coms[[s]][, i]
      ori_mat <- A[ori, ori]
      matO[t_ind[s, i],t_ind[s, i]] <- calculate_omega(ori_mat)
      
      if (s == 1) {
        matL[t_ind[s, i], t_ind[s, i]] <- NaN
        matR[t_ind[s, i], t_ind[s, i]] <- NaN
      }
      else {
        matL[t_ind[s, i], t_ind[s, i]] <- calculate_omega_overlap(ori_mat, ori_mat)
        matR[t_ind[s, i], t_ind[s, i]] <- calculate_omega_overlap(ori_mat, ori_mat)
      }
      
      for (k in (s + 1):num){
        targ2 <- sub_coms[[k]]
        targ3 <- extend_communities(ori, num, k-s)
        
        # currently works only for >= 2 species 
        if (s >= 2){
          for (j in (1:ncol(targ1))[-i]) {
            targ_mat <- A[targ1[, j], targ1[, j]]
            matL[t_ind[s, i], t_ind[s, j]] <- calculate_omega_overlap(ori_mat, targ_mat)
            matR[t_ind[s, i], t_ind[s, j]] <- calculate_omega_overlap(targ_mat, ori_mat)
            matO[t_ind[s, i], t_ind[s, j]] <- min(calculate_omega(ori_mat), calculate_omega(targ_mat))
            matS[t_ind[s, i], t_ind[s, j]] <- max(calculate_omega(ori_mat), calculate_omega(targ_mat))
          }
        }
        
        for (j in 1:ncol(targ2)) {
          #filter those that completely have ori as their subset
          if(is.vec_in_mat(targ2[, j], targ3)) {
            targ_mat <- A[targ2[, j], targ2[, j]]
            print(c(s,i,k,j))
            matL[t_ind[s, i], t_ind[k, j]] <- Omega_overlap_ext(ori_mat, targ_mat, ori, targ2[, j], 1)
            matR[t_ind[s, i], t_ind[k, j]] <- Omega_overlap_ext(ori_mat, targ_mat, ori, targ2[, j], 2)
            matO[t_ind[s, i], t_ind[k, j]] <- min(calculate_omega(ori_mat), calculate_omega(targ_mat))
            matS[t_ind[s, i], t_ind[s, j]] <- max(calculate_omega(ori_mat), calculate_omega(targ_mat))
          }
        }
      }
    }
}
  
  T1 <- (matL- matR)/((matL+matR)/2)
  T2 <- pmax(matL, matR)/matO
  
  T1 <- T1[!is.na(T1)]; T2 <- T2[!is.na(T2)]
  data1 <- data.frame(position = 1:length(T1), T1 = T1)
  data2 <- data.frame(position = 1:length(T2), T2 = T2)
  
  g1 <- g1 + 
    geom_point(data = data1,
               mapping = aes(x = position, y = T1),
               alpha = (abs(T1)/0.02+0.02),size = 2.5,
               position = "jitter")
  g2 <- g2 + 
    geom_point(data = data2,
               mapping = aes(x = position, y = T2),
               alpha = T2/max(2), size = 2.5,
               position = "jitter")
}

g31 <- g1+theme(axis.text=element_text(size=14),
         axis.title=element_text(size=16))
g32 <- g2+theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))

# nolint end