# R-code that generate a web of nodes and paths
#to illustrate possible path for the assembly of a given community
#in a probablistic approach
rm(list = ls())
source("toolbox.R")
library(deSolve)
library(igraph)

# ---- Initialize ----
num <- 4; stren <- 1; conne <- 1
A <- generate_Interaction_matrix(num, stren, conne)
#A <- as.matrix(read.table(paste("invivo.txt")))


# ---- Create the list of sub-communities ----
##>all subcmmunities of the same species number (defined as a group) is stored in one combn matrix
##>all combn matrix is stored in one list
sub_coms = list()
for (s in 1:num) {
  sub_coms[[s]] <- combn(num, s)
}

# ---- Compute transition possibilities via feasibility----
# Create the transition matrixes (T,D,H is forward, backward and bidirectional matrix respectively)
matT <- matrix(0, ncol = 2 ^ num, nrow = 2 ^ num)
matD <- matrix(0, ncol = 2 ^ num, nrow = 2 ^ num)
matH <- matrix(0, ncol = 2 ^ num, nrow = 2 ^ num)

##>For every node, compute all the possibilities (omega_overlap for different dim),
##>normalize again for each row?

# index transformation
t <- matrix(0, nrow = num, ncol = max(choose(num, 1:num)))
for (s in 1:num) {
  for (i in 1:choose(num, s)) {
    t[s,i] <- sum(choose(num,1:s-1)) + i
  }
}

##>just randomly set to 1/Z, for we have no a priori knowledge on them
matT[1, 2:(2 ^ num)] <- 1/(2 ^ num)
matD[2:(2 ^ num), 1] <- 1/(2 ^ num)

# Compute normalized omega for each node 
n_omega_node <- c(0.5, rep(0, 2^num-1))
for (s in 1:num){
  for (i in 1:choose(num, s)){
    ori <- sub_coms[[s]][, i]
    n_omega_node[t[s,i]] <- Omega(A[ori, ori])^(1/s)
  }
}

# Compute raw omega for each node 
omega_node <- c(0.5, rep(0,2^num-1))
for (s in 1:num){
  for (i in 1:choose(num, s)){
    ori <- sub_coms[[s]][, i]
    omega_node[t[s,i]] <- Omega(A[ori, ori])
  }
}

# Compute transitional probabilities
for (s in 1:(num - 1)) {
  targ1 <- sub_coms[[s]]
  
  for (i in 1:choose(num, s)) {
    ori <- sub_coms[[s]][, i]
    ori_mat <- A[ori, ori]
    
    for (k in (s + 1):num){
      targ2 <- sub_coms[[k]]
      targ3 <- extend_communities(ori, num, k-s)
      
      if (s >= 2){
        for (j in (1:ncol(targ1))[-i]) {
          targ_mat <- A[targ1[, j], targ1[, j]]
          matH[t[s, i], t[s, j]] <- Omega_overlap(ori_mat, targ_mat) / Omega(ori_mat)
        }
      }
      for (j in 1:ncol(targ2)) {
        #filter those that completely have ori as their subset
        if(is.vec_in_mat(targ2[, j], targ3)) {
          targ_mat <- A[targ2[, j], targ2[, j]]
          Omega_ij <- Omega_overlap_ext(ori_mat, targ_mat, ori, targ2[, j])
          matT[t[s, i], t[k, j]] <- Omega_ij / Omega(ori_mat)
          matD[t[k ,j], t[s, i]] <- Omega_ij / Omega(targ_mat)
        }
      }
    }
  }
}
# Generate the matrixes and their norm form
matH <- matH + matT + matD

matT_norm <- norm_row_sum(matT)
matD_norm <- norm_row_sum(matD)
matH_norm <- norm_row_sum(matH)
matT_norm <- (1-omega_node)*matT_norm
matD_norm <- (1-omega_node)*matD_norm
matH_norm <- (1-omega_node)*matH_norm
diag(matT_norm) <- omega_node
diag(matD_norm) <- omega_node
diag(matH_norm) <- omega_node


# ---- Visualize ----------
threshold <- 0
mat_disp <- matH_norm

# Specify row/col names for transition matrix
names <- c("0")
for (s in 1:num) {
  names <- c(names, c(apply(sub_coms[[s]], 2, convert2names)))
}
colnames(mat_disp) = rownames(mat_disp) = names

# Generate the adjacency matrix with the threshold
set_conne <- function(value, threshold){
  if(!is.na(value) && value>threshold){
    return(1)
  }else{
    return(0)
  }
}
mat_adj <- apply(mat_disp, c(1,2), set_conne, threshold)

# Get the value of transition possibility to display
value = c()
for (i in 1:nrow(mat_adj)){
  for (j in 1:ncol(mat_adj)){
    if(mat_adj[i,j] != 0){
      value = c(value, mat_disp[i,j])
    }
  }
}

# Set the grid to place nodes
grid <- matrix(0, nrow= 2^num, ncol = 2)
l <- 2
for (s in 1:num) {
  for (i in 1:choose(num,s)) {
    grid[l,1] <- s
    grid[l,2] <- (choose(num,s)-1)/2 -(i - 1)
    l <- l + 1
  }
}


# ---- Markov process and information analysis ----
mat_ent <- matH_norm
#>NA in Omega and thus in Trans needs to be fixed
mat_ent[is.na(mat_ent)] <- 0
#>use data from Hill2004 to check if the algorithm works well
# mat_ent <- read.csv("data_hill2004.csv", header = FALSE, sep=",") 
# mat_ent <- t(mat_ent)
# Stationary distribution
Stat_dist <- function(Trans){
  Trans <- t(Trans)
  vec_w <- eigen(Trans)$vectors[,1]
  vec_w <- vec_w/sum(vec_w)
  return(vec_w)
}
sd <- as.numeric(Stat_dist(mat_ent))


# Display the graph
diag(mat_adj) <- 0
# Display the graph
palf <- colorRampPalette(c("dodgerblue","darkorange"))
pal1 <- heat.colors(length(n_omega_node), alpha=1) 
pal2 <- palf(length(n_omega_node))
r <- rank(n_omega_node)
network <- graph_from_adjacency_matrix(mat_adj)
plot(network,
     layout = grid,
     edge.arrow.size = 0.4,
     edge.width = 3 * value / max(value),
     vertex.label.color = "black",
     vertex.frame.color = "black",
     vertex.color = pal2[r],
     #vertex.size = 30 * omega_node / max(omega_node)
     vertex.size = 30 * (sd) / max(sd)
)


# Entropy
Entropy <- function(Trans){
  Trans <- t(Trans)
  logsum <- function(x){
    sum <- 0
    for (i in 1:length(x)){
      if(x[i] > 0) sum <- sum + x[i]*log(x[i], base = exp(1))
    }
    return(sum)
  }
  vec_p <- apply(Trans, 2, logsum)
  vec_w <- eigen(Trans)$vectors[,1] #the dominated vector
  vec_w <- vec_w/(sum(vec_w))
  return(-sum(vec_p * vec_w))
}
Ent_a <- Entropy(mat_ent) #absolute entropy
Ent_r <- Entropy(mat_ent)/log(nrow(mat_ent)) #relative entropy

# One-species?
Pr <- c(rep(0, num))
for (i in 1:num){
  for (t in 1:2^num){
    #find the sub-community that contains spi and add the stationary distribution value
    if(is.element(i, convert2sets(names[t]))){
      Pr[i] <- Pr[i] + sd[t]
    }
  }
}
Pr

# ---- Least resiliance? Most probable path? ----