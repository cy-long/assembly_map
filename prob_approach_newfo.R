# nolint start
# Test if feasoveralp library works well in the assembly codes

rm(list = ls())
library(deSolve)
library(igraph)
source("toolbox.R") # feaoverlap pkg/source included
source("wrapper.R") # change defaults to raw omega if using feaoverlap

# ------ Initialize ------
num <- 4; stren <- 1; conne <- 0.9; order <- 1
set.seed(319)
A <- interaction_matrix_random(num, stren, conne)

# A <- as.matrix(read.table("data/Friedman_Matrix.csv",sep=","))
# num <- ncol(A)

# ------ Create the list of sub-communities ------
##>all subcmmunities of the same species number (defined as a group) is stored in one combn matrix
##>all combn matrix is stored in one list
sub_coms = list()
for (s in 1:num) {
  sub_coms[[s]] <- combn(num, s)
}
# index transformation
t_ind <- matrix(NA, nrow = num, ncol = max(choose(num, 1:num)))
for (s in 1:num) {
  for (i in 1:choose(num, s)) {
    t_ind[s,i] <- sum(choose(num,1:s-1)) + i
  }
}
l_ind <- matrix(0, 2^num, 2)
for (t in 2:nrow(l_ind)){
  l_ind[t,] <- which(t_ind == t, arr.ind = TRUE)
}


# ----- Compute transition possibilities via feasibility -----
# Create the transition matrixes (T,D,H is forward, backward and bidirectional matrix respectively)
matT <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
matD <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
matH <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)

##>For every node, compute all the possibilities (omega_overlap for different dim),
##>normalize again for each row?

# Compute normalized & raw omega for each node
n_omega_node <- c(0.5, rep(0, 2^num-1))
r_omega_node <- c(0.5, rep(0, 2^num-1))

for (s in 1:num){
  for (i in 1:choose(num, s)){
    ori <- sub_coms[[s]][, i]
    n_omega_node[t_ind[s,i]] <- calculate_omega(A[ori, ori], FALSE)
    r_omega_node[t_ind[s,i]] <- calculate_omega(A[ori, ori], TRUE)
  }
}

# Compute transitional probabilities
Overlap <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
overlap <- matrix(NA, ncol = 2 ^ num, nrow = 2 ^ num)
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
      targ_layer <- sub_coms[[p]]

      for (j in 1:ncol(targ_layer)) {
        targ <- targ_layer[, j]
        targ_mat <- A[targ, targ]

        #filter those that completely have ori as their subset
        if(is.vec_in_mat(targ, extend_communities(ori, num, p-s))) {
          Overlap_value <- Omega_overlap_ext(ori_mat, targ_mat, ori, targ)
          ti <- t_ind[s, i]
          tf <- t_ind[p, j]
          if(s == 0){
            matT[1, tf] <- Overlap_value
            matD[tf, 1] <- Overlap_value
            matH[1, tf] = matH[tf, 1] <- Overlap_value
            Overlap[1, tf] <- Overlap_value
            overlap[1, tf] <- Overlap_value^(1/length(targ))
          } else {
            matT[ti, tf] <- Overlap_value
            matD[tf, ti] <- Overlap_value
            matH[ti, tf] = matH[tf, ti] <- Overlap_value
            Overlap[ti, tf] <- Overlap_value
            overlap[ti, tf] <- Overlap_value^(1/length(targ))
          }
          #print(c(s,i,p,j))
        }
      }
    }
  }
}

prob_path(c(1,2,6,12),r_omega_node,Overlap,"r")
prob_path(c(1,2,6,12),r_omega_node,Overlap,"e")
prob_path(c(1,2,6,12),r_omega_node,Overlap,"s")




# ----- Operations on the matrices -----
# Generate the matrixes and their norm form
markov_norm <- function(mat, weights){
  mat_markov <- norm_row_sum(mat)
  mat <- (1-weights) * mat
  diag(mat) <- weights
  return(mat_markov)
}

matT_norm <- markov_norm(matT,r_omega_node)
matD_norm <- markov_norm(matD,r_omega_node)
matH_norm <- markov_norm(matH,r_omega_node)


# ------ Markov process and information analysis ------
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
  for (j in 1:2^num){
    #find the sub-community that contains spi and add the stationary distribution value
    if(is.element(i, convert2sets(names[j]))){
      Pr[i] <- Pr[i] + sd[j]
    }
  }
}
#Pr

# ----- Visualize -----
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
     #vertex.size = 30 * r_omega_node / max(r_omega_node)
     vertex.size = 30 * (sd) / max(sd)
)

# ------  pathwise probabilities ------
# find all paths from ti to tf
ti <- 1; tf <- 16
paths <- find_path(ti, tf)




# compute pathwise probability from Markov chain
entire_pr <- tibble(
  n_step = c(1),
  n_path = c(1),
  random_pr = c(1),
  environ_pr = c(1),
  species_pr = c(1)
)
entire_pr %>% add_row(n_step=1,n_path=2,random_pr=3,environ_pr=3,species_pr=4)


for (k in 1:(num-1)){
  if(k == 1){
    for (p in 1:length(paths[[k]])){
      prob <- prob_path(paths[[k]][p], mat_sto)
      entire_pr <- rbind(entire_pr, c(k,p,prob))
    }
  }
  else{
    for (p in 1:nrow(paths[[k]])){
      prob <- prob_path(paths[[k]][p,], mat_sto)
      entire_pr <- rbind(entire_pr, c(k,p,prob))
    }
  }
}

# filtering and visualization
positive_probs <- entire_pr[(entire_pr$v_prob >= 1e-12), ]
rownames(positive_probs) <- NULL
positive_probs

steps_log_probs <- list(); steps_nlog_probs <- list()

steps_log_probs[[as.character(0)]] <- log(positive_probs$v_prob[1])
for (k in 1:(num-1)){
  steps_log_probs[[as.character(k)]] <-
    log(positive_probs$v_prob[positive_probs$n_step %in% k])
  steps_nlog_probs[[as.character(k)]] <- (1/k)*steps_log_probs[[as.character(k)]]
}

# boxplot(steps_log_probs,
#         xlab = "middle_steps",
#         ylab = "log(probs)",
#         main = "Raw prob.dist. of assembly steps")

# boxplot(steps_nlog_probs,
#         xlab = "middle_steps",
#         ylab = "log(probs)",
#         main = "Normalized prob.dist. of assembly steps")

# nolint end