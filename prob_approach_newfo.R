# nolint start
# This code generate the assembly map of an ecological community
rm(list = ls())
library(deSolve)
library(igraph)
source("toolbox.R") # feaoverlap pkg/source included
source("wrapper.R") # change defaults to raw omega if using feaoverlap

# ------ Initialize ------
num <- 4; stren <- 1; conne <- 1; order <- 1
set.seed(354)
A <- interaction_matrix_random(num, stren, conne)

# A <- as.matrix(read.table("data/Gould_Matrix_median_12.csv",sep=","))[1:5,1:5]
# num <- ncol(A)

# ------ Create the list/indices of subcommunities ------
# Original list: location to composition
sub_coms = list()
for (s in 1:num) {
  sub_coms[[s]] <- combn(num, s)
}

# t_index: location to order
t_ind <- matrix(NA, nrow = num, ncol = max(choose(num, 1:num)))
for (s in 1:num) {
  for (i in 1:choose(num, s)) {
    t_ind[s,i] <- sum(choose(num,1:s-1)) + i
  }
}

# l_index: order to location
l_ind <- matrix(0, 2^num, 2)
for (t in 2:nrow(l_ind)){
  l_ind[t,] <- which(t_ind == t, arr.ind = TRUE)
}

# ----- Compute raw/norm/overlap omega values ------
# Compute raw/norm omega for each node
r_omega_node <- evaluate_nodes(A, TRUE)
n_omega_node <- evaluate_nodes(A, FALSE)

# Compute omega_overlaps between possible nodes
Overlap <- evaluate_overlaps(A, TRUE)
overlap <- evaluate_overlaps(A, FALSE)


# ------  Pathwise probability analysis ------
# find all paths from ti to tf
ti <- 1; tf <- 16
paths <- find_path(ti, tf)

entire_pr <- tibble(
  n_step = c(0), n_path = c(0),
  random_pr = c(0), environ_pr = c(0), species_pr = c(0)
)

#Notice: using normailized values now
for (k in 1:length(paths)){
  paths_k <- cbind(ti, paths[[k]], tf)
  for (p in 1:nrow(paths_k)){
    rp <- prob_path(paths_k[p,], n_omega_node, overlap, "r")
    ep <- prob_path(paths_k[p,], n_omega_node, overlap, "e")
    sp <- prob_path(paths_k[p,], n_omega_node, overlap, "s")
    # Do some filtering
    entire_pr <- entire_pr %>% 
      add_row(n_step=k, n_path=p, random_pr=rp, environ_pr=ep, species_pr=sp)
  }
}
view(entire_pr)


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
mat_disp <- Overlap

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
    # vertex.size = 30 * r_omega_node / max(r_omega_node)
    # vertex.size = 20 * (n_omega_node) / max(n_omega_node)
)



# nolint end