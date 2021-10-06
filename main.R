# R-code that generate a web of nodes and paths
# to illustrate possible path for the assembly of a given community
rm(list = ls())
source("toolbox.R")   #load the toolbox
library(deSolve)
library(igraph)

# ---- Initialize ----
num <- 3
stren <- 1
conne <- 1

A <- generate_Interaction_matrix(num, stren, conne) #Temporarily fix num = 3
r <- parameterization_feasible(A, num)


# ---- create matrixes that store f/s status for sub-communities----
mat_fea <- matrix(0, 3, 3)
mat_sta <- matrix(0, 3, 3)

# ---- create adjacency matrix for the network----
data <- matrix(rep(0, 2^ (2 * num)), nrow = 2^num)
colnames(data) <- c("0", "1", "2", "3", "12", "13", "23", "123")
rownames(data) <- c("0", "1", "2", "3", "12", "13", "23", "123")
data["0", "1"]  <- 1; data["0", "2"] <- 1; data["0", "3"] <- 1

# ---- check f/s states per node ----
## check f/s in 1-sp communities (no need to check)
mat_fea[ , 1] <- 1
mat_sta[ , 1] <- 1

## generate sub-communities and their parameters
## >this part should be conducted automatically in a mature code
sp1_list <- list(1, 2, 3)
sp1_extend_list <- list(c(2, 3), c(1, 3), c(1, 2))

int1_list <- list(A[1, 1], A[2, 2], A[3, 3])
gr1_list <- list(r[1], r[2], r[3])


## check f/s in 2-sp communities
### generate sub-communities and their parameters
sp2_list <- list(c(1, 2), c(1, 3), c(2, 3))
sp2_extend_list <- list(3, 2, 1)

int2_list <- list(A[c(1, 2), c(1, 2)], A[c(1, 3), c(1, 3)], A[c(2, 3), c(2, 3)])
gr2_list <- list(r[c(1, 2)], r[c(1, 3)], r[c(2, 3)])

for (i in 1:3) {
  mat_fea[i, 2] <- check_feasibility(int2_list[[i]], gr2_list[[i]])
  mat_sta[i, 2] <- check_stability(int2_list[[i]], gr2_list[[i]])
}

## check f/s in 3-sp communities
mat_fea[1, 3] <- check_feasibility(A, r)
mat_sta[1, 3] <- check_stability(A, r)

mat_node <- mat_fea * mat_sta

# ---- connect nodes with possible transition paths ----

## 2nd layer to 3rd layer
if (sum(mat_node[ , 2]) == 0) {
  print("No possible transitions")
} else if (mat_node[1, 3] == 1) {
  list_out <-  list()
  for (i in 1:3) {
    if (mat_node[i, 2] == 1) {
      Ni <- solve(-int2_list[[i]], gr2_list[[i]])
      N <- c(0, 0, 0)
      N[sp2_list[[i]]] <- Ni
      N[sp2_extend_list[[i]]] <- 0.01 * mean(Ni)
      ##set the new species abundance at 1% of the average (discount)
      ##>the discount factor may have a significant impact on the outcome
      list_out[[i]] <- LV(A, N, r)
      from_names <- list("12", "23", "13")
      a <- from_names[[i]]
      b <- convert2names(list_out[[i]])
      data[a, b] <- 1
    } else {
      list_out[[i]] <- "AN" #Already non-feasible or non-stable
    }
  }
} else{
  print("The final state is not feasible/stable")
}
list_out

## 1st layer to 2nd layer
list_out <- list()
for(i in 1:3) {
  list_out[[i]] <- list()
  if (mat_node[i, 1] == 1) { #always true
    Ni <- solve(-int1_list[[i]], gr1_list[[i]])
    N <- c(0, 0)
    N[1] <- Ni
    N[2] <- 0.01 * Ni
    for (j in 1:2) {
      tmp_list <- c(i, sp1_extend_list[[i]][[j]])
      A_sub <- A[tmp_list, tmp_list]
      r_sub <- r[tmp_list]
      list_out[[i]][[j]] <- tmp_list[LV(A_sub, N, r_sub)]
      from_names <- list("1", "2", "3")
      a <- from_names[[i]]
      b <- convert2names(tmp_list[LV(A_sub, N, r_sub)])
      data[a, b] <- 1
    }
  } else {
  }
}
list_out

# ---- compute omega for each f/s node ----
omega_node <- matrix(0, 3, 3)
for (i in 1:3) {
  omega_node[i, 1] <- 1/2
  omega_node[i, 2] <- (Omega(int2_list[[i]])) ^ (1/2)
}
omega_node[1, 3] <- (Omega(A)) ^ (1/3)
omega_node

omega <- as.vector(c(0, omega_node[1:7]))

# ---- depict all the potential paths on a network graph ----
## build the graph object
assemble_map <- graph_from_adjacency_matrix(data)

## plot it
grid <- matrix(
  c(0, 0, 1, 1, 1, 0, 1, -1, 2, 1, 2, 0, 2, -1, 3, 0),
  ncol = 2, byrow = TRUE)
plot(assemble_map,
     layout = grid,
     edge.arrow.size = 0.5,
     vertex.size = 20 * omega,
     vertex.label.dist = 1.5,
     vertex.label.color = "black",
     vertex.frame.color = NA)
