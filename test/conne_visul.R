# nolint start

# ------ Visualize ------
threshold <- 0.5
mat_disp <- Overlap

# Specify row/col names for transition matrix
names <- c("0")
for (s in 1:num) {
  names <- c(names, c(apply(sub_coms[[s]], 2, convert2names)))
}
colnames(mat_disp) = rownames(mat_disp) = names

mat_adj <- apply(mat_disp, c(1,2), set_conne, threshold)

# Get the value to display
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
# Generate the adjacency matrix with the threshold
set_conne <- function(value, threshold){
  if(!is.na(value) && value>threshold){
    return(1)
  }else{
    return(0)
  }
  
}
# Display the graph
diag(mat_adj) <- 0
palf <- colorRampPalette(c("dodgerblue","darkorange"))
pal1 <- heat.colors(length(n_omega), alpha=1) 
pal2 <- palf(length(n_omega))
r <- rank(n_omega)
network <- graph_from_adjacency_matrix(mat_adj)

plot(network,
     layout = grid,
     edge.arrow.size = 0.4,
     edge.width = 3 * value / max(value),
     vertex.label.color = "black",
     vertex.frame.color = "black",
     vertex.color = pal2[r],
    # vertex.size = 30 * r_omega / max(r_omega)
    # vertex.size = 20 * (n_omega) / max(n_omega)
)

# nolint end