# nolint start
# Test if the extened 2d code works well to compute overlap
source("overlap.R")

# Use brute-force method to compute omega_overlap
sampling_omega_overlap <- function(A,B,nsample = 10^5) {
  mat <- solve(B) %*% A
  abundance_all <- rmvnorm(n = nsample, mean = rep(0, nrow(A))) %>% 
    {abs(./sqrt(rowSums(.^2)))}
  get_feasibility <- function(N_A){
    N_B <- mat %*% matrix(N_A,ncol=1) %>% c()
    if_else(sum(N_B >= -1e-10) == length(N_B), 1, 0) 
  }
  percent <- 1:nsample %>% 
    map_dbl(~get_feasibility(abundance_all[.x,])) %>% 
    mean()

  return(calculate_omega(A) * percent^(1/nrow(A))) 
}

# Visualization of overlap in 2d cases
draw_omega_overlap <- function(A, B, canvas){
  span_A <- generate_span_vectors(A)
  span_B <- generate_span_vectors(B)
  canvas <- canvas + 
    geom_segment(aes(x=0,y=0, xend=span_A[1,],yend=span_A[2,],colour="A"),arrow = arrow()) +
    geom_segment(aes(x=0,y=0, xend=span_B[1,],yend=span_B[2,],colour="B"),arrow = arrow())
  canvas
}

# Analytical results for overlap in 2d
d2_omega <- function(A){
  if(nrow(A) != 2) {
    warning("Not computing 2d cases")
    return(0)
  }
  else {
    span_A <- generate_span_vectors(A)
    d2_O <- abs(atan(span_A[1,2]/span_A[1,1])-atan(span_A[2,2]/span_A[2,1])/(2*pi))
    return(d2_O)
  }
}

set.seed(2022)
A <- interaction_matrix_random(2,2,0.6)
B <- interaction_matrix_random(2,2,0.6)

sampling_omega_overlap(A,B)
calculate_omega_overlap(A,B)
#d2_omega(A); d2_omega(B)

canvas <- ggplot() + coord_fixed(xlim=c(-1,1),ylim=c(-1,1))
draw_omega_overlap(A,B,canvas)
# nolint end