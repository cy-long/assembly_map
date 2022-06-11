# nolint start
# Test if the extened 2d code works well to compute overlap
source("toolbox.R") # feaoverlap pkg/source included
source("wrapper.R") # change defaults to raw omega if using feaoverlap

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

  return(calculate_omega(A) * percent) 
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
    d2_O <- abs(acos(span_A[1,1]*span_A[1,2]+span_A[2,1]*span_A[2,2]))/(2*pi) #Always computes acute angles
    return(d2_O)
  }
}

x <- c(); y <- c(); z <- c()
for (i in 1:100){
  set.seed(89+64+i^2+i)
  A <- interaction_matrix_ill(2,1,1,10^-5)
  B <- interaction_matrix_random(2,1,1)
  tmp <- calculate_omega_overlap(A,B)
  if(abs(tmp)>10^(-5)){
    x <- c(x,sampling_omega_overlap(A,B))
    y <- c(y,tmp)
    AB <- vertex_detection(A,B)
    z <- c(z,d2_omega(AB))
  }
}

error1 <- cbind(z,t-z,1)
error2 <- cbind(z,y-z,2)
error <- as.data.frame(rbind(error1,error2))

colnames(error) <- c("Omega","Error","Method")
error$Method <- factor(error$Method)
g <- ggplot(data = error, mapping = aes(x = Omega, y = Error, colour = Method))+geom_point(size = 4)
g+theme(axis.text=element_text(size=14),
      axis.title=element_text(size=16))

# canvas <- ggplot() + coord_fixed(xlim=c(-1,1),ylim=c(-1,1))
# draw_omega_overlap(A,B,canvas)

#> The calculate results are closer analytical results than sampling are.
#> error(calculate) < error(sampling) in this seed
# nolint end