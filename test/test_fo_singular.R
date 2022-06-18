# nolint start
# This is a demonstration code showing how overlap could behave with nearly singular matrices

source("toolbox.R") # feaoverlap pkg/source included
source("wrapper.R") # change defaults to raw omega if using feaoverlap

interaction_matrix_ill <- function(num = 4, stren = 2, conne = 0.9, epsilon = 10^(-4), threshold = 0.001) {
  inte <- interaction_matrix_random(num, stren, conne)
  new_col <- floor(num/2) + 1
  inte[new_col,1] <- (max(abs(inte[new_col,1]), threshold)) #set a threshold to avoid zero entry
  #inte[new_col,1] <- (max(inte[new_col,1], threshold))
  factor <- (-1)/inte[new_col,1]
  inte[,new_col] <- factor*inte[,1] + epsilon * rnorm(num)
  inte[new_col,new_col] <- -1
  return(inte)
}

# Draw the plot of condition num. vs. eps.
eps_cond <- function(num = 5, stren = 2, conne = 1){
  relation <- tibble(ep=numeric(), cond=numeric())

  for (rep in 1:400) {
    eps <- 10^runif(1, -4, -1)
    A <- interaction_matrix_ill(num, stren, conne, eps)
    relation <- relation %>% add_row(ep = log(eps,10), cond = log(kappa(A,10)))
  }

  ggplot(relation, aes(ep, cond))+
    geom_point()+
    geom_smooth(method = lm)
}



set.seed(49)
eps_cond()

A <- interaction_matrix_ill()
B <- A[1:(num-1),1:(num-1)] %>% 
  cbind(c(rep(0,num-1))) %>% 
  rbind(c(rep(0,num-1),-1))

calculate_omega_overlap(A,B)
min(calculate_omega(A), calculate_omega(B))

ratio <- calculate_omega_overlap(A,B)/min(calculate_omega(A), calculate_omega(B))

# ratio of raw-omega-value: seems to be scaled with epsilon
c(ratio,ratio^num)

# nolint end