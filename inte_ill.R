# nolint start
# Test if feasoveralp library works well in the assembly codes
rm(list = ls())

source("toolbox.R") #Remember to exclude all the conflicting functions
source("overlap.R")

num <- 5; stren <- 2; conne <- 1
relation <- tibble(ep=numeric(), cond=numeric())

for (rep in 1:400) {
  eps <- 10^runif(1, -4, -1)
  A <- interaction_matrix_ill(num, stren, conne, eps)
  relation <- relation %>% add_row(ep = log(eps,10), cond = log(kappa(A,10)))
}

ggplot(relation, aes(ep, cond))+
  geom_point()+
  geom_smooth(method = lm)