# nolint start

source("toolbox.R")
source("wrapper.R")

num <- 5; stren <- 1; conne <- 0.5; order <- 1
set.seed(264)
A <- interaction_matrix_random(num, stren, conne)
ori <- c(1,2)
targ <- c(1,2,3,4,5)
ori_mat <- A[ori,ori]
targ_mat <- A[targ,targ]


sum(Omega_overlap_ext(ori_mat, targ_mat, ori, targ, 1))
calculate_omega(ori_mat)
calculate_omega(targ_mat)
sum(Omega_overlap_ext(ori_mat, targ_mat, ori, targ, 2))

# nolint end
