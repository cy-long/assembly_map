# nolint start

# ------ Initialize ------
rm(list = ls())
source("toolbox.R") # source feasoverlap pkg or overlap.R file
source("wrapper.R") # change defaults to raw omega when using feaoverlap

num <- 4; stren <- 1; conne <- 0.8;

# generate indices denoting the subcommunities
sub_coms <- generate_subcoms(num)
t_ind <- generate_t_index(num)
l_ind <- generate_l_index(num)
names <- generate_names(num, sub_coms)

# generate interaction matrix
set.seed(35)
A <- interaction_matrix_random(num, stren, conne)

# read interaction matrix
# A <- as.matrix(read.table("data/Gould_Matrix_median_12.csv",sep=","))[1:5,1:5]
# num <- ncol(A)


# ------ Compute raw/norm/overlap omega values ------
# Compute raw/norm omega for each node
r_omega <- evaluate_solo(A, TRUE)
n_omega <- evaluate_solo(A, FALSE)
# Compute omega_overlaps between possible nodes
Overlap <- evaluate_dual(A, TRUE, 2)
overlap <- evaluate_dual(A, FALSE, 2)
diag(Overlap) <- r_omega
diag(overlap) <- n_omega
# Compute extended omega for each node
s_omega <- evaluate_dual(A, TRUE, 1)
ns_omega <- evaluate_dual(A, FALSE, 1)


# ------  Pathwise probability analysis ------
# specify initial and final subcommunities
ti <- 1; tf <- 2^num
# find paths between i/f subcommunites
paths <- find_path(ti, tf)

# quantifying paths with RAW Omegas
# path_info <- quantify_paths(paths, r_omega, s_omega, Overlap)

# quantifying paths with NORMALIZED omegas
path_info <- quantify_paths(paths, n_omega, ns_omega, overlap)

path_pr <- select(path_info, random_pr, environ_pr, species_pr) %>% as.matrix()
view(path_pr)


# ------ Operations on matrices (Markov/Information) ------
mat_oper <- matrix(0, 2^num, 2^num)
mat_oper <- t(Markov_norm(Overlap, r_omega))

# Stationary Distribution
sd <- as.numeric(Stat_dist())
# Entropy
Ent_a <- Entropy(mat_oper) #absolute entropy
Ent_r <- Entropy(mat_oper)/log(nrow(mat_oper)) #relative entropy
# One-species
Pr <- Species_pr(num, sd, names)

# nolint end