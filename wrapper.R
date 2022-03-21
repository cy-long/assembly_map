# nolint start
# library(purrr)
# calculate_omega <- partial(feasoverlap::calculate_omega, raw = TRUE)
# calculate_omega_overlap <- partial(feasoverlap::calculate_omega_overlap, raw = TRUE)
calculate_omega <- function(vertex, raw = TRUE, nsamples = 100){
    if(is.null(ncol(vertex))) return(1/2)
    else {
        feasoverlap::calculate_omega(vertex, raw, nsamples)
    }
}
calculate_omega_overlap <- function(A, B, raw = TRUE, nsamples = 100){
    feasoverlap::calculate_omega_overlap(A, B, raw, nsamples)
}
# nolint end