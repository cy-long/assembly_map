# nolint start
library(purrr)
calculate_omega <- partial(feasoverlap::calculate_omega, raw = TRUE)
calculate_omega_overlap <- partial(feasoverlap::calculate_omega_overlap, raw = TRUE)
# nolint end