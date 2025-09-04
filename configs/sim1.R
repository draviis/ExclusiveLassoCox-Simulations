sim_config <- list(
  n_obs = 1000,
  n_vars = 500,
  group_sizes = c(250, 250),
  true_vars_per_group = c(2, 3),
  n_groups = 2,
  group_start = cumsum(c(1, head(c(250, 250), -1))),
  group_structure = rep(1:2, times = c(250, 250))
)
