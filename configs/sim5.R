sim_config <- list(
  n_obs = 1000,
  n_vars = 500,
  group_sizes = c(100, 100, 100, 100, 100),  # Custom group sizes
  true_vars_per_group = c(1, 1, 1, 1, 1),
  n_groups =  length(group_sizes),
  group_start = cumsum(c(1, head(group_sizes, -1))) ,
  group_structure =  rep(1:n_groups, times = group_sizes)
)
