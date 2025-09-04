sim_config <- list(
  n_obs = 1000,
  n_vars = 500,               # total original vars before transformation
  group_sizes = c(9, 1, 490), # note: group2 becomes categorical → 3 dummies
  true_vars_per_group = c(6, 4, 10),
  n_groups = 3,
  group_start = cumsum(c(1, head(c(9, 1, 490), -1)))
)
