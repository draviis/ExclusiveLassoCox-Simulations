generate_data_categorical <- function(sim_config) {
  with(sim_config, {
    # 1. Covariance + draw
    Sigma <- generate_cov_matrix(n_vars, group_sizes, group_start)
    X_all <- MASS::mvrnorm(n = n_obs, mu = rep(0, n_vars), Sigma = Sigma)
    colnames(X_all) <- paste0("V", 1:n_vars)

    # 2. Groups
    X_group1 <- X_all[, 1:group_sizes[1]]
    X_group2_raw <- X_all[, (group_sizes[1] + 1):(group_sizes[1] + group_sizes[2])]
    X_group3 <- X_all[, (group_sizes[1] + group_sizes[2] + 1):n_vars]

    # 3. Convert group2 (1 var) → categorical → 3 dummies
    breaks_4cat <- quantile(X_group2_raw, probs = c(0, 0.25, 0.5, 0.75, 1))
    X_cat <- cut(X_group2_raw, breaks = breaks_4cat, include.lowest = TRUE,
                 labels = c("cat1", "cat2", "cat3", "cat4"))
    X_cat_dummies <- model.matrix(~ X_cat)[, -1]  # 3 dummy vars
    colnames(X_cat_dummies) <- c("cat2", "cat3", "cat4")

    # 4. Final design matrix
    X_final <- cbind(X_group1, X_cat_dummies, X_group3)
    X <- scale(X_final)

    # 5. New group structure: group1 + dummies = group 1, group3 = group 2
    group_structure <- c(
      rep(1, ncol(X_group1)), 
      rep(1, ncol(X_cat_dummies)), 
      rep(2, ncol(X_group3))
    )

    # 6. Assign coefficients
    true_coefs <- numeric(ncol(X))
    true_indices <- c()
    group_sizes_adj <- c(ncol(X_group1) + ncol(X_cat_dummies), ncol(X_group3))

    # adjust to ensure correct sampling from dummy + continuous groups
    start_index <- 1
    for (i in seq_along(group_sizes_adj)) {
      end_index <- start_index + group_sizes_adj[i] - 1
      true_vars <- sample(start_index:end_index, true_vars_per_group[i])
      true_indices <- c(true_indices, true_vars)
      start_index <- end_index + 1
    }
    true_coefs[true_indices] <- runif(length(true_indices), 0.5, 1.5)

    # 7. Survival times
    survival_data <- generate_survival_data(X, true_coefs, n_obs)
    data <- data.frame(time = survival_data$time, status = survival_data$status, X)

    list(data = data, true_coefs = true_coefs, group_structure = group_structure)
  })
}
