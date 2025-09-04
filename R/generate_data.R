library(MASS)

# Covariance generator
generate_cov_matrix <- function(n_vars, group_sizes, group_start) {
  Sigma <- matrix(0, nrow = n_vars, ncol = n_vars)
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      group_i <- max(which(group_start <= i))
      group_j <- max(which(group_start <= j))
      if (group_i == group_j) {
        Sigma[i, j] <- .6^abs(i - j)  # within-group correlation
      } else {
        Sigma[i, j] <- .3^abs(i - j)  # between-group correlation
      }
    }
  }
  return(Sigma)
}

# Survival data generator
generate_survival_data <- function(X, true_coefs, n_obs) {
  lin_pred <- X %*% true_coefs
  baseline_hazard <- 0.08
  hazard <- baseline_hazard * exp(lin_pred)
  surv_time <- -log(runif(n_obs)) / hazard
  censoring_time <- rexp(n_obs, rate = 0.02)
  event <- as.numeric(surv_time <= censoring_time)
  observed_time <- pmin(surv_time, censoring_time)
  list(time = observed_time, status = event, y = survival::Surv(observed_time, event))
}

# Data wrapper
generate_data <- function(sim_config) {
  with(sim_config, {
    Sigma <- generate_cov_matrix(n_vars, group_sizes, group_start)
    X <- MASS::mvrnorm(n = n_obs, mu = rep(0, n_vars), Sigma = Sigma)

    true_coefs <- numeric(n_vars)
    true_indices <- c()
    for (i in 1:n_groups) {
      start_index <- group_start[i]
      end_index <- start_index + group_sizes[i] - 1
      true_vars <- sample(start_index:end_index, true_vars_per_group[i])
      true_indices <- c(true_indices, true_vars)
    }
    true_coefs[true_indices] <- runif(length(true_indices), 0.5, 1.5)

    survival_data <- generate_survival_data(X, true_coefs, n_obs)
    data <- data.frame(time = survival_data$time, status = survival_data$status, X)
    list(data = data, true_coefs = true_coefs)
  })
}
