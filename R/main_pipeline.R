source("R/generate_data.R")
source("R/evaluate_functions.R")
source("R/fit_models.R")

run_simulation <- function(sim_config, sim_name, n_iter = 50, categorical = FALSE) {
  all_results <- vector("list", n_iter)
  
  for (i in 1:n_iter) {
    cat("Iteration", i, "for", sim_name, "...\n")
    
    # Data generation: choose correct function
    gen <- if (categorical) {
      generate_data_categorical(sim_config)
    } else {
      generate_data(sim_config)
    }
    
    data <- gen$data
    true_coefs <- gen$true_coefs
    group_structure <- gen$group_structure %||% sim_config$group_structure  # fallback
    
    # Train/test split
    train_indices <- sample(1:nrow(data), size = sim_config$n_obs / 2)
    test_indices  <- setdiff(1:nrow(data), train_indices)
    train_data <- data[train_indices, ]
    test_data  <- data[test_indices, ]
    
    # Fit models
    res <- fit_models(train_data, test_data, true_coefs, group_structure)
    all_results[[i]] <- res
  }
  
  saveRDS(all_results, file = paste0("results/", sim_name, "_results.rds"))
}
