configs <- list.files("configs", full.names = TRUE)

for (cfg in configs) {
  source(cfg)  # loads sim_config
  sim_name <- gsub("\\.R$", "", basename(cfg))
  source("R/main_pipeline.R")
  
  # choose generator based on file name
  if (grepl("categorical", sim_name, ignore.case = TRUE)) {
    run_simulation(sim_config, sim_name, n_iter = 50, categorical = TRUE)
  } else {
    run_simulation(sim_config, sim_name, n_iter = 50, categorical = FALSE)
  }
}
