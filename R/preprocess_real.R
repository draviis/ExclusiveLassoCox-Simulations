# --- HNSC example ---
preprocess_hnsc <- function(path) {
  df_fulln <- read.csv(path, row.names = 1)
  df_fulln <- df_fulln[complete.cases(df_fulln$CLIN_overall_survival,
                                      df_fulln$CLIN_status), ]

  # Clinical subset
  df_clin <- df_fulln[, startsWith(colnames(df_fulln), "CLIN_"), drop = FALSE]
  df_mir  <- df_fulln[, startsWith(colnames(df_fulln), "MIR_"), drop = FALSE]
  df_rna  <- df_fulln[, startsWith(colnames(df_fulln), "RNA_"), drop = FALSE]
  df_mut  <- df_fulln[, startsWith(colnames(df_fulln), "MUT_"), drop = FALSE]

  # Fix clinical vars
  df_clin$CLIN_Tumor_purity <- as.numeric(df_clin$CLIN_Tumor_purity)
  df_clin$CLIN_years_to_birth <- as.numeric(df_clin$CLIN_years_to_birth)

  df_clin$time <- as.numeric(df_clin$CLIN_overall_survival)
  df_clin$status <- as.numeric(df_clin$CLIN_status)
  df_clin$CLIN_overall_survival <- NULL
  df_clin$CLIN_status <- NULL
  df_clin$time[df_clin$time == 0] <- 0.01

  # Convert gender to 0/1
  df_clin$CLIN_gender <- factor(df_clin$CLIN_gender)
  df_clin$CLIN_gender <- as.integer(df_clin$CLIN_gender == "male")

  # Scale omics
  df_mir <- as.data.frame(scale(df_mir))
  df_rna <- as.data.frame(scale(df_rna))
  df_mut <- as.data.frame(scale(df_mut))

  # Merge
  df_merged <- cbind(df_clin, df_mir, df_rna, df_mut)

  # Group structure
  colnames_X <- colnames(df_merged)[!colnames(df_merged) %in% c("time", "status")]
  g1 <- sum(grepl("^CLIN_", colnames_X))
  g2 <- sum(grepl("^MIR_", colnames_X))
  g3 <- sum(grepl("^RNA_", colnames_X))
  g4 <- sum(grepl("^MUT_", colnames_X))
  group_structure <- c(rep(1, g1), rep(2, g2), rep(3, g3), rep(4, g4))

  list(data = df_merged, group_structure = group_structure)
}

# --- GSE31684 example ---
preprocess_gse31684 <- function(gse_id) {
  gse <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  exprSet <- Biobase::exprs(gse[[1]])
  normalized_exprSet <- limma::normalizeBetweenArrays(exprSet, method = "quantile")

  # Top 5000 MAD
  gene_variability <- apply(normalized_exprSet, 1, mad)
  selected_genes <- names(sort(gene_variability, decreasing = TRUE))[1:5000]
  filtered_exprSet <- normalized_exprSet[selected_genes, ]
  filtered_exprSet <- t(impute::impute.knn(t(filtered_exprSet))$data)

  dfn <- as.data.frame(t(filtered_exprSet))
  pheno_data <- pData(gse[[1]])
  colData <- pheno_data[, c("survival.months:ch1","last known status:ch1",
                            "age at rc:ch1","rc_stage:ch1",
                            "nomogram score:ch1","smokingpack-years:ch1")]
  df_full <- cbind(dfn, colData)

  # Preprocess survival + clinical
  df_full$`last known status:ch1` <- ifelse(df_full$`last known status:ch1` == "DOD", 1, 0)
  df_fulln <- df_full[!is.na(df_full[["survival.months:ch1"]]), ]
  df_fulln$time <- as.numeric(df_fulln$`survival.months:ch1`)
  df_fulln$status <- as.numeric(df_fulln$`last known status:ch1`)
  df_fulln$time[df_fulln$time == 0] <- 0.01

  # Example: dummy encode stage
  dummies_rc_stage <- model.matrix(~ `rc_stage:ch1`, data = df_fulln)[, -1]
  df_fulln <- cbind(df_fulln, dummies_rc_stage)
  df_fulln$`rc_stage:ch1` <- NULL

  # Scale only continuous (skip dummies)
  X <- df_fulln[, !colnames(df_fulln) %in% c("time", "status")]
  is_dummy <- function(x) is.numeric(x) && all(na.omit(x) %in% c(0, 1))
  dummy_cols <- sapply(X, is_dummy)
  X_scaled <- X
  X_scaled[!dummy_cols] <- scale(X[!dummy_cols])
  df_merged <- cbind(X_scaled, df_fulln[, c("time", "status")])

  # Group structure (simplified: 5000 expr vs 7 clinical)
  g1 <- 5000
  g2 <- ncol(df_merged) - g1 - 2
  group_structure <- c(rep(1, g1), rep(2, g2))

  list(data = df_merged, group_structure = group_structure)
}

