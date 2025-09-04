fit_models <- function(train_data, test_data, true_coefs, group_structure) {
  X_train <- as.matrix(train_data[, -c(1, 2)])
  y_train <- survival::Surv(train_data$time, train_data$status)
  X_test  <- as.matrix(test_data[, -c(1, 2)])
  y_test  <- survival::Surv(test_data$time, test_data$status)

  results_list <- list()

  ############################################################
  # 1. Elastic Net (glmnet with alpha tuning)
  ############################################################
  library(glmnet)
  alpha_grid <- seq(0, 1, by = 0.1)
  k <- 10
  set.seed(1)
  foldid <- sample(rep(1:k, length.out = nrow(X_train)))

  cv_results <- lapply(alpha_grid, function(a) {
    cvfit <- cv.glmnet(X_train, y_train, family = "cox", alpha = a, foldid = foldid)
    list(cvm_min = min(cvfit$cvm),
         lambda_min = cvfit$lambda[which.min(cvfit$cvm)],
         cvfit = cvfit,
         alpha = a)
  })

  cvm_vals <- vapply(cv_results, `[[`, numeric(1), "cvm_min")
  best_idx <- which.min(cvm_vals)
  best_alpha <- cv_results[[best_idx]]$alpha
  best_lambda <- cv_results[[best_idx]]$lambda_min
  cvfit_best <- cv_results[[best_idx]]$cvfit

  fit_lasso <- glmnet(X_train, y_train, family = "cox",
                      alpha = best_alpha, lambda = best_lambda)

  pred_lasso_train <- predict(cvfit_best, newx = X_train, s = cvfit_best$lambda.min, type = "link")
  pred_lasso_test  <- predict(cvfit_best, newx = X_test,  s = cvfit_best$lambda.min, type = "link")
  selected_elasticnet <- coef(fit_lasso)

  results_list$ElasticNet <- tryCatch({
    evaluate_full_model(true_coefs, selected_elasticnet,
                        X_train, X_test,
                        y_train, y_test,
                        pred_lasso_train, pred_lasso_test,
                        train_data, test_data,
                        cvfit_best, cvfit_best$lambda.min,
                        eval_time = 5, tmax_auc = 10,
                        model = "glmnet")
  }, error = function(e) NULL)

  ############################################################
  # 2. Exclusive Lasso
  ############################################################
  library(ExclusiveLasso)
  cvfit_elasso <- cv.exclusive_lasso(
    X_train, y_train,
    groups = group_structure,
    family = "cox",
    type.measure = "cindex",
    skip_df = TRUE,
    intercept = FALSE,
    thresh = 1e-5
  )

  fit_elasso <- exclusive_lasso(
    X_train, y = y_train,
    lmd = cvfit_elasso$lmd.min,
    groups = group_structure,
    family = "cox",
    intercept = FALSE
  )

  pred_elasso_train <- predict(fit_elasso, newx = X_train, s = cvfit_elasso$lmd.min, type = "link")
  pred_elasso_test  <- predict(fit_elasso, newx = X_test,  s = cvfit_elasso$lmd.min, type = "link")
  selected_el <- coef(fit_elasso)[-1]

  results_list$ELasso <- tryCatch({
    evaluate_full_model(true_coefs, selected_el,
                        X_train, X_test,
                        y_train, y_test,
                        pred_elasso_train, pred_elasso_test,
                        train_data, test_data,
                        cvfit_elasso, cvfit_elasso$lmd.min,
                        eval_time = 5, tmax_auc = 10,
                        model = "glmnet")
  }, error = function(e) NULL)

  ############################################################
  # 3. Group Lasso
  ############################################################
  library(grpreg)
  cvfit_gl <- cv.grpsurv(X_train, y_train, group = group_structure, penalty = "grLasso", alpha = 1)

  fitgr <- grpsurv(X_train, y_train, group = group_structure,
                   penalty = "grLasso", alpha = 1, lambda = cvfit_gl$lambda.min)

  fit_gl1 <- grpsurv(X_train, y_train, group = group_structure,
                     penalty = "grLasso", alpha = 1)

  pred_gl_train <- predict(fit_gl1, X = X_train, type = "link", lambda = cvfit_gl$lambda.min)
  pred_gl_test  <- predict(fit_gl1, X = X_test,  type = "link", lambda = cvfit_gl$lambda.min)
  selected_gl <- coef(fitgr)

  results_list$GLasso <- tryCatch({
    evaluate_full_model(true_coefs, selected_gl,
                        X_train, X_test,
                        y_train, y_test,
                        matrix(pred_gl_train, ncol = 1), matrix(pred_gl_test, ncol = 1),
                        train_data, test_data,
                        cvfit_gl, cvfit_gl$lambda.min,
                        eval_time = 5, tmax_auc = 10,
                        model = "gl")
  }, error = function(e) NULL)

  ############################################################
  # 4. IPF-Lasso
  ############################################################
  library(ipflasso)

  group_sizes <- as.numeric(table(group_structure))
  group_ends <- cumsum(group_sizes)
  group_starts <- c(1, head(group_ends, -1) + 1)
  blocks <- mapply(function(start, end) start:end, group_starts, group_ends, SIMPLIFY = FALSE)

  # Combined Ridge for penalty factors
  cvfit_combined_ridge <- cv.glmnet(X_train, y_train, alpha = 0, family = "cox")
  coefs_combined_ridge <- abs(as.numeric(coef(cvfit_combined_ridge, s = "lambda.min")))
  split_indices <- rep(1:length(group_sizes), times = group_sizes)
  grouped_coefs <- split(coefs_combined_ridge, split_indices)
  means_combined_ridge <- sapply(grouped_coefs, mean)
  pflist <- c(1, means_combined_ridge[1] / means_combined_ridge[2])

  cvipf_fit <- cvr.ipflasso(
    X = X_train,
    Y = y_train,
    family = "cox",
    blocks = blocks,
    pf = pflist,
    nfolds = 5,
    ncv = 10,
    type.measure = "partial likelihood"
  )

  ipf_coefficients <- cvipf_fit$coeff[-1, cvipf_fit$ind.bestlambda]
  pred_ipf_train <- ipflasso.predict(object = cvipf_fit, Xtest = X_train)$linpredtest
  pred_ipf_test  <- ipflasso.predict(object = cvipf_fit, Xtest = X_test)$linpredtest
  selected_ipf <- cvipf_fit$coeff[-1, cvipf_fit$ind.bestlambda]

  results_list$IPF <- tryCatch({
    evaluate_full_model(true_coefs, selected_ipf,
                        X_train, X_test,
                        y_train, y_test,
                        pred_ipf_train, pred_ipf_test,
                        train_data, test_data,
                        cvipf_fit, selected_ipf,
                        eval_time = 5, tmax_auc = 10,
                        model = "ipf")
  }, error = function(e) NULL)

  return(results_list)
}
