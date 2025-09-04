evaluate_full_model <- function(
    true_coefs, selected,
    X_train, X_test,
    y_train, y_test,
    pred_train, pred_test,
    training_data, testing_data,
    cvfit_model, lambda.min,
    eval_time = 5, tmax_auc = 10,
    model = "glmnet"
) {
  # === Selection metrics ===
  true_coefs <- as.logical(true_coefs)
  selected <- as.logical(selected)
  TP <- sum(true_coefs & selected)
  FP <- sum(selected & !true_coefs)
  FN <- sum(true_coefs & !selected)
  TN <- sum(!true_coefs & !selected)

  precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  recall <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1_score <- if ((precision + recall) > 0) 2 * (precision * recall) / (precision + recall) else 0

  selection_results <- list(
    TP = TP, FP = FP, FN = FN, TN = TN,
    precision = precision, recall = recall, f1_score = f1_score
  )

  # === Survival metrics ===
  library(survival)
  library(risksetROC)
  library(riskRegression)
  library(survcomp)

  Cindex_train <- survcomp::concordance.index(pred_train[, 1], surv.time = y_train[, 1], surv.event = y_train[, 2])$c.index
  Cindex_test  <- survcomp::concordance.index(pred_test[, 1], surv.time = y_test[, 1], surv.event = y_test[, 2])$c.index

  AUC_train <- risksetAUC(Stime = y_train[, "time"], status = y_train[, "status"], marker = pred_train, tmax = max(y_train[, "time"]))
  AUC_test  <- risksetAUC(Stime = y_test[, "time"], status = y_test[, "status"], marker = pred_test, tmax = max(y_test[, "time"]))

  iAUC_train <- IntegrateAUC(AUC_train$AUC, AUC_train$utimes, AUC_train$St, tmax = tmax_auc)
  iAUC_test  <- IntegrateAUC(AUC_test$AUC, AUC_test$utimes, AUC_test$St, tmax = tmax_auc)

  survival_results <- list(
    Cindex = list(train = Cindex_train, test = Cindex_test),
    iAUC   = list(train = iAUC_train, test = iAUC_test)
  )

  return(list(selection = selection_results, survival = survival_results))
}
