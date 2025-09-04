############################################################
## Fair tuning under equal wall-clock budgets (full script)
############################################################
set.seed(1)
suppressPackageStartupMessages({
  library(ExclusiveLasso)
  library(survival)
  library(glmnet)
  library(grpreg)
  library(ipflasso)
  library(MASS)
  library(R.utils)  # withTimeout
})

## -----------------------------
## 0) Simulated data (your setup)
## -----------------------------
n_obs  <- 1000
n_vars <- 500
group_sizes <- c(100,100,100,100,100)
n_groups <- length(group_sizes)
true_vars_per_group <- rep(1, n_groups)

group_start <- cumsum(c(1, head(group_sizes, -1)))

generate_cov_matrix <- function(n_vars, group_sizes) {
  Sigma <- matrix(0, nrow = n_vars, ncol = n_vars)
  for (i in 1:n_vars) for (j in 1:n_vars) {
    gi <- max(which(group_start <= i))
    gj <- max(which(group_start <= j))
    Sigma[i, j] <- if (gi == gj) .6^abs(i - j) else .3^abs(i - j)
  }
  Sigma
}
generate_survival_data <- function(X, true_coefs, n_obs) {
  lp <- as.vector(X %*% true_coefs)
  hazard <- 0.08 * exp(lp)
  t_ev <- -log(runif(n_obs)) / hazard
  t_c  <- rexp(n_obs, rate = 0.02)
  status <- as.numeric(t_ev <= t_c)
  time <- pmin(t_ev, t_c)
  list(time = time, status = status, y = Surv(time, status))
}
split_data <- function(data, n_train = 500) {
  tr <- sample(seq_len(nrow(data)), n_train)
  list(train = data[tr, ], test = data[-tr, ])
}

Sigma <- generate_cov_matrix(n_vars, group_sizes)
X <- mvrnorm(n = n_obs, mu = rep(0, n_vars), Sigma = Sigma)

true_coefs <- numeric(n_vars)
for (g in seq_len(n_groups)) {
  s <- group_start[g]; e <- s + group_sizes[g] - 1
  idx <- sample(s:e, true_vars_per_group[g])
  true_coefs[idx] <- runif(length(idx), 0.5, 1.5)
}
sd <- generate_survival_data(X, true_coefs, n_obs)
dat <- data.frame(time = sd$time, status = sd$status, X)
sp  <- split_data(dat, n_train = 500)
train_data <- sp$train; test_data <- sp$test

X_train <- scale(as.matrix(train_data[, -c(1,2)]))
y_train <- Surv(train_data$time, train_data$status)
X_test  <- scale(as.matrix(test_data[, -c(1,2)]))
y_test  <- Surv(test_data$time, test_data$status)

group_structure <- rep(seq_len(n_groups), times = group_sizes)

## Fixed folds for methods that support it
set.seed(123)
k <- 10
foldid <- sample(rep(1:k, length.out = nrow(X_train)))

## ----------------------------------------
## 1) Utilities
## ----------------------------------------
now <- function() format(Sys.time(), "%H:%M:%S")
finite_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x_ok <- is.finite(x)
  if (!any(x_ok)) return(NA_real_)
  max(x[x_ok])
}
finite_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x_ok <- is.finite(x)
  if (!any(x_ok)) return(NA_real_)
  min(x[x_ok])
}

## Use glmnet once to anchor a sensible lambda range, then reuse
lambda_range_from_glmnet <- function(X, y) {
  fit <- glmnet::glmnet(X, y, family = "cox", alpha = 1, nlambda = 100)
  c(max(fit$lambda), min(fit$lambda))
}
exp_lambda_path_from_range <- function(lam_max, lam_min, length_out) {
  exp(seq(log(lam_max), log(lam_min), length.out = length_out))
}

## ----------------------------------------
## 2) Global wall-clock budget
## ----------------------------------------
time_budget_sec <- 2*60*60   # 2 hours (set 60 for a quick smoke test)

## ----------------------------------------
## 3) Elastic Net tuner (cycle + refine until time runs out)
## ----------------------------------------
tune_en_until <- function(time_budget_sec) {
  start <- Sys.time()
  best  <- list(score = -Inf, fit = NULL, meta = NULL)
  tried <- 0L; timings <- c()
  alphas <- seq(0, 1, by = 0.25)
  L_seq  <- c(50L, 100L, 200L)
  
  repeat {
    for (L in L_seq) {
      for (a in alphas) {
        elapsed <- as.numeric(difftime(Sys.time(), start, units="secs"))
        if (elapsed >= time_budget_sec) break
        per_call <- max(5, time_budget_sec - elapsed - 1)
        
        t1 <- Sys.time()
        res <- try(withTimeout({
          cvfit <- glmnet::cv.glmnet(
            X_train, y_train, family = "cox", alpha = a,
            foldid = foldid,  nlambda = L
          )
          sc <- finite_max(cvfit$cvm)
          if (is.na(sc)) NULL else list(
            score = sc,
            fit   = cvfit,
            meta  = list(alpha = a, nlambda = L, lambda.min = cvfit$lambda.min)
          )
        }, timeout = per_call, onTimeout = "silent"), silent = TRUE)
        t2 <- Sys.time()
        
        if (!inherits(res, "try-error") && !is.null(res)) {
          tried <- tried + 1L
          timings <- c(timings, as.numeric(difftime(t2, t1, units="secs")))
          if (res$score > best$score) best <- res
        }
      }
      if (as.numeric(difftime(Sys.time(), start, units="secs")) >= time_budget_sec) break
    }
    if (as.numeric(difftime(Sys.time(), start, units="secs")) >= time_budget_sec) break
    ## refine
    L_seq  <- unique(pmin(5000L, c(L_seq, L_seq*2L)))
    alphas <- sort(unique(c(alphas, (head(alphas,-1)+tail(alphas,-1))/2)))
  }
  
  list(
    best_fit = best$fit, best_meta = best$meta, best_score = best$score,
    models_tried = tried,
    mean_time_per_try = if (length(timings)) mean(timings) else NA_real_,
    total_time = as.numeric(difftime(Sys.time(), start, units="secs")),
    method = "Elastic Net"
  )
}

## ----------------------------------------
## 4) Exclusive Lasso tuner (Cox OK, no foldid) — cycle + refine
##    Your build uses nlmd / lmd.min.ratio in CV; final fit uses lmd.
## ----------------------------------------
tune_el_until <- function(time_budget_sec) {
  start <- Sys.time()
  best  <- list(score = -Inf, fit = NULL, meta = NULL)
  tried <- 0L; timings <- c()
  nlmds <- c(50L, 100L, 200L)
  lmrs  <- c(1e-2, 1e-3)
  
  repeat {
    for (nlmd in nlmds) {
      for (lmr in lmrs) {
        elapsed <- as.numeric(difftime(Sys.time(), start, units="secs"))
        if (elapsed >= time_budget_sec) break
        per_call <- max(5, time_budget_sec - elapsed - 1)
        
        t1 <- Sys.time()
        res <- try(withTimeout({
          set.seed(2025)  # reproducible internal fold split
          cvfit <- ExclusiveLasso::cv.exclusive_lasso(
            X = X_train, y = y_train, groups = group_structure,
            family = "cox", type.measure = "cindex",
            intercept = FALSE, skip_df = TRUE,
            nfolds = 10L,
            nlmd = nlmd,
            lmd.min.ratio = lmr
          )
          ## Try typical fields in different builds
          sc <- NA_real_
          if (!is.null(cvfit$cvm)) sc <- finite_max(cvfit$cvm)
          if (is.na(sc) && !is.null(cvfit$cve)) sc <- -finite_min(cvfit$cve)
          if (is.na(sc) && !is.null(cvfit$cv))  sc <- -finite_min(cvfit$cv)
          lmd_min <- if (!is.null(cvfit$lmd.min)) cvfit$lmd.min else
            if (!is.null(cvfit$lambda.min)) cvfit$lambda.min else NA_real_
          if (is.na(sc) || is.na(lmd_min)) NULL else list(
            score = sc,
            fit   = cvfit,
            meta  = list(nlmd = nlmd, lmd.min.ratio = lmr, lmd.min = lmd_min)
          )
        }, timeout = per_call, onTimeout = "silent"), silent = TRUE)
        t2 <- Sys.time()
        
        if (!inherits(res, "try-error") && !is.null(res)) {
          tried <- tried + 1L
          timings <- c(timings, as.numeric(difftime(t2, t1, units="secs")))
          if (res$score > best$score) best <- res
        }
      }
      if (as.numeric(difftime(Sys.time(), start, units="secs")) >= time_budget_sec) break
    }
    if (as.numeric(difftime(Sys.time(), start, units="secs")) >= time_budget_sec) break
    ## refine
    nlmds <- unique(pmin(5000L, c(nlmds, nlmds*2L)))
    lmrs  <- sort(unique(c(lmrs, lmrs/2)))
  }
  
  list(
    best_fit = best$fit, best_meta = best$meta, best_score = best$score,
    models_tried = tried,
    mean_time_per_try = if (length(timings)) mean(timings) else NA_real_,
    total_time = as.numeric(difftime(Sys.time(), start, units="secs")),
    method = "Exclusive Lasso"
  )
}

## ----------------------------------------
## 5) Group Lasso tuner (grpreg) — cycle + refine
## ----------------------------------------
tune_gl_until <- function(time_budget_sec) {
  start <- Sys.time()
  best  <- list(score = -Inf, fit = NULL, meta = NULL)
  tried <- 0L; timings <- c()
  
  ## Anchor a lambda range once (via glmnet) and reuse for grpreg paths
  lr <- lambda_range_from_glmnet(X_train, y_train)
  lam_max <- lr[1]; lam_min <- lr[2]
  L_seq <- c(50L, 100L, 200L, 400L)
  
  repeat {
    for (L in L_seq) {
      elapsed <- as.numeric(difftime(Sys.time(), start, units="secs"))
      if (elapsed >= time_budget_sec) break
      per_call <- max(5, time_budget_sec - elapsed - 1)
      
      lambda_path <- exp_lambda_path_from_range(lam_max, lam_min, L)
      
      t1 <- Sys.time()
      res <- try(withTimeout({
        cvfit <- grpreg::cv.grpsurv(
          X_train, y_train, group = group_structure,
          penalty = "grLasso",
          alpha = 1,
          lambda = lambda_path,
          fold = foldid
        )
        sc <- -finite_min(cvfit$cve)  # lower cve is better, so flip sign
        if (is.na(sc)) NULL else list(
          score = sc,
          fit   = cvfit,
          meta  = list(nlambda = L, lambda.min = cvfit$lambda.min)
        )
      }, timeout = per_call, onTimeout = "silent"), silent = TRUE)
      t2 <- Sys.time()
      
      if (!inherits(res, "try-error") && !is.null(res)) {
        tried <- tried + 1L
        timings <- c(timings, as.numeric(difftime(t2, t1, units="secs")))
        if (res$score > best$score) best <- res
      }
    }
    if (as.numeric(difftime(Sys.time(), start, units="secs")) >= time_budget_sec) break
    ## refine
    L_seq <- unique(pmin(5000L, c(L_seq, L_seq*2L)))
  }
  
  list(
    best_fit = best$fit, best_meta = best$meta, best_score = best$score,
    models_tried = tried,
    mean_time_per_try = if (length(timings)) mean(timings) else NA_real_,
    total_time = as.numeric(difftime(Sys.time(), start, units="secs")),
    method = "Group Lasso"
  )
}

## ----------------------------------------
## 6) IPF-Lasso tuner (ipflasso) — cycle + refine
## ----------------------------------------
# tune_ipf_until <- function(time_budget_sec) {
#   start <- Sys.time()
#   best  <- list(score = -Inf, fit = NULL, meta = NULL)
#   tried <- 0L; timings <- c()
#   
#   ## blocks
#   ends <- cumsum(group_sizes); starts <- c(1, head(ends, -1) + 1)
#   blocks <- mapply(function(s, e) s:e, starts, ends, SIMPLIFY = FALSE)
#   
#   ## priority factors via ridge pre-weighting (fixed folds)
#   cv_ridge <- glmnet::cv.glmnet(X_train, y_train, alpha = 0, family = "cox", foldid = foldid)
#   coefs_r  <- abs(as.numeric(coef(cv_ridge, s = "lambda.min")))
#   split_i  <- rep(seq_along(group_sizes), times = group_sizes)
#   grp_means <- sapply(split(coefs_r, split_i), function(v) mean(v, na.rm = TRUE))
#   pflist <- c(1,grp_means[1]/grp_means[2],
#                         grp_means[1]/grp_means[3],
#                         grp_means[1]/grp_means[4],
#                         grp_means[1]/grp_means[5])
#   
#   ## anchor lambda range via glmnet once (works fine for ipflasso too)
#   lr <- lambda_range_from_glmnet(X_train, y_train)
#   lam_max <- lr[1]; lam_min <- lr[2]
#   L_seq <- c(30L, 60L, 120L, 240L)
#   
#   repeat {
#     for (L in L_seq) {
#       elapsed <- as.numeric(difftime(Sys.time(), start, units="secs"))
#       if (elapsed >= time_budget_sec) break
#       per_call <- max(5, time_budget_sec - elapsed - 1)
#       
#       lambda_path <- exp_lambda_path_from_range(lam_max, lam_min, L)
#       
#       t1 <- Sys.time()
#       res <- try(withTimeout({
#         cvfit <- ipflasso::cvr.ipflasso(
#           X = X_train, Y = y_train, family = "cox",
#           blocks = blocks, pf = pflist, nfolds = 5, ncv = 5,
#           lambda = lambda_path,
#           type.measure = "partial likelihood"
#         )
#         sc <- NA_real_
#         if (!is.null(cvfit$cv)) sc <- -finite_min(cvfit$cv)  # lower is better
#         if (is.na(sc) && !is.null(cvfit$cve)) sc <- -finite_min(cvfit$cve)
#         ind <- cvfit$ind.bestlambda
#         if (is.na(sc) || is.null(ind)) NULL else list(
#           score = sc,
#           fit   = cvfit,
#           meta  = list(nlambda = L, ind = ind)
#         )
#       }, timeout = per_call, onTimeout = "silent"), silent = TRUE)
#       t2 <- Sys.time()
#       
#       if (!inherits(res, "try-error") && !is.null(res)) {
#         tried <- tried + 1L
#         timings <- c(timings, as.numeric(difftime(t2, t1, units="secs")))
#         if (res$score > best$score) best <- res
#       }
#     }
#     if (as.numeric(difftime(Sys.time(), start, units="secs")) >= time_budget_sec) break
#     ## refine
#     L_seq <- unique(pmin(5000L, c(L_seq, L_seq*2L)))
#   }
#   
#   list(
#     best_fit = best$fit, best_meta = best$meta, best_score = best$score,
#     models_tried = tried,
#     mean_time_per_try = if (length(timings)) mean(timings) else NA_real_,
#     total_time = as.numeric(difftime(Sys.time(), start, units="secs")),
#     method = "IPF-Lasso"
#   )
# }




tune_ipf_until <- function(time_budget_sec) {
  start <- Sys.time()
  best  <- list(score = -Inf, fit = NULL, meta = NULL)
  tried <- 0L; timings <- c()
  
  # blocks
  ends <- cumsum(group_sizes); starts <- c(1, head(ends, -1) + 1)
  blocks <- mapply(function(s, e) s:e, starts, ends, SIMPLIFY = FALSE)
  
  # pf via ridge Cox (fixed folds for fairness)
  cv_ridge <- glmnet::cv.glmnet(X_train, y_train, alpha = 0, family = "cox", foldid = foldid)
  coefs_r  <- abs(as.numeric(coef(cv_ridge, s = "lambda.min")))
  split_i  <- rep(seq_along(group_sizes), times = group_sizes)
  grp_means <- sapply(split(coefs_r, split_i), function(v) mean(v, na.rm = TRUE))
  
  # robust pf
  eps <- 1e-6
  grp_means <- pmax(grp_means, eps)
  pflist <- c(1, grp_means[1] / grp_means[-1])
  pflist <- pflist / mean(pflist)
  
  # lighter CV first
  nfolds_ipf <- 3; ncv_ipf <- 3
  
  repeat {
    # ---- seed once with ipflasso-chosen path ----
    if (tried == 0L) {
      elapsed <- as.numeric(difftime(Sys.time(), start, units="secs"))
      if (elapsed >= time_budget_sec) break
      per_call <- max(60, time_budget_sec - elapsed - 5)
      
      t1 <- Sys.time()
      res0 <- try(withTimeout({
        fit0 <- ipflasso::cvr.ipflasso(
          X = X_train, Y = y_train, family = "cox",
          blocks = blocks, pf = pflist, nfolds = nfolds_ipf, ncv = ncv_ipf,
          type.measure = "partial likelihood"
        )
        sc0 <- if (!is.null(fit0$cv)) -finite_min(fit0$cv) else if (!is.null(fit0$cve)) -finite_min(fit0$cve) else NA_real_
        if (is.finite(sc0) && !is.null(fit0$ind.bestlambda))
          list(score=sc0, fit=fit0, meta=list(nlambda=length(fit0$lambda), ind=fit0$ind.bestlambda))
        else NULL
      }, timeout = per_call, onTimeout = "silent"), silent = TRUE)
      t2 <- Sys.time()
      
      if (!inherits(res0, "try-error") && !is.null(res0)) {
        tried <- tried + 1L
        timings <- c(timings, as.numeric(difftime(t2, t1, units="secs")))
        if (res0$score > best$score) best <- res0
      }
    }
    
    # ---- refinement with custom paths ----
    for (L in L_seq) {
      elapsed <- as.numeric(difftime(Sys.time(), start, units="secs"))
      if (elapsed >= time_budget_sec) break
      per_call <- max(60, time_budget_sec - elapsed - 5)
      
      lambda_path <- exp_lambda_path_from_range(lam_max, lam_min, L)
      t1 <- Sys.time()
      res <- try(withTimeout({
        fit <- ipflasso::cvr.ipflasso(
          X = X_train, Y = y_train, family = "cox",
          blocks = blocks, pf = pflist, nfolds = nfolds_ipf, ncv = ncv_ipf,
          lambda = lambda_path,
          type.measure = "partial likelihood"
        )
        sc <- if (!is.null(fit$cv)) -finite_min(fit$cv) else if (!is.null(fit$cve)) -finite_min(fit$cve) else NA_real_
        if (is.finite(sc) && !is.null(fit$ind.bestlambda))
          list(score=sc, fit=fit, meta=list(nlambda=L, ind=fit$ind.bestlambda))
        else NULL
      }, timeout = per_call, onTimeout = "silent"), silent = TRUE)
      t2 <- Sys.time()
      
      if (!inherits(res, "try-error") && !is.null(res)) {
        tried <- tried + 1L
        timings <- c(timings, as.numeric(difftime(t2, t1, units="secs")))
        if (res$score > best$score) best <- res
      }
    }
    
    if (as.numeric(difftime(Sys.time(), start, units="secs")) >= time_budget_sec) break
    # refine
    L_seq <- unique(pmin(5000L, c(L_seq, L_seq*2L)))
    nfolds_ipf <- min(5, nfolds_ipf + 1)
    ncv_ipf    <- min(5, ncv_ipf + 1)
  }
  
  tune_ipf_until <- function(time_budget_sec) {
    start <- Sys.time()
    best  <- list(score = -Inf, fit = NULL, meta = NULL)
    tried <- 0L; timings <- c()
    
    # blocks
    ends <- cumsum(group_sizes); starts <- c(1, head(ends, -1) + 1)
    blocks <- mapply(function(s, e) s:e, starts, ends, SIMPLIFY = FALSE)
    
    # pf via ridge Cox (fixed folds for fairness)
    cv_ridge <- glmnet::cv.glmnet(X_train, y_train, alpha = 0, family = "cox", foldid = foldid)
    coefs_r  <- abs(as.numeric(coef(cv_ridge, s = "lambda.min")))
    split_i  <- rep(seq_along(group_sizes), times = group_sizes)
    grp_means <- sapply(split(coefs_r, split_i), function(v) mean(v, na.rm = TRUE))
    
    # robust pf
    eps <- 1e-6
    grp_means <- pmax(grp_means, eps)
    pflist <- c(1, grp_means[1] / grp_means[-1])
    pflist <- pflist / mean(pflist)
    
    # ---- MISSING INITIALIZATIONS (add these) ----
    # anchor lambda range for IPF; make it a bit wider than glmnet's
    lr <- lambda_range_from_glmnet(X_train, y_train)  # must exist
    lam_max <- as.numeric(lr[1]) * 5
    lam_min <- lam_max * 1e-4
    # initial lambda path lengths to try
    L_seq <- c(30L, 60L, 120L)   # will double later
    # ---------------------------------------------
    
    # lighter CV first
    nfolds_ipf <- 3; ncv_ipf <- 3
    
    repeat {
      # seed once with ipflasso-chosen path
      if (tried == 0L) {
        elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
        if (elapsed >= time_budget_sec) break
        per_call <- max(60, time_budget_sec - elapsed - 5)
        
        t1 <- Sys.time()
        res0 <- try(R.utils::withTimeout({
          fit0 <- ipflasso::cvr.ipflasso(
            X = X_train, Y = y_train, family = "cox",
            blocks = blocks, pf = pflist, nfolds = nfolds_ipf, ncv = ncv_ipf,
            type.measure = "partial likelihood"
          )
          sc0 <- if (!is.null(fit0$cv)) -finite_min(fit0$cv) else if (!is.null(fit0$cve)) -finite_min(fit0$cve) else NA_real_
          if (is.finite(sc0) && !is.null(fit0$ind.bestlambda))
            list(score = sc0, fit = fit0, meta = list(nlambda = length(fit0$lambda), ind = fit0$ind.bestlambda))
          else NULL
        }, timeout = per_call, onTimeout = "silent"), silent = TRUE)
        t2 <- Sys.time()
        
        if (!inherits(res0, "try-error") && !is.null(res0)) {
          tried <- tried + 1L
          timings <- c(timings, as.numeric(difftime(t2, t1, units = "secs")))
          if (res0$score > best$score) best <- res0
        }
      }
      
      # refinement with custom paths
      for (L in L_seq) {
        elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
        if (elapsed >= time_budget_sec) break
        per_call <- max(60, time_budget_sec - elapsed - 5)
        
        lambda_path <- exp_lambda_path_from_range(lam_max, lam_min, L)
        
        t1 <- Sys.time()
        res <- try(R.utils::withTimeout({
          fit <- ipflasso::cvr.ipflasso(
            X = X_train, Y = y_train, family = "cox",
            blocks = blocks, pf = pflist, nfolds = nfolds_ipf, ncv = ncv_ipf,
            lambda = lambda_path,
            type.measure = "partial likelihood"
          )
          sc <- if (!is.null(fit$cv)) -finite_min(fit$cv) else if (!is.null(fit$cve)) -finite_min(fit$cve) else NA_real_
          if (is.finite(sc) && !is.null(fit$ind.bestlambda))
            list(score = sc, fit = fit, meta = list(nlambda = L, ind = fit$ind.bestlambda))
          else NULL
        }, timeout = per_call, onTimeout = "silent"), silent = TRUE)
        t2 <- Sys.time()
        
        if (!inherits(res, "try-error") && !is.null(res)) {
          tried <- tried + 1L
          timings <- c(timings, as.numeric(difftime(t2, t1, units = "secs")))
          if (res$score > best$score) best <- res
        }
      }
      
      if (as.numeric(difftime(Sys.time(), start, units = "secs")) >= time_budget_sec) break
      # refine
      L_seq      <- unique(pmin(5000L, c(L_seq, L_seq * 2L)))
      nfolds_ipf <- min(5, nfolds_ipf + 1)
      ncv_ipf    <- min(5, ncv_ipf + 1)
    }
    
    list(
      best_fit = best$fit, best_meta = best$meta, best_score = best$score,
      models_tried = tried,
      mean_time_per_try = if (length(timings)) mean(timings) else NA_real_,
      total_time = as.numeric(difftime(Sys.time(), start, units = "secs")),
      method = "IPF-Lasso"
    )
  }
  

## ----------------------------------------
## 7) Run all tuners under the same budget
## ----------------------------------------
cat(sprintf("[%s] Tuning EN...\n", now()));  res_en  <- tune_en_until(time_budget_sec)
cat(sprintf("[%s] Tuning EL...\n", now()));  res_el  <- tune_el_until(time_budget_sec)
cat(sprintf("[%s] Tuning GL...\n", now()));  res_gl  <- tune_gl_until(time_budget_sec)
cat(sprintf("[%s] Tuning IPF...\n", now())); res_ipf <- tune_ipf_until(time_budget_sec)

## Fallbacks if a tuner couldn’t complete any valid grid (rare, but safe)
fallback_en <- function() {
  cv <- glmnet::cv.glmnet(X_train, y_train, family="cox", alpha=0.5, foldid=foldid, type.measure="C", nlambda=50)
  list(best_fit=cv, best_meta=list(alpha=0.5, nlambda=50, lambda.min=cv$lambda.min),
       best_score=finite_max(cv$cvm), models_tried=1L, mean_time_per_try=NA_real_, total_time=NA_real_,
       method="Elastic Net (fallback)")
}
fallback_el <- function() {
  set.seed(2025)
  cv <- ExclusiveLasso::cv.exclusive_lasso(
    X=X_train, y=y_train, groups=group_structure, family="cox", type.measure="cindex",
    intercept=FALSE, skip_df=TRUE, nfolds=10L, nlmd=100L, lmd.min.ratio=1e-3
  )
  sc <- if (!is.null(cv$cvm)) finite_max(cv$cvm) else if (!is.null(cv$cve)) -finite_min(cv$cve) else -finite_min(cv$cv)
  list(best_fit=cv, best_meta=list(nlmd=100L, lmd.min.ratio=1e-3, lmd.min=cv$lmd.min),
       best_score=sc, models_tried=1L, mean_time_per_try=NA_real_, total_time=NA_real_,
       method="Exclusive Lasso (fallback)")
}
fallback_gl <- function() {
  cv <- grpreg::cv.grpsurv(X_train, y_train, group=group_structure, penalty="grLasso", alpha=1, fold=foldid)
  list(best_fit=cv, best_meta=list(nlambda=length(cv$lambda), lambda.min=cv$lambda.min),
       best_score=-finite_min(cv$cve), models_tried=1L, mean_time_per_try=NA_real_, total_time=NA_real_,
       method="Group Lasso (fallback)")
}
fallback_ipf <- function() {
  ends <- cumsum(group_sizes); starts <- c(1, head(ends,-1)+1)
  blocks <- mapply(function(s,e) s:e, starts, ends, SIMPLIFY=FALSE)
  cv <- ipflasso::cvr.ipflasso(X=X_train, Y=y_train, family="cox",
                               blocks=blocks, pf=rep(1, length(blocks)), nfolds=5, ncv=5,
                               type.measure="partial likelihood")
  sc <- if (!is.null(cv$cv)) -finite_min(cv$cv) else -finite_min(cv$cve)
  list(best_fit=cv, best_meta=list(nlambda=length(cv$lambda), ind=cv$ind.bestlambda),
       best_score=sc, models_tried=1L, mean_time_per_try=NA_real_, total_time=NA_real_,
       method="IPF-Lasso (fallback)")
}

if (is.null(res_en$best_fit))  { message("EN fallback");  res_en  <- fallback_en() }
if (is.null(res_el$best_fit))  { message("EL fallback");  res_el  <- fallback_el() }
if (is.null(res_gl$best_fit))  { message("GL fallback");  res_gl  <- fallback_gl() }
if (is.null(res_ipf$best_fit)) { message("IPF fallback"); res_ipf <- fallback_ipf() }

## ----------------------------------------
## 8) Final refits / predictions using best configs
## ----------------------------------------
## EN
pred_en_train <- predict(res_en$best_fit, newx=X_train, s=res_en$best_meta$lambda.min, type="link")
pred_en_test  <- predict(res_en$best_fit, newx=X_test,  s=res_en$best_meta$lambda.min, type="link")
coef_en       <- coef(res_en$best_fit,   s=res_en$best_meta$lambda.min)

## EL (final non-CV fit requires the best lmd only)
fit_el <- ExclusiveLasso::exclusive_lasso(
  X=X_train, y=y_train, lmd=res_el$best_meta$lmd.min,
  groups=group_structure, family="cox", intercept=FALSE
)
pred_el_train <- predict(fit_el, newx=X_train, s=res_el$best_meta$lmd.min, type="link")
pred_el_test  <- predict(fit_el, newx=X_test,  s=res_el$best_meta$lmd.min, type="link")
coef_el       <- coef(fit_el)[-1]

## GL (refit at lambda.min)
fit_gl <- grpreg::grpsurv(
  X_train, y_train, group=group_structure, penalty="grLasso",
  alpha=1, lambda=res_gl$best_meta$lambda.min
) 
fit_gl1 <- grpreg::grpsurv(
  X_train, y_train, group=group_structure, penalty="grLasso",
  alpha=1
) 
pred_gl_train <- predict(fit_gl1, X=X_train, type="link", lambda=res_gl$best_meta$lambda.min)
pred_gl_test  <- predict(fit_gl1, X=X_test,  type="link", lambda=res_gl$best_meta$lambda.min)
coef_gl       <- coef(fit_gl)

## IPF
pred_ipf_train <- ipflasso::ipflasso.predict(object=res_ipf$best_fit, Xtest=X_train)$linpredtest
pred_ipf_test  <- ipflasso::ipflasso.predict(object=res_ipf$best_fit, Xtest=X_test)$linpredtest
coef_ipf       <- res_ipf$best_fit$coeff[-1, res_ipf$best_meta$ind]

## ----------------------------------------
## 9) Runtime summary (for the paper)
## ----------------------------------------
runtime_summary <- do.call(
  rbind,
  lapply(list(res_en, res_el, res_gl, res_ipf), function(z)
    data.frame(
      Method = z$method,
      Models_Tried = z$models_tried,
      Mean_Time_Per_Try_sec = round(z$mean_time_per_try, 3),
      Total_Time_sec = round(z$total_time, 2),
      Best_Tuning_Score = round(z$best_score, 6)
    )
  )
)
print(runtime_summary)

## At this point, compute C-index / iAUC / Brier on train/test from the
## four sets of predictions: pred_*_train / pred_*_test


library(survcomp)

Cindex <- function(pred, y) {
  survcomp::concordance.index(as.numeric(pred),
                              surv.time = y[,1],
                              surv.event = y[,2],
                              method = "noether")$c.index
}

test_cindex <- data.frame(
  Method = c("Elastic Net","Exclusive Lasso","Group Lasso","IPF-Lasso"),
  Cindex = c(
    Cindex(pred_en_test,  y_test),
    Cindex(pred_el_test,  y_test),
    Cindex(pred_gl_test,  y_test),
    Cindex(pred_ipf_test, y_test)
  )
)
print(test_cindex, row.names = FALSE)

# (Optional) how sparse are the fitted models?
nnz <- function(beta) sum(abs(as.numeric(beta)) > 0)
nz_summary <- data.frame(
  Method = c("Elastic Net","Exclusive Lasso","Group Lasso","IPF-Lasso"),
  Nonzero = c(
    nnz(coef_en),
    nnz(coef_el),
    nnz(coef_gl),
    nnz(coef_ipf)
  )
)
print(nz_summary, row.names = FALSE)



evaluate_selection <- function(true_coefs, selected) {
  # Ensure logical vectors
  true_coefs <- as.logical(true_coefs)
  selected <- as.logical(selected)

  # Confusion matrix components
  TP <- sum(true_coefs & selected)
  FP <- sum(selected & !true_coefs)
  FN <- sum(true_coefs & !selected)
  TN <- sum(!true_coefs & !selected)

  # Metrics
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
  false_negative_rate <- FN / (FN + TP)
  false_discovery_rate <- FP / (FP + TP)
  accuracy <- (TP + TN) / length(true_coefs)

  return(list(
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    false_negative_rate = false_negative_rate,
    false_discovery_rate = false_discovery_rate,
    accuracy = accuracy
  ))
}


results_elasticnet0 <- evaluate_selection(true_coefs, coef_en)
results_el0 <- evaluate_selection(true_coefs, coef_el)
results_gl0 <- evaluate_selection(true_coefs, coef_gl)
results_ipf0 <- evaluate_selection(true_coefs, coef_ipf)

results_elasticnet0$f1_score
results_el0$f1_score
results_gl0$f1_score
results_ipf0$f1_score


evaluate_survival_model <- function(
    X_train, X_test,
    y_train, y_test,
    pred_train, pred_test,
    training_data, testing_data,
    cvfit_model, lambda.min,
    eval_time = 5, tmax_auc = 10,
    model = "glmnet"
) {
  library(survival)
  library(risksetROC)
  library(riskRegression)
  library(survcomp)

  # Concordance Index
  Cindex_train <- survcomp::concordance.index(pred_train[, 1], surv.time = y_train[, 1], surv.event = y_train[, 2])$c.index
  Cindex_test  <- survcomp::concordance.index(pred_test[, 1], surv.time = y_test[, 1], surv.event = y_test[, 2])$c.index

  # AUC and iAUC
  AUC_train <- risksetAUC(Stime = y_train[, "time"], status = y_train[, "status"], marker = pred_train, tmax = max(y_train[, "time"]))
  AUC_test  <- risksetAUC(Stime = y_test[, "time"], status = y_test[, "status"], marker = pred_test, tmax = max(y_test[, "time"]))

  iAUC_train <- IntegrateAUC(AUC_train$AUC, AUC_train$utimes, AUC_train$St, tmax = tmax_auc)
  iAUC_test  <- IntegrateAUC(AUC_test$AUC, AUC_test$utimes, AUC_test$St, tmax = tmax_auc)

  # Linear predictor
  if (model == "ipf") {
    eta_train <- as.vector(ipflasso.predict(object = cvfit_model, Xtest = X_train)$linpredtest)
    eta_test  <- as.vector(ipflasso.predict(object = cvfit_model, Xtest = X_test)$linpredtest)
  } else if (model == "gl") {
    eta_train <- as.vector(predict(cvfit_model, X = X_train, type = "link", lambda = lambda.min))
    eta_test  <- as.vector(predict(cvfit_model, X = X_test, type = "link", lambda = lambda.min))
  } else {
    eta_train <- as.vector(predict(cvfit_model, newx = X_train, s = lambda.min, type = "link"))
    eta_test  <- as.vector(predict(cvfit_model, newx = X_test, s = lambda.min, type = "link"))
  }

  # Clean test data
  train_df <- training_data
  test_df  <- testing_data
  train_df$lp <- eta_train
  test_df$lp  <- eta_test

  keep_idx_test <- which(!is.na(eta_test))
  test_df <- test_df[keep_idx_test, ]
  eta_test <- eta_test[keep_idx_test]
  test_df$lp <- eta_test

  # Fit Cox models
  cox_model_train <- coxph(Surv(time, status) ~ lp, data = train_df, x = TRUE, y = TRUE)
  cox_model_test  <- coxph(Surv(time, status) ~ lp, data = test_df, x = TRUE, y = TRUE)

  # Integrated Brier Score (IBS)
  times_eval <- seq(365, 1825, by = 90)#c(365, 730, 1095, 1460, 1825)

  Brier_test <- Score(
    object = list("Model" = cox_model_test),
    formula = Surv(time, status) ~ 1,
    data = test_df,
    metrics = "brier",
    summary = "ibs",
    times = times_eval,
    conf.int = FALSE
  )

  Brier_train <- Score(
    object = list("Model" = cox_model_train),
    formula = Surv(time, status) ~ 1,
    data = train_df,
    metrics = "brier",
    summary = "ibs",
    times = times_eval,
    conf.int = FALSE
  )

  # Final Output: No plots, only numeric metrics
  list(
    Cindex = list(train = Cindex_train, test = Cindex_test),
    iAUC   = list(train = iAUC_train, test = iAUC_test),
    Brier  = list(
      train = Brier_train$Brier$score,
      test  = Brier_test$Brier$score
    )
  )
}

results_survival_elasticnet <- evaluate_full_model(X_train, X_test,
                    y_train, y_test,
                    pred_en_train, pred_en_test,
                    train_data, test_data,
                    cvfit_lasso, cvfit_lasso$lambda.min,
                    eval_time = 5, tmax_auc = 10,
                    model = "glmnet")














evaluate_survival_model <- function(
    X_train, X_test,
    y_train, y_test,
    training_data, testing_data,
    fit_obj,                  # cv.glmnet OR exclusive_lasso fit OR grpsurv fit OR cvr.ipflasso
    lambda_sel = NULL,        # s for glmnet, lambda for grpsurv, lmd for EL; ignored for IPF
    method = c("glmnet","el","gl","ipf"),
    tmax_auc = NULL,          # if NULL uses 90th percentile of train time
    times_eval =  seq(1, max(testing_data$time), by = 10)#seq(365, 1825, by = 90)# NULL 
    # if NULL picks 20 quantiles between 10% and 90% of test time
) {
  method <- match.arg(method)
  stopifnot(is.matrix(X_train), is.matrix(X_test))
  # 1) predictions (linear predictors)
  if (method == "ipf") {
    eta_train <- as.numeric(ipflasso::ipflasso.predict(object = fit_obj, Xtest = X_train)$linpredtest)
    eta_test  <- as.numeric(ipflasso::ipflasso.predict(object = fit_obj, Xtest = X_test)$linpredtest)
  } else if (method == "gl") {
    if (is.null(lambda_sel)) stop("Provide lambda_sel = lambda.min for grpsurv.")
    eta_train <- as.numeric(predict(fit_obj, X = X_train, type = "link", lambda = lambda_sel))
    eta_test  <- as.numeric(predict(fit_obj, X = X_test,  type = "link", lambda = lambda_sel))
  } else if (method == "el") {
    if (is.null(lambda_sel)) stop("Provide lambda_sel = best lmd for Exclusive Lasso.")
    eta_train <- as.numeric(predict(fit_obj, newx = X_train, s = lambda_sel, type = "link"))
    eta_test  <- as.numeric(predict(fit_obj,  newx = X_test,  s = lambda_sel, type = "link"))
  } else { # glmnet
    if (is.null(lambda_sel)) stop("Provide lambda_sel = lambda.min for glmnet.")
    eta_train <- as.numeric(predict(fit_obj, newx = X_train, s = lambda_sel, type = "link"))
    eta_test  <- as.numeric(predict(fit_obj,  newx = X_test,  s = lambda_sel, type = "link"))
  }
  
  # 2) package imports
  library(survival)
  library(risksetROC)
  library(riskRegression)
  library(survcomp)
  
  # 3) Concordance Index
  Cindex_train <- survcomp::concordance.index(eta_train, surv.time = y_train[, 1], surv.event = y_train[, 2])$c.index
  Cindex_test  <- survcomp::concordance.index(eta_test,  surv.time = y_test[, 1],  surv.event = y_test[, 2])$c.index
  
  # 4) AUC and iAUC
  if (is.null(tmax_auc)) tmax_auc <- as.numeric(quantile(y_train[, "time"], 0.9, na.rm = TRUE))
  AUC_train <- risksetAUC(Stime = y_train[, "time"], status = y_train[, "status"], marker = eta_train, tmax = tmax_auc)
  AUC_test  <- risksetAUC(Stime = y_test[,  "time"], status = y_test[,  "status"], marker = eta_test,  tmax = tmax_auc)
  
  iAUC_train <- IntegrateAUC(AUC_train$AUC, AUC_train$utimes, AUC_train$St, tmax = tmax_auc)
  iAUC_test  <- IntegrateAUC(AUC_test$AUC,  AUC_test$utimes,  AUC_test$St,  tmax = tmax_auc)
  
  # 5) Build data frames with the linear predictors
  train_df <- training_data
  test_df  <- testing_data
  train_df$lp <- eta_train
  test_df$lp  <- eta_test
  
  # 6) Cox models for Brier/IBS
  #    If your simulated times are not in days, do not use 365-based times.
  if (is.null(times_eval)) {
    q <- quantile(test_df$time, probs = seq(0.10, 0.90, length.out = 20), na.rm = TRUE)
    times_eval <- as.numeric(unique(sort(q)))
  }
  
  cox_model_train <- coxph(Surv(time, status) ~ lp, data = train_df, x = TRUE, y = TRUE)
  cox_model_test  <- coxph(Surv(time, status) ~ lp, data = test_df,  x = TRUE, y = TRUE)
  
  Brier_test <- Score(
    object  = list("Model" = cox_model_test),
    formula = Surv(time, status) ~ 1,
    data    = test_df,
    metrics = "brier",
    summary = "ibs",
    times   = times_eval,
    conf.int = FALSE
  )
  Brier_train <- Score(
    object  = list("Model" = cox_model_train),
    formula = Surv(time, status) ~ 1,
    data    = train_df,
    metrics = "brier",
    summary = "ibs",
    times   = times_eval,
    conf.int = FALSE
  )
  
  list(
    Cindex = list(train = Cindex_train, test = Cindex_test),
    iAUC   = list(train = iAUC_train,  test = iAUC_test),
    Brier  = list(
      train_curve = Brier_train$Brier$score,  # full curve if you want it
      test_curve  = Brier_test$Brier$score,
      IBS_train   = Brier_train$Brier$IBS,    # integrated summary
      IBS_test    = Brier_test$Brier$IBS
    )
  )
}


# EN
metrics_en <- evaluate_survival_model(
  X_train, X_test, y_train, y_test,                 # data
  training_data = train_data, testing_data = test_data,
  fit_obj = res_en$best_fit, lambda_sel = res_en$best_meta$lambda.min,
  method = "glmnet"
)

# EL (use the FINAL non-CV fit and best lmd)
metrics_el <- evaluate_survival_model(
  X_train, X_test, y_train, y_test,
  training_data = train_data, testing_data = test_data,
  fit_obj = fit_el, lambda_sel = res_el$best_meta$lmd.min,
  method = "el"
)

# GL (use the FINAL fit at lambda.min)
metrics_gl <- evaluate_survival_model(
  X_train, X_test, y_train, y_test,
  training_data = train_data, testing_data = test_data,
  fit_obj = fit_gl1, lambda_sel = res_gl$best_meta$lambda.min,
  method = "gl"
)

# IPF (use the CV object; lambda handled internally by ipflasso.predict)
metrics_ipf <- evaluate_survival_model(
  X_train, X_test, y_train, y_test,
  training_data = train_data, testing_data = test_data,
  fit_obj = res_ipf$best_fit, lambda_sel = NULL,
  method = "ipf"
)

# Save outputs for reproducibility
write.csv(runtime_summary,
          file = "results/time_compare/runtime_summary.csv", row.names = FALSE)
write.csv(test_cindex,
          file = "results/time_compare/cindex_summary.csv", row.names = FALSE)
write.csv(nz_summary,
          file = "results/time_compare/sparsity_summary.csv", row.names = FALSE)

saveRDS(list(
  EN = res_en, EL = res_el, GL = res_gl, IPF = res_ipf,
  metrics = list(EN = metrics_en, EL = metrics_el,
                 GL = metrics_gl, IPF = metrics_ipf)
), file = "results/time_compare/time_compare_results.rds")
