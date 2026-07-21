library(tidymodels)
library(tidyverse)
library(SuperLearner)
library(dplyr)
library(ggplot2)
library(GpGp)
library(Metrics)
library(yardstick)
library(glmnet)
library(kernlab)
library(parallel)



source('simulation/data_generation.R')
source('functions/funcs_kriging.R')
source('functions/loss_functions.R')
source('functions/loss_functions_RKHS.R')
source('functions/miscellaneous.r')

subset_element <- function(x, indices) {
  # Check if the element is a data frame or matrix
  if (is.data.frame(x) || is.matrix(x)) {
    return(x[indices, , drop = FALSE])
  } else {
    # Otherwise, treat it as a vector
    return(x[indices])
  }
}



get_smoothers <- function(data, data_test, outcome_mean_nonsp_model){

  ## for training data
  cf_predictions <- precompute_cf_predictions(data_obs = cbind(T = data$T, data$X),
                                              fit_outcome = outcome_mean_nonsp_model,
                                              # treatment_range = c(min(data$T) - 1, max(data$T) + 1),
                                              treatment_range = c(- 4, 4),
                                              treatment_step = 0.05,
                                              pred_args = list(onlySL = TRUE))
  cf_predictions <- isotomic_correction(treatment_grid = cf_predictions$treatment_values,
                                        pred_matrix_cf = cf_predictions$pred_matrix)
  smoothers <- make_prediction_smoothers(treatment_grid = cf_predictions$treatment_values,
                                         pred_matrix_cf = cf_predictions$pred_matrix,
                                         # smooth_method = "splinefun",
                                         smooth_method = "smooth.spline",
                                         spar = NULL)


  cumint_smoothers <- make_cumint_smoothers(treatment_grid = cf_predictions$treatment_values,
                                            pred_matrix_cf = cf_predictions$pred_matrix,
                                            # smooth_method = "splinefun",
                                            smooth_method = "smooth.spline",
                                            spar = NULL)

  ## for test data


  cf_predictions_test <- precompute_cf_predictions(data_obs = cbind(T = data_test$T, data_test$X),
                                                   fit_outcome = outcome_mean_nonsp_model,
                                                   treatment_range = c(-4, 4),
                                                   # treatment_range = c(min(data$T) - 1, max(data$T) + 1),
                                                   treatment_step = 0.02,
                                                   pred_args = list(onlySL = TRUE))

  cf_predictions_test <- isotomic_correction(treatment_grid = cf_predictions_test$treatment_values,
                                             pred_matrix_cf = cf_predictions_test$pred_matrix)

  smoothers_test <- make_prediction_smoothers(treatment_grid = cf_predictions_test$treatment_values,
                                              pred_matrix_cf = cf_predictions_test$pred_matrix,
                                              # smooth_method = "splinefun",
                                              smooth_method = "smooth.spline",
                                              spar = NULL)

  cumint_smoothers_test <- make_cumint_smoothers(treatment_grid = cf_predictions_test$treatment_values,
                                                 pred_matrix_cf = cf_predictions_test$pred_matrix,
                                                 # smooth_method = "splinefun",
                                                 smooth_method = "smooth.spline",
                                                 spar = NULL)
  return(list(smoothers = smoothers, cumint_smoothers = cumint_smoothers, smoothers_test = smoothers_test, cumint_smoothers_test = cumint_smoothers_test))

}



get_indirect_policy <- function(data, trt_bounds = c(-4,4), clip_epsilon = 30, threshold_val,
                                krige_values = NULL, gps_est = NULL, resids = NULL, kernel_bw = NULL,
                                smoothers, cumint_smoothers, loss_type = "or"){

  initial_value_set <-c(trt_bounds[1], mean(trt_bounds), trt_bounds[2])
  # for OR-based indirect, i.e., loss_type = "or" the kernel_bw, resids, gps_est does not make any effect, but we still need krige_values
  if(loss_type == "or"){
    kernel_bw <- 1.06 * sd(data$T) *nrow(data$X)^(-1/5)
    resids <- rep(0, nrow(data$X))
    gps_est <- rep(1, nrow(data$X))
  }

  if(is.null(krige_values)){
    krige_values <- rep(0, nrow(data$X))
  }

  indirect_policy <- rep(NA, nrow(data$X))

  for(i in 1:nrow(data$X)){
    # if(i%%100 == 0) print(i)
    # Use Inf for a more robust starting loss value
    tmp_loss <- Inf
    tmp_direct_policy <- NA # Default return value in case of errors

    # Inner loop to find the best initial point for the optimization
    for (initial_point in initial_value_set){

      # Use try() to catch potential errors during optimization

      direct_trt_result <- optim(par = initial_point,
                                 fn = compute_total_loss_smooth,
                                 gr = d_compute_total_loss_smooth,
                                 subject_idx = i,
                                 # Pass other necessary arguments using ...
                                 T = data$T,
                                 krige_adjust = krige_values,
                                 # krige_adjust = rep(0, nrow(data$X)), # setting kriging values equal to zero is the non-spatial case
                                 outcome_resid = resids,
                                 propensity_est = gps_est,
                                 smoothers =  smoothers,
                                 cumint_smoothers = cumint_smoothers,
                                 trt_bounds = trt_bounds,
                                 threshold_val = threshold_val,
                                 kernel_bw = kernel_bw,
                                 clip_epsilon = clip_epsilon ,
                                 surrogate_type = "Gaussian",
                                 loss_type = loss_type,
                                 method = "L-BFGS-B",
                                 # Set trace = 0 to prevent garbled console output from parallel workers
                                 control = list(maxit = 300, trace = 0))


      # Check if optim succeeded and if the new loss is an improvement
      if (direct_trt_result$value < tmp_loss) {
        tmp_loss <- direct_trt_result$value
        tmp_direct_policy <- direct_trt_result$par
      }
    }
    indirect_policy[i] <- tmp_direct_policy
  }

  return(indirect_policy)
}


get_kernel_design_matrix <- function(data_full, data, data_test, krige_values, krige_values_test){
  # 1. Define and prep the recipe on the TRAINING data
  rec <- recipe(~ ., data = data_full$X) %>%
    step_dummy(all_nominal_predictors()) %>%
    prep() # This "fits" the preprocessor

  # 2. Apply the *same* fitted recipe to both datasets, and add kriging values as additional
  kernel_design_matrix <- bake(rec, new_data = data$X)
  kernel_design_matrix <- cbind(kernel_design_matrix, U = krige_values)

  kernel_design_matrix_test <- bake(rec, new_data = data_test$X)
  kernel_design_matrix_test <- cbind(kernel_design_matrix_test, U = krige_values_test)


  # 3. standardize the matrix as mean equal to zero, sd equal to one
  train_means <- colMeans(kernel_design_matrix)
  train_sds <- apply(kernel_design_matrix, 2, sd)

  kernel_design_matrix <- scale(kernel_design_matrix)
  kernel_design_matrix_test <- scale(kernel_design_matrix_test, center = train_means, scale = train_sds)

  return(list(kernel_design_matrix = kernel_design_matrix, kernel_design_matrix_test = kernel_design_matrix_test))

}

# Assign k-fold membership for a training set of size n. Deterministic given seed_value, but
# uses a different RNG draw than the train/test split itself (offset seed) so fold membership
# isn't just a re-derivation of train_indices.
assign_folds <- function(n, n_folds = 5, seed_value){
  set.seed(seed_value + 100000)
  sample(rep(1:n_folds, length.out = n))
}

# Fit on one fold's training subset and evaluate on its held-out subset, given a Gram-matrix
# slice for each. K is expected to include the intercept column already (as produced by
# cv_hyperparam_group_for_m). This is the atomic per-fold fit.
fit_one_fold <- function(k_train, k_holdout, data_ft, T_holdout, Y_holdout,
                          krige_ft, gps_ft, resids_ft, smoothers_ft, cumint_smoothers_ft,
                          lambda, kernel_bw, threshold_val, trt_bounds = c(-4,4), clip_epsilon = 30,
                          maxit = 350){

  indirect_policy_DB_ft <- get_indirect_policy(data_ft, trt_bounds = trt_bounds, clip_epsilon = 30, threshold_val,
                                                krige_values = krige_ft, gps_est = gps_ft, resids = resids_ft, kernel_bw = kernel_bw,
                                                smoothers = smoothers_ft, cumint_smoothers = cumint_smoothers_ft, loss_type = "db")

  initial_glmnet <- glmnet(x = k_train[,-1], y = indirect_policy_DB_ft, alpha = 0, lambda = 0)
  params_initial_glmnet <- as.numeric(coef(initial_glmnet))

  fit <- optim(par = params_initial_glmnet,
               fn = compute_total_loss_smooth_RKHS,
               gr = d_compute_total_loss_smooth_RKHS,
               K = k_train,
               T = data_ft$T,
               krige_adjust = krige_ft,
               outcome_resid = resids_ft,
               propensity_est = pmax(gps_ft, 0.001),
               lambda = lambda,
               smoothers = smoothers_ft,
               cumint_smoothers = cumint_smoothers_ft,
               trt_bounds = trt_bounds,
               threshold_val = threshold_val,
               kernel_bw = kernel_bw,
               clip_epsilon = clip_epsilon,
               surrogate_type = "Gaussian",
               loss_type = "db",
               method = "L-BFGS-B",
               control = list(maxit = maxit, trace = 0))

  holdout_policy <- k_holdout %*% fit$par

  calculate_acc_mcc_two_sided_f1(T = T_holdout, Y = Y_holdout, policy = holdout_policy, threshold_val = threshold_val)
}

# k-fold cross-validate every (lambda, kernel_bw_mult) combination for ONE fixed m. Grouping
# by m matters: gamma (the RBF bandwidth) depends only on m, so the full n_train x n_train
# Gram matrix is built ONCE here and every fold's fitting matrix and held-out evaluation
# matrix are just row/column slices of it -- no need to recompute kernelMatrix() per fold or
# per (lambda, kernel_bw) combination.
cv_hyperparam_group_for_m <- function(m, lambda_lst, kernel_bw_mult_lst, data, kernel_design_matrix,
                                       kernel_bw_baseline, fold_id, n_folds = 5,
                                       threshold_val, trt_bounds = c(-4,4), clip_epsilon = 30,
                                       krige_values, gps_est, resids, smoothers, cumint_smoothers,
                                       progress_label = NULL, maxit = 350){

  gamma <- 2^(m) * median(fields::rdist(kernel_design_matrix))
  rbf <- rbfdot(sigma = 1/(gamma^2))
  k_raw <- kernelMatrix(rbf, kernel_design_matrix)  # n_train x n_train, no intercept column yet

  grid <- expand.grid(lambda = lambda_lst, kernel_bw_mult = kernel_bw_mult_lst)
  grid$acc <- NA; grid$mcc <- NA; grid$two_sided_f1 <- NA

  for(g in 1:nrow(grid)){

    lambda <- grid$lambda[g]
    kernel_bw <- grid$kernel_bw_mult[g] * kernel_bw_baseline

    fold_metrics <- vector("list", n_folds)

    for(f in 1:n_folds){
      fold_train_idx <- which(fold_id != f)
      fold_holdout_idx <- which(fold_id == f)

      k_train <- cbind(1, k_raw[fold_train_idx, fold_train_idx, drop = FALSE])
      k_holdout <- cbind(1, k_raw[fold_holdout_idx, fold_train_idx, drop = FALSE])

      data_ft <- lapply(data, subset_element, indices = fold_train_idx)

      fold_metrics[[f]] <- fit_one_fold(k_train = k_train, k_holdout = k_holdout, data_ft = data_ft,
                                        T_holdout = data$T[fold_holdout_idx], Y_holdout = data$Y[fold_holdout_idx],
                                        krige_ft = krige_values[fold_train_idx], gps_ft = gps_est[fold_train_idx],
                                        resids_ft = resids[fold_train_idx], smoothers_ft = smoothers[fold_train_idx],
                                        cumint_smoothers_ft = cumint_smoothers[fold_train_idx],
                                        lambda = lambda, kernel_bw = kernel_bw, threshold_val = threshold_val,
                                        trt_bounds = trt_bounds, clip_epsilon = clip_epsilon, maxit = maxit)
    }

    grid$acc[g] <- mean(vapply(fold_metrics, `[[`, numeric(1), "acc"))
    grid$mcc[g] <- mean(vapply(fold_metrics, `[[`, numeric(1), "mcc"))
    grid$two_sided_f1[g] <- mean(vapply(fold_metrics, `[[`, numeric(1), "two_sided_f1"))

    # Flag NA or suspiciously degenerate (exactly 0 or 1) metric values -- these usually
    # mean a fold ended up with a degenerate confusion matrix (e.g. every held-out
    # prediction landed in the same class), not a genuine, well-behaved CV result. Printed
    # via cat() rather than warning(), since warning() isn't reliably surfaced in the log
    # from inside a forked mclapply worker.
    metric_vals <- c(acc = grid$acc[g], mcc = grid$mcc[g], two_sided_f1 = grid$two_sided_f1[g])
    suspicious <- names(metric_vals)[is.na(metric_vals) | metric_vals %in% c(0, 1)]
    if (length(suspicious) > 0) {
      cat(sprintf("WARNING%s: suspicious metric value(s) for combo %d/%d (lambda=%.2f, kernel_bw_mult=%.2f): %s\n",
                  if (is.null(progress_label)) "" else paste0(" [", progress_label, "]"),
                  g, nrow(grid), lambda, grid$kernel_bw_mult[g],
                  paste(sprintf("%s=%s", suspicious, round(metric_vals[suspicious], 4)), collapse = ", ")))
    }

    if (!is.null(progress_label)) {
      cat(sprintf("%s: combo %d/%d done (lambda=%.2f, kernel_bw_mult=%.2f) -> acc=%.4f, mcc=%.4f, two_sided_f1=%.4f\n",
                  progress_label, g, nrow(grid), lambda, grid$kernel_bw_mult[g],
                  grid$acc[g], grid$mcc[g], grid$two_sided_f1[g]))
    }
  }

  grid$m <- m
  grid$kernel_bw <- grid$kernel_bw_mult * kernel_bw_baseline
  grid[, c("m","lambda","kernel_bw_mult","kernel_bw","acc","mcc","two_sided_f1")]
}

# Cross-validated hyperparameter search across the full (m, lambda, kernel_bw_mult) grid.
# Set mc.cores > 1 to run the m-groups in parallel (safe: each m-group is fully independent,
# since gamma -- and hence the Gram matrix -- is m-specific).
choose_hyperparameter_cv <- function(data, m_lst, lambda_lst, kernel_bw_mult_lst, kernel_bw_baseline,
                                      kernel_design_matrix, fold_id, n_folds = 5,
                                      threshold_val, trt_bounds = c(-4,4), clip_epsilon = 30,
                                      krige_values, gps_est, resids, smoothers, cumint_smoothers,
                                      mc.cores = 1){

  results_lst <- parallel::mclapply(m_lst, function(m){
    cv_hyperparam_group_for_m(m = m, lambda_lst = lambda_lst, kernel_bw_mult_lst = kernel_bw_mult_lst,
                               data = data, kernel_design_matrix = kernel_design_matrix, kernel_bw_baseline = kernel_bw_baseline,
                               fold_id = fold_id, n_folds = n_folds, threshold_val = threshold_val,
                               trt_bounds = trt_bounds, clip_epsilon = clip_epsilon,
                               krige_values = krige_values, gps_est = gps_est, resids = resids,
                               smoothers = smoothers, cumint_smoothers = cumint_smoothers)
  }, mc.cores = mc.cores)

  bind_rows(results_lst)
}
