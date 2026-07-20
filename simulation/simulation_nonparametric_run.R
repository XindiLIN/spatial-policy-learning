library(parallel)

source('simulation/simulation_nonparametric_funcs.R')

# DGP: nonlinear outcome_mean_nonlinear, spatial_dependence = TRUE and spatial_confounding =
# TRUE -- full spatial confounding, matching the paper's main setting.

seed_values <- 2000:2009
threshold_quantiles <- c(0.4, 0.5, 0.6)

# Hyperparameter grid: for each of m/lambda/kernel_bw, a too-small, a reasonable, and a
# too-high value. kernel_bw is expressed as a multiplier on the (seed-specific) Silverman
# baseline, since that baseline depends on sd(data$T) for that seed's training data.
m_lst <- c(-1, 0, 1)                    # gamma = 2^m * median pairwise distance; m=0 is the median-heuristic default
lambda_lst <- c(0, 0.25, 2)             # RKHS ridge penalty; 0 = unregularized, 2 = heavily shrunk
kernel_bw_mult_lst <- c(0.25, 1, 2)     # multiplier on the Silverman rule-of-thumb baseline
n_folds <- 5                            # k-fold CV for hyperparameter selection (not in-sample fit)

# This whole simulation is organized as four sequential phases, each run as a single flat
# parallel step across independent units of work -- rather than nesting mclapply calls
# inside mclapply calls, which would spawn far more processes than there are cores and
# cause thrashing instead of speedup. Phases run one after another, so each is free to use
# every available core:
#   A. per-seed setup       (10 independent units) -- data gen, nuisance models, smoothers, folds
#   B. per-(seed,threshold) (30 independent units) -- optimal/indirect policies
#   C. per-(seed,threshold,m) (90 independent units) -- k-fold CV hyperparameter search,
#      grouped by m because gamma (and hence the whole Gram matrix) only depends on m: each
#      task builds the Gram matrix once and reuses row/column slices of it across all
#      (lambda, kernel_bw) combinations and all folds
#   D. per-(seed,threshold) (30 independent units) -- pick winner by CV mcc, final refit on
#      the full training set
# On a laptop with ~8-10 cores, phase A/B/D's 10-30-way parallelism already saturates
# available cores; phase C's 90-way granularity is what lets this scale further on a cluster.

# On a SLURM cluster (e.g. Yale's Bouchet), detectCores() reports the whole node's core
# count, not what the job was actually allocated -- read SLURM's own env var when present,
# and only fall back to detectCores() - 1 for local/laptop runs outside SLURM.
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
cat(sprintf("Using up to %d cores per phase\n", n_cores))

## ---- Phase A: per-seed setup ----

get_seed_context <- function(seed_value){

  cat('seed value: ', seed_value, '\n')
  
  # data generation
  data <- data_generation_general(d = 50,
                                  treatment_mean,
                                  treatment_fun,
                                  outcome_mean_nonlinear,
                                  outcome_fun,
                                  gp_params_U = c(0.25,0.2,0.5,0),
                                  spatial_dependence = TRUE,
                                  spatial_confounding = TRUE,
                                  seed = seed_value)

  set.seed(seed = seed_value)
  
  # split data into training (70%) and test set (30%)
  train_indices <- sample(nrow(data$X), nrow(data$X) * 0.7)
  data_full <- data

  data <- lapply(data_full, subset_element, indices = train_indices)
  data_test <- lapply(data_full, subset_element, indices = -train_indices)
  
  # nuisance estimation
  ## fit non-spatial outcome regression. Step 1 in OR estimation in manuscript
  outcome_mean_nonsp_model <- SuperLearner(
    Y = data$Y,
    X = cbind(T = data$T, data$X),
    family = gaussian(),
    SL.library = c("SL.gam")
  )
  
  ## fit GP covariance matrix. Step 2 in OR estimation in manuscript
  outcome_resid <- data$Y - outcome_mean_nonsp_model$Z[,1]
  resid_Gp <- GpGp::fit_model(y = outcome_resid, locs = data$grd, silent = TRUE,convtol = 1e-03)
  
  # leave-one-out kriging. Step 3 in OR estimation in manuscript
  krige_values <- leave_one_out_kriging(locs = data$grd,
                                        y_obs = outcome_resid,
                                        gp_model = resid_Gp$covfun_name,
                                        gp_params = resid_Gp$covparms,
                                        order = c("coordinate"))
  
  ## kriging in test set
  krige_values_test <- GpGp::predictions(fit = resid_Gp, locs_pred = data_test$grd, X_pred = rep(1,nrow(data_test$grd)))
  
  ## propensity estimation 
  treatment_mean_model <- SuperLearner(
    Y = data$T,
    X = cbind(data$X, U = krige_values),
    family = gaussian(),
    SL.library = c("SL.lm")
  )

  treatment_sd = sd(data$T - treatment_mean_model$Z)
  ### generalized propensity in training set
  gps_est <- dnorm(data$T, mean = treatment_mean_model$SL.predict, sd = treatment_sd)
  
  
  # calculate the main component of loss function
  
  ## get smoothers of the non-spatial outcome regression model
  smoothers_lst <- get_smoothers(data = data, data_test = data_test, outcome_mean_nonsp_model = outcome_mean_nonsp_model)
  
  ## build the design matrix of the kernel svm
  kernel_design_matrix_lst <- get_kernel_design_matrix(data_full, data, data_test, krige_values, krige_values_test)

  ## fitted residuals of outcome regression as the numerator of the IPW term
  resids <- data$Y - outcome_mean_nonsp_model$SL.predict - krige_values

  ## baseline bandwidth of the surrogate kernel
  kernel_bw_baseline <- 1.06 * sd(data$T) * nrow(data$X)^(-1/5)
  
  ## fold id of each sample in training in cross-validation hyper-parameter selction
  fold_id <- assign_folds(nrow(data$X), n_folds = n_folds, seed_value = seed_value)

  list(seed_value = seed_value,
       data = data, data_test = data_test, data_full = data_full, train_indices = train_indices,
       krige_values = krige_values, krige_values_test = krige_values_test,
       gps_est = gps_est,
       smoothers = smoothers_lst$smoothers, cumint_smoothers = smoothers_lst$cumint_smoothers,
       smoothers_test = smoothers_lst$smoothers_test, cumint_smoothers_test = smoothers_lst$cumint_smoothers_test,
       kernel_design_matrix = kernel_design_matrix_lst$kernel_design_matrix,
       kernel_design_matrix_test = kernel_design_matrix_lst$kernel_design_matrix_test,
       resids = resids,
       kernel_bw_baseline = kernel_bw_baseline,
       fold_id = fold_id)
}

cat("Phase A: per-seed setup\n")
seed_ctxs <- parallel::mclapply(seed_values, get_seed_context, mc.cores = min(length(seed_values), n_cores))
names(seed_ctxs) <- as.character(seed_values)

failed <- vapply(seed_ctxs, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) stop(sprintf("Phase A failed for seed(s): %s", paste(seed_values[failed], collapse = ", ")))

## ---- Phase B: per-(seed, threshold) optimal/indirect policies ----

thr_grid <- expand.grid(seed_idx = seq_along(seed_values), threshold_quantile = threshold_quantiles)

get_threshold_context <- function(seed_idx, threshold_quantile){

  ctx <- seed_ctxs[[seed_idx]]
  cat(sprintf("seed=%d | threshold_quantile=%.2f: computing optimal/indirect policies\n", ctx$seed_value, threshold_quantile))

  threshold_val <- as.numeric(quantile(ctx$data_full$Y, threshold_quantile))

  data_test_lst_entry <- data.frame(
    threshold_val = threshold_val, T = ctx$data_test$T, Y = ctx$data_test$Y,
    seed = ctx$seed_value, threshold_quantile = threshold_quantile,
    test_obs_id = rep(1:nrow(ctx$data_test$X)))

  optimal_policy <- sapply(1:nrow(ctx$data_full$X), function(i) {
    res <- optim(
      par = 0,
      fn = function(t) {
        (outcome_mean_nonlinear(ctx$data_full$X[i, ], t, ctx$data_full$U[i]) - threshold_val)^2
      },
      gr = function(t) {
        2 * outcome_mean_nonlinear_derivative(ctx$data_full$X[i, ], t, ctx$data_full$U[i]) *
          (outcome_mean_nonlinear(ctx$data_full$X[i, ], t, ctx$data_full$U[i]) - threshold_val)
      },
      method = "L-BFGS-B"
    )
    res$par
  })

  indirect_nonsp_policy_test <- get_indirect_policy(data = ctx$data_test, trt_bounds = c(-4,4), krige_values = NULL, clip_epsilon = 30, loss_type = "or",
                                                  threshold_val = threshold_val, smoothers = ctx$smoothers_test, cumint_smoothers = ctx$cumint_smoothers_test)

  indirect_policy_test <- get_indirect_policy(data = ctx$data_test, trt_bounds = c(-4,4), krige_values = ctx$krige_values_test, clip_epsilon = 30, loss_type = "or",
                                                       threshold_val = threshold_val, smoothers = ctx$smoothers_test, cumint_smoothers = ctx$cumint_smoothers_test)

  c(ctx, list(threshold_quantile = threshold_quantile, threshold_val = threshold_val,
              data_test_lst_entry = data_test_lst_entry,
              optimal_policy = optimal_policy,
              indirect_nonsp_policy_test = indirect_nonsp_policy_test,
              indirect_policy_test = indirect_policy_test))
}

cat("Phase B: per-(seed, threshold) optimal/indirect policies\n")
thr_ctxs <- parallel::mclapply(1:nrow(thr_grid), function(j){
  get_threshold_context(thr_grid$seed_idx[j], thr_grid$threshold_quantile[j])
}, mc.cores = min(nrow(thr_grid), n_cores))

failed <- vapply(thr_ctxs, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) stop(sprintf("Phase B failed for task(s): %s", paste(which(failed), collapse = ", ")))

## ---- Phase C: per-(seed, threshold, m) k-fold CV hyperparameter search ----

cv_task_grid <- expand.grid(thr_idx = seq_len(nrow(thr_grid)), m = m_lst)

run_one_cv_group <- function(thr_idx, m){
  ctx <- thr_ctxs[[thr_idx]]
  label <- sprintf("seed=%d | threshold_quantile=%.2f | m=%.2f", ctx$seed_value, ctx$threshold_quantile, m)
  n_combos <- length(lambda_lst) * length(kernel_bw_mult_lst)
  cat(sprintf("%s: starting (%d combos x %d folds = %d fits)\n", label, n_combos, n_folds, n_combos * n_folds))

  result <- cv_hyperparam_group_for_m(m = m, lambda_lst = lambda_lst, kernel_bw_mult_lst = kernel_bw_mult_lst,
                                       data = ctx$data, kernel_design_matrix = ctx$kernel_design_matrix,
                                       kernel_bw_baseline = ctx$kernel_bw_baseline,
                                       fold_id = ctx$fold_id, n_folds = n_folds,
                                       threshold_val = ctx$threshold_val, trt_bounds = c(-4,4), clip_epsilon = 30,
                                       krige_values = ctx$krige_values, gps_est = ctx$gps_est, resids = ctx$resids,
                                       smoothers = ctx$smoothers, cumint_smoothers = ctx$cumint_smoothers,
                                       progress_label = label)
  result$thr_idx <- thr_idx

  cat(sprintf("%s: finished\n", label))
  result
}

cat(sprintf("Phase C: %d-fold CV hyperparameter search (%d seed x threshold x m tasks)\n", n_folds, nrow(cv_task_grid)))
hp_results <- parallel::mclapply(1:nrow(cv_task_grid), function(k){
  run_one_cv_group(cv_task_grid$thr_idx[k], cv_task_grid$m[k])
}, mc.cores = n_cores)

failed <- vapply(hp_results, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) stop(sprintf("Phase C failed for task(s): %s", paste(which(failed), collapse = ", ")))

hp_results_all <- bind_rows(hp_results)

# Attach seed/threshold_quantile (rather than just the internal thr_idx) so the saved
# hyperparameter-search results are self-describing, and mark the winning row per
# (seed, threshold) so Phase D and the saved output agree on exactly which combo was chosen.
hp_results_all$seed <- vapply(hp_results_all$thr_idx, function(i) thr_ctxs[[i]]$seed_value, numeric(1))
hp_results_all$threshold_quantile <- vapply(hp_results_all$thr_idx, function(i) thr_ctxs[[i]]$threshold_quantile, numeric(1))
hp_results_all$selected <- FALSE
for (idx in unique(hp_results_all$thr_idx)) {
  rows <- which(hp_results_all$thr_idx == idx)
  best <- rows[which.max(hp_results_all$mcc[rows])]
  hp_results_all$selected[best] <- TRUE
}

## ---- Phase D: pick winner per (seed, threshold), final refit ----

run_final_fit <- function(thr_idx){

  ctx <- thr_ctxs[[thr_idx]]
  winner <- hp_results_all[hp_results_all$thr_idx == thr_idx & hp_results_all$selected, ]

  m <- winner$m[1]
  lambda <- winner$lambda[1]
  kernel_bw <- winner$kernel_bw[1]

  cat(sprintf("seed=%d | threshold_quantile=%.2f: final refit with m=%.2f, lambda=%.2f, kernel_bw=%.3f\n",
              ctx$seed_value, ctx$threshold_quantile, m, lambda, kernel_bw))

  indirect_policy_DB <- get_indirect_policy(ctx$data, trt_bounds = c(-4,4), clip_epsilon = 30, ctx$threshold_val,
                                            krige_values = ctx$krige_values, gps_est = ctx$gps_est, resids = ctx$resids, kernel_bw = kernel_bw,
                                            smoothers = ctx$smoothers, cumint_smoothers = ctx$cumint_smoothers, loss_type = "db")

  gamma <- 2^(m) * median(fields::rdist(ctx$kernel_design_matrix))
  rbf <- rbfdot(sigma = 1/(gamma^2))

  k_matrix <- kernelMatrix(rbf, ctx$kernel_design_matrix)
  k_matrix <- cbind(rep(1,nrow(k_matrix)), k_matrix)

  k_matrix_test <- kernelMatrix(rbf, ctx$kernel_design_matrix_test, ctx$kernel_design_matrix)
  k_matrix_test <- cbind(rep(1,nrow(k_matrix_test)), k_matrix_test)

  initial_glmnet <- glmnet(x = k_matrix[,-1], y = indirect_policy_DB, alpha = 0, lambda = 0)
  params_initial_glmnet <- coef(initial_glmnet)
  params_initial_glmnet <- as.numeric(params_initial_glmnet)

  kernel_optim_trt_DB <- optim(par = params_initial_glmnet,
                               fn = compute_total_loss_smooth_RKHS,
                               gr = d_compute_total_loss_smooth_RKHS,
                               K = k_matrix,
                               T = ctx$data$T,
                               krige_adjust = ctx$krige_values,
                               outcome_resid = ctx$resids,
                               propensity_est = pmax(ctx$gps_est, 0.01),
                               lambda = lambda,
                               smoothers = ctx$smoothers,
                               cumint_smoothers = ctx$cumint_smoothers,
                               trt_bounds = c(-4,4),
                               threshold_val = ctx$threshold_val,
                               kernel_bw = kernel_bw,
                               clip_epsilon = 30,
                               surrogate_type = "Gaussian",
                               loss_type = "db",
                               method = "L-BFGS-B",
                               control = list(maxit = 600, trace = 0))

  direct_policy_test <- k_matrix_test %*% kernel_optim_trt_DB$par

  iteration_df <- bind_rows(
    data.frame(method = "direct", value = direct_policy_test),
    data.frame(method = "indirect", value = ctx$indirect_policy_test),
    data.frame(method = "indirect_nonsp", value = ctx$indirect_nonsp_policy_test),
    data.frame(method = "optimal", value = ctx$optimal_policy[-ctx$train_indices])
  ) %>%
    mutate(
      seed = ctx$seed_value,
      threshold_quantile = ctx$threshold_quantile,
      test_obs_id = rep(1:nrow(ctx$data_test$X), 4)
    )

  cat(sprintf("seed=%d | threshold_quantile=%.2f: finished\n", ctx$seed_value, ctx$threshold_quantile))
  cat('The Performance of Indirect Nonspatial Method: \n')
  cat('Bias: ', Metrics::bias(ctx$optimal_policy[-ctx$train_indices], ctx$indirect_nonsp_policy_test),'\n')
  cat('rMSE: ', Metrics::rmse(ctx$optimal_policy[-ctx$train_indices], ctx$indirect_nonsp_policy_test),'\n')
  cat('MDAE: ', Metrics::mdae(ctx$optimal_policy[-ctx$train_indices], ctx$indirect_nonsp_policy_test),'\n')

  cat('The Performance of Indirect Method: \n')
  cat('Bias: ', Metrics::bias(ctx$optimal_policy[-ctx$train_indices], ctx$indirect_policy_test),'\n')
  cat('rMSE: ', Metrics::rmse(ctx$optimal_policy[-ctx$train_indices], ctx$indirect_policy_test),'\n')
  cat('MDAE: ', Metrics::mdae(ctx$optimal_policy[-ctx$train_indices], ctx$indirect_policy_test),'\n')

  cat('The Performance of Direct Method: \n')
  cat('Bias: ', Metrics::bias(ctx$optimal_policy[-ctx$train_indices], direct_policy_test),'\n')
  cat('rMSE: ', Metrics::rmse(ctx$optimal_policy[-ctx$train_indices], direct_policy_test),'\n')
  cat('MDAE: ', Metrics::mdae(ctx$optimal_policy[-ctx$train_indices], direct_policy_test),'\n')

  list(data_test_lst_entry = ctx$data_test_lst_entry, iteration_df = iteration_df)
}

cat("Phase D: final refit per (seed, threshold)\n")
final_results <- parallel::mclapply(seq_len(nrow(thr_grid)), run_final_fit, mc.cores = min(nrow(thr_grid), n_cores))

failed <- vapply(final_results, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) stop(sprintf("Phase D failed for task(s): %s", paste(which(failed), collapse = ", ")))

data_test_lst <- lapply(final_results, `[[`, "data_test_lst_entry")
estimated_policy_lst <- lapply(final_results, `[[`, "iteration_df")

dir.create('simulation/output', showWarnings = FALSE)
saveRDS(estimated_policy_lst, file = 'simulation/output/simulation_nonparametric_estimated_policy.rds')
saveRDS(data_test_lst, file = 'simulation/output/simulation_nonparametric_data_test.rds')
saveRDS(hp_results_all, file = 'simulation/output/simulation_nonparametric_hp_results.rds')
