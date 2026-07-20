library(parallel)

source('simulation/simulation_parametric_funcs.R')

# DGP: linear outcome_mean (parametric case), spatial_dependence = TRUE (default) so U is a
# spatially-correlated factor affecting the outcome, but spatial_confounding = FALSE (default)
# so U does not bias treatment assignment -- matches the original parametric_simulation.Rmd setup.

seed_values <- 2000:2009
threshold_quantiles <- c(0.4, 0.5, 0.6)

# Two sequential phases, each run as a single flat parallel step (no nested mclapply, to
# avoid spawning more processes than there are cores):
#   A. per-seed setup       (10 independent units) -- data gen, nuisance models, smoothers
#   B. per-(seed,threshold) (30 independent units) -- optimal/indirect/direct policies
# Unlike the nonparametric (RKHS) script, the parametric direct method fits a single fixed
# linear model with no hyperparameter grid search, so there's no third phase here.

# On a SLURM cluster (e.g. Yale's Bouchet), detectCores() reports the whole node's core
# count, not what the job was actually allocated -- read SLURM's own env var when present,
# and only fall back to detectCores() - 1 for local/laptop runs outside SLURM.
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
cat(sprintf("Using up to %d cores per phase\n", n_cores))

## ---- Phase A: per-seed setup ----

get_seed_context <- function(seed_value){

  cat('seed value: ', seed_value, '\n')

  data <- data_generation_general(d = 50,
                                  treatment_mean,
                                  treatment_fun,
                                  outcome_mean,
                                  outcome_fun,
                                  gp_params_U = c(0.25,0.2,0.5,0),
                                  seed = seed_value)

  set.seed(seed = seed_value)
  train_indices <- sample(nrow(data$X), nrow(data$X) * 0.7)
  data_full <- data

  data <- lapply(data_full, subset_element, indices = train_indices)
  data_test <- lapply(data_full, subset_element, indices = -train_indices)

  outcome_mean_nonsp_model <- SuperLearner(
    Y = data$Y,
    X = cbind(T = data$T, data$X),
    family = gaussian(),
    SL.library = c("SL.gam")
  )

  outcome_resid <- data$Y - outcome_mean_nonsp_model$Z[,1]
  resid_Gp <- GpGp::fit_model(y = outcome_resid, locs = data$grd, silent = TRUE, convtol = 1e-03)
  krige_values <- leave_one_out_kriging(locs = data$grd,
                                        y_obs = outcome_resid,
                                        gp_model = resid_Gp$covfun_name,
                                        gp_params = resid_Gp$covparms,
                                        order = c("coordinate"))

  krige_values_test <- GpGp::predictions(fit = resid_Gp, locs_pred = data_test$grd, X_pred = rep(1,nrow(data_test$grd)))

  treatment_mean_model <- SuperLearner(
    Y = data$T,
    X = data$X,
    family = gaussian(),
    SL.library = c("SL.gam")
  )

  treatment_sd = sd(data$T - treatment_mean_model$Z)
  gps_est <- dnorm(data$T, mean = treatment_mean_model$SL.predict, sd = treatment_sd)

  smoothers_lst <- get_smoothers(data = data, data_test = data_test, outcome_mean_nonsp_model = outcome_mean_nonsp_model)

  kernel_bw <- 1.06 * sd(data$T) * nrow(data$X)^(-1/5)
  resids <- data$Y - outcome_mean_nonsp_model$SL.predict - krige_values

  list(seed_value = seed_value,
       data = data, data_test = data_test, data_full = data_full, train_indices = train_indices,
       krige_values = krige_values, krige_values_test = krige_values_test,
       gps_est = gps_est,
       smoothers = smoothers_lst$smoothers, cumint_smoothers = smoothers_lst$cumint_smoothers,
       smoothers_test = smoothers_lst$smoothers_test, cumint_smoothers_test = smoothers_lst$cumint_smoothers_test,
       kernel_bw = kernel_bw, resids = resids)
}

cat("Phase A: per-seed setup\n")
seed_ctxs <- parallel::mclapply(seed_values, get_seed_context, mc.cores = min(length(seed_values), n_cores))
names(seed_ctxs) <- as.character(seed_values)

failed <- vapply(seed_ctxs, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) stop(sprintf("Phase A failed for seed(s): %s", paste(seed_values[failed], collapse = ", ")))

## ---- Phase B: per-(seed, threshold) -- optimal/indirect/direct policies ----

thr_grid <- expand.grid(seed_idx = seq_along(seed_values), threshold_quantile = threshold_quantiles)

run_one_threshold <- function(seed_idx, threshold_quantile){

  ctx <- seed_ctxs[[seed_idx]]
  cat(sprintf("seed=%d | threshold_quantile=%.2f: starting\n", ctx$seed_value, threshold_quantile))

  threshold_val <- as.numeric(quantile(ctx$data_full$Y, threshold_quantile))

  data_test_lst_entry <- data.frame(
    threshold_val = threshold_val, T = ctx$data_test$T, Y = ctx$data_test$Y,
    seed = ctx$seed_value, threshold_quantile = threshold_quantile,
    test_obs_id = rep(1:nrow(ctx$data_test$X)))

  # get optimal policy (using the true linear outcome_mean/outcome_mean_derivative)

  optimal_policy <- sapply(1:nrow(ctx$data_full$X), function(i) {
    res <- optim(
      par = 0,
      fn = function(t) {
        (outcome_mean(ctx$data_full$X[i, ], t, ctx$data_full$U[i]) - threshold_val)^2
      },
      gr = function(t) {
        2 * outcome_mean_derivative(ctx$data_full$X[i, ], t, ctx$data_full$U[i]) *
          (outcome_mean(ctx$data_full$X[i, ], t, ctx$data_full$U[i]) - threshold_val)
      },
      method = "L-BFGS-B"
    )
    res$par
  })

  # get indirect policy

  indirect_nonsp_policy_test <- get_indirect_policy(data = ctx$data_test, trt_bounds = c(-4,4), krige_values = NULL, clip_epsilon = 30, loss_type = "or",
                                                  threshold_val = threshold_val, smoothers = ctx$smoothers_test, cumint_smoothers = ctx$cumint_smoothers_test)

  indirect_policy_test <- get_indirect_policy(data = ctx$data_test, trt_bounds = c(-4,4), krige_values = ctx$krige_values_test, clip_epsilon = 30, loss_type = "or",
                                                       threshold_val = threshold_val, smoothers = ctx$smoothers_test, cumint_smoothers = ctx$cumint_smoothers_test)

  # get direct policy using the parametric (linear index) policy class, doubly-robust loss

  direct_result <- get_parametric_direct_policy(data = ctx$data, data_test = ctx$data_test,
                                                krige_values = ctx$krige_values, krige_values_test = ctx$krige_values_test,
                                                resids = ctx$resids, gps_est = ctx$gps_est,
                                                smoothers = ctx$smoothers, cumint_smoothers = ctx$cumint_smoothers,
                                                threshold_val = threshold_val, trt_bounds = c(-4,4),
                                                kernel_bw = ctx$kernel_bw, clip_epsilon = 1)

  direct_policy_test <- direct_result$direct_policy_test

  iteration_df <- bind_rows(
    data.frame(method = "direct", value = direct_policy_test),
    data.frame(method = "indirect", value = indirect_policy_test),
    data.frame(method = "indirect_nonsp", value = indirect_nonsp_policy_test),
    data.frame(method = "optimal", value = optimal_policy[-ctx$train_indices])
  ) %>%
    mutate(
      seed = ctx$seed_value,
      threshold_quantile = threshold_quantile,
      test_obs_id = rep(1:nrow(ctx$data_test$X), 4)
    )

  cat(sprintf("seed=%d | threshold_quantile=%.2f: finished\n", ctx$seed_value, threshold_quantile))

  cat('The Performance of Indirect Nonspatial Method: \n')
  cat('Bias: ', Metrics::bias(optimal_policy[-ctx$train_indices], indirect_nonsp_policy_test),'\n')
  cat('rMSE: ', Metrics::rmse(optimal_policy[-ctx$train_indices], indirect_nonsp_policy_test),'\n')
  cat('MDAE: ', Metrics::mdae(optimal_policy[-ctx$train_indices], indirect_nonsp_policy_test),'\n')

  cat('The Performance of Indirect Method: \n')
  cat('Bias: ', Metrics::bias(optimal_policy[-ctx$train_indices], indirect_policy_test),'\n')
  cat('rMSE: ', Metrics::rmse(optimal_policy[-ctx$train_indices], indirect_policy_test),'\n')
  cat('MDAE: ', Metrics::mdae(optimal_policy[-ctx$train_indices], indirect_policy_test),'\n')

  cat('The Performance of Direct Method: \n')
  cat('Bias: ', Metrics::bias(optimal_policy[-ctx$train_indices], direct_policy_test),'\n')
  cat('rMSE: ', Metrics::rmse(optimal_policy[-ctx$train_indices], direct_policy_test),'\n')
  cat('MDAE: ', Metrics::mdae(optimal_policy[-ctx$train_indices], direct_policy_test),'\n')

  list(data_test_lst_entry = data_test_lst_entry, iteration_df = iteration_df)
}

cat("Phase B: per-(seed, threshold) optimal/indirect/direct policies\n")
final_results <- parallel::mclapply(1:nrow(thr_grid), function(j){
  run_one_threshold(thr_grid$seed_idx[j], thr_grid$threshold_quantile[j])
}, mc.cores = min(nrow(thr_grid), n_cores))

failed <- vapply(final_results, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) stop(sprintf("Phase B failed for task(s): %s", paste(which(failed), collapse = ", ")))

data_test_lst <- lapply(final_results, `[[`, "data_test_lst_entry")
estimated_policy_lst <- lapply(final_results, `[[`, "iteration_df")

dir.create('simulation/output', showWarnings = FALSE)
saveRDS(estimated_policy_lst, file = 'simulation/output/simulation_parametric_estimated_policy.rds')
saveRDS(data_test_lst, file = 'simulation/output/simulation_parametric_data_test.rds')
