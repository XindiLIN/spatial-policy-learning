library(parallel)

source('simulation/simulation_nonparametric_funcs.R')

# Phase D: pick winner per (seed, threshold) from phase C's CV search, final refit on the
# full training set. Split into its own script/job because this final refit uses a
# full-sample kernel matrix per task -- much more memory-hungry than phase C's CV folds --
# and needs its own SLURM resource request (more memory, fewer concurrent workers) to avoid
# OOM kills. Run simulation_nonparametric_run_phaseABC.R first; it saves the two files this
# script reads below.

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
cat(sprintf("Using up to %d cores for phase D\n", n_cores))

thr_ctxs_path <- 'simulation/output/simulation_nonparametric_thr_ctxs.rds'
hp_results_path <- 'simulation/output/simulation_nonparametric_hp_results.rds'

if (!file.exists(thr_ctxs_path) || !file.exists(hp_results_path)) {
  stop(sprintf("Missing phase A-C output. Run simulation_nonparametric_run_phaseABC.R first (expected %s and %s).",
               thr_ctxs_path, hp_results_path))
}

thr_ctxs <- readRDS(thr_ctxs_path)
hp_results_all <- readRDS(hp_results_path)

## ---- Phase D: pick winner per (seed, threshold), final refit ----

phase_d_partial_dir <- 'simulation/output/phase_d_partial'
dir.create(phase_d_partial_dir, recursive = TRUE, showWarnings = FALSE)

phase_d_partial_path <- function(thr_idx) {
  file.path(phase_d_partial_dir, sprintf("thr%02d.rds", thr_idx))
}

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

  result <- list(data_test_lst_entry = ctx$data_test_lst_entry, iteration_df = iteration_df)

  # Incremental checkpoint: each of phase D's 30 tasks writes its own result to disk the
  # moment it finishes, rather than only after the whole mclapply call returns -- mirrors
  # phase C's hp_partial/ checkpointing, since phase D's full-sample refit is exactly the
  # step that OOM'd in the original run.
  saveRDS(result, file = phase_d_partial_path(thr_idx))

  result
}

# Skip any thr_idx that already has a checkpoint file from a previous run, and retry
# whatever's still missing in a loop -- each retry only tackles the shrinking set of missing
# tasks, re-checking phase_d_partial/ after every pass, so mc.cores naturally scales down
# with it and later retries face less concurrent memory pressure than the first attempt did.
max_phase_d_attempts <- 5
all_thr_idx <- seq_along(thr_ctxs)

cat("Phase D: final refit per (seed, threshold)\n")
for (attempt in seq_len(max_phase_d_attempts)) {
  already_done <- file.exists(phase_d_partial_path(all_thr_idx))
  todo_idx <- all_thr_idx[!already_done]

  if (length(todo_idx) == 0) {
    if (attempt == 1) {
      cat(sprintf("Phase D: all %d tasks already done, nothing to run\n", length(all_thr_idx)))
    }
    break
  }

  cat(sprintf("Phase D attempt %d/%d: %d of %d tasks remaining\n",
              attempt, max_phase_d_attempts, length(todo_idx), length(all_thr_idx)))

  parallel::mclapply(todo_idx, run_final_fit, mc.cores = min(length(todo_idx), n_cores))
}

# Reassemble the full task-list result set from whatever's on disk now, rather than relying
# on mclapply's return values directly, so tasks a previous run or earlier retry already
# finished aren't lost even if a later retry's mclapply call fails partway through.
final_results <- lapply(all_thr_idx, function(thr_idx) {
  path <- phase_d_partial_path(thr_idx)
  if (file.exists(path)) readRDS(path) else NULL
})

missing <- vapply(final_results, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
if (any(missing)) {
  cat(sprintf("WARNING: %d of %d phase D tasks are still missing after %d attempt(s): thr_idx = %s\n",
              sum(missing), length(all_thr_idx), max_phase_d_attempts, paste(all_thr_idx[missing], collapse = ", ")))
  cat("Re-run this script again to retry just the missing tasks (already-completed ones will be skipped).\n")
}

data_test_lst <- lapply(final_results[!missing], `[[`, "data_test_lst_entry")
estimated_policy_lst <- lapply(final_results[!missing], `[[`, "iteration_df")

dir.create('simulation/output', recursive = TRUE, showWarnings = FALSE)
saveRDS(estimated_policy_lst, file = 'simulation/output/simulation_nonparametric_estimated_policy.rds')
saveRDS(data_test_lst, file = 'simulation/output/simulation_nonparametric_data_test.rds')

if (any(missing)) {
  cat(sprintf("Phase D incomplete: %d of %d tasks missing. Re-run this script to fill in the rest.\n",
              sum(missing), length(all_thr_idx)))
} else {
  cat("Phase D complete.\n")
}
