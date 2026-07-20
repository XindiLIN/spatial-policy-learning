# Sanity check for the 5-fold CV hyperparameter search (choose_hyperparameter_cv /
# cv_hyperparam_group_for_m) in simulation_nonparametric_funcs.R.
# Run for a single seed/threshold_quantile to confirm the CV search converges and produces
# sensible, non-degenerate held-out acc/mcc/two_sided_f1 values before committing to the
# full seed x threshold_quantile loop in simulation_nonparametric_run.R.
#
# Run with working directory set to the repo root.

source('simulation/simulation_nonparametric_funcs.R')

seed_value <- 2000
threshold_quantile <- 0.5

cat('seed value: ', seed_value, '\n')

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
resid_Gp <- GpGp::fit_model(y = outcome_resid, locs = data$grd, silent = TRUE,convtol = 1e-03)
krige_values <- leave_one_out_kriging(locs = data$grd,
                                      y_obs = outcome_resid,
                                      gp_model = resid_Gp$covfun_name,
                                      gp_params = resid_Gp$covparms,
                                      order = c("coordinate"))

krige_values_test <- GpGp::predictions(fit = resid_Gp, locs_pred = data_test$grd, X_pred = rep(1,nrow(data_test$grd)))

treatment_mean_model <- SuperLearner(
  Y = data$T,
  X = cbind(data$X, U = krige_values),
  family = gaussian(),
  SL.library = c("SL.lm")
)

treatment_sd = sd(data$T - treatment_mean_model$Z)
gps_est <- dnorm(data$T, mean = treatment_mean_model$SL.predict, sd = treatment_sd)

smoothers_lst <- get_smoothers(data = data, data_test = data_test, outcome_mean_nonsp_model = outcome_mean_nonsp_model)
smoothers <- smoothers_lst$smoothers
cumint_smoothers <- smoothers_lst$cumint_smoothers

threshold_val <- as.numeric(quantile(data_full$Y, threshold_quantile))

kernel_design_matrix_lst <- get_kernel_design_matrix(data_full, data, data_test, krige_values, krige_values_test)

kernel_bw_baseline <- 1.06 * sd(data$T) * nrow(data$X)^(-1/5)
m_lst <- c(-1, 0, 1)
lambda_lst <- c(0, 0.25, 2)
kernel_bw_mult_lst <- c(0.25, 1, 2)
n_folds <- 5
clip_epsilon <- 30

resids <- data$Y -  outcome_mean_nonsp_model$SL.predict - krige_values

fold_id <- assign_folds(nrow(data$X), n_folds = n_folds, seed_value = seed_value)

cat('\n\n=== calling choose_hyperparameter_cv with 5-fold CV over 27-combo grid ===\n\n')

t0 <- Sys.time()
params_grid <- choose_hyperparameter_cv(data = data,
                                    m_lst = m_lst, lambda_lst = lambda_lst, kernel_bw_mult_lst = kernel_bw_mult_lst,
                                    kernel_bw_baseline = kernel_bw_baseline,
                                    kernel_design_matrix = kernel_design_matrix_lst$kernel_design_matrix,
                                    fold_id = fold_id, n_folds = n_folds,
                                    clip_epsilon = 30,
                                    threshold_val = threshold_val,
                                    trt_bounds = c(-4,4),
                                    krige_values = krige_values,
                                    gps_est = gps_est,
                                    resids = resids,
                                    smoothers = smoothers,
                                    cumint_smoothers = cumint_smoothers)
t1 <- Sys.time()

cat('\n\n=== RESULT params_grid (CV-averaged over', n_folds, 'folds) ===\n\n')
print(params_grid[order(-params_grid$mcc), ])

cat('\nBest row by CV mcc:\n')
print(params_grid[which.max(params_grid$mcc), ])

cat('\nElapsed time for 27-combination x', n_folds, '-fold CV grid (1 seed, 1 threshold):', as.numeric(difftime(t1, t0, units = "secs")), "seconds\n")

dir.create('simulation/output', showWarnings = FALSE)
saveRDS(params_grid, 'simulation/output/test_nonparametric_hyperparam_grid.rds')
