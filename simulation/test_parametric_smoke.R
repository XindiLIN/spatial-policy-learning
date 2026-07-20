# Smoke test for the parametric-policy-class simulation functions
# (simulation_parametric_funcs.R): confirms get_indirect_policy and
# get_parametric_direct_policy both converge and produce sensible, non-degenerate
# bias/rMSE/MDAE and classification metrics for a single seed/threshold_quantile,
# before committing to the full seed x threshold_quantile loop in
# simulation_parametric_run.R.
#
# Run with working directory set to the repo root.

source('simulation/simulation_parametric_funcs.R')

seed_value <- 2000
threshold_quantile <- 0.5

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
resid_Gp <- GpGp::fit_model(y = outcome_resid, locs = data$grd, silent = TRUE,convtol = 1e-03)
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
smoothers <- smoothers_lst$smoothers
cumint_smoothers <- smoothers_lst$cumint_smoothers
smoothers_test <- smoothers_lst$smoothers_test
cumint_smoothers_test <- smoothers_lst$cumint_smoothers_test

kernel_bw <- 1.06 * sd(data$T) * nrow(data$X)^(-1/5)
resids <- data$Y - outcome_mean_nonsp_model$SL.predict - krige_values

threshold_val <- as.numeric(quantile(data_full$Y, threshold_quantile))

cat('\n\n=== optimal policy ===\n\n')
optimal_policy <- sapply(1:nrow(data_full$X), function(i) {
  res <- optim(
    par = 0,
    fn = function(t) {
      (outcome_mean(data_full$X[i, ], t, data_full$U[i]) - threshold_val)^2
    },
    gr = function(t) {
      2 * outcome_mean_derivative(data_full$X[i, ], t, data_full$U[i]) *
        (outcome_mean(data_full$X[i, ], t, data_full$U[i]) - threshold_val)
    },
    method = "L-BFGS-B"
  )
  res$par
})

cat('\n\n=== indirect methods ===\n\n')
indirect_nonsp_policy_test <- get_indirect_policy(data = data_test, trt_bounds = c(-4,4), krige_values = NULL, clip_epsilon = 30, loss_type = "or",
                                                threshold_val = threshold_val, smoothers = smoothers_test, cumint_smoothers = cumint_smoothers_test)

indirect_policy_test <- get_indirect_policy(data = data_test, trt_bounds = c(-4,4), krige_values = krige_values_test, clip_epsilon = 30, loss_type = "or",
                                                     threshold_val = threshold_val, smoothers = smoothers_test, cumint_smoothers = cumint_smoothers_test)

cat('\n\n=== direct method (parametric DB) ===\n\n')
direct_result <- get_parametric_direct_policy(data = data, data_test = data_test,
                                              krige_values = krige_values, krige_values_test = krige_values_test,
                                              resids = resids, gps_est = gps_est,
                                              smoothers = smoothers, cumint_smoothers = cumint_smoothers,
                                              threshold_val = threshold_val, trt_bounds = c(-4,4),
                                              kernel_bw = kernel_bw, clip_epsilon = 1)

direct_policy_test <- direct_result$direct_policy_test

cat('\n\n=== RESULTS ===\n\n')
cat('coefs:\n'); print(direct_result$coefs)

cat('\nIndirect Nonspatial:\n')
cat('Bias: ', Metrics::bias(optimal_policy[-train_indices], indirect_nonsp_policy_test),'\n')
cat('rMSE: ', Metrics::rmse(optimal_policy[-train_indices], indirect_nonsp_policy_test),'\n')
cat('MDAE: ', Metrics::mdae(optimal_policy[-train_indices], indirect_nonsp_policy_test),'\n')

cat('\nIndirect (spatial):\n')
cat('Bias: ', Metrics::bias(optimal_policy[-train_indices], indirect_policy_test),'\n')
cat('rMSE: ', Metrics::rmse(optimal_policy[-train_indices], indirect_policy_test),'\n')
cat('MDAE: ', Metrics::mdae(optimal_policy[-train_indices], indirect_policy_test),'\n')

cat('\nDirect (parametric DB):\n')
cat('Bias: ', Metrics::bias(optimal_policy[-train_indices], direct_policy_test),'\n')
cat('rMSE: ', Metrics::rmse(optimal_policy[-train_indices], direct_policy_test),'\n')
cat('MDAE: ', Metrics::mdae(optimal_policy[-train_indices], direct_policy_test),'\n')

cat('\n\n=== Binary classification metrics (test) ===\n\n')
cat('Indirect nonsp:\n')
print(calculate_acc_mcc_two_sided_f1(T = data_test$T, Y = data_test$Y, policy = indirect_nonsp_policy_test, threshold_val = threshold_val))
cat('Indirect:\n')
print(calculate_acc_mcc_two_sided_f1(T = data_test$T, Y = data_test$Y, policy = indirect_policy_test, threshold_val = threshold_val))
cat('Direct:\n')
print(calculate_acc_mcc_two_sided_f1(T = data_test$T, Y = data_test$Y, policy = direct_policy_test, threshold_val = threshold_val))
