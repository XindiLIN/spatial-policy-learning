library(tidymodels)
library(tidyverse)
library(SuperLearner)
library(dplyr)
library(ggplot2)
library(GpGp)
library(Metrics)
library(yardstick)

source('simulation/data_generation.R')
source('functions/funcs_kriging.R')
source('functions/loss_functions.R')
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
                                              treatment_range = c(- 4, 4),
                                              treatment_step = 0.05,
                                              pred_args = list(onlySL = TRUE))
  cf_predictions <- isotomic_correction(treatment_grid = cf_predictions$treatment_values,
                                        pred_matrix_cf = cf_predictions$pred_matrix)
  smoothers <- make_prediction_smoothers(treatment_grid = cf_predictions$treatment_values,
                                         pred_matrix_cf = cf_predictions$pred_matrix,
                                         smooth_method = "smooth.spline",
                                         spar = NULL)


  cumint_smoothers <- make_cumint_smoothers(treatment_grid = cf_predictions$treatment_values,
                                            pred_matrix_cf = cf_predictions$pred_matrix,
                                            smooth_method = "smooth.spline",
                                            spar = NULL)

  ## for test data

  cf_predictions_test <- precompute_cf_predictions(data_obs = cbind(T = data_test$T, data_test$X),
                                                   fit_outcome = outcome_mean_nonsp_model,
                                                   treatment_range = c(-4, 4),
                                                   treatment_step = 0.02,
                                                   pred_args = list(onlySL = TRUE))

  cf_predictions_test <- isotomic_correction(treatment_grid = cf_predictions_test$treatment_values,
                                             pred_matrix_cf = cf_predictions_test$pred_matrix)

  smoothers_test <- make_prediction_smoothers(treatment_grid = cf_predictions_test$treatment_values,
                                              pred_matrix_cf = cf_predictions_test$pred_matrix,
                                              smooth_method = "smooth.spline",
                                              spar = NULL)

  cumint_smoothers_test <- make_cumint_smoothers(treatment_grid = cf_predictions_test$treatment_values,
                                                 pred_matrix_cf = cf_predictions_test$pred_matrix,
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
    tmp_loss <- Inf
    tmp_direct_policy <- NA

    for (initial_point in initial_value_set){

      direct_trt_result <- optim(par = initial_point,
                                 fn = compute_total_loss_smooth,
                                 gr = d_compute_total_loss_smooth,
                                 subject_idx = i,
                                 T = data$T,
                                 krige_adjust = krige_values,
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
                                 control = list(maxit = 300, trace = 0))


      if (direct_trt_result$value < tmp_loss) {
        tmp_loss <- direct_trt_result$value
        tmp_direct_policy <- direct_trt_result$par
      }
    }
    indirect_policy[i] <- tmp_direct_policy
  }

  return(indirect_policy)
}



# Direct method for the parametric (linear index) policy class: theta'X is optimized
# against the doubly-robust "Integral" loss (compute_total_loss_db_sum_smooth), which is
# well-behaved enough for direct L-BFGS-B optimization -- no DC-algorithm decomposition
# needed here, unlike the higher-dimensional RKHS policy class.
get_parametric_direct_policy <- function(data, data_test, krige_values, krige_values_test,
                                          resids, gps_est, smoothers, cumint_smoothers,
                                          threshold_val, trt_bounds = c(-4,4),
                                          kernel_bw, clip_epsilon = 1){

  design_matrix <- model.matrix(~ X1 + X2 + X3 + X4 + U, cbind(data$X, U = krige_values))
  design_matrix_test <- model.matrix(~ X1 + X2 + X3 + X4 + U, cbind(data_test$X, U = krige_values_test))

  fit <- optim(par = runif(ncol(design_matrix), -1, 1),
               fn = function(coefs){
                 assigned_trt <- design_matrix %*% coefs
                 compute_total_loss_db_sum_smooth(assigned_trt = assigned_trt,
                                                   subject_idx = 1:nrow(data$X),
                                                   T = data$T,
                                                   krige_adjust = krige_values,
                                                   outcome_resid = resids,
                                                   propensity_est = gps_est,
                                                   cumint_smoothers = cumint_smoothers,
                                                   trt_bounds = trt_bounds,
                                                   threshold_val = threshold_val,
                                                   kernel_bw = kernel_bw,
                                                   clip_epsilon = clip_epsilon,
                                                   surrogate_type = "Gaussian",
                                                   loss_type = "Integral")
               },
               gr = function(coefs){
                 assigned_trt <- design_matrix %*% coefs
                 as.numeric(t(design_matrix) %*%
                   d_compute_total_loss_db_smooth(assigned_trt,
                                                   subject_idx = 1:nrow(data$X),
                                                   T = data$T,
                                                   krige_adjust = krige_values,
                                                   outcome_resid = resids,
                                                   propensity_est = gps_est,
                                                   smoothers = smoothers,
                                                   trt_bounds = trt_bounds,
                                                   threshold_val = threshold_val,
                                                   kernel_bw = kernel_bw,
                                                   clip_epsilon = clip_epsilon,
                                                   surrogate_type = "Gaussian"))
               },
               method = "L-BFGS-B",
               control = list(maxit = 600))

  list(coefs = fit$par,
       direct_policy = as.numeric(design_matrix %*% fit$par),
       direct_policy_test = as.numeric(design_matrix_test %*% fit$par))
}
