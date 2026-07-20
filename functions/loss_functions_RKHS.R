

# x is only within the treatment range
gaussian_kernel_plus <- function(x, mu = 0, sigma = 1, a = 0.4289) {
  (1 / (sigma * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sigma)^2) +  (1/(sqrt(pi) * sigma ** 2)) * a * x
}

# x is only within the treatment range
gaussian_kernel_minus <- function(x, mu = 0, sigma = 1, a = 0.4289) {
  (1/(sqrt(pi) * sigma ** 2)) * a * x
}


integral_gaussian_kernel_plus <- function(x, mu = 0, sigma = 1, a = 0.4289, trt_bounds) {
  pnorm(x, mean = mu, sd = sigma) - pnorm(trt_bounds[1], mean = mu, sd = sigma) + 
    (1/(2 * sqrt(pi) * sigma ** 2)) * a * x ** 2 - (1/(2 * sqrt(pi) * sigma ** 2)) * a * trt_bounds[1] ** 2
}


integral_gaussian_kernel_minus <- function(x, mu = 0, sigma = 1, a = 0.4289, trt_bounds) {
  (1/(2 * sqrt(pi) * sigma ** 2)) * a * x ** 2 - (1/(2 * sqrt(pi) * sigma ** 2)) * a * trt_bounds[1] ** 2
}


# we would re-write this function

# compute_total_loss_smooth_plus <- function(assigned_trt,
#                                           subject_idx,
#                                           T,
#                                           krige_adjust,
#                                           outcome_resid,
#                                           propensity_est,
#                                           cumint_smoothers,
#                                           trt_bounds,
#                                           threshold_val,
#                                           kernel_bw,
#                                           clip_epsilon,
#                                           clip_epsilon_bar,
#                                           surrogate_type = c("Gaussian", "Box"),
#                                           loss_type = c("db","ipw","or")) {
#   
#   surrogate_type <- match.arg(surrogate_type)
#   loss_type <- match.arg(loss_type)
#   n <- length(subject_idx)
#   loss_vec <- rep(NA, n)
#   
#   # Identify clipped observations
#   index_left_clip  <- assigned_trt < trt_bounds[1]
#   index_right_clip <- assigned_trt > trt_bounds[2]
#   index_no_clip    <- !index_left_clip & !index_right_clip
#   
#   # Right-clipped treatment values used in integration
#   trt_right_clipped <- pmin(assigned_trt, trt_bounds[2])
#   trt_clipped <- pmax(trt_right_clipped, trt_bounds[1])
#   
#   
#   
#   # ---------------------- uncapped and right-clipped loss ----------------------
#   
#   # --- Term 1: Threshold integral ---
#   term_threshold <- threshold_val * (trt_clipped - trt_bounds[1])
#   
#   # --- Term 2: Doubly-robust kernel term or IPW term---
#   if (loss_type == "db" | loss_type == "ipw"){
#     
#     outcome_resid_scaled <- (outcome_resid[subject_idx] / propensity_est[subject_idx])
#     outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
#     outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
#     
#     
#     if(surrogate_type == "Box"){
#       kernel_integral_plus <- kernel_integral_plus(trt = trt_clipped,
#                                                    trt_obs = T[subject_idx],
#                                                    bandwidth = kernel_bw)
#       kernel_integral_minus <- kernel_integral_minus(trt = trt_clipped,
#                                                      trt_obs = T[subject_idx],
#                                                      bandwidth = kernel_bw)
#     } else if (surrogate_type == "Gaussian"){
#       kernel_integral_plus <- integral_gaussian_kernel_plus(trt_clipped,
#                                                             T[subject_idx],
#                                                             sigma = kernel_bw,
#                                                             trt_bounds = trt_bounds)
#       kernel_integral_minus <- integral_gaussian_kernel_minus(trt_clipped,
#                                                               T[subject_idx],
#                                                               sigma = kernel_bw,
#                                                               trt_bounds = trt_bounds)
#     }
#     term_db <-  - outcome_resid_scaled_plus * kernel_integral_minus - outcome_resid_scaled_minus * kernel_integral_plus
#   } else {
#     term_db = rep(0, length(subject_idx))
#   }
#   
#   
#   # --- Term 3: Outcome regression integral using prediction_smoothers ---
#   if(loss_type == "db" | loss_type == "or"){
#     term_outcome <- mapply(
#       FUN        = integrate_prediction_range,
#       lower      = trt_bounds[1],
#       upper      = trt_clipped,
#       obs_index  = subject_idx,
#       krige_value= krige_adjust[subject_idx],
#       MoreArgs   = list(cumint_smoothers = cumint_smoothers)
#     )
#   } else {
#     term_outcome <- rep(0, length(subject_idx))
#   }
#   
#   
#   # --- Combine terms ---
#   loss_vec <- term_threshold - term_db - term_outcome
#   
#   # if (any(!index_left_clip)) {
#   #   
#   #   idx <- !index_left_clip 
#   #   
#   #   # --- Term 1: Threshold integral ---
#   #   term_threshold <- threshold_val * (trt_right_clipped[idx] - trt_bounds[1])
#   #   
#   #   # --- Term 2: Doubly-robust kernel term or IPW term---
#   #   if (loss_type == "db" | loss_type == "ipw"){
#   #     
#   #     outcome_resid_scaled <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]])
#   #     outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
#   #     outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
#   #     
#   #     
#   #     if(surrogate_type == "Box"){
#   #       kernel_integral_plus <- kernel_integral_plus(trt = trt_right_clipped[idx],
#   #                                                     trt_obs = T[subject_idx[idx]],
#   #                                                     bandwidth = kernel_bw)
#   #       kernel_integral_minus <- kernel_integral_minus(trt = trt_right_clipped[idx],
#   #                                                     trt_obs = T[subject_idx[idx]],
#   #                                                     bandwidth = kernel_bw)
#   #     } else if (surrogate_type == "Gaussian"){
#   #       kernel_integral_plus <- integral_gaussian_kernel_plus(trt_right_clipped[idx],
#   #                                                       T[subject_idx[idx]],
#   #                                                       sigma = kernel_bw)
#   #       kernel_integral_minus <- integral_gaussian_kernel_minus(trt_right_clipped[idx],
#   #                                                             T[subject_idx[idx]],
#   #                                                             sigma = kernel_bw)
#   #     }
#   #     term_db <-  - outcome_resid_scaled_plus * kernel_integral_minus - outcome_resid_scaled_minus * kernel_integral_plus
#   #   } else {
#   #     term_db = rep(0, sum(idx))
#   #   }
#   #   
#   #   
#   #   # --- Term 3: Outcome regression integral using prediction_smoothers ---
#   #   if(loss_type == "db" | loss_type == "or"){
#   #     term_outcome <- mapply(
#   #       FUN        = integrate_prediction_range,
#   #       lower      = trt_bounds[1],
#   #       upper      = trt_right_clipped[idx],
#   #       obs_index  = subject_idx[idx],
#   #       krige_value= krige_adjust[subject_idx[idx]],
#   #       MoreArgs   = list(cumint_smoothers = cumint_smoothers)
#   #     )
#   #   } else {
#   #     term_outcome <- rep(0, sum(idx))
#   #   }
#   #   
#   #   
#   #   # --- Combine terms ---
#   #   loss_vec[idx] <- term_threshold - term_db - term_outcome
#   # }
#   
#   # ---------------------- left-clipped loss ----------------------
#   if (any(index_left_clip)) {
#     # when we revise loss outside boundary, we also need to revise loss in dc decomposition
#     
#     loss_vec[index_left_clip] <- loss_vec[index_left_clip] + clip_epsilon * (assigned_trt[index_left_clip] - trt_bounds[1])^2 + 
#       clip_epsilon_bar * (trt_bounds[1] - assigned_trt[index_left_clip])
#   }
#   
#   # ---------------------- right-clipped additional loss ----------------------
#   
#   if (any(index_right_clip)) {
#     # loss_vec[index_right_clip] <- loss_vec[index_right_clip] +  
#     #   (clip_epsilon_bar + 2 * clip_epsilon) * (assigned_trt[index_right_clip] - trt_bounds[2])
#     
#     # when we revise loss outside boundary, we also need to revise loss in dc decomposition
#     loss_vec[index_right_clip] <- loss_vec[index_right_clip] +  
#       clip_epsilon * (assigned_trt[index_right_clip] - trt_bounds[2])^2 + clip_epsilon_bar * (assigned_trt[index_right_clip] - trt_bounds[2])
#   }
#   
#   if (is.null(loss_vec) || any(!is.finite(loss_vec))) {
#     browser()
#   }
#   
#   return(loss_vec)
# }


# the re-written smoothed version 
compute_total_loss_smooth_plus <- function(assigned_trt,
                                           subject_idx,
                                           T,
                                           krige_adjust,
                                           outcome_resid,
                                           propensity_est,
                                           cumint_smoothers,
                                           trt_bounds,
                                           threshold_val,
                                           kernel_bw,
                                           clip_epsilon,
                                           clip_epsilon_bar,
                                           surrogate_type = c("Gaussian", "Box"),
                                           loss_type = c("db","ipw","or")) {
  
  surrogate_type <- match.arg(surrogate_type)
  loss_type <- match.arg(loss_type)
  n <- length(subject_idx)
  loss_vec <- rep(NA, n)
  
  
  # Identify clipped observations
  index_left_clip  <- assigned_trt < trt_bounds[1]
  index_right_clip <- assigned_trt > trt_bounds[2]
  index_no_clip    <- !index_left_clip & !index_right_clip
  
  # Clipped treatment values used 
  trt_clipped <- pmax(pmin(assigned_trt, trt_bounds[2]),trt_bounds[1])
  
  d_loss_minus_clipped_vec <- d_compute_total_loss_smooth_minus(assigned_trt = trt_clipped,
                                                        subject_idx = subject_idx,
                                                        T = T,
                                                        krige_adjust = krige_adjust,
                                                        outcome_resid = outcome_resid,
                                                        propensity_est = propensity_est,
                                                        trt_bounds = trt_bounds,
                                                        threshold_val = threshold_val,
                                                        kernel_bw = kernel_bw,
                                                        clip_epsilon_bar = clip_epsilon_bar,
                                                        surrogate_type = surrogate_type,
                                                        loss_type = loss_type)
  
  
  # ---------------------- uncapped and right-clipped loss ----------------------
  
  # --- Term 1: Threshold integral ---
  term_threshold <- threshold_val * (trt_clipped - trt_bounds[1])
  
  # --- Term 2: Doubly-robust kernel term or IPW term---
  if (loss_type == "db" | loss_type == "ipw"){
    
    outcome_resid_scaled <- (outcome_resid[subject_idx] / propensity_est[subject_idx])
    outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
    outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
    
    
    if(surrogate_type == "Box"){
      kernel_integral_plus <- kernel_integral_plus(trt = trt_clipped,
                                                   trt_obs = T[subject_idx],
                                                   bandwidth = kernel_bw)
      kernel_integral_minus <- kernel_integral_minus(trt = trt_clipped,
                                                     trt_obs = T[subject_idx],
                                                     bandwidth = kernel_bw)
    } else if (surrogate_type == "Gaussian"){
      kernel_integral_plus <- integral_gaussian_kernel_plus(trt_clipped,
                                                            T[subject_idx],
                                                            sigma = kernel_bw,
                                                            trt_bounds = trt_bounds)
      kernel_integral_minus <- integral_gaussian_kernel_minus(trt_clipped,
                                                              T[subject_idx],
                                                              sigma = kernel_bw,
                                                              trt_bounds = trt_bounds)
    }
    term_db <-  - outcome_resid_scaled_plus * kernel_integral_minus - outcome_resid_scaled_minus * kernel_integral_plus
  } else {
    term_db = rep(0, length(subject_idx))
  }
  
  
  # --- Term 3: Outcome regression integral using prediction_smoothers ---
  if(loss_type == "db" | loss_type == "or"){
    term_outcome <- mapply(
      FUN        = integrate_prediction_range,
      lower      = trt_bounds[1],
      upper      = trt_clipped,
      obs_index  = subject_idx,
      krige_value= krige_adjust[subject_idx],
      MoreArgs   = list(cumint_smoothers = cumint_smoothers)
    )
  } else {
    term_outcome <- rep(0, length(subject_idx))
  }
  
  
  # --- Combine terms ---
  loss_vec <- term_threshold - term_db - term_outcome
  
  # ---------------------- left-clipped loss ----------------------
  if (any(index_left_clip)) {
    # first add boundary loss in initial loss
    # loss_vec[index_left_clip] <- loss_vec[index_left_clip] + clip_epsilon * (assigned_trt[index_left_clip] - trt_bounds[1])^2
    loss_vec[index_left_clip] <- loss_vec[index_left_clip] + clip_epsilon - clip_epsilon / log( exp(1) + trt_bounds[1] - assigned_trt[index_left_clip])
    # then we add the boundary loss in loss_minus
    loss_vec[index_left_clip] <- loss_vec[index_left_clip] + d_loss_minus_clipped_vec[index_left_clip] * (assigned_trt[index_left_clip] - trt_bounds[1])
    loss_vec[index_left_clip] <- loss_vec[index_left_clip] + clip_epsilon_bar * (assigned_trt[index_left_clip] - trt_bounds[1])^2
  }
  
  # ---------------------- right-clipped additional loss ----------------------
  
  if (any(index_right_clip)) {
  
    # loss_vec[index_right_clip] <- loss_vec[index_right_clip] +  clip_epsilon * (assigned_trt[index_right_clip] - trt_bounds[2])^2 
    loss_vec[index_right_clip] <- loss_vec[index_right_clip] + clip_epsilon - clip_epsilon / log(exp(1)+ assigned_trt[index_right_clip] - trt_bounds[2])
    loss_vec[index_right_clip] <- loss_vec[index_right_clip] + d_loss_minus_clipped_vec[index_right_clip] * (assigned_trt[index_right_clip] - trt_bounds[2])
    loss_vec[index_right_clip] <- loss_vec[index_right_clip] + clip_epsilon_bar * (assigned_trt[index_right_clip] - trt_bounds[2])^2
  }
  
  if (is.null(loss_vec) || any(!is.finite(loss_vec))) {
    browser()
  }
  
  return(loss_vec)
}

## we would re-write this function to a smooth version

# compute_total_loss_smooth_minus <- function(assigned_trt,
#                                            subject_idx,
#                                            T,
#                                            krige_adjust,
#                                            outcome_resid,
#                                            propensity_est,
#                                            cumint_smoothers,
#                                            trt_bounds,
#                                            threshold_val,
#                                            kernel_bw,
#                                            clip_epsilon,
#                                            clip_epsilon_bar,
#                                            surrogate_type = c("Gaussian", "Box"),
#                                            loss_type = c("db","ipw","or")) {
#   
#   surrogate_type <- match.arg(surrogate_type)
#   loss_type <- match.arg(loss_type)
#   n <- length(subject_idx)
#   loss_vec <- rep(NA, n)
#   
#   # Identify clipped observations
#   index_left_clip  <- assigned_trt < trt_bounds[1]
#   index_right_clip <- assigned_trt > trt_bounds[2]
#   index_no_clip    <- !index_left_clip & !index_right_clip
#   
#   # Right-clipped treatment values used in integration
#   trt_right_clipped <- pmin(assigned_trt, trt_bounds[2])
#   trt_clipped <- pmax(trt_right_clipped, trt_bounds[1])
#   
#   # ---------------------- left-clipped loss ----------------------
#   # if (any(index_left_clip)) {
#   #   loss_vec[index_left_clip] <- - 2 * clip_epsilon * (assigned_trt[index_left_clip] - trt_bounds[1]) -
#   #     clip_epsilon + clip_epsilon * exp(assigned_trt[index_left_clip] - trt_bounds[1])
#   # }
# 
#   # ---------------------- uncapped and right-clipped loss ----------------------
#   # if (any(!index_left_clip)) {
#   # 
#   #   idx <- !index_left_clip
#   # 
#   #   # --- Term 2: Doubly-robust kernel term or IPW term---
#   #   if (loss_type == "db" | loss_type == "ipw"){
#   # 
#   #     outcome_resid_scaled <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]])
#   #     outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
#   #     outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
#   # 
#   # 
#   #     if(surrogate_type == "Box"){
#   #       kernel_integral_plus <- kernel_integral_plus(trt = trt_right_clipped[idx],
#   #                                                    trt_obs = T[subject_idx[idx]],
#   #                                                    bandwidth = kernel_bw)
#   #       kernel_integral_minus <- kernel_integral_minus(trt = trt_right_clipped[idx],
#   #                                                      trt_obs = T[subject_idx[idx]],
#   #                                                      bandwidth = kernel_bw)
#   #     } else if (surrogate_type == "Gaussian"){
#   #       kernel_integral_plus <- integral_gaussian_kernel_plus(trt_right_clipped[idx],
#   #                                                             T[subject_idx[idx]],
#   #                                                             sigma = kernel_bw)
#   #       kernel_integral_minus <- integral_gaussian_kernel_minus(trt_right_clipped[idx],
#   #                                                               T[subject_idx[idx]],
#   #                                                               sigma = kernel_bw)
#   #     }
#   #     term_db <-  - outcome_resid_scaled_plus * kernel_integral_plus - outcome_resid_scaled_minus * kernel_integral_minus
#   #   } else {
#   #     term_db = rep(0, sum(idx))
#   #   }
#   # 
#   #   # --- Combine terms ---
#   #   loss_vec[idx] <- - term_db
#   # }
#     
#   
#   # --- Term 2: Doubly-robust kernel term or IPW term---
#   if (loss_type == "db" | loss_type == "ipw"){
#     
#     outcome_resid_scaled <- (outcome_resid[subject_idx] / propensity_est[subject_idx])
#     outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
#     outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
#     
#     
#     if(surrogate_type == "Box"){
#       kernel_integral_plus <- kernel_integral_plus(trt = trt_clipped,
#                                                    trt_obs = T[subject_idx],
#                                                    bandwidth = kernel_bw)
#       kernel_integral_minus <- kernel_integral_minus(trt = trt_clipped,
#                                                      trt_obs = T[subject_idx],
#                                                      bandwidth = kernel_bw)
#     } else if (surrogate_type == "Gaussian"){
#       kernel_integral_plus <- integral_gaussian_kernel_plus(trt_clipped,
#                                                             T[subject_idx],
#                                                             sigma = kernel_bw,
#                                                             trt_bounds = trt_bounds)
#       kernel_integral_minus <- integral_gaussian_kernel_minus(trt_clipped,
#                                                               T[subject_idx],
#                                                               sigma = kernel_bw,
#                                                               trt_bounds = trt_bounds)
#     }
#     term_db <-  - outcome_resid_scaled_plus * kernel_integral_plus - outcome_resid_scaled_minus * kernel_integral_minus
#   } else {
#     term_db = rep(0, length(subject_idx))
#   }
#   
#   # --- Combine terms ---
#   loss_vec <- - term_db
#     
#   # ---------------------- left-clipped additional loss ----------------------
#   if (any(index_left_clip)) {
#     # loss_vec[index_left_clip] <- loss_vec[index_left_clip] - clip_epsilon_bar * (assigned_trt[index_left_clip] - trt_bounds[1]) - 
#     #   clip_epsilon + clip_epsilon * exp(assigned_trt[index_left_clip] - trt_bounds[1])
#     
#     # when we revise the loss function outside the boundary, we need to revise the loss function in DC decomposition
#     loss_vec[index_left_clip] <- loss_vec[index_left_clip] + clip_epsilon_bar * (trt_bounds[1] - assigned_trt[index_left_clip]) 
#   }
#   
#   # ---------------------- right-clipped additional loss ----------------------
#   if (any(index_right_clip)) {
#     # loss_vec[index_right_clip] <- loss_vec[index_right_clip] + 
#     #   (clip_epsilon_bar + 2 * clip_epsilon) * (assigned_trt[index_right_clip] - trt_bounds[2]) - 
#     #   clip_epsilon + clip_epsilon * exp(- assigned_trt[index_right_clip] + trt_bounds[2])
#     
#     # when we revise the loss function outside the boundary, we need to revise the loss function in DC decomposition
#     
#     loss_vec[index_right_clip] <- loss_vec[index_right_clip] + clip_epsilon_bar * (assigned_trt[index_right_clip] - trt_bounds[2])
#   }
#   
#   if (is.null(loss_vec) || any(!is.finite(loss_vec))) {
#     browser()
#   }
#   
#   return(loss_vec)
# }


# the re-written version
compute_total_loss_smooth_minus <- function(assigned_trt,
                                            subject_idx,
                                            T,
                                            krige_adjust,
                                            outcome_resid,
                                            propensity_est,
                                            # cumint_smoothers,
                                            trt_bounds,
                                            threshold_val,
                                            kernel_bw,
                                            # clip_epsilon,
                                            clip_epsilon_bar,
                                            surrogate_type = c("Gaussian", "Box"),
                                            loss_type = c("db","ipw","or")) {
  
  surrogate_type <- match.arg(surrogate_type)
  loss_type <- match.arg(loss_type)
  n <- length(subject_idx)
  loss_vec <- rep(NA, n)
  
  
  
  # Identify clipped observations
  index_left_clip  <- assigned_trt < trt_bounds[1]
  index_right_clip <- assigned_trt > trt_bounds[2]
  index_no_clip    <- !index_left_clip & !index_right_clip
  
  # Right-clipped treatment values used in integration
  trt_right_clipped <- pmin(assigned_trt, trt_bounds[2])
  trt_clipped <- pmax(trt_right_clipped, trt_bounds[1])
  
  # calculate derivative
  
  d_loss_minus_clipped_vec <- d_compute_total_loss_smooth_minus(assigned_trt = trt_clipped,
                                                  subject_idx = subject_idx,
                                                  T = T,
                                                  krige_adjust = krige_adjust,
                                                  outcome_resid = outcome_resid,
                                                  propensity_est = propensity_est,
                                                  trt_bounds = trt_bounds,
                                                  threshold_val = threshold_val,
                                                  kernel_bw = kernel_bw,
                                                  clip_epsilon_bar = clip_epsilon_bar,
                                                  surrogate_type = surrogate_type,
                                                  loss_type = loss_type)

  
  # --- Term 2: Doubly-robust kernel term or IPW term---
  if (loss_type == "db" | loss_type == "ipw"){
    
    outcome_resid_scaled <- (outcome_resid[subject_idx] / propensity_est[subject_idx])
    outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
    outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
    
    
    if(surrogate_type == "Box"){
      kernel_integral_plus <- kernel_integral_plus(trt = trt_clipped,
                                                   trt_obs = T[subject_idx],
                                                   bandwidth = kernel_bw)
      kernel_integral_minus <- kernel_integral_minus(trt = trt_clipped,
                                                     trt_obs = T[subject_idx],
                                                     bandwidth = kernel_bw)
    } else if (surrogate_type == "Gaussian"){
      kernel_integral_plus <- integral_gaussian_kernel_plus(trt_clipped,
                                                            T[subject_idx],
                                                            sigma = kernel_bw,
                                                            trt_bounds = trt_bounds)
      kernel_integral_minus <- integral_gaussian_kernel_minus(trt_clipped,
                                                              T[subject_idx],
                                                              sigma = kernel_bw,
                                                              trt_bounds = trt_bounds)
    }
    term_db <-  - outcome_resid_scaled_plus * kernel_integral_plus - outcome_resid_scaled_minus * kernel_integral_minus
  } else {
    term_db = rep(0, length(subject_idx))
  }
  
  # --- Combine terms ---
  loss_vec <- - term_db
  
  # ---------------------- left-clipped additional loss ----------------------
  if (any(index_left_clip)) {
    # we use the derivative at boundary to expand a linear term
    
    loss_vec[index_left_clip] <- loss_vec[index_left_clip] + d_loss_minus_clipped_vec[index_left_clip] * (assigned_trt[index_left_clip] - trt_bounds[1]) 
    loss_vec[index_left_clip] <- loss_vec[index_left_clip] + clip_epsilon_bar * (assigned_trt[index_left_clip] - trt_bounds[1])^2
  }
  
  # ---------------------- right-clipped additional loss ----------------------
  if (any(index_right_clip)) {
    # we use the derivative at boundary to expand a linear term
    
    loss_vec[index_right_clip] <- loss_vec[index_right_clip] + d_loss_minus_clipped_vec[index_right_clip] * (assigned_trt[index_right_clip] - trt_bounds[2])
    loss_vec[index_right_clip] <- loss_vec[index_right_clip] + clip_epsilon_bar * (assigned_trt[index_right_clip] - trt_bounds[2])^2
  }
  
  if (is.null(loss_vec) || any(!is.finite(loss_vec))) {
    browser()
  }
  
  return(loss_vec)
}



# ## derivative of compute_total_loss_smooth_plus, we would rewrite-this function
# d_compute_total_loss_smooth_plus <- function(assigned_trt,
#                                            subject_idx,
#                                            T,
#                                            krige_adjust,
#                                            outcome_resid,
#                                            propensity_est,
#                                            smoothers,
#                                            trt_bounds,
#                                            threshold_val,
#                                            kernel_bw,
#                                            clip_epsilon,
#                                            clip_epsilon_bar,
#                                            surrogate_type = c("Gaussian", "Box"),
#                                            loss_type = c("db","ipw","or"),
#                                            delta = 0.001) {
#   
#   surrogate_type <- match.arg(surrogate_type)
#   n <- length(subject_idx)
#   d_loss_vec <- rep(NA, n)
#   
#   # delta <- 0.001 # boundary smoother
#   # delta <- 0 # boundary smoother
#   
#   # Identify clipped observations
#   index_left_clip  <- assigned_trt < (trt_bounds[1] - delta)
#   index_right_clip <- assigned_trt > (trt_bounds[2] + delta)
#   index_no_clip    <- !index_left_clip & !index_right_clip
#   
#   # Right-clipped treatment values used in integration
#   trt_right_clipped <- pmin(assigned_trt, trt_bounds[2])
#   
#   # ---------------------- left-clipped loss ----------------------
#   if (any(index_left_clip)) {
#     d_loss_vec[index_left_clip] <- - clip_epsilon_bar 
#     # change when boundary loss changed
#     d_loss_vec[index_left_clip] <- 2 * (clip_epsilon) * (assigned_trt[index_left_clip] - trt_bounds[1]) - clip_epsilon_bar
#   }
#   
#   # ---------------------- uncapped and right-clipped loss ----------------------
#   if (any(!index_left_clip)) {
#     
#     idx <- !index_left_clip 
#     
#     # --- Term 1: Derivative of Threshold integral ---
#     # term_threshold <- rep(threshold_val, sum(idx)) 
#     term_threshold <- rep(threshold_val, sum(idx))
#     
#     # --- Term 2: Derivative of Doubly-robust kernel term ---
#     
#     if (loss_type == "db"| loss_type == "ipw"){
#     
#       outcome_resid_scaled <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]])
#       outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
#       outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
#       
#       
#       if(surrogate_type == "Box"){
#         kernel_plus <- kernel_integral_plus_derivative(trt = trt_right_clipped[idx],
#                                                      trt_obs = T[subject_idx[idx]],
#                                                      bandwidth = kernel_bw)
#         kernel_minus <- kernel_integral_minus_derivative(trt = trt_right_clipped[idx],
#                                                        trt_obs = T[subject_idx[idx]],
#                                                        bandwidth = kernel_bw)
#       } else if (surrogate_type == "Gaussian"){
#         kernel_plus <- gaussian_kernel_plus(trt_right_clipped[idx],
#                                                               T[subject_idx[idx]],
#                                                               sigma = kernel_bw)
#         kernel_minus <- gaussian_kernel_minus(trt_right_clipped[idx],
#                                                                 T[subject_idx[idx]],
#                                                                 sigma = kernel_bw)
#       }
#       term_db <-  - outcome_resid_scaled_plus * kernel_minus - outcome_resid_scaled_minus * kernel_plus
#     } else {
#       term_db <- rep(0, sum(idx))
#     }
#     
#     # --- Term 3: Outcome regression integral using prediction_smoothers ---
#     
#     if(loss_type == "db" | loss_type == "or"){
#       term_outcome <- sapply(which(idx), 
#                              function(i){ smoothers[[subject_idx[i]]](trt_right_clipped[i]) + 
#                                  krige_adjust[subject_idx[i]]})
#     } else {
#       term_outcome <- rep(0, sum(idx))
#     }
#     
#     
#     # so we can use predict.smooth.spline() to get the derivative 
#     
#     # --- Combine terms ---
#     # is it idx not subject[idx], yes, it is idx
#     d_loss_vec[idx] <- term_threshold - term_db - term_outcome
#   }
#   
#   # ---------------------- right-clipped additional loss ----------------------
#   if (any(index_right_clip)) {
#     # d_loss_vec[index_right_clip] <- 2 * clip_epsilon + clip_epsilon_bar
#     
#     # change when boundary loss changed
#     d_loss_vec[index_right_clip] <- 2 * (clip_epsilon) * (assigned_trt[index_right_clip] - trt_bounds[2]) + clip_epsilon_bar
#   }
#   
#   if (is.null(d_loss_vec) || any(!is.finite(d_loss_vec))) {
#     browser()
#   }
#   
#   return(d_loss_vec)
# }

# the re-written smoothed version

d_compute_total_loss_smooth_plus <- function(assigned_trt,
                                             subject_idx,
                                             T,
                                             krige_adjust,
                                             outcome_resid,
                                             propensity_est,
                                             smoothers,
                                             trt_bounds,
                                             threshold_val,
                                             kernel_bw,
                                             clip_epsilon,
                                             clip_epsilon_bar,
                                             surrogate_type = c("Gaussian", "Box"),
                                             loss_type = c("db","ipw","or")) {
  
  surrogate_type <- match.arg(surrogate_type)
  n <- length(subject_idx)
  d_loss_vec <- rep(NA, n)
  
  # Identify clipped observations
  index_left_clip  <- assigned_trt < (trt_bounds[1])
  index_right_clip <- assigned_trt > (trt_bounds[2])
  index_no_clip    <- !index_left_clip & !index_right_clip
  
  # Right-clipped treatment values used in integration
  trt_right_clipped <- pmin(assigned_trt, trt_bounds[2])
  trt_clipped <- pmax(trt_right_clipped, trt_bounds[1])
  
  d_loss_minus_clipped_vec <- d_compute_total_loss_smooth_minus(assigned_trt = trt_clipped,
                                                        subject_idx = subject_idx,
                                                        T = T,
                                                        krige_adjust = krige_adjust,
                                                        outcome_resid = outcome_resid,
                                                        propensity_est = propensity_est,
                                                        trt_bounds = trt_bounds,
                                                        threshold_val = threshold_val,
                                                        kernel_bw = kernel_bw,
                                                        clip_epsilon_bar = clip_epsilon_bar,
                                                        surrogate_type = surrogate_type,
                                                        loss_type = loss_type)
  
  # ---------------------- left-clipped loss ----------------------
  if (any(index_left_clip)) {
    
    # first add the outside-boundary d_loss of the original function, # then add the d_loss_minus
    # d_loss_vec[index_left_clip] <- 2 * (clip_epsilon) * (assigned_trt[index_left_clip] - trt_bounds[1]) + d_loss_minus_vec[index_left_clip]
    
    d_loss_vec[index_left_clip] <-  - clip_epsilon /(log(exp(1) + trt_bounds[1] - assigned_trt[index_left_clip])^2 * (exp(1) + trt_bounds[1] - assigned_trt[index_left_clip])) + 
      d_loss_minus_clipped_vec[index_left_clip] + 2 * clip_epsilon_bar * (assigned_trt[index_left_clip] - trt_bounds[1])
    
  }
  
  # ---------------------- uncapped and right-clipped loss ----------------------
  if (any(!index_left_clip)) {
    
    idx <- !index_left_clip 
    
    # --- Term 1: Derivative of Threshold integral ---
    # term_threshold <- rep(threshold_val, sum(idx)) 
    term_threshold <- rep(threshold_val, sum(idx))
    
    # --- Term 2: Derivative of Doubly-robust kernel term ---
    
    if (loss_type == "db"| loss_type == "ipw"){
      
      outcome_resid_scaled <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]])
      outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
      outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
      
      
      if(surrogate_type == "Box"){
        kernel_plus <- kernel_integral_plus_derivative(trt = trt_right_clipped[idx],
                                                       trt_obs = T[subject_idx[idx]],
                                                       bandwidth = kernel_bw)
        kernel_minus <- kernel_integral_minus_derivative(trt = trt_right_clipped[idx],
                                                         trt_obs = T[subject_idx[idx]],
                                                         bandwidth = kernel_bw)
      } else if (surrogate_type == "Gaussian"){
        kernel_plus <- gaussian_kernel_plus(trt_right_clipped[idx],
                                            T[subject_idx[idx]],
                                            sigma = kernel_bw)
        kernel_minus <- gaussian_kernel_minus(trt_right_clipped[idx],
                                              T[subject_idx[idx]],
                                              sigma = kernel_bw)
      }
      term_db <-  - outcome_resid_scaled_plus * kernel_minus - outcome_resid_scaled_minus * kernel_plus
    } else {
      term_db <- rep(0, sum(idx))
    }
    
    # --- Term 3: Outcome regression integral using prediction_smoothers ---
    
    if(loss_type == "db" | loss_type == "or"){
      term_outcome <- sapply(which(idx), 
                             function(i){ smoothers[[subject_idx[i]]](trt_right_clipped[i]) + 
                                 krige_adjust[subject_idx[i]]})
    } else {
      term_outcome <- rep(0, sum(idx))
    }
    
    
    # so we can use predict.smooth.spline() to get the derivative 
    
    # --- Combine terms ---
    # is it idx not subject[idx], yes, it is idx
    d_loss_vec[idx] <- term_threshold - term_db - term_outcome
  }
  
  # ---------------------- right-clipped additional loss ----------------------
  if (any(index_right_clip)) {
    
    # first add the outside-boundary d_loss of the original function
    # d_loss_vec[index_right_clip] <- 2 * (clip_epsilon) * (assigned_trt[index_right_clip] - trt_bounds[2]) 
    d_loss_vec[index_right_clip] <- clip_epsilon /(log(exp(1) + assigned_trt[index_right_clip] - trt_bounds[2])^2 * (exp(1) + assigned_trt[index_right_clip] - trt_bounds[2]))
    d_loss_vec[index_right_clip] <- d_loss_vec[index_right_clip] + d_loss_minus_clipped_vec[index_right_clip] + 2 * clip_epsilon_bar * (assigned_trt[index_right_clip] - trt_bounds[2])
  }
  
  if (is.null(d_loss_vec) || any(!is.finite(d_loss_vec))) {
    browser()
  }
  
  return(d_loss_vec)
}


# we would rewrite this function
# d_compute_total_loss_smooth_minus <- function(assigned_trt,
#                                              subject_idx,
#                                              T,
#                                              krige_adjust,
#                                              outcome_resid,
#                                              propensity_est,
#                                              smoothers,
#                                              trt_bounds,
#                                              threshold_val,
#                                              kernel_bw,
#                                              clip_epsilon,
#                                              clip_epsilon_bar,
#                                              surrogate_type = c("Gaussian", "Box"),
#                                              loss_type = c("db","ipw","or"),
#                                              delta = 0.001) {
# 
#   surrogate_type <- match.arg(surrogate_type)
#   n <- length(subject_idx)
#   d_loss_vec <- rep(NA, n)
# 
#   # delta <- 0.001 # boundary smoother
#   # delta <- 0 # boundary smoother
# 
#   # Identify clipped observations
# 
#   index_left_clip  <- assigned_trt < trt_bounds[1] - delta
#   index_right_clip <- assigned_trt > trt_bounds[2] + delta
#   # index_no_clip    <- !index_left_clip & !index_right_clip
#   index_no_clip    <- !index_left_clip & !index_right_clip 
# 
#   # Right-clipped treatment values used in integration
#   # trt_right_clipped <- pmin(assigned_trt, trt_bounds[2])
#   trt_clipped <- pmax(pmin(assigned_trt, trt_bounds[2]), trt_bounds[1])
# 
# 
#   # ---------------------- uncapped and right-clipped loss ----------------------
#   # if (any(!index_left_clip)) {
#   if (any(index_no_clip)) {
# 
#     idx <- !index_left_clip
#     idx <- index_left_clip
# 
#     # --- Term 2: Derivative of Doubly-robust kernel term ---
# 
#     if (loss_type == "db"| loss_type == "ipw"){
# 
#       outcome_resid_scaled <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]])
#       outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
#       outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
# 
# 
#       if(surrogate_type == "Box"){
#         kernel_plus <- kernel_integral_plus_derivative(trt = trt_right_clipped[idx],
#                                                        trt_obs = T[subject_idx[idx]],
#                                                        bandwidth = kernel_bw)
#         kernel_minus <- kernel_integral_minus_derivative(trt = trt_right_clipped[idx],
#                                                          trt_obs = T[subject_idx[idx]],
#                                                          bandwidth = kernel_bw)
#       } else if (surrogate_type == "Gaussian"){
#         kernel_plus <- gaussian_kernel_plus(trt_right_clipped[idx],
#                                             T[subject_idx[idx]],
#                                             sigma = kernel_bw)
#         kernel_minus <- gaussian_kernel_minus(trt_right_clipped[idx],
#                                               T[subject_idx[idx]],
#                                               sigma = kernel_bw)
#       }
#       term_db <-  - outcome_resid_scaled_plus * kernel_plus - outcome_resid_scaled_minus * kernel_minus
#     } else {
#       term_db <- rep(0, sum(idx))
#     }
# 
#     # --- Term 3: Outcome regression integral using prediction_smoothers ---
# 
# 
# 
#     # so we can use predict.smooth.spline() to get the derivative
# 
#     # --- Combine terms ---
#     # is it idx not subject[idx], yes, it is idx
#     d_loss_vec[idx] <- - term_db
#   }
# 
# 
#   # ---------------------- left-clipped loss ----------------------
#   if (any(index_left_clip)) {
#     # d_loss_vec[index_left_clip] <- - clip_epsilon_bar + clip_epsilon * exp(assigned_trt[index_left_clip] - trt_bounds[1])
# 
#     # when we revise the boundary loss, we need to revise derivative of dc decomposition
# 
#     d_loss_vec[index_left_clip] <- - clip_epsilon_bar
#     
#     # how about if we directly use the boundary value as the derivative to make it smooth
#     
# 
# 
#   }
# 
# 
# 
#   # ---------------------- right-clipped additional loss ----------------------
#   if (any(index_right_clip)) {
#     # d_loss_vec[index_right_clip] <- (clip_epsilon_bar + 2 * clip_epsilon) -
#     #   clip_epsilon * exp( - assigned_trt[index_right_clip] + trt_bounds[2])
# 
#     # when # when we revise the boundary loss, we need to revise derivative of dc decomposition
# 
#     d_loss_vec[index_right_clip] <- clip_epsilon_bar
# 
#   }
# 
#   if (is.null(d_loss_vec) || any(!is.finite(d_loss_vec))) {
#     browser()
#   }
# 
#   return(d_loss_vec)
# }
# this is not the boudary version, this is a smooth version
d_compute_total_loss_smooth_minus <- function(assigned_trt,
                                              subject_idx,
                                              T,
                                              krige_adjust,
                                              outcome_resid,
                                              propensity_est,
                                              trt_bounds,
                                              threshold_val,
                                              kernel_bw,
                                              clip_epsilon_bar,
                                              surrogate_type = c("Gaussian", "Box"),
                                              loss_type = c("db","ipw","or")) {
  
  surrogate_type <- match.arg(surrogate_type)
  n <- length(subject_idx)
  d_loss_vec <- rep(NA, n)
  
  
  
  # Identify clipped observations
  
  index_left_clip  <- assigned_trt < trt_bounds[1] 
  index_right_clip <- assigned_trt > trt_bounds[2] 
  index_no_clip    <- !index_left_clip & !index_right_clip 
  
  # clipped treatment values
  trt_clipped <- pmax(pmin(assigned_trt, trt_bounds[2]), trt_bounds[1])
  
  # The derivative of the derivative outside the boundary is equal to its value at the boundary point, in this case, it is linear (convex)
  
  idx <- 1:n
  # --- Term 2: Derivative of Doubly-robust kernel term ---
  
  if (loss_type == "db"| loss_type == "ipw"){
    
    outcome_resid_scaled <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]])
    outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
    outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
    
    
    if(surrogate_type == "Box"){
      kernel_plus <- kernel_integral_plus_derivative(trt = trt_clipped[idx],
                                                     trt_obs = T[subject_idx[idx]],
                                                     bandwidth = kernel_bw)
      kernel_minus <- kernel_integral_minus_derivative(trt = trt_clipped[idx],
                                                       trt_obs = T[subject_idx[idx]],
                                                       bandwidth = kernel_bw)
    } else if (surrogate_type == "Gaussian"){
      kernel_plus <- gaussian_kernel_plus(trt_clipped[idx],
                                          T[subject_idx[idx]],
                                          sigma = kernel_bw)
      kernel_minus <- gaussian_kernel_minus(trt_clipped[idx],
                                            T[subject_idx[idx]],
                                            sigma = kernel_bw)
    }
    term_db <-  - outcome_resid_scaled_plus * kernel_plus - outcome_resid_scaled_minus * kernel_minus
  } else {
    term_db <- rep(0, sum(idx))
  }
  
  
  d_loss_vec[idx] <- - term_db
  
  if(any(index_left_clip)){
    d_loss_vec[index_left_clip] <- d_loss_vec[index_left_clip] + 2 * clip_epsilon_bar * (assigned_trt[index_left_clip] - trt_bounds[1])
  }
  
  if(any(index_right_clip)){
    d_loss_vec[index_right_clip] <- d_loss_vec[index_right_clip] + 2 * clip_epsilon_bar * (assigned_trt[index_right_clip] - trt_bounds[2])
  }
  
  if (is.null(d_loss_vec) || any(!is.finite(d_loss_vec))) {
    browser()
  }
  
  return(d_loss_vec)
}

## This is the re-write version
## it turns out a smoothing version of average as sub-gradient is not usefull
# d_compute_total_loss_smooth_minus <- function(assigned_trt,
#                                               subject_idx,
#                                               T,
#                                               krige_adjust,
#                                               outcome_resid,
#                                               propensity_est,
#                                               smoothers,
#                                               trt_bounds,
#                                               threshold_val,
#                                               kernel_bw,
#                                               clip_epsilon,
#                                               clip_epsilon_bar,
#                                               surrogate_type = c("Gaussian", "Box"),
#                                               loss_type = c("db","ipw","or")) {
#   
#   surrogate_type <- match.arg(surrogate_type)
#   n <- length(subject_idx)
#   d_loss_vec <- rep(NA, n)
#   
#   delta_1 <- 0.01 # boundary smoother
#   delta_2 <- - 0.001 # boundary smoother
#   # when specify delta_1 and delta_2, the boundary of left is (trt_bounds[1] - delta_1,trt_bounds[1] + delta_2)
#   # the boundary of right is (trt_bounds[2] - delta_2,trt_bounds[2] + delta_1)
#   # delta <- 0 # boundary smoother
#   
#   # Identify clipped observations
#   
#   index_left_clip  <- assigned_trt < trt_bounds[1] - delta_1
#   index_left_boundary <- (assigned_trt >= (trt_bounds[1] - delta_1)) & (assigned_trt <= (trt_bounds[1] + delta_2))
#   index_right_clip <- assigned_trt > trt_bounds[2] + delta_1 
#   index_right_boundary <- (assigned_trt >= (trt_bounds[2] - delta_2)) & (assigned_trt <= (trt_bounds[2] + delta_1))
#   index_no_clip    <- !index_left_clip & !index_right_clip
#   index_no_boundary    <- !index_left_clip & !index_right_clip & !index_left_boundary & !index_right_boundary
#   
#   # Right-clipped treatment values used in integration
#   # trt_right_clipped <- pmin(assigned_trt, trt_bounds[2])
#   
#   # boundary-clipped treatment values used in pre-calculation
#   trt_no_boundary <-pmax(pmin(assigned_trt, trt_bounds[2] - delta_2), trt_bounds[1] + delta_2)
#   
#   
#   
#   # ---------------------- non-boundary loss and pre-compute boundary loss ----------------------
#   if (any(index_no_clip)) {
#     
#     idx <- index_no_clip
#     
#     # --- Term 2: Derivative of Doubly-robust kernel term ---
#     
#     if (loss_type == "db"| loss_type == "ipw"){
#       
#       outcome_resid_scaled <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]])
#       outcome_resid_scaled_plus <- pmax(0, outcome_resid_scaled)
#       outcome_resid_scaled_minus <- - pmin(0, outcome_resid_scaled)
#       
#       
#       if(surrogate_type == "Box"){
#         kernel_plus <- kernel_integral_plus_derivative(trt = trt_no_boundary[idx],
#                                                        trt_obs = T[subject_idx[idx]],
#                                                        bandwidth = kernel_bw)
#         kernel_minus <- kernel_integral_minus_derivative(trt = trt_no_boundary[idx],
#                                                          trt_obs = T[subject_idx[idx]],
#                                                          bandwidth = kernel_bw)
#       } else if (surrogate_type == "Gaussian"){
#         kernel_plus <- gaussian_kernel_plus(trt_no_boundary[idx],
#                                             T[subject_idx[idx]],
#                                             sigma = kernel_bw)
#         kernel_minus <- gaussian_kernel_minus(trt_no_boundary[idx],
#                                               T[subject_idx[idx]],
#                                               sigma = kernel_bw)
#       }
#       term_db <-  - outcome_resid_scaled_plus * kernel_plus - outcome_resid_scaled_minus * kernel_minus
#     } else {
#       term_db <- rep(0, sum(idx))
#     }
#     
#     # --- Term 3: Outcome regression integral using prediction_smoothers ---
#     
#     
#     
#     # so we can use predict.smooth.spline() to get the derivative 
#     
#     # --- Combine terms ---
#     # is it idx not subject[idx], yes, it is idx
#     d_loss_vec[idx] <- - term_db 
#   }
#   
#   # ---------------------- left-boundary loss ----------------------
#   if (any(index_left_boundary)){
#     idx <- index_left_boundary
#     # d_loss_vec[idx] <- (assigned_trt[idx] - (trt_bounds[1] - delta))/(2*delta) * d_loss_vec[idx]  + (trt_bounds[1] + delta - assigned_trt[idx])/(2*delta) * (-clip_epsilon_bar) 
#     d_loss_vec[idx] <- (assigned_trt[idx] - (trt_bounds[1] - delta_1))/(delta_1 + delta_2) * d_loss_vec[idx]  + (trt_bounds[1] + delta_2 - assigned_trt[idx])/(delta_1 + delta_2) * (-clip_epsilon_bar) 
#   }
#   
#   
#   
#   # ---------------------- left-clipped loss ----------------------
#   if (any(index_left_clip)) {
#     # d_loss_vec[index_left_clip] <- - clip_epsilon_bar + clip_epsilon * exp(assigned_trt[index_left_clip] - trt_bounds[1])
#     
#     # when we revise the boundary loss, we need to revise derivative of dc decomposition
#     
#     d_loss_vec[index_left_clip] <- - clip_epsilon_bar 
#     
#     
#   }
#   
#   # ---------------------- right-boundary loss ----------------------
#   if (any(index_right_boundary)){
#     idx <- index_right_boundary
#     # d_loss_vec[idx] <- (assigned_trt[idx] - (trt_bounds[2] - delta))/(2*delta) * (clip_epsilon_bar)   + (trt_bounds[2] + delta - assigned_trt[idx])/(2*delta) * d_loss_vec[idx]
#     d_loss_vec[idx] <- (assigned_trt[idx] - (trt_bounds[2] - delta_2))/(delta_1 + delta_2) * (clip_epsilon_bar)   + (trt_bounds[2] + delta_1 - assigned_trt[idx])/(delta_1 + delta_2) * d_loss_vec[idx]
#   }
#   
#   
#   
#   # ---------------------- right-clipped  loss ----------------------
#   if (any(index_right_clip)) {
#     # d_loss_vec[index_right_clip] <- (clip_epsilon_bar + 2 * clip_epsilon) - 
#     #   clip_epsilon * exp( - assigned_trt[index_right_clip] + trt_bounds[2])
#     
#     # when # when we revise the boundary loss, we need to revise derivative of dc decomposition
#     
#     d_loss_vec[index_right_clip] <- clip_epsilon_bar
#     
#   }
#   
#   if (is.null(d_loss_vec) || any(!is.finite(d_loss_vec))) {
#     browser()
#   }
#   
#   return(d_loss_vec)
# }
# 


compute_total_loss_smooth_dc_approx <- function(assigned_trt,
                                                assigned_trt_old,
                                                subject_idx,
                                                T,
                                                krige_adjust,
                                                outcome_resid,
                                                propensity_est,
                                                cumint_smoothers,
                                                smoothers,
                                                trt_bounds,
                                                threshold_val,
                                                kernel_bw,
                                                clip_epsilon,
                                                clip_epsilon_bar,
                                                surrogate_type = c("Gaussian", "Box"),
                                                loss_type = c("db","ipw","or")){
  
  surrogate_type <- match.arg(surrogate_type)
  loss_type <- match.arg(loss_type)
  
  total_loss_plus <- compute_total_loss_smooth_plus(assigned_trt = assigned_trt,
                                                    subject_idx = subject_idx,
                                                    T = T,
                                                    krige_adjust = krige_adjust,
                                                    outcome_resid = outcome_resid,
                                                    propensity_est = propensity_est,
                                                    cumint_smoothers = cumint_smoothers,
                                                    trt_bounds = trt_bounds,
                                                    threshold_val = threshold_val,
                                                    kernel_bw = kernel_bw,
                                                    clip_epsilon = clip_epsilon,
                                                    clip_epsilon_bar = clip_epsilon_bar,
                                                    surrogate_type = surrogate_type,
                                                    loss_type = loss_type) 
  # print('total_loss_plus')
  # print(total_loss_plus)
  
  total_loss_minus_old <- compute_total_loss_smooth_minus(assigned_trt = assigned_trt_old,
                                                    subject_idx = subject_idx,
                                                    T = T,
                                                    krige_adjust = krige_adjust,
                                                    outcome_resid = outcome_resid,
                                                    propensity_est = propensity_est,
                                                    # cumint_smoothers = cumint_smoothers,
                                                    trt_bounds = trt_bounds,
                                                    threshold_val = threshold_val,
                                                    kernel_bw = kernel_bw,
                                                    # clip_epsilon = clip_epsilon,
                                                    clip_epsilon_bar = clip_epsilon_bar,
                                                    surrogate_type = surrogate_type,
                                                    loss_type = loss_type) 
  # print('total_loss_minus_old')
  # print(total_loss_minus_old)
  
  d_total_loss_minus_old <- d_compute_total_loss_smooth_minus(assigned_trt = assigned_trt_old,
                                                              subject_idx = subject_idx,
                                                              T = T,
                                                              krige_adjust = krige_adjust,
                                                              outcome_resid = outcome_resid,
                                                              propensity_est = propensity_est,
                                                              # smoothers = smoothers,
                                                              trt_bounds = trt_bounds,
                                                              threshold_val = threshold_val,
                                                              kernel_bw = kernel_bw,
                                                              # clip_epsilon = clip_epsilon,
                                                              clip_epsilon_bar = clip_epsilon_bar,
                                                              surrogate_type = surrogate_type,
                                                              loss_type = loss_type)
  # print('d_total_loss_minus_old')
  # print(d_total_loss_minus_old)
  
  
  total_loss_dc_approx <- total_loss_plus - (total_loss_minus_old + d_total_loss_minus_old * (assigned_trt - assigned_trt_old))
  
  if (sum(!is.finite(total_loss_dc_approx)) | sum(is.na(total_loss_dc_approx))) {
    cat("NA in total_loss_dc_approx")
    browser()
  }
  
  return(total_loss_dc_approx)
  
}

d_compute_total_loss_smooth_dc_approx <- function(assigned_trt,
                                                  assigned_trt_old,
                                                  subject_idx,
                                                  T,
                                                  krige_adjust,
                                                  outcome_resid,
                                                  propensity_est,
                                                  cumint_smoothers,
                                                  smoothers,
                                                  trt_bounds,
                                                  threshold_val,
                                                  kernel_bw,
                                                  clip_epsilon,
                                                  clip_epsilon_bar,
                                                  surrogate_type = c("Gaussian", "Box"),
                                                  loss_type = c("db","ipw","or")){
  
  surrogate_type <- match.arg(surrogate_type)
  loss_type <- match.arg(loss_type)
  
  d_total_loss_plus <- d_compute_total_loss_smooth_plus(assigned_trt = assigned_trt,
                                                        subject_idx = subject_idx,
                                                        T = T,
                                                        krige_adjust = krige_adjust,
                                                        outcome_resid = outcome_resid,
                                                        propensity_est = propensity_est,
                                                        smoothers = smoothers,
                                                        trt_bounds = trt_bounds,
                                                        threshold_val = threshold_val,
                                                        kernel_bw = kernel_bw,
                                                        clip_epsilon = clip_epsilon,
                                                        clip_epsilon_bar = clip_epsilon_bar,
                                                        surrogate_type = surrogate_type,
                                                        loss_type = loss_type) 
  
  d_total_loss_minus_old <- d_compute_total_loss_smooth_minus(assigned_trt = assigned_trt_old,
                                                              subject_idx = subject_idx,
                                                              T = T,
                                                              krige_adjust = krige_adjust,
                                                              outcome_resid = outcome_resid,
                                                              propensity_est = propensity_est,
                                                              # smoothers = smoothers,
                                                              trt_bounds = trt_bounds,
                                                              threshold_val = threshold_val,
                                                              kernel_bw = kernel_bw,
                                                              # clip_epsilon = clip_epsilon,
                                                              clip_epsilon_bar = clip_epsilon_bar,
                                                              surrogate_type = surrogate_type,
                                                              loss_type = loss_type)
  
  
  d_total_loss_dc_approx <- d_total_loss_plus - d_total_loss_minus_old
  
  
  
  if (sum(!is.finite(d_total_loss_dc_approx)) | sum(is.na(d_total_loss_dc_approx))) {
    cat("NA in d_total_loss_dc_approx")
    browser()
  }
  
  return(d_total_loss_dc_approx)
  
}

compute_total_loss_smooth_dc_approx_RKHS <-function(params, 
                                                    params_old, 
                                                    K, 
                                                    T, 
                                                    krige_adjust, 
                                                    outcome_resid, 
                                                    propensity_est, 
                                                    lambda, 
                                                    smoothers, 
                                                    cumint_smoothers, 
                                                    kernel_bw, 
                                                    trt_bounds,
                                                    threshold_val, 
                                                    clip_epsilon,
                                                    clip_epsilon_bar,
                                                    surrogate_type = c("Gaussian", "Box"),
                                                    loss_type = c("db","ipw","or")){
      
      surrogate_type <- match.arg(surrogate_type)
      loss_type <- match.arg(loss_type)
      
      n <- nrow(K)
      
      assigned_trt <- as.vector(K %*% params)
      
      assigned_trt_old <- as.vector(K %*% params_old)
      
      total_loss_dc_approx <- compute_total_loss_smooth_dc_approx(assigned_trt = as.numeric(K %*% params),
                                                                  assigned_trt_old = assigned_trt_old,
                                                                  subject_idx = 1:n,
                                                                  T = T,
                                                                  krige_adjust = krige_adjust,
                                                                  outcome_resid = outcome_resid,
                                                                  propensity_est = propensity_est,
                                                                  cumint_smoothers = cumint_smoothers,
                                                                  smoothers = smoothers,
                                                                  trt_bounds = trt_bounds,
                                                                  threshold_val = threshold_val,
                                                                  kernel_bw = kernel_bw,
                                                                  clip_epsilon = clip_epsilon,
                                                                  clip_epsilon_bar = clip_epsilon_bar,
                                                                  surrogate_type = surrogate_type,
                                                                  loss_type = loss_type)
      
      
      if(lambda == 0){
        penalty <- 0
      } else {
        if(ncol(K) < nrow(K)){
          penalty <- (lambda / 2) *(t(params[-1]) %*% params[-1])
        } else {
          penalty <- (lambda / 2) * (t(params[-1]) %*% K[,-1] %*% params[-1])
        }
      }
      
      if (sum(!is.finite(total_loss_dc_approx)) | sum(is.na(total_loss_dc_approx))) {
        cat("NA in total_loss_dc_approx")
        browser()
      }
      
      return(sum(total_loss_dc_approx) + penalty + 15000)
      
    }


d_compute_total_loss_smooth_dc_approx_RKHS <- function(params, 
                                                      params_old, 
                                                      K, 
                                                      T, 
                                                      krige_adjust, 
                                                      outcome_resid, 
                                                      propensity_est, 
                                                      lambda, 
                                                      smoothers, 
                                                      cumint_smoothers, 
                                                      kernel_bw, 
                                                      trt_bounds,
                                                      threshold_val, 
                                                      clip_epsilon,
                                                      clip_epsilon_bar,
                                                      surrogate_type = c("Gaussian", "Box"),
                                                      loss_type = c("db","ipw","or")){
      
  surrogate_type <- match.arg(surrogate_type)
  loss_type <- match.arg(loss_type)
  
  n <- nrow(K)
  
  assigned_trt <- as.vector(K %*% params)
  
  assigned_trt_old <- as.vector(K %*% params_old)

  d_total_loss_dc_approx <- d_compute_total_loss_smooth_dc_approx(assigned_trt = assigned_trt,
                                                                  assigned_trt_old = assigned_trt_old,
                                                                  subject_idx = 1:n,
                                                                  T = T,
                                                                  krige_adjust = krige_adjust,
                                                                  outcome_resid = outcome_resid,
                                                                  propensity_est = propensity_est,
                                                                  cumint_smoothers = cumint_smoothers,
                                                                  smoothers = smoothers,
                                                                  trt_bounds = trt_bounds,
                                                                  threshold_val = threshold_val,
                                                                  kernel_bw = kernel_bw,
                                                                  clip_epsilon = clip_epsilon,
                                                                  clip_epsilon_bar = clip_epsilon_bar,
                                                                  surrogate_type = surrogate_type,
                                                                  loss_type = loss_type)
  
  return(as.numeric((t(K) %*% d_total_loss_dc_approx) +  lambda * c(0, (K[,-1] %*% params[-1]))))
}




compute_total_loss_smooth_RKHS <- function(params, 
                                           K, 
                                           T,
                                           krige_adjust,
                                           outcome_resid,
                                           propensity_est,
                                           lambda, 
                                           smoothers, 
                                           cumint_smoothers, 
                                           kernel_bw, 
                                           trt_bounds,
                                           threshold_val, 
                                           clip_epsilon,
                                           surrogate_type = c("Gaussian", "Box"),
                                           loss_type = c("db","ipw","or")){
  surrogate_type <- match.arg(surrogate_type)
  loss_type <- match.arg(loss_type)
  
  n <- nrow(K)
  
  assigned_trt <- as.vector(K %*% params)

  total_loss <- compute_total_loss_sum_smooth(assigned_trt = as.numeric(K %*% params),
                                          subject_idx = 1:n,
                                          T = T,
                                          krige_adjust = krige_adjust,
                                          outcome_resid = outcome_resid,
                                          propensity_est = propensity_est,
                                          cumint_smoothers = cumint_smoothers,
                                          trt_bounds = trt_bounds,
                                          threshold_val = threshold_val,
                                          kernel_bw = kernel_bw,
                                          clip_epsilon = clip_epsilon,
                                          surrogate_type = surrogate_type,
                                          loss_type = loss_type)
  if(lambda == 0){
    penalty <- 0
  } else {
    if(ncol(K) < nrow(K)){
      penalty <-  (lambda / 2) *(t(params[-1]) %*% params[-1])
    } else {
      penalty <-  (lambda / 2) * (t(params[-1]) %*% K[,-1] %*% params[-1])
    }
  }
  
  return(sum(total_loss) + penalty + 17000)
  
}


# the derivative of loss of ipw, or, and db loss

d_compute_total_loss_smooth_RKHS <-  function(params, 
                                              K, 
                                              T,
                                              krige_adjust,
                                              outcome_resid,
                                              propensity_est,
                                              smoothers, 
                                              cumint_smoothers, 
                                              lambda,
                                              trt_bounds,
                                              threshold_val,
                                              kernel_bw,
                                              clip_epsilon,
                                              surrogate_type = c("Gaussian", "Box"),
                                              loss_type = c("db","ipw","or")){
  surrogate_type <- match.arg(surrogate_type)
  loss_type <- match.arg(loss_type)
  
  n <- nrow(K)
  
  assigned_trt <- as.vector(K %*% params)
  
  
  d_total_loss <- d_compute_total_loss_smooth(assigned_trt = assigned_trt,
                                              subject_idx = 1:n,
                                              T = T,
                                              krige_adjust = krige_adjust,
                                              outcome_resid = outcome_resid,
                                              propensity_est = propensity_est,
                                              smoothers = smoothers,
                                              trt_bounds = trt_bounds,
                                              threshold_val = threshold_val,
                                              kernel_bw = kernel_bw,
                                              clip_epsilon = clip_epsilon,
                                              surrogate_type = surrogate_type,
                                              loss_type = loss_type)
  if(lambda==0){
    penalty <- rep(0, length(params))
  } else {
    penalty <- lambda * c(0, (K[,-1] %*% params[-1]))
  }
  
  return(as.numeric((t(K) %*% d_total_loss) + penalty))
  
}

