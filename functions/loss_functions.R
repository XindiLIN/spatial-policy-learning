library(parallel)
library(stats)

# ---------------------------------------------------------------------------- 
# Surrogate Kernels
# ---------------------------------------------------------------------------- 



############################################################
# Gaussian Kernel Functions
#
# These functions define the Gaussian kernel (PDF), its 
# integral (CDF), and a treatment-based shifted integral
# similar to the Box kernel's "plus" / "minus" versions.
#
# Notation:
#   x        - point(s) at which to evaluate
#   mu       - mean of the Gaussian distribution
#   sigma    - standard deviation of the Gaussian distribution
#   trt      - treatment value(s)
#   trt_obs  - observed treatment value(s)
#   bandwidth- kernel bandwidth (scalar)
#   side     - "plus" or "minus" for shifted integration
############################################################

#' Gaussian kernel (PDF)
#'
#' Calculates the value of the Gaussian kernel at given points.
#'
#' @param x Numeric vector or scalar: points at which to evaluate.
#' @param mu Mean of the Gaussian kernel (default 0).
#' @param sigma Standard deviation of the Gaussian kernel (default 1).
#' @return Numeric vector: PDF values of the Gaussian kernel at \code{x}.
#' @examples
#' gaussian_kernel(0)             # PDF at 0 for N(0,1)
#' gaussian_kernel(c(-1, 0, 1))
gaussian_kernel <- function(x, mu = 0, sigma = 1) {
  (1 / (sigma * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sigma)^2)
}

#' Integral of the Gaussian kernel (CDF)
#'
#' Computes the cumulative probability from -Inf to x.
#'
#' @param x Numeric vector or scalar: upper limit(s) of integration.
#' @param mu Mean of the Gaussian kernel (default 0).
#' @param sigma Standard deviation of the Gaussian kernel (default 1).
#' @return Numeric vector: CDF values of the Gaussian kernel at \code{x}.
#' @examples
#' integral_gaussian_kernel(0)    # CDF at 0 for N(0,1)
#' integral_gaussian_kernel(c(-1, 0, 1))
integral_gaussian_kernel <- function(x, mu = 0, sigma = 1) {
  pnorm(x, mean = mu, sd = sigma)
}

#' Shifted Gaussian kernel integral (plus / minus)
#'
#' Integrates the Gaussian kernel from -Inf up to a shifted treatment point,
#' where the shift depends on the "side" argument:
#'   - "plus"  : shift = trt_obs - bandwidth
#'   - "minus" : shift = trt_obs + bandwidth
#'
#' @param trt Numeric: treatment values at which to integrate.
#' @param trt_obs Numeric: observed treatment values.
#' @param bandwidth Positive numeric: kernel bandwidth.
#' @param side Character: "plus" or "minus".
#' @return Numeric vector: shifted CDF values.
#' @examples
#' gaussian_kernel_integral(0.5, 0.3, 0.1, "plus")
#' gaussian_kernel_integral(0.5, 0.3, 0.1, "minus")
gaussian_kernel_integral <- function(trt, trt_obs, bandwidth, side = c("plus", "minus")) {
  side <- match.arg(side)
  
  shift <- if (side == "plus") {
    trt_obs - bandwidth
  } else {
    trt_obs + bandwidth
  }
  
  pnorm(trt, mean = shift, sd = bandwidth)
}

# The above integral is not correct. we want to split to subtraction of convex function


############################################################
# Box Kernel Integral and Derivative (Plus / Minus variants)
#
# This implementation covers the "plus" and "minus" cases 
# in a single function each, using a 'side' argument.
# 
# Notation:
#   trt      - treatment value (can be scalar or vector)
#   trt_obs  - observed treatment value (scalar or vector)
#   bandwidth- kernel bandwidth (scalar)
#   side     - "plus" or "minus" kernel type
#
# Box kernel shape:
#   K(u) = 0.5 * I(|u| <= 1)
#   Here we integrate K from a boundary to trt, shifting by ±bandwidth
#
# The "plus" integral starts from (trt_obs - bandwidth).
# The "minus" integral starts from (trt_obs + bandwidth).
############################################################

# ---- Kernel integral function ----
kernel_integral <- function(trt, trt_obs, bandwidth, side = c("plus", "minus")) {
  side <- match.arg(side)   # Ensure side is either "plus" or "minus"
  
  inv_bw <- 1 / (2 * bandwidth)  # Precompute 1/(2h), constant factor in Box kernel
  # Determine the lower bound for integration depending on the side
  shift <- if (side == "plus") {
    trt_obs - bandwidth
  } else {
    trt_obs + bandwidth
  }
  
  # Integral value: zero if trt < shift; otherwise inv_bw * (trt - shift)
  # pmax(0, ...) enforces the non-negative region without explicit if/else
  pmax(0, inv_bw * (trt - shift))
}

# ---- Derivative of the kernel integral ----
kernel_integral_derivative <- function(trt, trt_obs, bandwidth, side = c("plus", "minus")) {
  side <- match.arg(side)
  
  # The derivative is constant (1 / (2h)) when trt >= shift, otherwise 0
  shift <- if (side == "plus") {
    trt_obs - bandwidth
  } else {
    trt_obs + bandwidth
  }
  
  as.numeric(trt >= shift) * (1 / (2 * bandwidth))
}

############################################################
# Example usage:
# kernel_integral(0.5, 0.3, 0.1, "plus")
# kernel_integral(0.5, 0.3, 0.1, "minus")
# kernel_integral_derivative(0.5, 0.3, 0.1, "plus")
# kernel_integral_derivative(0.5, 0.3, 0.1, "minus")
############################################################

# ---------------------------------------------------------------------------- 

# ---------------------------------------------------------------------------- #
# Outcome Regression Integrals Approximation
# ---------------------------------------------------------------------------- 

############################################################
# Precompute Counterfactual Prediction Matrix
#
# This function generates a matrix of predicted outcomes
# for each observation in the dataset, across a grid of
# possible treatment values. This is useful for causal
# inference tasks where we need to integrate or average
# predicted outcomes over the treatment distribution.
#
# Arguments:
#   data_obs        - Data frame of observed covariates and treatment.
#   fit_outcome     - Fitted outcome regression model (lm, gam, SuperLearner, etc.).
#   treatment_range - Length-2 numeric vector: min and max treatment values.
#   treatment_step  - Numeric: step size between treatment values.
#   pred_args       - Named list of extra arguments to pass to predict().
#                     Examples:
#                       GLM: list(type = "response")
#                       SuperLearner: list(type = "response")
#
# Returns:
#   A list with:
#     pred_matrix       - Matrix of predictions:
#                           rows = treatment values
#                           cols = original observations
#     treatment_values  - Numeric vector of treatment values (row index)
#
# Example:
#   fit <- lm(mpg ~ wt + hp + trt, data = df)
#   result <- precompute_cf_predictions(df, fit, c(0, 1), 0.05,
#                                       pred_args = list(type = "response"))
#   result$pred_matrix  # predicted outcomes
#   result$treatment_values  # grid of treatment values
############################################################

precompute_cf_predictions <- function(data_obs, fit_outcome,
                                      treatment_range = c(0, 1),
                                      treatment_step = 0.01,
                                      pred_args = list(),
                                      treatment_name = 'T') {
  
  # ---- 1. Create treatment grid ----
  treatment_values <- seq(treatment_range[1], treatment_range[2], by = treatment_step)
  
  # ---- 2. Expand original dataset for all treatment values ----
  # Each original observation is replicated for each treatment value
  data_expanded <- data_obs[rep(seq_len(nrow(data_obs)), each = length(treatment_values)), ]
  data_expanded[,treatment_name] <- rep(treatment_values, times = nrow(data_obs))
  
  # ---- 3. Generate predictions ----
  # Combine model + expanded data + user-supplied predict() args
  preds <- do.call(predict, c(list(object = fit_outcome, newdata = data_expanded), pred_args))
  if(is.list(preds)){
    preds <- preds$pred
  }
  
  # If the model returns a matrix or data.frame (e.g., SuperLearner), take the first column
  if (is.matrix(preds) || is.data.frame(preds)) {
    preds <- preds[, 1]
  }
  
  # ---- 4. Reshape into a matrix ----
  # Rows = treatment values, Columns = original observations
  pred_matrix <- matrix(preds, nrow = length(treatment_values), ncol = nrow(data_obs))
  
  # ---- 5. Return results ----
  list(pred_matrix = pred_matrix, treatment_values = treatment_values)
}


# ------------------------------------------------------------------------------
# make_prediction_smoothers
# ------------------------------------------------------------------------------
# Purpose:
#   Given a treatment grid and predicted outcomes for each observation
#   (e.g., output from precompute_cf_predictions), this function builds
#   a list of smooth functions mapping treatment -> predicted outcome.
#
# Arguments:
#   treatment_grid : numeric vector of treatment values (length m)
#   pred_matrix_cf : numeric matrix of predicted outcomes (m x n),
#                    rows = treatment values, columns = observations
#   smooth_method  : smoothing approach, one of:
#                      "smooth.spline" - cubic smoothing spline
#                      "splinefun"     - exact interpolation spline
#                      "loess"         - local polynomial regression
#   spar           : smoothing parameter for smooth.spline (NULL = auto)
#
# Details:
#   - Uses local() to ensure each smoother function captures its *own* fit object
#     instead of referencing the last fit in the loop.
#   - Each function in the output takes a numeric vector x (treatment values)
#     and returns the smoothed predicted outcome(s).
#
# Returns:
#   List of length n (number of observations). Each element is a function f(x)
#   giving smoothed predictions at new treatment x.
# ------------------------------------------------------------------------------


isotomic_correction <- function(treatment_grid, pred_matrix_cf){
  n_subjects <- ncol(pred_matrix_cf)       # number of observations
  
  # --- Helper for decreasing isotonic regression ---
  # isoreg() only fits increasing functions, so we use the trick:
  # decreasing_fit(y) = -increasing_fit(-y)
  iso_decrease <- function(x, y) {
    # Fit a monotonically INCREASING function to the NEGATED values
    iso_fit <- isoreg(x, -y)
    # Return the NEGATED fitted values to make them decreasing
    return(-iso_fit$yf)
  }
  
  
  # Loop over each observation
  for (j in seq_len(n_subjects)) {
    
    # Extract predicted values for this observation across treatment grid
    pred_values <- pred_matrix_cf[, j]
    pred_values <- iso_decrease(x = treatment_grid, y = pred_values)
    pred_matrix_cf[, j] <- pred_values
  }
  return(list(treatment_values = treatment_grid, pred_matrix_cf = pred_matrix_cf))
}
make_prediction_smoothers <- function(treatment_grid,
                                      pred_matrix_cf,
                                      smooth_method = c("smooth.spline", "splinefun", "loess", "gam"),
                                      spar = NULL) {
  
  # Ensure chosen smoothing method is valid
  smooth_method <- match.arg(smooth_method)
  
  n_subjects <- ncol(pred_matrix_cf)       # number of observations
  smoothers <- vector("list", n_subjects)  # allocate list for smooth functions
  
  # Loop over each observation
  for (j in seq_len(n_subjects)) {
    
    # Extract predicted values for this observation across treatment grid
    pred_values <- pred_matrix_cf[, j]
    
    # ---------------------- smoothing methods ----------------------
    if (smooth_method == "gam"){
      # this is now not usable, see pcls() and mono.con() in mgcv
      df <- data.frame(pred_values, treatment_grid)
      fit <- mgcv::gam(pred_values ~ s(treatment_grid, bs="cr", m=2), data=df)
      
      smoother_fn <- local({
        fit_copy <- fit
        function(x) predict(fit_copy, newdata = data.frame(treatment_grid = x))
      })
      
    } else if (smooth_method == "smooth.spline") {
      # ---- Option 1: Smoothing spline ----
      # Fit a cubic smoothing spline (controlled by 'spar' or automatic)
      fit <- smooth.spline(treatment_grid, pred_values, spar = spar)
      
      # Wrap in a closure to preserve the fit object
      smoother_fn <- local({
        fit_copy <- fit
        function(x) predict(fit_copy, x)$y
      })
      
    } else if (smooth_method == "splinefun") {
      # ---- Option 2: Interpolation spline ----
      # Exact fit through points; "fmm" = cubic spline
      smoother_fn <- local({
        f_raw_copy <- splinefun(treatment_grid, pred_values, method = "monoH.FC")
        function(x) f_raw_copy(x)
      })
      
    } else { 
      # ---- Option 3: Loess smoothing ----
      # Local polynomial regression (degree 2, span 0.5)
      fit <- loess(pred_values ~ treatment_grid, span = 0.5, degree = 2)
      smoother_fn <- local({
        fit_copy <- fit
        function(x) predict(fit_copy, newdata = data.frame(treatment_grid = x))
      })
    }
    
    # Store the smoother function for this observation
    smoothers[[j]] <- smoother_fn
  }
  
  return(smoothers)
}


# ------------------------------------------------------------------------------
# make_cumint_smoothers
# ------------------------------------------------------------------------------
# Purpose:
#   Given a treatment grid and predicted outcomes for each observation
#   (matrix of size m x n), compute *cumulative integrals* of predictions
#   using the trapezoidal rule, then smooth those cumulative integrals.
#
# Arguments:
#   treatment_grid : numeric vector of treatment values (length m)
#   pred_matrix_cf : numeric matrix of predicted outcomes (m x n),
#                    rows = treatment values, columns = observations
#   smooth_method  : smoothing approach, one of:
#                      "smooth.spline" - cubic smoothing spline
#                      "splinefun"     - exact interpolation spline
#   spar           : smoothing parameter for smooth.spline (NULL = auto)
#
# Details:
#   - Uses the trapezoidal rule to approximate cumulative integrals of
#     predicted outcomes across the treatment grid.
#   - Uses local() to ensure each smoother function captures its *own* fit
#     object instead of referencing the last fit in the loop.
#   - Each function in the output takes a numeric vector x (treatment values)
#     and returns the smoothed cumulative integral at x.
#
# Returns:
#   List of length n (number of observations). Each element is a function f(x)
#   giving smoothed cumulative integrals at new treatment x.
# ------------------------------------------------------------------------------

make_cumint_smoothers <- function(treatment_grid,
                                   pred_matrix_cf,
                                   smooth_method = c("smooth.spline", "splinefun", "gam"),
                                   spar = NULL) {
  
  # Ensure chosen smoothing method is valid
  smooth_method <- match.arg(smooth_method)
  
  n_subjects <- ncol(pred_matrix_cf)         # number of observations
  cumint_smoothers <- vector("list", n_subjects)  # allocate list
  
  # Pre-compute dt for trapezoidal rule
  dt <- diff(treatment_grid)
  
  # Loop over each observation
  for (j in seq_len(n_subjects)) {
    
    # Predicted values for this subject
    pred_values <- pred_matrix_cf[, j]
    
    # Trapezoidal cumulative integral
    trapz_vals <- dt * (pred_values[-length(pred_values)] + pred_values[-1]) / 2
    cum_int <- c(0, cumsum(trapz_vals))  # same length as treatment_grid
    
    # ---------------------- smoothing methods ----------------------
    if (smooth_method == "gam"){
      # this is now not usable, see pcls() and mono.con() in mgcv
      df <- data.frame(cum_int, treatment_grid)
      fit <- mgcv::gam(cum_int ~ s(treatment_grid, bs="cr", m=2), data=df)
      
      smoother_fn <- local({
        fit_copy <- fit
        function(x) predict(fit_copy, newdata = data.frame(treatment_grid = x))
      })
      
    } else if (smooth_method == "smooth.spline") {
      # ---- Option 1: Smoothing spline ----
      fit <- smooth.spline(treatment_grid, cum_int, spar = spar)
      smoother_fn <- local({
        fit_copy <- fit
        function(x) predict(fit_copy, x)$y
      })
      
    } else { 
      # ---- Option 2: Interpolation spline ----
      smoother_fn <- local({
        f_raw_copy <- splinefun(treatment_grid, cum_int, method = "fmm")
        function(x) f_raw_copy(x)
      })
    }
    
    # Store smoother for this subject
    cumint_smoothers[[j]] <- smoother_fn
  }
  
  return(cumint_smoothers)
}

# calculate the derivative of cumint_smoothers
make_d_cumint_smoothers <- function(treatment_grid,
                                    pred_matrix_cf,
                                    smooth_method = c("smooth.spline", "splinefun", "gam"),
                                    spar = NULL) {
  
  # Ensure chosen smoothing method is valid
  smooth_method <- match.arg(smooth_method)
  
  n_subjects <- ncol(pred_matrix_cf)         # number of observations
  cumint_smoothers <- vector("list", n_subjects)  # allocate list
  
  # Pre-compute dt for trapezoidal rule
  dt <- diff(treatment_grid)
  
  # Loop over each observation
  for (j in seq_len(n_subjects)) {
    
    # Predicted values for this subject
    pred_values <- pred_matrix_cf[, j]
    
    # Trapezoidal cumulative integral
    trapz_vals <- dt * (pred_values[-length(pred_values)] + pred_values[-1]) / 2
    cum_int <- c(0, cumsum(trapz_vals))  # same length as treatment_grid
    
    # ---------------------- smoothing methods ----------------------
   if (smooth_method == "smooth.spline") {
      # ---- Option 1: Smoothing spline ----
      fit <- smooth.spline(treatment_grid, cum_int, spar = spar)
      smoother_fn <- local({
        fit_copy <- fit
        function(x) predict(fit_copy, x, deriv = 1)$y
      })
      
    }
    # Store smoother for this subject
    cumint_smoothers[[j]] <- smoother_fn
  }
  
  return(cumint_smoothers)
}




# ------------------------------------------------------------------------------
# integrate_prediction_range
# ------------------------------------------------------------------------------
# Purpose:
#   Numerically integrate a prediction smoother function over a treatment interval [lower, upper]
#   for a single observation, and add a linear adjustment term.
#
# Arguments:
#   lower, upper     : numeric, integration bounds
#   obs_index        : integer, index of the observation (column in pred_smoothers)
#   krige_value      : numeric, linear adjustment to add over the interval
#   pred_smoothers   : list of prediction smoother functions (from make_prediction_smoothers)
#
# Returns:
#   Numeric: integral of predicted outcome over [lower, upper] + krige_value adjustment
# ------------------------------------------------------------------------------

integrate_prediction_range <- function(lower, upper,
                                       obs_index,
                                       krige_value,
                                       cumint_smoothers) {
  
  # Validate inputs
  if (is.null(cumint_smoothers)) {
    stop("cumint_smoothers cannot be NULL")
  }
  
  if (upper <= lower) return(0)  # Return 0 for invalid or empty intervals
  
  # Extract the smoother function for the observation
  f_fn <- cumint_smoothers[[obs_index]]
  if (is.null(f_fn)) stop("cumint_smoothers[[obs_index]] is NULL")
  
  # Add linear adjustment
  integral_val <- f_fn(upper) - f_fn(lower) + krige_value * (upper - lower)
  
  return(integral_val)
}


# ------------------------------------------------------------------------------ 
# Loss Function: Doubly Robust 
# ------------------------------------------------------------------------------ 

# ------------------------------------------------------------------------------
# compute_total_loss_db_smooth
# ------------------------------------------------------------------------------
# Purpose:
#   Compute the "total loss" for a vector of assigned treatments, including:
#     1. Threshold term (linear)
#     2. Doubly-robust (DB) kernel term (Gaussian or Box)
#     3. Integrated predicted outcome (from prediction_smoothers)
#
# Arguments:
#   assigned_trt         : numeric vector of assigned treatment values
#   subject_idx          : integer vector of observation indices
#   data                 : data.frame with original observations (must contain "trt")
#   krige_adjust         : numeric vector of krige adjustment values
#   outcome_resid        : numeric vector of residuals from outcome regression
#   propensity_est       : numeric vector of propensity score estimates
#   prediction_smoothers : list of prediction smoother functions (from make_prediction_smoothers)
#   trt_bounds           : numeric vector c(min, max) defining valid treatment range
#   threshold_val        : numeric, threshold weight for linear term
#   kernel_bw            : numeric, bandwidth for DB kernel
#   clip_epsilon         : numeric, small constant for clipped treatment loss
#   surrogate_type       : character, "Gaussian" or "Box" for DB term
#
# Returns:
#   Numeric vector: total loss for each observation
# ------------------------------------------------------------------------------
compute_total_loss_db_smooth <- function(assigned_trt,
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
                                         surrogate_type = c("Gaussian", "Box")) {
  
  surrogate_type <- match.arg(surrogate_type)
  n <- length(subject_idx)
  loss_vec <- rep(NA, n)
  
  # Identify clipped observations
  index_left_clip  <- assigned_trt < trt_bounds[1]
  index_right_clip <- assigned_trt > trt_bounds[2]
  index_no_clip    <- !index_left_clip & !index_right_clip
  
  # Right-clipped treatment values used in integration
  trt_right_clipped <- pmin(assigned_trt, trt_bounds[2])
  
  # ---------------------- left-clipped loss ----------------------
  if (any(index_left_clip)) {
    loss_vec[index_left_clip] <- clip_epsilon - clip_epsilon *
      exp(assigned_trt[index_left_clip] - trt_bounds[1])
  }
  
  # ---------------------- uncapped and right-clipped loss ----------------------
  if (any(!index_left_clip)) {
    
    idx <- !index_left_clip 
    
    # --- Term 1: Threshold integral ---
    term_threshold <- threshold_val * (trt_right_clipped[idx] - trt_bounds[1])
    
    # --- Term 2: Doubly-robust kernel term ---
    if (surrogate_type == "Box") {
      term_db <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]]) *
        kernel_integral_vectorized(trt = trt_right_clipped[idx],
                                   trt_obs = T[subject_idx[idx]],
                                   bandwidth = kernel_bw)
    } else if (surrogate_type == "Gaussian") {
      term_db <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]]) *
        integral_gaussian_kernel(trt_right_clipped[idx],
                                 T[subject_idx[idx]],
                                 sigma = kernel_bw)
    }
    
    # --- Term 3: Outcome regression integral using prediction_smoothers ---
    term_outcome <- mapply(
      FUN        = integrate_prediction_range,
      lower      = trt_bounds[1],
      upper      = trt_right_clipped[idx],
      obs_index  = subject_idx[idx],
      krige_value= krige_adjust[subject_idx[idx]],
      MoreArgs   = list(cumint_smoothers = cumint_smoothers)
    )
    
    # --- Combine terms ---
    loss_vec[idx] <- term_threshold - term_db - term_outcome
  }
  
  # ---------------------- right-clipped additional loss ----------------------
  if (any(index_right_clip)) {
    loss_vec[index_right_clip] <- loss_vec[index_right_clip] +
      clip_epsilon - clip_epsilon *
      exp(-assigned_trt[index_right_clip] + trt_bounds[2])
  }
  
  if (is.null(loss_vec) || any(!is.finite(loss_vec))) {
    browser()
  }
  
  return(loss_vec)
}


compute_total_loss_smooth <- function(assigned_trt,
                                         subject_idx,
                                         T,
                                         krige_adjust,
                                         outcome_resid,
                                         propensity_est,
                                         smoothers,
                                         cumint_smoothers,
                                         trt_bounds,
                                         threshold_val,
                                         kernel_bw,
                                         clip_epsilon,
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
  
  # ---------------------- left-clipped loss ----------------------
  if (any(index_left_clip)) {
    # loss_vec[index_left_clip] <- clip_epsilon - clip_epsilon *
    #   exp(assigned_trt[index_left_clip] - trt_bounds[1]) # for this one, the derivative goes to zero as the treatment gets further from to boundary
    # loss_vec[index_left_clip] <- clip_epsilon * (trt_bounds[1] - assigned_trt[index_left_clip]) # this one it not continuous at the boundary
    # loss_vec[index_left_clip] <- clip_epsilon * (trt_bounds[1] - assigned_trt[index_left_clip])^2 # so we try this more smooth one
    loss_vec[index_left_clip] <- clip_epsilon - clip_epsilon / log(exp(1) + trt_bounds[1] - assigned_trt[index_left_clip]) # the loss function has to be bounded
  }
  
  # ---------------------- uncapped and right-clipped loss ----------------------
  if (any(!index_left_clip)) {
    
    idx <- !index_left_clip 
    
    # --- Term 1: Threshold integral ---
    term_threshold <- threshold_val * (trt_right_clipped[idx] - trt_bounds[1])
    
    # --- Term 2: Doubly-robust kernel term or IPW term---
    if (loss_type == "db"| loss_type == "ipw"){
      if (surrogate_type == "Box") {
        term_db <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]]) *
          kernel_integral_vectorized(trt = trt_right_clipped[idx],
                                     trt_obs = T[subject_idx[idx]],
                                     bandwidth = kernel_bw)
      } else if (surrogate_type == "Gaussian") {
        term_db <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]]) *
          integral_gaussian_kernel(trt_right_clipped[idx],
                                   T[subject_idx[idx]],
                                   sigma = kernel_bw)
      }
    } else {
      term_db = rep(0, sum(idx))
    }
    
    
    # --- Term 3: Outcome regression integral using prediction_smoothers ---
    if(loss_type == "db" | loss_type == "or"){
      term_outcome <- mapply(
        FUN        = integrate_prediction_range,
        lower      = trt_bounds[1],
        upper      = trt_right_clipped[idx],
        obs_index  = subject_idx[idx],
        krige_value= krige_adjust[subject_idx[idx]],
        MoreArgs   = list(cumint_smoothers = cumint_smoothers)
      )
    } else {
      term_outcome <- rep(0, sum(idx))
    }
    
    
    # --- Combine terms ---
    loss_vec[idx] <- term_threshold - term_db - term_outcome
  }
  
  # ---------------------- right-clipped additional loss ----------------------
  if (any(index_right_clip)) {
    # loss_vec[index_right_clip] <- loss_vec[index_right_clip] +
    #   clip_epsilon - clip_epsilon *
    #   exp(-assigned_trt[index_right_clip] + trt_bounds[2])
    # loss_vec[index_right_clip] <- loss_vec[index_right_clip] + clip_epsilon * (assigned_trt[index_right_clip] - trt_bounds[2])
    # loss_vec[index_right_clip] <- loss_vec[index_right_clip] + clip_epsilon * (assigned_trt[index_right_clip] - trt_bounds[2])^2
    loss_vec[index_right_clip] <- loss_vec[index_right_clip] + clip_epsilon - clip_epsilon / log(exp(1)+ assigned_trt[index_right_clip] - trt_bounds[2])
  }
  
  if (is.null(loss_vec) || any(!is.finite(loss_vec))) {
    cat('null or infinite in loss_vec')
    browser()
  }
  
  return(loss_vec)
}






## derivative of compute_total_loss_db_sum_smooth
d_compute_total_loss_db_smooth <- function(assigned_trt,
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
                                          surrogate_type = c("Gaussian", "Box")) {
  
  surrogate_type <- match.arg(surrogate_type)
  n <- length(subject_idx)
  d_loss_vec <- rep(NA, n)
  
  # Identify clipped observations
  index_left_clip  <- assigned_trt < trt_bounds[1]
  index_right_clip <- assigned_trt > trt_bounds[2]
  index_no_clip    <- !index_left_clip & !index_right_clip
  
  # Right-clipped treatment values used in integration
  trt_right_clipped <- pmin(assigned_trt, trt_bounds[2])
  
  # ---------------------- left-clipped loss ----------------------
  if (any(index_left_clip)) {
    d_loss_vec[index_left_clip] <- clip_epsilon *
      exp(assigned_trt[index_left_clip] - trt_bounds[1])
  }
  
  # ---------------------- uncapped and right-clipped loss ----------------------
  if (any(!index_left_clip)) {
    
    idx <- !index_left_clip 
    
    # --- Term 1: Derivative of Threshold integral ---
    # term_threshold <- rep(threshold_val, sum(idx)) 
    term_threshold <- threshold_val
    
    # --- Term 2: Derivative of Doubly-robust kernel term ---
   
    
    term_db <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]]) *
                gaussian_kernel(trt_right_clipped[idx],
                               T[subject_idx[idx]],
                               sigma = kernel_bw)
    
    # --- Term 3: Outcome regression integral using prediction_smoothers ---
    
    
    # I am not sure if this is caused by the difference of smoothers and cumint_smoothers so
    # the derivative is not correct
    term_outcome <- sapply(which(idx), function(i){ smoothers[[subject_idx[i]]](trt_right_clipped[i]) + krige_adjust[subject_idx[i]]})
    # so we can use predict.smooth.spline() to get the derivative 
    
    # --- Combine terms ---
    # is it idx not subject[idx], yes, it is idx
    d_loss_vec[idx] <- term_threshold - term_db - term_outcome
  }
  
  # ---------------------- right-clipped additional loss ----------------------
  if (any(index_right_clip)) {
    d_loss_vec[index_right_clip] <- d_loss_vec[index_right_clip] + clip_epsilon *
      exp(-assigned_trt[index_right_clip] + trt_bounds[2])
  }
  
  if (is.null(d_loss_vec) || any(!is.finite(d_loss_vec))) {
    browser()
  }
  
  return(d_loss_vec)
}

# the derivative of loss of ipw, or, and db loss
d_compute_total_loss_smooth <- function(assigned_trt,
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
                                       surrogate_type = c("Gaussian", "Box"),
                                       loss_type = c("db","ipw","or")) {
  
  surrogate_type <- match.arg(surrogate_type)
  loss_type <- match.arg(loss_type)
  n <- length(subject_idx)
  d_loss_vec <- rep(NA, n)
  
  # Identify clipped observations
  index_left_clip  <- assigned_trt < trt_bounds[1]
  index_right_clip <- assigned_trt > trt_bounds[2]
  index_no_clip    <- !index_left_clip & !index_right_clip
  
  # Right-clipped treatment values used in integration
  trt_right_clipped <- pmin(assigned_trt, trt_bounds[2])
  
  # ---------------------- left-clipped loss ----------------------
  if (any(index_left_clip)) {
    # d_loss_vec[index_left_clip] <-  - clip_epsilon *
    #   exp(assigned_trt[index_left_clip] - trt_bounds[1])
    # d_loss_vec[index_left_clip] <- - clip_epsilon
    # d_loss_vec[index_left_clip] <- 2 * clip_epsilon * (assigned_trt[index_left_clip] - trt_bounds[1]) 
    
    d_loss_vec[index_left_clip] <-  - clip_epsilon /(log(exp(1) + trt_bounds[1] - assigned_trt[index_left_clip])^2 * (exp(1) + trt_bounds[1] - assigned_trt[index_left_clip])) # make it upper bounded 
  }
  
  # ---------------------- uncapped and right-clipped loss ----------------------
  if (any(!index_left_clip)) {
    
    idx <- !index_left_clip 
    
    # --- Term 1: Derivative of Threshold integral ---
    # term_threshold <- rep(threshold_val, sum(idx)) 
    term_threshold <- rep(threshold_val, sum(idx))
    
    # --- Term 2: Derivative of Doubly-robust kernel term ---
    
    if(loss_type == "db"| loss_type == "ipw"){
      term_db <- (outcome_resid[subject_idx[idx]] / propensity_est[subject_idx[idx]]) *
        gaussian_kernel(trt_right_clipped[idx],
                        T[subject_idx[idx]],
                        sigma = kernel_bw)
    } else {
      term_db <- rep(0, sum(idx))
    }
    
    
    # --- Term 3: Outcome regression integral using prediction_smoothers ---
    
    
    # I am not sure if this is caused by the difference of smoothers and cumint_smoothers so
    # the derivative is not correct
    if(loss_type == "db"| loss_type == "or"){
      term_outcome <- sapply(which(idx), function(i){ smoothers[[subject_idx[i]]](trt_right_clipped[i]) + krige_adjust[subject_idx[i]]})
    # so we can use predict.smooth.spline() to get the derivative 
    } else {
      term_outcome <- rep(0, sum(idx))
    }
    # --- Combine terms ---
    # is it idx not subject[idx], yes, it is idx
    d_loss_vec[idx] <- term_threshold - term_db - term_outcome
  }
  
  # ---------------------- right-clipped additional loss ----------------------
  if (any(index_right_clip)) {
    # d_loss_vec[index_right_clip] <- clip_epsilon * exp(-assigned_trt[index_right_clip] + trt_bounds[2])
    # d_loss_vec[index_right_clip] <- clip_epsilon
    # d_loss_vec[index_right_clip] <- 2 * clip_epsilon * (assigned_trt[index_right_clip] - trt_bounds[2])
    d_loss_vec[index_right_clip] <- clip_epsilon /(log(exp(1) + assigned_trt[index_right_clip] - trt_bounds[2])^2 * (exp(1) + assigned_trt[index_right_clip] - trt_bounds[2])) # make it upper bounded 
  }
  
  if (is.null(d_loss_vec) || any(!is.finite(d_loss_vec))) {
    browser()
  }
  
  return(d_loss_vec)
}

compute_total_least_square_loss_db_smooth <- function(assigned_trt,
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
                                             surrogate_type = c("Gaussian", "Box")) {
  
  surrogate_type <- match.arg(surrogate_type)
  n <- length(subject_idx)
  loss_vec <- rep(NA, n)
  
  # Identify clipped observations
  index_left_clip  <- assigned_trt < trt_bounds[1]
  index_right_clip <- assigned_trt > trt_bounds[2]
  index_no_clip    <- !index_left_clip & !index_right_clip
  
  # Left and right-clipped treatment values used in integration
  trt_clipped <- pmin(assigned_trt, trt_bounds[2])
  trt_clipped <- pmax(assigned_trt, trt_bounds[1])
  
  
  # ---------------------- clipped loss -------------------------
  
  term_or <- smoothers[[i]](trt_clipped[i])
  
 
  term_db <- (outcome_resid[subject_idx] / propensity_est[subject_idx]) *
              gaussian_kernel(trt_clipped,
                             T[subject_idx],
                             sigma = kernel_bw)

  
  (threshold_val - term_or - term_db)^2
  
  # ---------------------- left-clipped additional loss ----------------------
  if (any(index_left_clip)) {
    loss_vec[index_left_clip] <- loss_vec[index_left_clip] + clip_epsilon - clip_epsilon *
      exp(assigned_trt[index_left_clip] - trt_bounds[1])
  }
  
  # ---------------------- right-clipped additional loss ----------------------
  if (any(index_right_clip)) {
    loss_vec[index_right_clip] <- loss_vec[index_right_clip] +
      clip_epsilon - clip_epsilon *
      exp(-assigned_trt[index_right_clip] + trt_bounds[2])
  }
  
  return(loss_vec)
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# compute_total_loss_db_sum_smooth
# ------------------------------------------------------------------------------
# Purpose:
#   Sum the total loss (including DB term) over all observations
# ------------------------------------------------------------------------------
compute_total_loss_db_sum_smooth <- function(assigned_trt,
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
                                             surrogate_type = c("Gaussian", "Box"),
                                             loss_type =c ("LS","Integral")) {
  surrogate_type <- match.arg(surrogate_type)
  loss_type <- match.arg(loss_type)
  if (loss_type == "Integral"){
    return(sum(
    compute_total_loss_db_smooth(assigned_trt,
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
                                 surrogate_type)
  ))}
  if (loss_type == "LS"){
    return(sum(
      compute_total_least_square_loss_db_smooth(assigned_trt,
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
                                   surrogate_type)
    ))
  }
    
}

compute_total_loss_sum_smooth <- function(assigned_trt,
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
                                           surrogate_type = c("Gaussian", "Box"),
                                           loss_type =c("db","or","ipw")) {
  surrogate_type <- match.arg(surrogate_type)
  loss_type <- match.arg(loss_type)
  
  return(sum(
    compute_total_loss_smooth(assigned_trt = assigned_trt,
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
                             surrogate_type = surrogate_type,
                             loss_type = loss_type)
  ))
  
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Loss Function: Parametric Method
# ---------------------------------------------------------------------------- #

# parametric loss minimization for outcome regression
total_loss_true_sum_parametric_or_storage_smooth =  function(coefs, param_model_initial, data, cum_smoothers, threshold, krige_values, trt_range, epsilon){
  # this is the predicted assigned treatment levels
  trt_assigned = c(model.matrix(param_model_initial) %*% coefs)
  total_loss_sum_or_storage_smooth(trt_assigned = trt_assigned, data = data, cum_smoothers = cum_smoothers, threshold = threshold, krige_values = krige_values, trt_range = trt_range, epsilon = epsilon)
}


# parametric loss minimization for doubly robust regression
total_loss_sum_parametric_db_storage_smooth =  function(coefs, param_model_initial,data,krige_values, resids, prop_estimation,cum_smoothers,trt_range, threshold, bandwidth = bandwidth){
  # this is the predicted assigned treatment levels
  trt_assigned = model.matrix(param_model_initial) %*% coefs
  total_loss_db_sum_storage_smooth(trt_assigned = trt_assigned, col_idx = 1:nrow(data), data = data, krige_values = krige_values, resids = resids, prop_estimation = prop_estimation, cum_smoothers = cum_smoothers, trt_range = trt_range, threshold = threshold, bandwidth = bandwidth)
}

