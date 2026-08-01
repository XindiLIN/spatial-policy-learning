# ============================================================
# region_pipeline_funcs.R
# ============================================================
# Per-region building blocks for the minimum-well-depth policy
# pipeline: the direct RKHS/DC method and the non-spatial indirect
# method (the spatial-indirect variant has been dropped). Every
# function below operates on ONE region at a time and takes that
# region's data/objects as arguments -- no region name is hardcoded
# anywhere in this file. A separate driver script is expected to
# loop over regions, call these in order, and combine the resulting
# per-region sf objects into the final map.
#
# Pipeline order for one region:
#   (i)   fit_region_base()               -> outcome_reg, weights_cbps, smoothers
#   (ii)  get_indirect_nonspatial_policy() -> indirect (non-spatial) policy, train + test
#   (v)   choose_direct_hyperparameters_cv() -> CV-selected m, lambda, kernel_bw
#   (iii) prep_direct_method_inputs()     -> kernel design matrices, K/K_test/K_plss, DC warm start
#   (iv)  fit_direct_dc_policy()          -> DC algorithm fit, train/test metrics, PLSS map layer
# ============================================================

# Resolve this file's own location so its dependencies can be sourced
# regardless of the caller's working directory -- this file is always
# meant to be source()d (never Rscript'd directly), and different callers
# use different cwd conventions (Rscript from repo root, RStudio's
# interactive-chunk execution defaulting to the calling .Rmd's own
# directory, knitting with a custom root.dir, etc). Falls back to
# assuming cwd is already the repo root if the file wasn't reached via
# source() (e.g. pasted directly into the console).
.region_pipeline_funcs_path <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    if (identical(sys.function(i), base::source)) {
      ofile <- evalq(ofile, envir = sys.frame(i))
      if (!is.null(ofile)) return(normalizePath(ofile))
    }
  }
  NULL
})
.repo_root <- if (!is.null(.region_pipeline_funcs_path)) {
  dirname(dirname(.region_pipeline_funcs_path))
} else {
  getwd()
}

source(file.path(.repo_root, "functions", "miscellaneous.r"))
source(file.path(.repo_root, "functions", "loss_functions_RKHS.R"))

library(dplyr)
library(sf)
library(recipes)
library(WeightIt)
library(glmnet)
library(kernlab)
library(parallel)

RKHS_COVARIATE_NAMES_DEFAULT <- c("StaticLevel", "crop_type_combine", "drainagecl",
                                   "precipitation", "cafolog")


# ============================================================
# (i) Outcome regression + propensity weights + RKHS smoothers
# ============================================================
# Fits (or loads cached) the SVM+kriging outcome-regression model, the
# CBPS propensity weights used by the direct method, and the RKHS
# counterfactual smoothers used by both methods. Cached per region to
# `output_dir` so a rerun skips the slow SVM/kriging and smoother
# fitting steps entirely.
fit_region_base <- function(data_split, area_key, output_dir,
                             RKHS_covariate_names = RKHS_COVARIATE_NAMES_DEFAULT,
                             use_cache = TRUE) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  outcome_reg_path <- file.path(output_dir, sprintf("outcome_regression_%s.rds", area_key))
  if (use_cache && file.exists(outcome_reg_path)) {
    cat("Loading cached outcome regression...\n")
    outcome_reg <- readRDS(outcome_reg_path)
  } else {
    cat("Fitting outcome regression...\n")
    outcome_reg <- outcome_regression_SVM(
      data      = data_split$data,
      data_test = data_split$data_test,
      tunning   = FALSE
    )
    saveRDS(outcome_reg, outcome_reg_path)
  }

  cat(sprintf(
    "RMSE  train: %.4f   test: %.4f\n",
    Metrics::rmse(data_split$data$logconcentration_plus_median,
                  outcome_reg$pred + outcome_reg$krige_values),
    Metrics::rmse(data_split$data_test$logconcentration_plus_median,
                  outcome_reg$pred_test + outcome_reg$krige_values_test)
  ))

  depth_range <- c(min(data_split$data$logWellDepth), max(data_split$data$logWellDepth))

  weights_cbps <- weightit(
    logWellDepth ~ StaticLevel + crop_type_combine + drainagecl + precipitation + cafolog,
    data = data_split$data, method = "cbps"
  )

  # Silverman's rule-of-thumb bandwidth -- used as the direct method's
  # DC-algorithm warm-start bandwidth in prep_direct_method_inputs().
  # This is NOT the RKHS Gram-matrix bandwidth (gamma, governed by m)
  # and NOT the loss function's kernel_bw used in the final DC fit --
  # both of those come from choose_direct_hyperparameters_cv() instead.
  kernel_bw_silverman <- 1.06 * sd(data_split$data$logWellDepth) *
    nrow(data_split$data)^(-1 / 5)

  smoothers_path <- file.path(output_dir, sprintf("smoothers_%s.rds", area_key))
  if (use_cache && file.exists(smoothers_path)) {
    cat("Loading cached smoothers...\n")
    smoothers <- readRDS(smoothers_path)
  } else {
    cat("Computing smoothers (slow)...\n")
    smoothers <- get_smoothers_RKHS(
      design_matrix      = outcome_reg$design_matrix,
      design_matrix_test = outcome_reg$design_matrix_test,
      svm_auto           = outcome_reg$svm,
      depth_range        = depth_range
    )
    saveRDS(smoothers, smoothers_path)
  }

  list(
    outcome_reg          = outcome_reg,
    weights_cbps         = weights_cbps,
    smoothers            = smoothers,
    depth_range          = depth_range,
    kernel_bw_silverman  = kernel_bw_silverman
  )
}


# ============================================================
# (ii) Indirect non-spatial policy (train + test)
# ============================================================
# spatial = FALSE deliberately omits the kriging/spatial-correction
# term from the plug-in outcome model -- this is the "non-spatial
# indirect" comparison method. The spatial-indirect variant is not
# computed here (dropped per current analysis plan).
get_indirect_nonspatial_policy <- function(region_base, data_split, threshold_val) {

  smoothers <- region_base$smoothers

  policy_train <- compute_indirect_policy(
    smoothers     = smoothers$smoothers_RKHS,
    trt_bounds    = region_base$depth_range,
    threshold_val = threshold_val,
    spatial       = FALSE
  )
  policy_test <- compute_indirect_policy(
    smoothers     = smoothers$smoothers_test_RKHS,
    trt_bounds    = region_base$depth_range,
    threshold_val = threshold_val,
    spatial       = FALSE
  )

  metrics_train <- calculate_acc_mcc_two_sided_f1(
    T = data_split$data$logWellDepth, Y = data_split$data$logconcentration_plus_median,
    policy = policy_train, threshold_val = threshold_val
  )
  metrics_test <- calculate_acc_mcc_two_sided_f1(
    T = data_split$data_test$logWellDepth, Y = data_split$data_test$logconcentration_plus_median,
    policy = policy_test, threshold_val = threshold_val
  )

  cat(sprintf("Indirect non-spatial -- Train: acc=%.4f mcc=%.4f f1=%.4f\n",
              metrics_train$acc, metrics_train$mcc, metrics_train$two_sided_f1))
  cat(sprintf("Indirect non-spatial -- Test:  acc=%.4f mcc=%.4f f1=%.4f\n",
              metrics_test$acc, metrics_test$mcc, metrics_test$two_sided_f1))

  list(policy_train = policy_train, policy_test = policy_test,
       metrics_train = metrics_train, metrics_test = metrics_test)
}


# ============================================================
# Shared plumbing (used by both (iii) and (v), not one of the 5
# numbered pieces): fit the RKHS recipe on training covariates and
# bake+scale the train/test kernel design matrices. Independent of
# the kernel bandwidth exponent m, so it doesn't belong to either the
# CV loop or the final DC prep alone -- both need the exact same
# recipe/scaling to stay consistent.
# ============================================================
build_kernel_design_matrix <- function(data_split, region_base,
                                        RKHS_covariate_names = RKHS_COVARIATE_NAMES_DEFAULT) {

  # step_dummy() below learns its dummy columns from TRAINING levels only.
  # A category present in training but absent from test is harmless (the
  # dummy column still exists, just constant zero in kdm_test). The
  # reverse is not: a category in test that step_dummy() never saw during
  # training has no learned dummy column, so bake() silently fills that
  # row with NA for every dummy derived from that variable (with only a
  # warning, not an error) -- and that NA then propagates silently through
  # scale()/kernelMatrix() into K_test. Check for this up front and fail
  # loudly instead.
  nominal_vars <- RKHS_covariate_names[
    !vapply(data_split$data[, RKHS_covariate_names], is.numeric, logical(1))
  ]
  for (v in nominal_vars) {
    train_levels <- unique(as.character(data_split$data[[v]]))
    test_levels  <- unique(as.character(data_split$data_test[[v]]))
    novel <- setdiff(test_levels, train_levels)
    if (length(novel) > 0) {
      stop(sprintf(
        "build_kernel_design_matrix(): test set has categor%s in `%s` not seen in training: %s. step_dummy() would silently NA those test rows (see ?recipes::step_novel). Re-check the train/test split for this region.",
        if (length(novel) > 1) "ies" else "y", v, paste(novel, collapse = ", ")
      ))
    }
  }

  rec <- recipe(~ ., data = data_split$data[, RKHS_covariate_names]) %>%
    step_dummy(all_nominal_predictors()) %>%
    prep()

  kdm      <- bake(rec, new_data = data_split$data[, RKHS_covariate_names])
  kdm      <- cbind(kdm, U = region_base$outcome_reg$krige_values)
  kdm_test <- bake(rec, new_data = data_split$data_test[, RKHS_covariate_names])
  kdm_test <- cbind(kdm_test, U = region_base$outcome_reg$krige_values_test)

  train_means <- colMeans(kdm)
  train_sds   <- apply(kdm, 2, sd)

  list(
    rec         = rec,
    kdm         = scale(kdm),
    kdm_test    = scale(kdm_test, center = train_means, scale = train_sds),
    train_means = train_means,
    train_sds   = train_sds
  )
}

# PLSS-only counterfactual smoothers -- mirrors get_smoothers_RKHS()'s
# PLSS branch without recomputing the (already-cached) train/test
# smoothers. Needs the OR design matrix purely for its column order.
get_smoothers_plss_RKHS <- function(design_matrix_plss, design_matrix_train, svm_auto,
                                     depth_range, treatment_step = 0.02) {

  design_matrix_plss$logWellDepth <- 0
  design_matrix_plss <- design_matrix_plss[, colnames(design_matrix_train)]

  cf_predictions_plss <- precompute_cf_predictions(
    data_obs = design_matrix_plss, fit_outcome = svm_auto,
    treatment_range = depth_range, treatment_step = treatment_step,
    treatment_name = "logWellDepth"
  )
  cf_predictions_plss <- isotomic_correction(
    treatment_grid = cf_predictions_plss$treatment_values,
    pred_matrix_cf = cf_predictions_plss$pred_matrix
  )
  make_prediction_smoothers(
    treatment_grid = cf_predictions_plss$treatment_values,
    pred_matrix_cf = cf_predictions_plss$pred_matrix_cf,
    smooth_method  = "smooth.spline",
    spar           = NULL
  )
}


# ============================================================
# (iii) Preparations before the direct method's DC algorithm
# ============================================================
# Builds the Gram matrices (K, K_test, K_plss) at the chosen kernel
# bandwidth exponent `m`, subsets/scores the PLSS prediction grid for
# this region, and computes the DC algorithm's warm-start parameters
# (per-observation direct-policy initialization, projected into the
# RKHS via ridge regression). `m` is expected to come from
# choose_direct_hyperparameters_cv() (piece v).
prep_direct_method_inputs <- function(data_split, region_base, kdm_info, m, kernel_bw,
                                       plss_covariates_all, area, area_key, output_dir,
                                       plss_county_override = NULL,
                                       plss_crop_exclude = NULL,
                                       RKHS_covariate_names = RKHS_COVARIATE_NAMES_DEFAULT,
                                       threshold_val, clip_epsilon_init = 20,
                                       mc.cores = 1, use_cache = TRUE) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  outcome_reg <- region_base$outcome_reg
  depth_range <- region_base$depth_range
  kdm      <- kdm_info$kdm
  kdm_test <- kdm_info$kdm_test
  rec      <- kdm_info$rec

  gamma <- 2^m * median(fields::rdist(kdm))
  rbf   <- rbfdot(sigma = 1 / (gamma^2))

  K      <- cbind(1, kernelMatrix(rbf, kdm))
  K_test <- cbind(1, kernelMatrix(rbf, kdm_test, kdm))

  # ---- PLSS covariates for this area ----
  if (!is.null(plss_county_override)) {
    plss <- plss_covariates_all[plss_covariates_all$County %in% plss_county_override, ]
  } else {
    plss <- plss_covariates_all[plss_covariates_all$County %in% unique(data_split$data$County), ]
  }
  plss <- na.omit(plss)

  if (!is.null(plss_crop_exclude)) {
    plss <- plss[!plss$crop_type_combine %in% plss_crop_exclude, ]
  } else {
    plss <- plss[plss$crop_type_combine %in% unique(data_split$data$crop_type_combine), ]
  }

  # Spatial random effect (kriging) at PLSS locations
  plss_krige <- GpGp::predictions(
    fit       = outcome_reg$gpfit,
    locs_pred = plss[, c("longitude", "latitude")],
    X_pred    = rep(1, nrow(plss))
  )
  plss$U <- plss_krige

  plss_kdm <- bake(rec, new_data = plss[, RKHS_covariate_names])
  plss_kdm <- cbind(plss_kdm, U = plss_krige)
  plss_kdm <- scale(plss_kdm, center = kdm_info$train_means, scale = kdm_info$train_sds)
  K_plss   <- cbind(1, kernelMatrix(rbf, plss_kdm, kdm))

  smoothers_plss_RKHS <- get_smoothers_plss_RKHS(
    design_matrix_plss  = bake(rec, new_data = plss[, RKHS_covariate_names]),
    design_matrix_train = outcome_reg$design_matrix,
    svm_auto            = outcome_reg$svm,
    depth_range         = depth_range
  )

  # ---- DC algorithm warm start: per-observation direct policy, then
  # projected into the RKHS via ridge regression (alpha = 0, lambda = 0).
  # find_best_policy()'s kernel_bw is the SAME kind of smoothing bandwidth
  # as the DC/CV loss's own kernel_bw (both smooth the same Gaussian
  # surrogate loss), so the warm start is computed at the CV-chosen
  # kernel_bw here rather than an unrelated fixed heuristic. It still does
  # NOT depend on m or lambda (find_best_policy optimizes per observation
  # directly in treatment space, not the RKHS coefficient space), so it's
  # cached per (region, threshold, kernel_bw) rather than per full grid
  # point -- and since this function only runs once per FINAL chosen
  # hyperparameter set (not once per CV grid combo), that's still cheap.
  threshold_label   <- paste0("log", round(exp(threshold_val), 0))
  direct_init_path  <- file.path(output_dir,
    sprintf("direct_init_%s_%s_bw%.4f.rds", area_key, threshold_label, kernel_bw))
  if (use_cache && file.exists(direct_init_path)) {
    cat("Loading cached direct_init warm start...\n")
    direct_init <- readRDS(direct_init_path)
  } else {
    cat("  Initializing direct policy...\n")
    init_vals <- c(depth_range[1] - 1, mean(depth_range), depth_range[2] + 1)
    direct_init <- mcmapply(
      FUN      = find_best_policy,
      i        = seq_len(nrow(data_split$data)),
      MoreArgs = list(
        initial_value_set         = init_vals,
        kernel_bw                 = kernel_bw,
        outcome_regression_object = outcome_reg,
        data                      = data_split$data,
        weights                   = region_base$weights_cbps$weights,
        smoothers                 = region_base$smoothers,
        depth_range               = depth_range,
        threshold_val             = threshold_val,
        clip_epsilon              = clip_epsilon_init
      ),
      SIMPLIFY = TRUE,
      mc.cores = mc.cores
    )
    saveRDS(direct_init, direct_init_path)
  }
  cat("  Init acc/mcc: ")
  print_acc_mcc(data_split$data$logWellDepth, data_split$data$logconcentration_plus_median,
                direct_init, threshold_val)

  init_glmnet  <- glmnet(x = K[, -1], y = direct_init, alpha = 0, lambda = 0)
  params_init  <- as.numeric(coef(init_glmnet))

  # Warm-start-only (pre-DC) train/test policy -- a cheap diagnostic for how
  # much the (slow) DC algorithm actually improves over the warm start alone.
  warmstart_policy_train <- as.numeric(K %*% params_init)
  warmstart_policy_test  <- as.numeric(K_test %*% params_init)
  warmstart_metrics_train <- calculate_acc_mcc_two_sided_f1(
    T = data_split$data$logWellDepth, Y = data_split$data$logconcentration_plus_median,
    policy = warmstart_policy_train, threshold_val = threshold_val
  )
  warmstart_metrics_test <- calculate_acc_mcc_two_sided_f1(
    T = data_split$data_test$logWellDepth, Y = data_split$data_test$logconcentration_plus_median,
    policy = warmstart_policy_test, threshold_val = threshold_val
  )
  cat(sprintf("  Warm start (pre-DC) -- Train: acc=%.4f mcc=%.4f f1=%.4f\n",
              warmstart_metrics_train$acc, warmstart_metrics_train$mcc, warmstart_metrics_train$two_sided_f1))
  cat(sprintf("  Warm start (pre-DC) -- Test:  acc=%.4f mcc=%.4f f1=%.4f\n",
              warmstart_metrics_test$acc, warmstart_metrics_test$mcc, warmstart_metrics_test$two_sided_f1))

  list(
    K = K, K_test = K_test, K_plss = K_plss,
    rbf = rbf, gamma = gamma,
    plss = plss, smoothers_plss_RKHS = smoothers_plss_RKHS,
    direct_init = direct_init, params_init = params_init,
    warmstart_policy_train = warmstart_policy_train, warmstart_policy_test = warmstart_policy_test,
    warmstart_metrics_train = warmstart_metrics_train, warmstart_metrics_test = warmstart_metrics_test
  )
}


# ============================================================
# (iv) Run the DC algorithm for the direct method
# ============================================================
# Fits the direct policy via the DC algorithm (warm-started from
# prep$params_init), evaluates it on train/test, and scores it on the
# PLSS grid to produce the sf object used for mapping.
fit_direct_dc_policy <- function(data_split, region_base, prep,
                                  lambda, kernel_bw, threshold_val,
                                  clip_epsilon = 30, clip_epsilon_bar = 50,
                                  tol = 0.001, max_iter = 100,
                                  optim_maxit = 300, optim_trace = 0) {

  outcome_reg   <- region_base$outcome_reg
  outcome_resid <- data_split$data$logconcentration_plus_median -
    outcome_reg$pred - outcome_reg$krige_values

  # optim_trace defaults to 0 (silent) here, unlike run_dc_algorithm()'s own
  # default of 3 -- the outer "Iter: k, Change: ..., Objective: ..." line
  # printed once per DC iteration (below) is enough signal for routine runs
  # and Bouchet logs; bump optim_trace back up to 3 for a one-off deep dive
  # into a single inner L-BFGS-B solve.
  dc_result <- run_dc_algorithm(
    params_initial   = prep$params_init,
    k_matrix         = prep$K,
    data             = data_split$data,
    T                = data_split$data$logWellDepth,
    krige_adjust     = outcome_reg$krige_values,
    outcome_resid    = outcome_resid,
    propensity_est   = 1 / region_base$weights_cbps$weights,
    lambda           = lambda,
    smoothers        = region_base$smoothers$smoothers_RKHS,
    cumint_smoothers = region_base$smoothers$cumint_smoothers_RKHS,
    trt_bounds       = region_base$depth_range,
    threshold_val    = threshold_val,
    kernel_bw        = kernel_bw,
    clip_epsilon     = clip_epsilon,
    clip_epsilon_bar = clip_epsilon_bar,
    tol              = tol,
    max_iter         = max_iter,
    optim_maxit      = optim_maxit,
    optim_trace      = optim_trace
  )

  metrics_train <- calculate_acc_mcc_two_sided_f1(
    T = data_split$data$logWellDepth, Y = data_split$data$logconcentration_plus_median,
    policy = as.numeric(prep$K %*% dc_result$par), threshold_val = threshold_val
  )
  metrics_test <- calculate_acc_mcc_two_sided_f1(
    T = data_split$data_test$logWellDepth, Y = data_split$data_test$logconcentration_plus_median,
    policy = as.numeric(prep$K_test %*% dc_result$par), threshold_val = threshold_val
  )
  cat(sprintf("Direct DC -- Train: acc=%.4f mcc=%.4f f1=%.4f\n",
              metrics_train$acc, metrics_train$mcc, metrics_train$two_sided_f1))
  cat(sprintf("Direct DC -- Test:  acc=%.4f mcc=%.4f f1=%.4f\n",
              metrics_test$acc, metrics_test$mcc, metrics_test$two_sided_f1))

  plss_dir        <- prep$plss
  plss_dir$policy <- as.numeric(prep$K_plss %*% dc_result$par)
  plss_sf         <- st_as_sf(plss_dir, coords = c("longitude", "latitude"), crs = 4326)

  list(dc_result = dc_result, metrics_train = metrics_train, metrics_test = metrics_test,
       plss_sf = plss_sf)
}


# ============================================================
# (v) Hyperparameter selection via 5-fold cross-validation
# ============================================================
# Structure follows simulation/simulation_nonparametric_funcs.R's
# choose_hyperparameter_cv()/cv_hyperparam_group_for_m(): gamma (and
# hence the full Gram matrix) depends only on m, so it's built ONCE
# per m and every fold's train/holdout matrix is just a slice of it --
# no repeated kernelMatrix() calls per (lambda, kernel_bw) combo.
#
# Deviation from the simulation reference: the per-fold DC warm start
# uses compute_indirect_policy() (closed-form root-finding) rather
# than a per-observation optim() -- the same speed tradeoff already
# made in the existing nitrate_data_analysis/cv_hyperparameter.R,
# which is necessary here since CV runs this many times per region.

assign_folds <- function(n, n_folds = 5, seed_value = 42) {
  set.seed(seed_value)
  sample(rep(1:n_folds, length.out = n))
}

fit_one_fold_direct <- function(k_train, k_holdout, data_fold_train, T_holdout, Y_holdout,
                                 krige_ft, weights_ft, resid_ft, smoothers_ft, cumint_smoothers_ft,
                                 lambda, kernel_bw, threshold_val, depth_range,
                                 clip_epsilon = 30, maxit = 350) {

  indirect_init <- compute_indirect_policy(
    smoothers = smoothers_ft, trt_bounds = depth_range,
    threshold_val = threshold_val, krige_adjust = krige_ft, spatial = TRUE
  )
  init_glmnet <- glmnet(x = k_train[, -1], y = indirect_init, alpha = 0, lambda = 0)
  params_init <- as.numeric(coef(init_glmnet))

  fit <- optim(
    par = params_init, fn = compute_total_loss_smooth_RKHS, gr = d_compute_total_loss_smooth_RKHS,
    K = k_train, T = data_fold_train$logWellDepth, krige_adjust = krige_ft,
    outcome_resid = resid_ft, propensity_est = 1 / weights_ft,
    lambda = lambda, smoothers = smoothers_ft, cumint_smoothers = cumint_smoothers_ft,
    trt_bounds = depth_range, threshold_val = threshold_val, kernel_bw = kernel_bw,
    clip_epsilon = clip_epsilon, surrogate_type = "Gaussian", loss_type = "db",
    method = "L-BFGS-B", control = list(maxit = maxit, trace = 0)
  )

  holdout_policy <- as.numeric(k_holdout %*% fit$par)
  calculate_acc_mcc_two_sided_f1(T = T_holdout, Y = Y_holdout, policy = holdout_policy,
                                  threshold_val = threshold_val)
}

cv_hyperparam_group_for_m <- function(m, lambda_lst, kernel_bw_mult_lst, data_train, kernel_design_matrix,
                                       kernel_bw_baseline, fold_id, n_folds = 5,
                                       threshold_val, depth_range, clip_epsilon = 30,
                                       krige_values, weights, outcome_resid,
                                       smoothers_RKHS, cumint_smoothers_RKHS,
                                       progress_label = NULL, maxit = 350) {

  gamma <- 2^m * median(fields::rdist(kernel_design_matrix))
  rbf   <- rbfdot(sigma = 1 / (gamma^2))
  k_raw <- kernelMatrix(rbf, kernel_design_matrix)  # n_train x n_train, no intercept yet

  grid <- expand.grid(lambda = lambda_lst, kernel_bw_mult = kernel_bw_mult_lst)
  grid$acc <- NA_real_; grid$mcc <- NA_real_; grid$two_sided_f1 <- NA_real_

  for (g in seq_len(nrow(grid))) {

    lambda    <- grid$lambda[g]
    kernel_bw <- grid$kernel_bw_mult[g] * kernel_bw_baseline

    fold_metrics <- vector("list", n_folds)
    for (f in seq_len(n_folds)) {
      train_idx   <- which(fold_id != f)
      holdout_idx <- which(fold_id == f)

      k_train   <- cbind(1, k_raw[train_idx, train_idx, drop = FALSE])
      k_holdout <- cbind(1, k_raw[holdout_idx, train_idx, drop = FALSE])

      fold_metrics[[f]] <- fit_one_fold_direct(
        k_train = k_train, k_holdout = k_holdout,
        data_fold_train = data_train[train_idx, ],
        T_holdout = data_train$logWellDepth[holdout_idx],
        Y_holdout = data_train$logconcentration_plus_median[holdout_idx],
        krige_ft = krige_values[train_idx], weights_ft = weights[train_idx],
        resid_ft = outcome_resid[train_idx],
        smoothers_ft = smoothers_RKHS[train_idx], cumint_smoothers_ft = cumint_smoothers_RKHS[train_idx],
        lambda = lambda, kernel_bw = kernel_bw, threshold_val = threshold_val,
        depth_range = depth_range, clip_epsilon = clip_epsilon, maxit = maxit
      )
    }

    # na.rm = TRUE: mcc (and, more rarely, acc/f1) can be genuinely undefined
    # for a single fold if that fold's held-out predictions collapse to one
    # class (e.g. severe class imbalance with a small fold) -- see
    # calculate_acc_mcc_two_sided_f1()'s explicit factor(levels=c(FALSE,TRUE))
    # guard, which now returns NA for that fold instead of erroring. One
    # degenerate fold shouldn't zero out this whole combo's score.
    fold_mcc <- vapply(fold_metrics, `[[`, numeric(1), "mcc")
    grid$acc[g]          <- mean(vapply(fold_metrics, `[[`, numeric(1), "acc"), na.rm = TRUE)
    grid$mcc[g]           <- mean(fold_mcc, na.rm = TRUE)
    grid$two_sided_f1[g] <- mean(vapply(fold_metrics, `[[`, numeric(1), "two_sided_f1"), na.rm = TRUE)

    if (!is.null(progress_label)) {
      n_na_folds <- sum(is.na(fold_mcc))
      na_note <- if (n_na_folds > 0) sprintf("  [%d/%d folds had undefined mcc]", n_na_folds, n_folds) else ""
      cat(sprintf("%s: combo %d/%d done (lambda=%.2f, kernel_bw_mult=%.2f) -> acc=%.4f mcc=%.4f f1=%.4f%s\n",
                  progress_label, g, nrow(grid), lambda, grid$kernel_bw_mult[g],
                  grid$acc[g], grid$mcc[g], grid$two_sided_f1[g], na_note))
    }
  }

  grid$m         <- m
  grid$kernel_bw <- grid$kernel_bw_mult * kernel_bw_baseline
  grid[, c("m", "lambda", "kernel_bw_mult", "kernel_bw", "acc", "mcc", "two_sided_f1")]
}

# Top-level CV driver for one region. `region_base` and `kdm_info` come
# from fit_region_base() and build_kernel_design_matrix() respectively.
# Set mc.cores > 1 to run m-groups in parallel (safe: each m-group's
# Gram matrix is independent of the others).
choose_direct_hyperparameters_cv <- function(data_split, region_base, kdm_info,
                                              m_lst, lambda_lst, kernel_bw_mult_lst,
                                              threshold_val, n_folds = 5,
                                              seed_value = 42, mc.cores = 1, maxit = 350) {

  data_train <- data_split$data
  fold_id    <- assign_folds(nrow(data_train), n_folds = n_folds, seed_value = seed_value)
  outcome_resid <- data_train$logconcentration_plus_median -
    region_base$outcome_reg$pred - region_base$outcome_reg$krige_values

  results_lst <- parallel::mclapply(m_lst, function(m) {
    cv_hyperparam_group_for_m(
      m = m, lambda_lst = lambda_lst, kernel_bw_mult_lst = kernel_bw_mult_lst,
      data_train = data_train, kernel_design_matrix = kdm_info$kdm,
      kernel_bw_baseline = region_base$kernel_bw_silverman,
      fold_id = fold_id, n_folds = n_folds,
      threshold_val = threshold_val, depth_range = region_base$depth_range,
      krige_values = region_base$outcome_reg$krige_values,
      weights = region_base$weights_cbps$weights,
      outcome_resid = outcome_resid,
      smoothers_RKHS = region_base$smoothers$smoothers_RKHS,
      cumint_smoothers_RKHS = region_base$smoothers$cumint_smoothers_RKHS,
      progress_label = sprintf("m=%.1f", m), maxit = maxit
    )
  }, mc.cores = mc.cores)

  cv_results <- dplyr::bind_rows(results_lst)

  # Guard against a silent downstream crash: if every combo's mcc ended up
  # NA/NaN (e.g. severe class imbalance collapsing every fold, at every
  # combo, to a single predicted class), which.max(cv_results$mcc) in the
  # caller returns integer(0), silently producing an EMPTY `best` row --
  # best$m/best$kernel_bw become numeric(0), and that propagates into a
  # cryptic kernlab error ("replacement has length zero") deep inside
  # prep_direct_method_inputs(), far from the actual cause. Fail loudly
  # here instead, where the cause is still clear.
  if (all(is.na(cv_results$mcc) | is.nan(cv_results$mcc))) {
    stop(sprintf(
      "choose_direct_hyperparameters_cv(): every combo's mcc is NA/undefined (%d combos x %d folds). ",
      nrow(cv_results), n_folds
    ), "Likely severe class imbalance in Y > threshold_val relative to n_folds, causing every fold's ",
    "held-out predictions to collapse to a single class. Try fewer folds, a larger training set, ",
    "or check table(data_split$data$logconcentration_plus_median > threshold_val).")
  }

  cv_results
}
