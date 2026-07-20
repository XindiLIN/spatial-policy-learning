# ============================================================
# cv_hyperparameter.R
# ============================================================
# K-fold cross-validation to select optimal hyperparameters
# for a fixed area and threshold:
#
#   m              : kernel bandwidth exponent  (gamma = 2^m * median_dist)
#   lambda         : RKHS regularisation penalty
#   kernel_bw_mult : loss-function bandwidth multiplier,
#                    kernel_bw = bw_mult / n  (proportional to 1/n)
#
# Speed tips (see bottom of file for full discussion):
#   - Reduce cv_folds to 3, or switch to a single hold-out split
#   - Reduce optim_maxit_cv (e.g. 200)
#   - Narrow the grid once you have a rough best region
#   - Only tune m first (fix lambda and kernel_bw_mult), then tune the rest
# ============================================================


# ============================================================
# CONFIGURATION  — edit only this section
# ============================================================

cv_area          <- "Northwest"    # one area at a time
cv_threshold     <- log(5)         # e.g. log(2), log(5), log(10)
cv_seed          <- 42             # fixed seed for fold assignment
cv_folds         <- 5              # number of CV folds
optim_maxit_cv   <- 300            # maxit for optim inside each CV fold
                                   # (lower = faster, less precise)
output_dir       <- "."
plss_path        <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/plss_covariates.csv"
year_filter      <- 2020           # set NULL to use all years

# Hyperparameter grid
m_grid          <- c(-3, -2, -1.5, -1, 0, 0.5)
lambda_grid     <- c(0.1, 0.25, 0.5, 1.0)
# kernel_bw = bw_mult / n_fold_train  (proportional to 1/n)
bw_mult_grid    <- c(0.5, 1, 2, 5, 10)


# ============================================================
# SOURCES AND PACKAGES
# ============================================================

source("data_generation.R")
source("funcs_kriging.R")
source("loss_functions.R")
source("loss_functions_RKHS.R")
source("miscellaneous.r")


# ============================================================
# DATA AND SHARED OBJECTS
# ============================================================

data_all <- load_nitrate_data(zero_inflated = FALSE)
if (!is.null(year_filter))
  data_all <- data_all[data_all$SampleYear >= year_filter, ]

area_key        <- gsub(" ", "_", tolower(cv_area))
data_area       <- data_all[data_all$area == cv_area, ]
data_area_split <- split_nitrate_data(data_area)
data_train      <- data_area_split$data
n_train         <- nrow(data_train)

RKHS_covariate_names <- c("StaticLevel", "crop_type_combine", "drainagecl",
                           "precipitation", "cafolog")

# ---- Outcome regression (load cache if available) ----
outcome_reg_path <- file.path(output_dir, sprintf("outcome_regression_%s.rds", area_key))
if (file.exists(outcome_reg_path)) {
  cat("Loading cached outcome regression...\n")
  outcome_reg <- readRDS(outcome_reg_path)
} else {
  cat("Fitting outcome regression...\n")
  outcome_reg <- outcome_regression_SVM(
    data      = data_train,
    data_test = data_area_split$data_test,
    tunning   = FALSE
  )
  saveRDS(outcome_reg, outcome_reg_path)
}

depth_range <- c(min(data_train$logWellDepth), max(data_train$logWellDepth))

# ---- Smoothers (load cache if available) ----
smoothers_path <- file.path(output_dir, sprintf("smoothers_%s.rds", area_key))
if (file.exists(smoothers_path)) {
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

# ---- CBPS propensity weights ----
weights_cbps <- weightit(
  logWellDepth ~ StaticLevel + crop_type_combine + drainagecl +
                 precipitation + cafolog,
  data = data_train, method = "cbps"
)

# ---- Recipe fitted on full training data (reused across folds) ----
rec <- recipe(~ ., data = data_train[, RKHS_covariate_names]) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep()

# ---- Residuals for loss function ----
outcome_resid_full <- data_train$logconcentration_plus_median -
                      outcome_reg$pred -
                      outcome_reg$krige_values


# ============================================================
# CV FOLD ASSIGNMENT
# ============================================================

set.seed(cv_seed)
fold_ids <- sample(rep(1:cv_folds, length.out = n_train))


# ============================================================
# HYPERPARAMETER GRID
# ============================================================

params_grid        <- expand.grid(m = m_grid, lambda = lambda_grid, bw_mult = bw_mult_grid)
params_grid$cv_mcc <- NA_real_
params_grid$cv_acc <- NA_real_

n_combos <- nrow(params_grid)
cat(sprintf(
  "\nCV for area='%s', threshold=log(%g)\n%d combos x %d folds = %d optim runs\n\n",
  cv_area, round(exp(cv_threshold)), n_combos, cv_folds, n_combos * cv_folds
))


# ============================================================
# MAIN CV LOOP
# ============================================================

for (g in seq_len(n_combos)) {

  m_val      <- params_grid$m[g]
  lambda_val <- params_grid$lambda[g]
  bw_mult    <- params_grid$bw_mult[g]

  fold_mcc <- numeric(cv_folds)
  fold_acc <- numeric(cv_folds)

  for (fold in seq_len(cv_folds)) {

    val_idx   <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    n_fold    <- length(train_idx)

    # kernel_bw proportional to 1/n
    kernel_bw_fold <- bw_mult / n_fold

    data_fold_tr  <- data_train[train_idx, ]
    data_fold_val <- data_train[val_idx, ]

    # ---- Kernel design matrices for this fold ----
    kdm_tr  <- bake(rec, new_data = data_fold_tr[,  RKHS_covariate_names])
    kdm_tr  <- cbind(kdm_tr,  U = outcome_reg$krige_values[train_idx])
    kdm_val <- bake(rec, new_data = data_fold_val[, RKHS_covariate_names])
    kdm_val <- cbind(kdm_val, U = outcome_reg$krige_values[val_idx])

    fold_means <- colMeans(kdm_tr)
    fold_sds   <- apply(kdm_tr, 2, sd)
    kdm_tr     <- scale(kdm_tr)
    kdm_val    <- scale(kdm_val, center = fold_means, scale = fold_sds)

    gamma  <- 2^(m_val) * median(fields::rdist(kdm_tr))
    rbf    <- rbfdot(sigma = 1 / (gamma^2))
    K_tr   <- cbind(1, kernelMatrix(rbf, kdm_tr))
    K_val  <- cbind(1, kernelMatrix(rbf, kdm_val, kdm_tr))

    # ---- Initialization via indirect policy (faster than find_best_policy) ----
    indirect_init <- compute_indirect_policy(
      smoothers     = smoothers$smoothers_RKHS[train_idx],
      trt_bounds    = depth_range,
      threshold_val = cv_threshold,
      krige_adjust  = outcome_reg$krige_values[train_idx],
      spatial       = TRUE
    )
    init_glmnet <- glmnet(x = K_tr[, -1], y = indirect_init, alpha = 0, lambda = 0)
    params_init <- as.numeric(coef(init_glmnet))

    # ---- Optimisation on training fold ----
    fold_optim <- tryCatch(
      optim(
        par              = params_init,
        fn               = compute_total_loss_smooth_RKHS,
        gr               = d_compute_total_loss_smooth_RKHS,
        K                = K_tr,
        T                = data_fold_tr$logWellDepth,
        krige_adjust     = outcome_reg$krige_values[train_idx],
        outcome_resid    = outcome_resid_full[train_idx],
        propensity_est   = 1 / weights_cbps$weights[train_idx],
        lambda           = lambda_val,
        smoothers        = smoothers$smoothers_RKHS[train_idx],
        cumint_smoothers = smoothers$cumint_smoothers_RKHS[train_idx],
        trt_bounds       = depth_range,
        threshold_val    = cv_threshold,
        kernel_bw        = kernel_bw_fold,
        clip_epsilon     = 30,
        surrogate_type   = "Gaussian",
        loss_type        = "db",
        method           = "L-BFGS-B",
        control          = list(maxit = optim_maxit_cv, trace = 0)
      ),
      error = function(e) {
        message(sprintf("  optim failed at combo %d fold %d: %s", g, fold, e$message))
        NULL
      }
    )

    if (is.null(fold_optim)) {
      fold_mcc[fold] <- NA_real_
      fold_acc[fold] <- NA_real_
      next
    }

    # ---- Evaluate on validation fold ----
    val_policy <- as.numeric(K_val %*% fold_optim$par)
    result <- calculate_acc_mcc(
      T             = data_fold_val$logWellDepth,
      Y             = data_fold_val$logconcentration_plus_median,
      policy        = val_policy,
      threshold_val = cv_threshold
    )
    fold_mcc[fold] <- result$mcc
    fold_acc[fold] <- result$acc

  }  # end fold loop

  params_grid$cv_mcc[g] <- mean(fold_mcc, na.rm = TRUE)
  params_grid$cv_acc[g] <- mean(fold_acc, na.rm = TRUE)

  cat(sprintf(
    "Combo %3d/%d | m=%5.1f  lambda=%.2f  bw_mult=%4.1f | CV MCC=%6.4f  ACC=%6.4f\n",
    g, n_combos, m_val, lambda_val, bw_mult,
    params_grid$cv_mcc[g], params_grid$cv_acc[g]
  ))

}  # end grid loop


# ============================================================
# RESULTS
# ============================================================

threshold_label <- paste0("log", round(exp(cv_threshold), 0))
cv_rds_path <- file.path(output_dir,
  sprintf("cv_results_%s_%s.rds", area_key, threshold_label))
saveRDS(params_grid, cv_rds_path)
cat(sprintf("\nFull CV grid saved to: %s\n", cv_rds_path))

# ---- Best by CV MCC ----
best_idx    <- which.max(params_grid$cv_mcc)
best_params <- params_grid[best_idx, ]

cat("\n========== Best Hyperparameters (max CV MCC) ==========\n")
cat(sprintf("  m              = %.1f\n",  best_params$m))
cat(sprintf("  lambda         = %.2f\n",  best_params$lambda))
cat(sprintf("  bw_mult        = %.1f  =>  kernel_bw = %.1f / n\n",
            best_params$bw_mult, best_params$bw_mult))
cat(sprintf("  CV MCC         = %.4f\n",  best_params$cv_mcc))
cat(sprintf("  CV Accuracy    = %.4f\n",  best_params$cv_acc))
cat("=======================================================\n")


# ============================================================
# VISUALISATION
# ============================================================

# 1. Heatmap of CV MCC for (m, lambda), averaged over bw_mult
cv_ml <- params_grid %>%
  group_by(m, lambda) %>%
  summarise(mean_mcc = mean(cv_mcc, na.rm = TRUE), .groups = "drop")

p1 <- ggplot(cv_ml, aes(x = factor(lambda), y = factor(m), fill = mean_mcc)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_mcc, 3)), size = 3) +
  scale_fill_viridis_c(option = "plasma", name = "Mean CV MCC") +
  labs(
    title = sprintf("CV MCC — %s, threshold = log(%g)", cv_area, round(exp(cv_threshold))),
    subtitle = "Averaged over bw_mult",
    x = "lambda", y = "m"
  ) +
  theme_minimal()
print(p1)

# 2. MCC vs bw_mult for the best (m, lambda)
cv_bw <- params_grid %>%
  filter(m == best_params$m, lambda == best_params$lambda)

p2 <- ggplot(cv_bw, aes(x = bw_mult, y = cv_mcc)) +
  geom_line() +
  geom_point(size = 2) +
  geom_point(data = best_params, color = "red", size = 4) +
  labs(
    title   = sprintf("CV MCC vs kernel_bw multiplier  (m=%.1f, lambda=%.2f)",
                      best_params$m, best_params$lambda),
    x       = "bw_mult  (kernel_bw = bw_mult / n)",
    y       = "CV MCC"
  ) +
  theme_minimal()
print(p2)


# ============================================================
# IF CV IS TOO SLOW — OPTIONS (in order of impact)
# ============================================================
#
# 1. SWITCH TO A SINGLE HOLD-OUT SPLIT (biggest speedup)
#    Instead of cv_folds = 5 (5x optim runs per combo), use
#    cv_folds = 1 with a fixed 80/20 split.
#    Change:  fold_ids <- rep(1, n_train)
#             fold_ids[sample(n_train, floor(0.2 * n_train))] <- 2
#    This gives a 5x speedup at the cost of higher variance in the estimate.
#
# 2. REDUCE optim_maxit_cv
#    For hyperparameter search you only need a rough solution.
#    optim_maxit_cv = 100-200 is usually sufficient.
#    The final model uses maxit = 1000.
#
# 3. SEQUENTIAL TUNING — tune one parameter at a time
#    Step 1: fix lambda=0.25, bw_mult=1, sweep over m.
#    Step 2: fix best m, bw_mult=1, sweep over lambda.
#    Step 3: fix best m and lambda, sweep over bw_mult.
#    Reduces grid from |m| x |lambda| x |bw_mult| runs
#    to |m| + |lambda| + |bw_mult| runs.
#
# 4. PARALLELIZE ACROSS COMBOS
#    Replace the outer for-loop with:
#      library(parallel)
#      results <- mclapply(seq_len(n_combos), function(g) { ... },
#                          mc.cores = detectCores() - 1)
#    Each combo is independent so this is embarrassingly parallel.
#
# 5. REUSE SMOOTHERS AND OUTCOME REGRESSION ACROSS FOLDS
#    Already implemented above — both are precomputed on the full
#    training set and sub-indexed per fold, avoiding the most
#    expensive recomputations.
#
# 6. RANDOM SEARCH INSTEAD OF GRID SEARCH
#    Sample n_random = 20-30 random combos from the grid
#    rather than exhaustively evaluating all combinations:
#      set.seed(cv_seed)
#      grid_idx <- sample(n_combos, min(30, n_combos))
#      params_grid <- params_grid[grid_idx, ]
