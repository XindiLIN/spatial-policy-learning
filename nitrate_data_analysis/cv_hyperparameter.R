# ============================================================
# cv_hyperparameter.R
# ============================================================
# Selects (m, lambda, kernel_bw_mult) via k-fold CV, one or more regions
# at a time. Thin wrapper around region_pipeline_funcs.R's fit_region_base()
# (piece i) and cv_hyperparam_group_for_m() (the per-m building block of
# piece v) -- the interactive/exploratory counterpart to
# run_all_regions.R's automatic per-region CV.
#
# Checkpointed and retryable, same pattern as
# simulation/simulation_nonparametric_run_phaseABC.R: each (area, m) task
# writes its own result to nitrate_data_analysis/output/cv_partial/ the
# moment it finishes, so a killed/OOM'd run can be re-run and pick up where
# it left off instead of redoing already-completed (area, m) work. Tasks
# still missing after all retry attempts are reported explicitly -- not
# silently dropped -- so an incomplete run is never mistaken for a
# complete one.
#
# cv_areas is a vector so this scales from "just North East, to see what
# happens" (its current default) up to all 9 regions later without any
# other changes -- the (area, m) task grid, checkpointing, and per-area
# results/plots below already handle either case.
#
# Speed tips (see bottom of file for full discussion):
#   - Reduce cv_folds to 3
#   - Reduce optim_maxit_cv (e.g. 200)
#   - Narrow the grid once you have a rough best region
#   - Only tune m first (fix lambda and kernel_bw_mult), then tune the rest
# ============================================================

source("nitrate_data_analysis/region_pipeline_funcs.R")
library(ggplot2)


# ============================================================
# CONFIGURATION  — edit only this section
# ============================================================

cv_areas       <- c("North East")  # start with one; add more areas to this vector later
cv_threshold   <- log(5)           # e.g. log(2), log(5), log(10)
cv_seed        <- 42               # fixed seed for fold assignment
cv_folds       <- 5                # number of CV folds
optim_maxit_cv <- 300              # maxit for optim inside each CV fold
                                    # (lower = faster, less precise)
output_dir     <- "nitrate_data_analysis/output"
year_filter    <- 2020             # set NULL to use all years
max_cv_attempts <- 5               # retry passes over still-missing (area, m) tasks

# On a SLURM cluster (e.g. Yale's Bouchet), detectCores() reports the whole
# node's core count, not what the job was actually allocated -- read SLURM's
# own env var when present, and only fall back to detectCores() - 1 for
# local/laptop runs outside SLURM. Parallelizes across (area, m) tasks.
mc.cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(mc.cores)) mc.cores <- max(1, parallel::detectCores() - 1)

# Hyperparameter grid
m_grid       <- c(-3, -2, -1.5, -1, 0, 0.5)
lambda_grid  <- c(0.1, 0.25, 0.5, 1.0)
bw_mult_grid <- c(0.5, 1, 2, 5, 10)

# Must match run_all_regions.R's area_key_map -- "Northwest" has no space
# in `data$area` but its cache files use "north_west", so this can't be
# derived as gsub(" ", "_", tolower(cv_area)) for every area.
area_key_map <- c(
  "Central"       = "central",
  "Northwest"     = "north_west",
  "North East"    = "north_east",
  "North Central" = "north_central",
  "West Central"  = "west_central",
  "East Central"  = "east_central",
  "South West"    = "south_west",
  "South Central" = "south_central",
  "South East"    = "south_east"
)

threshold_label <- paste0("log", round(exp(cv_threshold), 0))


# ============================================================
# DATA + PER-AREA SHARED OBJECTS (piece i)
# ============================================================

data_all <- load_nitrate_data(file_path = "data/data_Nitrate_with_covar.csv", zero_inflated = FALSE)
if (!is.null(year_filter)) data_all <- data_all[data_all$SampleYear >= year_filter, ]

# fit_region_base()/build_kernel_design_matrix() are already cached per
# area_key by output_dir, so re-running this script for an area already
# set up here is fast regardless of the checkpointing below.
area_ctxs <- lapply(cv_areas, function(cv_area) {
  area_key   <- area_key_map[[cv_area]]
  data_area  <- data_all[data_all$area == cv_area, ]
  data_split <- split_nitrate_data(data_area)

  region_base <- fit_region_base(data_split = data_split, area_key = area_key, output_dir = output_dir)
  kdm_info    <- build_kernel_design_matrix(data_split = data_split, region_base = region_base)
  fold_id     <- assign_folds(nrow(data_split$data), n_folds = cv_folds, seed_value = cv_seed)
  outcome_resid <- data_split$data$logconcentration_plus_median -
    region_base$outcome_reg$pred - region_base$outcome_reg$krige_values

  list(cv_area = cv_area, area_key = area_key, data_split = data_split,
       region_base = region_base, kdm_info = kdm_info,
       fold_id = fold_id, outcome_resid = outcome_resid)
})
names(area_ctxs) <- cv_areas


# ============================================================
# CHECKPOINTED, RETRYABLE CV — one task per (area, m)
# ============================================================

n_combos <- length(lambda_grid) * length(bw_mult_grid)
cat(sprintf(
  "\nCV for %d area(s) x %d m-values, threshold=log(%g)\n%d combos x %d folds = %d optim runs per (area, m)\n\n",
  length(cv_areas), length(m_grid), round(exp(cv_threshold)), n_combos, cv_folds, n_combos * cv_folds
))

cv_partial_dir <- file.path(output_dir, "cv_partial")
dir.create(cv_partial_dir, recursive = TRUE, showWarnings = FALSE)

cv_partial_path <- function(area_key, m) {
  file.path(cv_partial_dir, sprintf("%s_%s_m%s.rds", area_key, threshold_label, gsub("-", "neg", sprintf("%.2f", m))))
}

cv_task_grid <- expand.grid(cv_area = cv_areas, m = m_grid, stringsAsFactors = FALSE)

run_one_cv_group <- function(cv_area, m) {
  ctx   <- area_ctxs[[cv_area]]
  label <- sprintf("area=%s | m=%.2f", cv_area, m)
  cat(sprintf("%s: starting (%d combos x %d folds = %d fits)\n", label, n_combos, cv_folds, n_combos * cv_folds))

  result <- cv_hyperparam_group_for_m(
    m = m, lambda_lst = lambda_grid, kernel_bw_mult_lst = bw_mult_grid,
    data_train = ctx$data_split$data, kernel_design_matrix = ctx$kdm_info$kdm,
    kernel_bw_baseline = ctx$region_base$kernel_bw_silverman,
    fold_id = ctx$fold_id, n_folds = cv_folds,
    threshold_val = cv_threshold, depth_range = ctx$region_base$depth_range,
    krige_values = ctx$region_base$outcome_reg$krige_values,
    weights = ctx$region_base$weights_cbps$weights,
    outcome_resid = ctx$outcome_resid,
    smoothers_RKHS = ctx$region_base$smoothers$smoothers_RKHS,
    cumint_smoothers_RKHS = ctx$region_base$smoothers$cumint_smoothers_RKHS,
    progress_label = label, maxit = optim_maxit_cv
  )
  result$cv_area  <- cv_area
  result$area_key <- ctx$area_key

  # Incremental checkpoint: each (area, m) task writes its own result to
  # disk the moment it finishes, rather than only after the whole mclapply
  # call returns -- if a later task crashes/OOMs, this one's result is
  # already safely on disk.
  saveRDS(result, file = cv_partial_path(ctx$area_key, m))
  cat(sprintf("%s: finished\n", label))
  result
}

# Skip any (area, m) task that already has a checkpoint file from a
# previous run (e.g. one that got OOM-killed partway through), and retry
# whatever's still missing in a loop -- each retry only tackles the
# shrinking set of missing tasks, re-checking cv_partial/ after every pass.
for (attempt in seq_len(max_cv_attempts)) {
  area_keys_for_task <- area_key_map[cv_task_grid$cv_area]
  already_done <- file.exists(cv_partial_path(area_keys_for_task, cv_task_grid$m))
  todo_grid <- cv_task_grid[!already_done, , drop = FALSE]

  if (nrow(todo_grid) == 0) {
    if (attempt == 1) cat(sprintf("CV: all %d task(s) already done, nothing to run\n", nrow(cv_task_grid)))
    break
  }

  cat(sprintf("CV attempt %d/%d: %d of %d task(s) remaining\n",
              attempt, max_cv_attempts, nrow(todo_grid), nrow(cv_task_grid)))

  # mclapply silently returns NULL (not a "try-error") for a task whose
  # forked worker gets OOM-killed by the OS -- the missing-task check
  # below catches that case, not just R-level errors caught via try().
  parallel::mclapply(seq_len(nrow(todo_grid)), function(k) {
    run_one_cv_group(todo_grid$cv_area[k], todo_grid$m[k])
  }, mc.cores = min(nrow(todo_grid), mc.cores))
}

# Reassemble from disk (not from mclapply's return values), so tasks
# completed in an earlier attempt or a previous run of this script aren't
# lost even if a later retry's mclapply call itself fails partway through
# on the remaining tasks.
area_keys_for_task <- area_key_map[cv_task_grid$cv_area]
cv_task_results <- lapply(seq_len(nrow(cv_task_grid)), function(k) {
  path <- cv_partial_path(area_keys_for_task[k], cv_task_grid$m[k])
  if (file.exists(path)) readRDS(path) else NULL
})

missing <- vapply(cv_task_results, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
if (any(missing)) {
  missing_desc <- sprintf("area=%s, m=%.2f", cv_task_grid$cv_area[missing], cv_task_grid$m[missing])
  cat(sprintf("\nWARNING: %d of %d CV task(s) are still missing after %d attempt(s):\n%s\n",
              sum(missing), nrow(cv_task_grid), max_cv_attempts, paste(missing_desc, collapse = "\n")))
  cat("Re-run this script to retry just the missing tasks (already-completed ones will be skipped).\n\n")
}

cv_results_all <- dplyr::bind_rows(cv_task_results[!missing])


# ============================================================
# RESULTS + VISUALISATION — per area
# ============================================================

dir.create("nitrate_data_analysis/figures", recursive = TRUE, showWarnings = FALSE)

for (cv_area in cv_areas) {
  area_key   <- area_key_map[[cv_area]]
  cv_results <- cv_results_all[cv_results_all$cv_area == cv_area, ]

  if (nrow(cv_results) == 0) {
    cat(sprintf("\nNo completed results for %s (all its tasks are still missing) -- skipping results/plots.\n", cv_area))
    next
  }
  if (nrow(cv_results) < length(m_grid)) {
    cat(sprintf("\nNote: %s has results for only %d of %d m-values -- best-hyperparameter selection below is based on partial results.\n",
                cv_area, nrow(cv_results), length(m_grid)))
  }

  cv_rds_path <- file.path(output_dir, sprintf("cv_results_%s_%s.rds", area_key, threshold_label))
  saveRDS(cv_results, cv_rds_path)
  cat(sprintf("\nFull CV grid for %s saved to: %s\n", cv_area, cv_rds_path))

  # ---- Best by CV MCC ----
  best_idx    <- which.max(cv_results$mcc)
  best_params <- cv_results[best_idx, ]

  cat(sprintf("\n========== Best Hyperparameters for %s (max CV MCC) ==========\n", cv_area))
  cat(sprintf("  m               = %.1f\n",  best_params$m))
  cat(sprintf("  lambda          = %.2f\n",  best_params$lambda))
  cat(sprintf("  kernel_bw_mult  = %.1f  =>  kernel_bw = %.4f\n",
              best_params$kernel_bw_mult, best_params$kernel_bw))
  cat(sprintf("  CV MCC          = %.4f\n",  best_params$mcc))
  cat(sprintf("  CV Accuracy     = %.4f\n",  best_params$acc))
  cat(sprintf("  CV Two-sided F1 = %.4f\n",  best_params$two_sided_f1))
  cat("=======================================================\n")

  # 1. Heatmap of CV MCC for (m, lambda), averaged over kernel_bw_mult
  cv_ml <- cv_results %>%
    group_by(m, lambda) %>%
    summarise(mean_mcc = mean(mcc, na.rm = TRUE), .groups = "drop")

  p1 <- ggplot(cv_ml, aes(x = factor(lambda), y = factor(m), fill = mean_mcc)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(mean_mcc, 3)), size = 3) +
    scale_fill_viridis_c(option = "plasma", name = "Mean CV MCC") +
    labs(
      title = sprintf("CV MCC — %s, threshold = log(%g)", cv_area, round(exp(cv_threshold))),
      subtitle = "Averaged over kernel_bw_mult",
      x = "lambda", y = "m"
    ) +
    theme_minimal()
  print(p1)

  p1_path <- file.path("nitrate_data_analysis/figures", sprintf("cv_heatmap_%s_%s.png", area_key, threshold_label))
  ggsave(p1_path, p1, width = 7, height = 6, dpi = 150, bg = "white")
  cat(sprintf("Heatmap saved to: %s\n", p1_path))

  # 2. MCC vs kernel_bw_mult for the best (m, lambda)
  cv_bw <- cv_results %>%
    filter(m == best_params$m, lambda == best_params$lambda)

  p2 <- ggplot(cv_bw, aes(x = kernel_bw_mult, y = mcc)) +
    geom_line() +
    geom_point(size = 2) +
    geom_point(data = best_params, color = "red", size = 4) +
    labs(
      title   = sprintf("CV MCC vs kernel_bw multiplier — %s (m=%.1f, lambda=%.2f)",
                        cv_area, best_params$m, best_params$lambda),
      x       = "kernel_bw_mult",
      y       = "CV MCC"
    ) +
    theme_minimal()
  print(p2)

  p2_path <- file.path("nitrate_data_analysis/figures", sprintf("cv_bwmult_%s_%s.png", area_key, threshold_label))
  ggsave(p2_path, p2, width = 7, height = 6, dpi = 150, bg = "white")
  cat(sprintf("kernel_bw_mult plot saved to: %s\n", p2_path))
}


# ============================================================
# IF CV IS TOO SLOW — OPTIONS (in order of impact)
# ============================================================
#
# 1. SWITCH TO A SINGLE HOLD-OUT SPLIT (biggest speedup)
#    Instead of cv_folds = 5 (5x optim runs per combo), use cv_folds = 1
#    with a fixed 80/20 split -- would need a small change to
#    assign_folds() in region_pipeline_funcs.R to support a single split.
#
# 2. REDUCE optim_maxit_cv
#    For hyperparameter search you only need a rough solution.
#    optim_maxit_cv = 100-200 is usually sufficient.
#    The final full-data DC fit (fit_direct_dc_policy()) uses its own
#    max_iter, separately from this.
#
# 3. SEQUENTIAL TUNING — tune one parameter at a time
#    Step 1: fix lambda=0.25, bw_mult=1, sweep over m.
#    Step 2: fix best m, bw_mult=1, sweep over lambda.
#    Step 3: fix best m and lambda, sweep over bw_mult.
#    Reduces the grid from |m| x |lambda| x |bw_mult| runs
#    to |m| + |lambda| + |bw_mult| runs.
#
# 4. PARALLELIZE ACROSS (area, m) TASKS
#    Already implemented above: the retry loop's mclapply runs over the
#    full (area, m) task grid, not just m within one area, so adding more
#    areas to cv_areas scales the available parallelism up with it. Set
#    mc.cores in the CONFIGURATION section above (bump this on Bouchet).
#
# 5. REUSE SMOOTHERS AND OUTCOME REGRESSION ACROSS FOLDS
#    Already implemented: fit_region_base() computes/caches these ONCE
#    per area; cv_hyperparam_group_for_m() (in region_pipeline_funcs.R)
#    also builds the Gram matrix ONCE per m and slices it per fold, rather
#    than recomputing kernelMatrix() per (lambda, kernel_bw_mult) combo.
#
# 6. RANDOM SEARCH INSTEAD OF GRID SEARCH
#    Sample n_random = 20-30 random combos rather than exhaustively
#    evaluating the full grid, e.g. by subsetting m_grid/lambda_grid/
#    bw_mult_grid before running the CV loop above.
