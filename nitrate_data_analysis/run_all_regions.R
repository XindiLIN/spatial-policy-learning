# ============================================================
# run_all_regions.R
# ============================================================
# Driver: runs the per-region pipeline (region_pipeline_funcs.R) for
# every Wisconsin region and combines the results into a single map
# (direct method + non-spatial indirect method).
#
# Two things live here:
#   run_region_pipeline(area, ...) -- runs pieces (i)-(iv) for ONE region (in
#     order, loading piece (v)'s hyperparameters from cv_hyperparameter.R's
#     output rather than computing them here -- see that function's cv_path
#     comment) and returns everything needed for the map.
#   The main loop below applies it to every region and rbinds the
#     per-region sf outputs into the combined map layers.
#
# Each region's run is itself expensive (SVM+kriging, per-observation DC warm
# start), so:
#   - fit_region_base() caches its (slow) results to `output_dir` (see
#     region_pipeline_funcs.R); hyperparameter CV must already be cached
#     there too, by cv_hyperparameter.R, before this script is run
#   - this driver additionally caches each region's FULL result to
#     `output_dir/area_results/region_result_<area_key>.rds`, so a rerun
#     after an interruption skips regions that already finished
# ============================================================

source("nitrate_data_analysis/region_pipeline_funcs.R")
library(tigris)
library(ggplot2)


# ============================================================
# CONFIGURATION
# ============================================================

areas <- c("Central", "Northwest", "North East", "North Central",
           "West Central", "East Central", "South West", "South Central", "South East")

# area_key is used for cache file names. NOT derivable as
# gsub(" ", "_", tolower(area)) for every area: "Northwest" has no
# space in `data$area` but the region's cache files (and the original
# per-region scripts) use "north_west". Kept explicit to avoid a
# silent cache-miss -> full, unnecessary Northwest refit.
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

threshold_vals <- c(log(2), log(5), log(10))  # any number of thresholds; all get run
# Must match cv_hyperparameter.R's own threshold_label_for() so this driver
# looks up the right cv_results_<area_key>_<threshold_label>.rds, and so this
# run's own outputs (region_result_*, plss_sf_combine_*,
# region_pipeline_summary, sanity maps) don't collide across thresholds.
threshold_label_for <- function(threshold_val) paste0("log", round(exp(threshold_val), 0))
year_filter   <- 2020

# Explicit county lists for PLSS sub-setting; NULL means use the counties
# actually observed in that region's training wells
plss_county_override <- list(
  "South West"    = c("Vernon", "Crawford", "Grant", "Richland",
                      "Sauk", "Iowa", "Lafayette"),
  "South Central" = c("Columbia", "Dodge", "Dane", "Jefferson",
                      "Green", "Rock")
)

# Crop types to exclude from PLSS per area. NOTE: fixed from the original
# draw_minimum_depth_map.R, which had "Vegetables" (no effect -- the
# crop_type_combine level is actually spelled "Vegtables") and
# "Other Crop" (also no effect -- the PLSS data's crop_type_combine has no
# "Other Crop" level at all, only "Undefined", which is already dropped by
# the default county/crop filtering below since it never matches any
# well-data crop type). Left in place in case region-specific exclusions
# are still wanted, but currently neither entry does anything.
plss_crop_exclude <- list(
  "South West"    = "Vegtables",
  "North Central" = "Other Crop"
)

output_dir <- "nitrate_data_analysis/output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Final per-region results and the combined cross-region outputs go in their
# own subfolder, separate from output_dir's CV artifacts (cv_results_*.rds,
# cv_partial/) and the piece (i)/(iii) caches (outcome_regression_*.rds,
# smoothers_*.rds, direct_init_*.rds) -- those all stay directly under
# output_dir since fit_region_base()/prep_direct_method_inputs() and the CV
# results lookup below must keep reading from the same location
# cv_hyperparameter.R already wrote them to.
region_results_dir <- file.path(output_dir, "area_results")
dir.create(region_results_dir, recursive = TRUE, showWarnings = FALSE)

# On a SLURM cluster (e.g. Yale's Bouchet), detectCores() reports the whole
# node's core count, not what the job was actually allocated -- read SLURM's
# own env var when present, and only fall back to detectCores() - 1 for
# local/laptop runs outside SLURM (matches cv_hyperparameter.R's mc.cores
# pattern). This is the TOTAL budget, split across two levels below: across
# (area, threshold) tasks (the main loop), and within each task's own piece
# (iii) per-observation DC warm start. Piece (iv), the DC algorithm itself,
# is single-threaded -- running regions sequentially left most of a 16-core
# allocation idle for that entire serial phase, every region, every
# threshold (measured: 6.4% CPU efficiency over a 4-hour run). Running many
# (area, threshold) tasks concurrently instead fills that idle time with
# OTHER tasks' work, since all 9 areas x however many thresholds are fully
# independent (separate data, separate area_key/threshold-scoped caches).
mc.cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(mc.cores)) mc.cores <- max(1, parallel::detectCores() - 1)

max_region_attempts <- 5  # retry passes over still-missing (area, threshold) tasks

use_cache <- TRUE  # reuse cached outcome_reg / smoothers / CV / region results


# ============================================================
# SHARED DATA (loaded once, reused across all regions)
# ============================================================

data_all <- load_nitrate_data(file_path = "data/data_Nitrate_with_covar.csv", zero_inflated = FALSE)
if (!is.null(year_filter)) data_all <- data_all[data_all$SampleYear >= year_filter, ]

plss_covariates_all <- read.csv("data/plss_covariates.csv")


# ============================================================
# PER-REGION PIPELINE
# ============================================================
# Runs pieces (i)-(v) from region_pipeline_funcs.R, in order, for one
# region, and scores both methods on that region's PLSS grid.
run_region_pipeline <- function(area,
                                 data_all, plss_covariates_all, output_dir,
                                 threshold_val,
                                 plss_county_override = NULL, plss_crop_exclude = NULL,
                                 RKHS_covariate_names = RKHS_COVARIATE_NAMES_DEFAULT,
                                 init_mc.cores = 1,
                                 use_cache = TRUE) {

  area_key <- area_key_map[[area]]
  cat(sprintf("\n%s\nRegion: %s (%s)\n%s\n", strrep("=", 60), area, area_key, strrep("=", 60)))
  flush(stdout())

  data_area  <- data_all[data_all$area == area, ]
  data_split <- split_nitrate_data(data_area)

  # (i) outcome regression + weights_cbps + smoothers
  region_base <- fit_region_base(
    data_split = data_split, area_key = area_key, output_dir = output_dir,
    RKHS_covariate_names = RKHS_covariate_names, use_cache = use_cache
  )

  # (ii) indirect non-spatial policy, train + test
  indirect_result <- get_indirect_nonspatial_policy(
    region_base = region_base, data_split = data_split, threshold_val = threshold_val
  )

  # shared kernel design matrix (needed by both CV and the final DC prep)
  kdm_info <- build_kernel_design_matrix(
    data_split = data_split, region_base = region_base,
    RKHS_covariate_names = RKHS_covariate_names
  )

  # (v) hyperparameter selection -- REQUIRES cv_hyperparameter.R to already
  # have been run for this area. There is deliberately no in-line fallback
  # to compute CV here: an earlier version of this file built the cache path
  # without the threshold-label suffix that cv_hyperparameter.R actually
  # uses (cv_results_<area_key>.rds vs. the real cv_results_<area_key>_<label>.rds),
  # so file.exists() always missed and silently recomputed CV from scratch
  # with the older, unparallelized choose_direct_hyperparameters_cv() path --
  # slow, and a second CV implementation that can drift out of sync with the
  # first. Failing loudly here instead keeps there being exactly one way to
  # compute CV results.
  threshold_label <- threshold_label_for(threshold_val)
  cv_path <- file.path(output_dir, "hyper_parameter_cv", sprintf("cv_results_%s_%s.rds", area_key, threshold_label))
  if (!file.exists(cv_path)) {
    stop(sprintf(
      "run_region_pipeline(): no CV results found for %s at %s. Run cv_hyperparameter.R for this area first (add \"%s\" to its cv_areas).",
      area, cv_path, area
    ))
  }
  cat("Loading CV results...\n")
  flush(stdout())
  cv_results <- readRDS(cv_path)
  if (all(is.na(cv_results$mcc) | is.nan(cv_results$mcc))) {
    stop(sprintf("run_region_pipeline(): every combo's mcc in %s is NA/undefined.", cv_path))
  }
  best <- cv_results[which.max(cv_results$mcc), ]
  cat(sprintf("Best CV hyperparameters: m=%.1f  lambda=%.2f  kernel_bw_mult=%.1f  (CV mcc=%.4f)\n",
              best$m, best$lambda, best$kernel_bw_mult, best$mcc))
  flush(stdout())

  # (iii) preparations before the DC algorithm, at the CV-chosen m/kernel_bw
  prep <- prep_direct_method_inputs(
    data_split = data_split, region_base = region_base, kdm_info = kdm_info,
    m = best$m, kernel_bw = best$kernel_bw,
    plss_covariates_all = plss_covariates_all, area = area, area_key = area_key, output_dir = output_dir,
    plss_county_override = plss_county_override[[area]],
    plss_crop_exclude    = plss_crop_exclude[[area]],
    RKHS_covariate_names = RKHS_covariate_names, threshold_val = threshold_val,
    mc.cores = init_mc.cores, use_cache = use_cache
  )

  # (iv) run the DC algorithm at the CV-chosen lambda/kernel_bw
  direct_result <- fit_direct_dc_policy(
    data_split = data_split, region_base = region_base, prep = prep,
    lambda = best$lambda, kernel_bw = best$kernel_bw, threshold_val = threshold_val
  )

  # Score the indirect non-spatial policy on the same PLSS grid, reusing
  # the PLSS smoothers already built in prep() -- gives the indirect
  # method's map layer alongside the direct method's.
  plss_indirect        <- prep$plss
  plss_indirect$policy <- compute_indirect_policy(
    smoothers = prep$smoothers_plss_RKHS, trt_bounds = region_base$depth_range,
    threshold_val = threshold_val, spatial = FALSE
  )
  plss_sf_indirect <- st_as_sf(plss_indirect, coords = c("longitude", "latitude"), crs = 4326)

  list(
    area = area, area_key = area_key,
    indirect_metrics_train = indirect_result$metrics_train,
    indirect_metrics_test  = indirect_result$metrics_test,
    direct_metrics_train   = direct_result$metrics_train,
    direct_metrics_test    = direct_result$metrics_test,
    best_hyperparams  = best,
    plss_sf_direct    = direct_result$plss_sf,
    plss_sf_indirect  = plss_sf_indirect
  )
}


# ============================================================
# CHECKPOINTED, RETRYABLE MAIN LOOP -- one task per (area, threshold)
# ============================================================
# Parallelized across (area, threshold) pairs -- fully independent, so
# running many concurrently fills the idle time piece (iv)'s single-threaded
# DC algorithm otherwise leaves on the table (see mc.cores comment above).
# Each concurrently-running task gets a share of mc.cores for its OWN piece
# (iii) internal parallelism, sized so the two levels together don't
# oversubscribe the allocation.
#
# Checkpointed and retryable the same way as cv_hyperparameter.R: each task
# writes its own region_result_<area_key>_<threshold_label>.rds the moment
# it finishes, so a killed/OOM'd run resumes without redoing already-
# completed (area, threshold) work, and adding a new threshold to
# threshold_vals later only computes the new pairs, not the old ones.

region_threshold_grid <- expand.grid(area = areas, threshold_val = threshold_vals, stringsAsFactors = FALSE)

n_concurrent_tasks <- min(nrow(region_threshold_grid), mc.cores)
# Cores left over for each concurrently-running task's own piece (iii)
# parallelism, after the outer (area, threshold) level's share.
init_mc.cores_per_task <- max(1, floor(mc.cores / n_concurrent_tasks))

region_result_path <- function(area_key, threshold_label) {
  file.path(region_results_dir, sprintf("region_result_%s_%s.rds", area_key, threshold_label))
}

run_one_region_threshold <- function(area, threshold_val) {
  area_key        <- area_key_map[[area]]
  threshold_label <- threshold_label_for(threshold_val)
  label <- sprintf("area=%s | threshold=%s", area, threshold_label)

  result <- run_region_pipeline(
    area = area,
    data_all = data_all, plss_covariates_all = plss_covariates_all, output_dir = output_dir,
    threshold_val = threshold_val,
    plss_county_override = plss_county_override, plss_crop_exclude = plss_crop_exclude,
    init_mc.cores = init_mc.cores_per_task,
    use_cache = use_cache
  )
  saveRDS(result, region_result_path(area_key, threshold_label))
  cat(sprintf("%s: done\n", label))
  flush(stdout())
  result
}

cat(sprintf("\nrun_all_regions: %d area(s) x %d threshold(s) = %d (area, threshold) tasks, %d concurrent, %d core(s) each for piece iii\n\n",
            length(areas), length(threshold_vals), nrow(region_threshold_grid), n_concurrent_tasks, init_mc.cores_per_task))
flush(stdout())

for (attempt in seq_len(max_region_attempts)) {
  area_keys_for_task <- area_key_map[region_threshold_grid$area]
  threshold_labels_for_task <- vapply(region_threshold_grid$threshold_val, threshold_label_for, character(1))
  already_done <- use_cache & file.exists(region_result_path(area_keys_for_task, threshold_labels_for_task))
  todo_grid <- region_threshold_grid[!already_done, , drop = FALSE]

  if (nrow(todo_grid) == 0) {
    if (attempt == 1) cat(sprintf("run_all_regions: all %d task(s) already done, nothing to run\n", nrow(region_threshold_grid)))
    break
  }

  cat(sprintf("run_all_regions attempt %d/%d: %d of %d task(s) remaining\n",
              attempt, max_region_attempts, nrow(todo_grid), nrow(region_threshold_grid)))
  flush(stdout())

  # mclapply silently returns NULL (not a "try-error") for a task whose
  # forked worker gets OOM-killed by the OS -- the missing-task check
  # below catches that case, not just R-level errors caught via try().
  parallel::mclapply(seq_len(nrow(todo_grid)), function(k) {
    run_one_region_threshold(todo_grid$area[k], todo_grid$threshold_val[k])
  }, mc.cores = n_concurrent_tasks)
}


# ============================================================
# RESULTS + VISUALISATION -- per threshold, reassembled from disk
# ============================================================
# Reads from disk (not from mclapply's return values), so tasks completed in
# an earlier attempt or a previous run of this script aren't lost even if a
# later retry's mclapply call itself fails partway through on the
# remaining tasks. wi_counties.rds is threshold-independent (just Wisconsin
# county boundaries), so it's loaded/cached once outside the threshold loop.

wi_counties_path <- file.path(output_dir, "wi_counties.rds")
if (file.exists(wi_counties_path)) {
  counties_wi <- readRDS(wi_counties_path)
} else {
  counties_wi <- tigris::counties(state = "WI", cb = TRUE, year = 2022)
  counties_wi <- counties_wi[, c("NAME", "geometry")]
  colnames(counties_wi) <- c("County", "geometry")
  saveRDS(counties_wi, wi_counties_path)
}

make_sanity_map <- function(plss_sf, title_str) {
  ggplot() +
    geom_sf(data = plss_sf, aes(color = policy), size = 0.2) +
    geom_sf(data = counties_wi, color = "black", fill = NA) +
    scale_color_viridis_c(option = "plasma", name = "Log Min. Depth",
                          labels = function(b) round(exp(b), 0)) +
    labs(title = title_str, x = "Longitude", y = "Latitude") +
    theme_void()
}

for (threshold_val in threshold_vals) {
  threshold_label <- threshold_label_for(threshold_val)
  cat(sprintf("\n%s\n%s THRESHOLD = %s %s\n%s\n", strrep("#", 60), strrep("#", 10), threshold_label, strrep("#", 10), strrep("#", 60)))
  flush(stdout())

  region_results <- lapply(areas, function(area) {
    path <- region_result_path(area_key_map[[area]], threshold_label)
    if (file.exists(path)) readRDS(path) else NULL
  })
  names(region_results) <- areas

  missing <- vapply(region_results, is.null, logical(1))
  if (any(missing)) {
    cat(sprintf("\nWARNING: %d of %d region(s) still missing for threshold %s after %d attempt(s): %s\n",
                sum(missing), length(areas), threshold_label, max_region_attempts, paste(areas[missing], collapse = ", ")))
    cat("Re-run this script to retry just the missing (area, threshold) tasks.\n\n")
  }
  region_results <- region_results[!missing]
  if (length(region_results) == 0) {
    cat(sprintf("No completed regions for threshold %s -- skipping combine/summary/sanity map.\n", threshold_label))
    next
  }


  # ---- Combine results across regions, for this threshold ----

  plss_sf_combine_direct <- do.call(rbind, lapply(region_results, `[[`, "plss_sf_direct"))
  plss_sf_combine_indirect_nonspatial <- do.call(rbind, lapply(region_results, `[[`, "plss_sf_indirect"))

  saveRDS(plss_sf_combine_direct, file.path(region_results_dir, sprintf("plss_sf_combine_direct_%s.rds", threshold_label)))
  saveRDS(plss_sf_combine_indirect_nonspatial,
          file.path(region_results_dir, sprintf("plss_sf_combine_indirect_nonspatial_%s.rds", threshold_label)))

  # Summary table: chosen hyperparameters + train/test metrics per region,
  # for a quick sanity check before handing off to policy_visualization.r
  summary_df <- do.call(rbind, lapply(region_results, function(r) {
    data.frame(
      area = r$area,
      m = r$best_hyperparams$m, lambda = r$best_hyperparams$lambda,
      kernel_bw_mult = r$best_hyperparams$kernel_bw_mult, cv_mcc = r$best_hyperparams$mcc,
      direct_train_mcc = r$direct_metrics_train$mcc, direct_test_mcc = r$direct_metrics_test$mcc,
      indirect_train_mcc = r$indirect_metrics_train$mcc, indirect_test_mcc = r$indirect_metrics_test$mcc
    )
  }))
  rownames(summary_df) <- NULL
  saveRDS(summary_df, file.path(region_results_dir, sprintf("region_pipeline_summary_%s.rds", threshold_label)))
  cat(sprintf("\n=== Summary for threshold %s ===\n", threshold_label))
  print(summary_df)


  # ---- Sanity-check maps for this threshold (quick look only -- see
  # policy_visualization.r for the polished, paper-ready figures) ----

  ggsave(file.path("nitrate_data_analysis/figures", sprintf("sanity_map_direct_%s.png", threshold_label)),
         make_sanity_map(plss_sf_combine_direct, "Direct Method"),
         width = 8, height = 7, dpi = 200, bg = "white")
  ggsave(file.path("nitrate_data_analysis/figures", sprintf("sanity_map_indirect_nonspatial_%s.png", threshold_label)),
         make_sanity_map(plss_sf_combine_indirect_nonspatial, "Non-Spatial Indirect Method"),
         width = 8, height = 7, dpi = 200, bg = "white")
}
