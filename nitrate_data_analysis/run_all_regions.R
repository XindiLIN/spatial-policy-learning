# ============================================================
# run_all_regions.R
# ============================================================
# Driver: runs the per-region pipeline (region_pipeline_funcs.R) for
# every Wisconsin region and combines the results into a single map
# (direct method + non-spatial indirect method).
#
# Two things live here:
#   run_region_pipeline(area, ...) -- runs pieces (i)-(v) for ONE
#     region, in order, and returns everything needed for the map.
#   The main loop below applies it to every region and rbinds the
#     per-region sf outputs into the combined map layers.
#
# Each region's run is itself expensive (SVM+kriging, hyperparameter
# CV, per-observation DC warm start), so:
#   - fit_region_base() and choose_direct_hyperparameters_cv() cache
#     their (slow) results to `output_dir` (see region_pipeline_funcs.R)
#   - this driver additionally caches each region's FULL result to
#     `output_dir/region_result_<area_key>.rds`, so a rerun after an
#     interruption skips regions that already finished
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

threshold_val <- log(5)
year_filter   <- 2020

# Hyperparameter CV grid (matches nitrate_data_analysis/cv_hyperparameter.R's
# established defaults)
m_lst              <- c(-3, -2, -1.5, -1, 0, 0.5)
lambda_lst         <- c(0.1, 0.25, 0.5, 1.0)
kernel_bw_mult_lst <- c(0.5, 1, 2, 5, 10)
n_folds            <- 5
cv_seed            <- 42

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

# Bump these when running on Bouchet; kept conservative for a local run.
cv_mc.cores   <- 1  # parallelizes across `m` values inside CV (piece v)
init_mc.cores <- 1  # parallelizes the per-observation DC warm start (piece iii)

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
                                 m_lst, lambda_lst, kernel_bw_mult_lst,
                                 plss_county_override = NULL, plss_crop_exclude = NULL,
                                 RKHS_covariate_names = RKHS_COVARIATE_NAMES_DEFAULT,
                                 n_folds = 5, cv_seed = 42,
                                 cv_mc.cores = 1, init_mc.cores = 1,
                                 use_cache = TRUE) {

  area_key <- area_key_map[[area]]
  cat(sprintf("\n%s\nRegion: %s (%s)\n%s\n", strrep("=", 60), area, area_key, strrep("=", 60)))

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

  # (v) hyperparameter selection via 5-fold CV
  cv_path <- file.path(output_dir, sprintf("cv_results_%s.rds", area_key))
  if (use_cache && file.exists(cv_path)) {
    cat("Loading cached CV results...\n")
    cv_results <- readRDS(cv_path)
  } else {
    cat("Running hyperparameter CV...\n")
    cv_results <- choose_direct_hyperparameters_cv(
      data_split = data_split, region_base = region_base, kdm_info = kdm_info,
      m_lst = m_lst, lambda_lst = lambda_lst, kernel_bw_mult_lst = kernel_bw_mult_lst,
      threshold_val = threshold_val, n_folds = n_folds, seed_value = cv_seed,
      mc.cores = cv_mc.cores
    )
    saveRDS(cv_results, cv_path)
  }
  best <- cv_results[which.max(cv_results$mcc), ]
  cat(sprintf("Best CV hyperparameters: m=%.1f  lambda=%.2f  kernel_bw_mult=%.1f  (CV mcc=%.4f)\n",
              best$m, best$lambda, best$kernel_bw_mult, best$mcc))

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
# MAIN LOOP -- apply the pipeline to every region
# ============================================================
# Sequential across regions (each region already parallelizes
# internally via cv_mc.cores/init_mc.cores); nesting mclapply at both
# levels would oversubscribe cores. Each region's full result is
# cached, so re-running this script after an interruption picks up
# where it left off.

region_results <- list()

for (area in areas) {
  area_key <- area_key_map[[area]]
  result_path <- file.path(output_dir, sprintf("region_result_%s.rds", area_key))

  if (use_cache && file.exists(result_path)) {
    cat(sprintf("Loading cached region result for %s...\n", area))
    region_results[[area]] <- readRDS(result_path)
    next
  }

  region_results[[area]] <- run_region_pipeline(
    area = area,
    data_all = data_all, plss_covariates_all = plss_covariates_all, output_dir = output_dir,
    threshold_val = threshold_val,
    m_lst = m_lst, lambda_lst = lambda_lst, kernel_bw_mult_lst = kernel_bw_mult_lst,
    plss_county_override = plss_county_override, plss_crop_exclude = plss_crop_exclude,
    n_folds = n_folds, cv_seed = cv_seed,
    cv_mc.cores = cv_mc.cores, init_mc.cores = init_mc.cores,
    use_cache = use_cache
  )
  saveRDS(region_results[[area]], result_path)
}


# ============================================================
# COMBINE RESULTS ACROSS REGIONS
# ============================================================

plss_sf_combine_direct <- do.call(rbind, lapply(region_results, `[[`, "plss_sf_direct"))
plss_sf_combine_indirect_nonspatial <- do.call(rbind, lapply(region_results, `[[`, "plss_sf_indirect"))

saveRDS(plss_sf_combine_direct, file.path(output_dir, "plss_sf_combine_direct.rds"))
saveRDS(plss_sf_combine_indirect_nonspatial,
        file.path(output_dir, "plss_sf_combine_indirect_nonspatial.rds"))

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
saveRDS(summary_df, file.path(output_dir, "region_pipeline_summary.rds"))
print(summary_df)


# ============================================================
# SANITY-CHECK MAP (quick look only -- see policy_visualization.r
# for the polished, paper-ready figures)
# ============================================================

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

ggsave(file.path("nitrate_data_analysis/figures", "sanity_map_direct.png"),
       make_sanity_map(plss_sf_combine_direct, "Direct Method"),
       width = 8, height = 7, dpi = 200, bg = "white")
ggsave(file.path("nitrate_data_analysis/figures", "sanity_map_indirect_nonspatial.png"),
       make_sanity_map(plss_sf_combine_indirect_nonspatial, "Non-Spatial Indirect Method"),
       width = 8, height = 7, dpi = 200, bg = "white")
