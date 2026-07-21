# ============================================================
# draw_minimum_depth_map.R
# ============================================================
# Unified script to train and visualize minimum well depth
# policy maps across Wisconsin regions.
#
# Configuration (edit section below):
#   areas_to_run   : regions to process
#   threshold_vals : nitrate thresholds, e.g. c(log(2.5), log(5.5), log(10.5))
#   methods_to_run : any subset of
#                    "indirect_nonspatial", "indirect_spatial", "direct"
#   map_mode       : "individual"  -> one plot per area
#                    "combined"    -> all areas on one plot
# ============================================================


# ============================================================
# CONFIGURATION
# ============================================================

areas_to_run <- c(
  "Central", "Northwest", "North East", "North Central",
  "South West", "South Central", "East Central", "South East",
  "West Central"
)

areas_to_run <- c("Central")

# threshold_vals <- c(log(5.5))
threshold_vals <- c(log(10.5))
# e.g. c(log(2.5), log(5.5), log(10.5))

# methods_to_run <- c("indirect_spatial", "direct")
methods_to_run <- c("indirect_nonspatial","indirect_spatial", "direct")
# Options: "indirect_nonspatial", "indirect_spatial", "direct"

map_mode    <- "combined"   # "combined" or "individual"
year_filter <- 2020         # set NULL to use all years

# output_dir <- "."
output_dir <- "/Users/xindilin/Desktop/policy learning/simulation/clear_version"
# plss_path  <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/plss_covariates.csv"
plss_path <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/plss_covariates_static_no_na.csv" # this file has cafolog, StaticLevel and omit NAs


# ============================================================
# AREA-SPECIFIC TUNING PARAMETERS  (from cross-validation)
# ============================================================

# Chosen kernel bandwidth exponent m per area
m_by_area <- c(
  "Central"       = -2.0,
  "Northwest"     = -1.0,
  "North East"    = -2.0,
  "North Central" = -1.5,
  "South West"    = -3.0,
  "South Central" = -1.5,
  "East Central"  = -2.0,
  "South East"    = -1.5,
  "West Central"  = -2.0   # update after cross-validation
)

# Explicit county lists for PLSS sub-setting; NULL means use training data counties
plss_county_override <- list(
  "South West"    = c("Vernon", "Crawford", "Grant", "Richland",
                      "Sauk", "Iowa", "Lafayette"),
  "South Central" = c("Columbia", "Dodge", "Dane", "Jefferson",
                      "Green", "Rock")
)

# Crop types to EXCLUDE from PLSS per area
# For all other areas, PLSS is filtered to match training crop types
plss_crop_exclude <- list(
  "South West"    = "Vegetables",
  "North Central" = "Other Crop"
)


# ============================================================
# SOURCES AND PACKAGES
# ============================================================

source("data_generation.R")
source("funcs_kriging.R")
source("loss_functions.R")
source("loss_functions_RKHS.R")
source("miscellaneous.r")
library(rnaturalearthdata)
library(rnaturalearth)
library(WeightIt)
library(glmnet)
library(kernlab)
library(sf)
library(tigris)
library(ggplot2)
library(Metrics)
library(parallel)



# ============================================================
# SHARED SETUP
# ============================================================

command <- paste("brctl download", shQuote(plss_path))
system(command)
plss_covariates_all         <- read.csv(plss_path)
# plss_covariates_all$cafolog <- log(plss_covariates_all$cafo + 1)

data_all <- load_nitrate_data(zero_inflated = FALSE)
if (!is.null(year_filter))
  data_all <- data_all[data_all$SampleYear >= year_filter, ]

counties_wi           <- counties(state = "WI", cb = TRUE, year = 2022)
counties_wi           <- counties_wi[, c("NAME", "geometry")]
colnames(counties_wi) <- c("County", "geometry")

RKHS_covariate_names <- c("StaticLevel", "crop_type_combine", "drainagecl",
                           "precipitation", "cafolog")

# Results stored as: sf_store[[area_key]][[threshold_label]][[method]]
sf_store <- list()


# ============================================================
# MAIN LOOP OVER AREAS
# ============================================================

area <- areas_to_run[1]
for (area in areas_to_run) {

  cat(sprintf("\n%s\nArea: %s\n%s\n", strrep("=", 60), area, strrep("=", 60)))
  area_key <- gsub(" ", "_", tolower(area))

  # ---- Data ----
  data_area       <- data_all[data_all$area == area, ]
  data_area_split <- split_nitrate_data(data_area)

  # ---- Outcome regression ----
  outcome_reg_path <- file.path(output_dir, sprintf("outcome_regression_%s.rds", area_key))
  if (file.exists(outcome_reg_path)) {
    cat("Loading cached outcome regression...\n")
    outcome_reg <- readRDS(outcome_reg_path)
  } else {
    cat("Fitting outcome regression...\n")
    outcome_reg <- outcome_regression_SVM(
      data      = data_area_split$data,
      data_test = data_area_split$data_test,
      tunning   = FALSE
    )
    saveRDS(outcome_reg, outcome_reg_path)
  }

  cat(sprintf("RMSE  train: %.4f   test: %.4f\n",
    Metrics::rmse(data_area_split$data$logconcentration_plus_median,
                  outcome_reg$pred + outcome_reg$krige_values),
    Metrics::rmse(data_area_split$data_test$logconcentration_plus_median,
                  outcome_reg$pred_test + outcome_reg$krige_values_test)))

  depth_range <- c(min(data_area_split$data$logWellDepth),
                   max(data_area_split$data$logWellDepth))


  # ---- Propensity weights and kernel bandwidth (direct method only) ----
  if ("direct" %in% methods_to_run) {
    weights_cbps <- weightit(
      logWellDepth ~ StaticLevel + crop_type_combine + drainagecl +
                     precipitation + cafolog,
      data = data_area_split$data, method = "cbps"
    )
    kernel_bw <- 1.06 * sd(data_area_split$data$logWellDepth) *
                 nrow(data_area_split$data)^(-1/5)
  }

  # ---- Kernel design matrices ----
  rec <- recipe(~ ., data = data_area_split$data[, RKHS_covariate_names]) %>%
    step_dummy(all_nominal_predictors()) %>%
    prep()

  kdm      <- bake(rec, new_data = data_area_split$data[, RKHS_covariate_names])
  kdm      <- cbind(kdm, U = outcome_reg$krige_values)
  kdm_test <- bake(rec, new_data = data_area_split$data_test[, RKHS_covariate_names])
  kdm_test <- cbind(kdm_test, U = outcome_reg$krige_values_test)

  train_means <- colMeans(kdm)
  train_sds   <- apply(kdm, 2, sd)
  kdm         <- scale(kdm)
  kdm_test    <- scale(kdm_test, center = train_means, scale = train_sds)

  m_area <- m_by_area[area]
  gamma  <- 2^(m_area) * median(fields::rdist(kdm))
  rbf    <- rbfdot(sigma = 1 / (gamma^2))

  K      <- cbind(1, kernelMatrix(rbf, kdm))
  K_test <- cbind(1, kernelMatrix(rbf, kdm_test, kdm))

  # ---- PLSS covariates for this area ----
  if (area %in% names(plss_county_override)) {
    plss <- plss_covariates_all[
      plss_covariates_all$County %in% plss_county_override[[area]], ]
  } else {
    plss <- plss_covariates_all[
      plss_covariates_all$County %in% unique(data_area$County), ]
  }
  plss <- na.omit(plss)

  if (area %in% names(plss_crop_exclude)) {
    plss <- plss[!plss$crop_type_combine %in% plss_crop_exclude[[area]], ]
  } else {
    plss <- plss[plss$crop_type_combine %in% unique(data_area$crop_type_combine), ]
  }

  # Spatial random effect for PLSS
  plss_krige <- GpGp::predictions(
    fit       = outcome_reg$gpfit,
    locs_pred = plss[, c("longitude", "latitude")],
    X_pred    = rep(1, nrow(plss))
  )
  plss$U <- plss_krige

  # PLSS kernel matrix
  plss_kdm <- bake(rec, new_data = plss[, RKHS_covariate_names])
  plss_kdm <- cbind(plss_kdm, U = plss_krige)
  plss_kdm <- scale(plss_kdm, center = train_means, scale = train_sds)
  K_plss   <- cbind(1, kernelMatrix(rbf, plss_kdm, kdm))
  
  
  # ---- Smoothers ----
  # Needed for all methods: indirect methods use them directly;
  # the direct method uses them inside find_best_policy (initialization).
  
  smoothers_path <- file.path(output_dir, sprintf("smoothers_%s.rds", area_key))
  if(file.exists(smoothers_path)) {
    cat("Loading cached smoothers...\n")
  } else {
    cat("Computing smoothers (slow)...\n")
    smoothers <- get_smoothers_RKHS(
      design_matrix      = outcome_reg$design_matrix,
      design_matrix_test = outcome_reg$design_matrix_test,
      design_matrix_plss = bake(rec, new_data = plss[, RKHS_covariate_names]),
      svm_auto           = outcome_reg$svm,
      depth_range        = depth_range
    )
    saveRDS(smoothers,
            file.path(output_dir, sprintf("smoothers_%s.rds", area_key)))
  }

  sf_store[[area_key]] <- list()
  
  
  


  # ==========================================================
  # THRESHOLD LOOP
  # ==========================================================
  threshold_val <- threshold_vals[1]
  for (threshold_val in threshold_vals) {

    threshold_label <- paste0("log", round(exp(threshold_val), 0))
    cat(sprintf("\n  Threshold: %s\n", threshold_label))
    sf_store[[area_key]][[threshold_label]] <- list()

    # ---- Helper: evaluate indirect policy and save sf ----
    run_indirect <- function(spatial) {
      pol_train <- compute_indirect_policy(
        smoothers     = smoothers$smoothers_RKHS,
        trt_bounds    = depth_range,
        threshold_val = threshold_val,
        krige_adjust  = outcome_reg$krige_values,
        spatial       = spatial
      )
      pol_test <- compute_indirect_policy(
        smoothers     = smoothers$smoothers_test_RKHS,
        trt_bounds    = depth_range,
        threshold_val = threshold_val,
        krige_adjust  = outcome_reg$krige_values_test,
        spatial       = spatial
      )
      cat("    Training: "); print_acc_mcc(data_area_split$data$logWellDepth,
        data_area_split$data$logconcentration_plus_median, pol_train, threshold_val)
      cat("    Test:     "); print_acc_mcc(data_area_split$data_test$logWellDepth,
        data_area_split$data_test$logconcentration_plus_median, pol_test, threshold_val)

      
      pol_plss <- compute_indirect_policy(
        smoothers     = smoothers$smoothers_plss_RKHS,
        trt_bounds    = depth_range,
        threshold_val = threshold_val,
        krige_adjust  = plss$U,
        spatial       = spatial
      )
      return(list(pol_train = pol_train,
                  pol_test = pol_test,
                  pol_plss = pol_plss))
    }

    # ---- Indirect nonspatial ----
    if ("indirect_nonspatial" %in% methods_to_run) {
      cat("  Method: indirect_nonspatial\n")
      indirect_nonsp_pol_lst <- run_indirect(spatial = FALSE)
      sf_store[[area_key]][[threshold_label]][["indirect_nonspatial"]] <- plss_sf
      saveRDS(plss_sf, file.path(output_dir,
        sprintf("plss_sf_%s_%s_indirect_nonspatial.rds", area_key, threshold_label)))
    }

    # ---- Indirect spatial ----
    if ("indirect_spatial" %in% methods_to_run) {
      cat("  Method: indirect_spatial\n")
      indirect_pol_lst <- run_indirect(spatial = TRUE)
      sf_store[[area_key]][[threshold_label]][["indirect_spatial"]] <- plss_sf
      saveRDS(plss_sf, file.path(output_dir,
        sprintf("plss_sf_%s_%s_indirect_spatial.rds", area_key, threshold_label)))
    }

    # ---- Direct method ----
    if ("direct" %in% methods_to_run) {
      cat("  Method: direct\n")

      # Initialization via per-observation optimization
      cat("    Initializing direct policy...\n")
      init_vals   <- c(depth_range[1] - 1, mean(depth_range), depth_range[2] + 1)
      direct_init <- mcmapply(
        FUN      = find_best_policy,
        i        = 1:nrow(data_area_split$data),
        MoreArgs = list(
          initial_value_set         = init_vals,
          kernel_bw                 = kernel_bw,
          outcome_regression_object = outcome_reg,
          data                      = data_area_split$data,
          weights                   = weights_cbps$weights,
          smoothers                 = smoothers,
          depth_range               = depth_range,
          threshold_val             = threshold_val,
          clip_epsilon              = 20
        ),
        SIMPLIFY = TRUE
      )
      saveRDS(direct_init, file.path(output_dir,
        sprintf("direct_init_%s_%s.rds", area_key, threshold_label)))
      cat("    Init acc/mcc: "); print_acc_mcc(data_area_split$data$logWellDepth,
        data_area_split$data$logconcentration_plus_median, direct_init, threshold_val)

      # Project initialization onto kernel space
      init_glmnet <- glmnet(x = K[, -1], y = direct_init, alpha = 0, lambda = 0)
      params_init <- as.numeric(coef(init_glmnet))

      # RKHS optimization
      cat("    Optimizing...\n")
      kernel_optim <- optim(
        par              = params_init,
        fn               = compute_total_loss_smooth_RKHS,
        gr               = d_compute_total_loss_smooth_RKHS,
        K                = K,
        T                = data_area_split$data$logWellDepth,
        krige_adjust     = outcome_reg$krige_values,
        outcome_resid    = data_area_split$data$logconcentration_plus_median -
                           outcome_reg$pred - outcome_reg$krige_values,
        propensity_est   = 1 / weights_cbps$weights,
        lambda           = 0.25,
        smoothers        = smoothers$smoothers_RKHS,
        cumint_smoothers = smoothers$cumint_smoothers_RKHS,
        trt_bounds       = depth_range,
        threshold_val    = threshold_val,
        kernel_bw        = kernel_bw,
        clip_epsilon     = 30,
        surrogate_type   = "Gaussian",
        loss_type        = "db",
        method           = "L-BFGS-B",
        control          = list(maxit = 1000, trace = 3)
      )
      saveRDS(kernel_optim, file.path(output_dir,
        sprintf("kernel_optim_%s_%s.rds", area_key, threshold_label)))

      cat("    Non-DC Training: "); print_acc_mcc(data_area_split$data$logWellDepth,
        data_area_split$data$logconcentration_plus_median,
        K      %*% kernel_optim$par, threshold_val)
      cat("    Non-DC Test:     "); print_acc_mcc(data_area_split$data_test$logWellDepth,
        data_area_split$data_test$logconcentration_plus_median,
        K_test %*% kernel_optim$par, threshold_val)

      # DC algorithm (warm-started from non-DC result)
      cat("    Running DC algorithm...\n")
      dc_result <- run_dc_algorithm(
        params_initial   = kernel_optim$par,
        k_matrix         = K,
        data             = data_area_split$data,
        T                = data_area_split$data$logWellDepth,
        krige_adjust     = outcome_reg$krige_values,
        outcome_resid    = data_area_split$data$logconcentration_plus_median -
                           outcome_reg$pred - outcome_reg$krige_values,
        propensity_est   = 1 / weights_cbps$weights,
        lambda           = 0.25,
        smoothers        = smoothers$smoothers_RKHS,
        cumint_smoothers = smoothers$cumint_smoothers_RKHS,
        trt_bounds       = depth_range,
        threshold_val    = threshold_val,
        kernel_bw        = kernel_bw,
        clip_epsilon     = 30,
        clip_epsilon_bar = 50,
        tol              = 0.001,
        max_iter         = 100
      )
      saveRDS(dc_result, file.path(output_dir,
        sprintf("dc_result_%s_%s.rds", area_key, threshold_label)))

      cat("    DC Training: "); print_acc_mcc(data_area_split$data$logWellDepth,
        data_area_split$data$logconcentration_plus_median,
        K      %*% dc_result$par, threshold_val)
      cat("    DC Test:     "); print_acc_mcc(data_area_split$data_test$logWellDepth,
        data_area_split$data_test$logconcentration_plus_median,
        K_test %*% dc_result$par, threshold_val)

      plss_dir        <- plss
      plss_dir$policy <- as.numeric(K_plss %*% dc_result$par)
      plss_sf         <- st_as_sf(plss_dir, coords = c("longitude", "latitude"), crs = 4326)
      sf_store[[area_key]][[threshold_label]][["direct"]] <- plss_sf
      saveRDS(plss_sf, file.path(output_dir,
        sprintf("plss_sf_%s_%s_direct.rds", area_key, threshold_label)))
    }

  }  # end threshold loop

}  # end area loop


# ============================================================
# MAP VISUALIZATION
# ============================================================

for (threshold_val in threshold_vals) {
  threshold_label <- paste0("log", round(exp(threshold_val), 0))

  for (method in methods_to_run) {

    # Collect sf objects for this (threshold, method) across all areas
    sf_list <- Filter(Negate(is.null), lapply(names(sf_store), function(ak) {
      sf_store[[ak]][[threshold_label]][[method]]
    }))
    if (length(sf_list) == 0) next

    method_label <- gsub("_", " ", method)
    title_str    <- sprintf("%g mg/L Threshold — %s", exp(threshold_val), method_label)

    if (map_mode == "combined") {

      p <- ggplot()
      for (sf_obj in sf_list)
        p <- p + geom_sf(data = sf_obj, aes(color = policy), size = 0.2)
      p <- p +
        geom_sf(data = counties_wi, color = "black") +
        scale_color_viridis_c(
          option = "plasma",
          name   = "Minimum Depth (Feet)",
          labels = function(b) round(exp(b), 0)
        ) +
        labs(title = title_str, x = "Longitude", y = "Latitude") +
        theme_void()
      print(p)

    } else {  # individual

      for (i in seq_along(sf_list)) {
        ak <- names(sf_store)[i]
        p  <- ggplot() +
          geom_sf(data = sf_list[[i]], aes(color = policy), size = 0.2) +
          geom_sf(data = counties_wi, color = "black") +
          scale_color_viridis_c(
            option = "plasma",
            name   = "Minimum Depth (Feet)",
            labels = function(b) round(exp(b), 0)
          ) +
          labs(title = paste(title_str, "-", gsub("_", " ", ak)),
               x = "Longitude", y = "Latitude") +
          theme_void()
        print(p)
      }
    }

  }
}
