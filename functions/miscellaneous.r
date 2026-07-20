library(caret)
library(recipes)
library(e1071)
source('functions/funcs_kriging.R')
source('functions/loss_functions.R')

## I think this one using GP estimation (krige_values) from the outcome regression modeling is problematic

# get_generalized_propensity <- function(GPS_covariate_names = c("StaticLevel","crop_type_combine","drainagecl","precipitation","cafolog"),
#                                        data,krige_values){
#   
#   treatment_mean_model <- SuperLearner(
#     Y = data$logWellDepth,
#     X = cbind(data[,GPS_covariate_names], 
#               U = krige_values), # notice that we might add U
#     family = gaussian(),
#     SL.library = c("SL.ranger")
#   )
#   
#   treatment_sd = sd(data$logWellDepth - treatment_mean_model$Z)
#   gps_est <- dnorm(treatment_mean_model$Z, mean = treatment_mean_model$SL.predict, sd = treatment_sd)
#   gps_est <- pmax(gps_est, 0.02)
#   return(list(treatment_mean_model = treatment_mean_model,
#               treatment_sd = treatment_sd,
#               gps_est = gps_est))
# }


# x is the design matrix, y is the response
get_SVM_cv_residuals <- function(x,y,seed_value = 1998, fold = 5){
  # Setup Cross-Validation
  k <- fold  # Number of folds
  n <- nrow(x)
  set.seed(seed_value) # For reproducibility
  
  # Create random folds
  # 'folds' is a list where each element contains the indices for that fold
  folds <- split(sample(1:n), rep(1:k, length.out = n))
  
  # Initialize a vector to store CV predictions (same length as data)
  cv_preds <- numeric(n)
  
  # 3. Run the Loop
  for(i in 1:k) {
    test_indices <- folds[[i]]
    
    # Split data
    x_train <- x[-test_indices, ]
    x_test <- x[test_indices, ]
    
    y_train <- y[-test_indices]
    y_test <- y[test_indices]
    
    # Train SVM on training set
    # type="eps-regression" is standard for continuous outcomes
    # model <- svm(mpg ~ ., data = train_data, type = "eps-regression")
    model <- svm(x = x_train, y = y_train, type='eps-regression')
    
    # Predict on held-out test set
    cv_preds[test_indices] <- predict(model, newdata = x_test)
  }
  
  # 4. Calculate Cross-Validated Residuals
  # Residual = Observed - Predicted
  cv_residuals <- y - cv_preds
}


# using spatial mixed effect model to estimate propensity is much better (in terms of the distribution of the residual, looks more i.i.d.)
# Also, it make more sense and align with two stage estimator for geo-referenced paper
propensity_SVM_regresion <- function(GPS_covariate_names = c("StaticLevel","crop_type_combine","drainagecl","precipitation","cafolog"),
                                     data, trim_value = 0.01, tunning = FALSE){
  # 1. Define and prep the recipe on the TRAINING data
  rec <- recipe(~ ., data = data[, GPS_covariate_names]) %>%
    step_dummy(all_nominal_predictors()) %>%
    prep() # This "fits" the preprocessor
  
  # 2. Apply the *same* fitted recipe to both datasets
  design_matrix <- bake(rec, new_data = data[, GPS_covariate_names])
  
  if(tunning){
    gamma_default <- dim(design_matrix)[2]
    svm_auto <- best.svm(x = design_matrix, y = data$logWellDepth, type='nu-regression', gamma = gamma_default*c(0.25,0.5,1,1.5,2,4,8), tunecontrol = tune.control(cross = 5))
  } else {
    svm_auto <- svm(x = design_matrix, y = data$logWellDepth, type='nu-regression')
  }
  
  
  cv_residuals <- get_SVM_cv_residuals(x = design_matrix, y = data$logWellDepth)
  
  gpfit_depth = GpGp::fit_model(cv_residuals,
                               locs = data[,c("longitude", "latitude")],
                               covfun_name = "matern_sphere")
  
  krige_values_depth <- leave_one_out_kriging(locs = data[,c("longitude", "latitude")], 
                                              y_obs = cv_residuals, 
                                              gp_model = gpfit_depth$covfun_name,
                                              gp_params = gpfit_depth$covparms,
                                              order = c("coordinate")) 
  
  # if we directly use the residuals, we might under-estimate the gps_sigma
  gps_sigma <- sd(cv_residuals - krige_values_depth)
  # gps_est <- dnorm(data$logWellDepth, mean = predict(svm_auto,new_data = design_matrix) + krige_values_depth, sd = gps_sigma)
  gps_est <- dnorm(cv_residuals - krige_values_depth, sd = gps_sigma)
  gps_est <- pmax(gps_est, trim_value)
  
  return(list(svm_auto = svm_auto, 
              gpfit_depth = gpfit_depth,
              krige_values_depth = krige_values_depth, 
              gps_sigma = gps_sigma, 
              gps_est = gps_est))
  
}

# get_generalized_propensity <- function(GPS_covariate_names = c("StaticLevel","crop_type_combine","drainagecl","precipitation","cafolog"),
#                                       data){
#   
#   # # according to our assumption, we should use svm rather than random forest, so we use the function propensity_SVM_regresion()
#   gps_mean_model <- SuperLearner(
#     Y = data$logWellDepth,
#     X = data[, GPS_covariate_names],
#     family = gaussian(),
#     SL.library = c("SL.ranger")
#   )
#   
#   gpfit_gps <- GpGp::fit_model(data$logWellDepth - gps_mean_model$Z,
#                                locs = data[,c("longitude", "latitude")],
#                                covfun_name = "matern_sphere")
#   
#   krige_values_depth <- leave_one_out_kriging(locs = data[,c("longitude", "latitude")], 
#                                               y_obs = data$logWellDepth - gps_mean_model$Z, 
#                                               gp_model = gpfit_gps$covfun_name,
#                                               gp_params = gpfit_gps$covparms,
#                                               order = c("coordinate")) 
#   
#   gps_sigma <- sd(data$logWellDepth - gps_mean_model$Z - krige_values_depth)
#   gps_est <- dnorm(data$logWellDepth, mean = gps_mean_model$SL.predict, sd = gps_sigma)
#   
#   return(list(gps_mean_model = gps_mean_model, 
#               krige_values_depth = krige_values_depth, 
#               gps_sigma = gps_sigma, 
#               gps_est = gps_est))
#   
# }



# This function contains the work for a single iteration of your original loop
find_best_policy <- function(i, initial_value_set, kernel_bw, outcome_regression_object, data, weights, smoothers, depth_range, threshold_val, clip_epsilon) {
  
  resids <- data$logconcentration_plus_median - outcome_regression_object$pred - outcome_regression_object$krige_values
  
  # Use Inf for a more robust starting loss value
  tmp_loss <- Inf
  tmp_direct_policy <- NA # Default return value in case of errors
  
  # Inner loop to find the best initial point for the optimization
  for (initial_point in initial_value_set) {
    
    # Use try() to catch potential errors during optimization
    
    direct_trt_result <- optim(par = initial_point,
                               fn = compute_total_loss_smooth,
                               gr = d_compute_total_loss_smooth,
                               subject_idx = i,
                               # Pass other necessary arguments using ...
                               T = data$logWellDepth,
                               krige_adjust = outcome_regression_object$krige_values,
                               outcome_resid = resids,
                               propensity_est = 1/weights,
                               smoothers =  smoothers$smoothers_RKHS,
                               cumint_smoothers = smoothers$cumint_smoothers_RKHS, 
                               trt_bounds = depth_range,
                               threshold_val = threshold_val, 
                               kernel_bw = kernel_bw, 
                               clip_epsilon = clip_epsilon ,
                               surrogate_type = "Gaussian",
                               loss_type = "db",
                               method = "L-BFGS-B",
                               # Set trace = 0 to prevent garbled console output from parallel workers
                               control = list(maxit = 300, trace = 0)) 
    
    
    # Check if optim succeeded and if the new loss is an improvement
    if (direct_trt_result$value < tmp_loss) {
      tmp_loss <- direct_trt_result$value
      tmp_direct_policy <- direct_trt_result$par
    }
  }
  
  return(tmp_direct_policy)
}



get_smoothers_RKHS <- function(design_matrix,
                          design_matrix_test = NULL,
                          design_matrix_plss = NULL,
                          svm_auto,
                          depth_range,
                          treatment_step = 0.02){
  # xgboost use the design matrix
  cf_predictions_RKHS <- precompute_cf_predictions(data_obs = design_matrix, 
                                                   fit_outcome = svm_auto,
                                                   treatment_range = depth_range,
                                                   treatment_step = treatment_step,
                                                   treatment_name = "logWellDepth")

  cf_predictions_RKHS <- isotomic_correction(treatment_grid = cf_predictions_RKHS$treatment_values,
                                             pred_matrix_cf = cf_predictions_RKHS$pred_matrix)
  
  smoothers_RKHS <- make_prediction_smoothers(treatment_grid = cf_predictions_RKHS$treatment_values,
                                              pred_matrix_cf = cf_predictions_RKHS$pred_matrix,
                                              # smooth_method = "gam",
                                              smooth_method = "smooth.spline",
                                              # smooth_method = "splinefun",
                                              spar = NULL)
  cumint_smoothers_RKHS <- make_cumint_smoothers(treatment_grid = cf_predictions_RKHS$treatment_values,
                                                 pred_matrix_cf = cf_predictions_RKHS$pred_matrix,
                                                 # smooth_method = "splinefun",
                                                 # smooth_method = "gam",
                                                 smooth_method = "smooth.spline",
                                                 spar = NULL)
  if(is.null(design_matrix_test)){
    smoothers_test_RKHS <- NULL
    cumint_smoothers_test_RKHS <- NULL
  } else {
    cf_predictions_test_RKHS <- precompute_cf_predictions(data_obs = design_matrix_test, 
                                                          fit_outcome = svm_auto,
                                                          treatment_range = depth_range,
                                                          treatment_step = treatment_step,
                                                          treatment_name = "logWellDepth")
    
    cf_predictions_test_RKHS <- isotomic_correction(treatment_grid = cf_predictions_test_RKHS$treatment_values,
                                                    pred_matrix_cf = cf_predictions_test_RKHS$pred_matrix)
    
    smoothers_test_RKHS <- make_prediction_smoothers(treatment_grid = cf_predictions_test_RKHS$treatment_values,
                                                     pred_matrix_cf = cf_predictions_test_RKHS$pred_matrix,
                                                     smooth_method = "smooth.spline",
                                                     spar = NULL)
    
    cumint_smoothers_test_RKHS <- make_cumint_smoothers(treatment_grid = cf_predictions_test_RKHS$treatment_values,
                                                        pred_matrix_cf = cf_predictions_test_RKHS$pred_matrix,
                                                        # smooth_method = "splinefun",
                                                        smooth_method = "smooth.spline",
                                                        spar = NULL)
  }
  
  if(is.null(design_matrix_plss)){
    smoothers_plss_RKHS <- NULL
  } else {
    # we do this because the svm_auto requires same order of the columns
    design_matrix_plss$logWellDepth <- 0
    design_matrix_plss <- design_matrix_plss[,colnames(design_matrix)]
    cf_predictions_plss_RKHS <- precompute_cf_predictions(data_obs = design_matrix_plss, 
                                                          fit_outcome = svm_auto,
                                                          treatment_range = depth_range,
                                                          treatment_step = treatment_step,
                                                          treatment_name = "logWellDepth")
    
    cf_predictions_plss_RKHS <- isotomic_correction(treatment_grid = cf_predictions_plss_RKHS$treatment_values,
                                                    pred_matrix_cf = cf_predictions_plss_RKHS$pred_matrix)
    
    smoothers_plss_RKHS <- make_prediction_smoothers(treatment_grid = cf_predictions_plss_RKHS$treatment_values,
                                                     pred_matrix_cf = cf_predictions_plss_RKHS$pred_matrix,
                                                     smooth_method = "smooth.spline",
                                                     spar = NULL)
    
  }
  
  return(list(smoothers_RKHS = smoothers_RKHS,
         cumint_smoothers_RKHS = cumint_smoothers_RKHS,
         smoothers_test_RKHS = smoothers_test_RKHS,
         cumint_smoothers_test_RKHS = cumint_smoothers_test_RKHS,
         smoothers_plss_RKHS = smoothers_plss_RKHS))
  
}

outcome_regression_SVM <- function(data,data_test = NULL, tunning = FALSE,
                                   OR_covariate_names = c("logWellDepth","crop_type_combine","drainagecl","precipitation","cafolog","StaticLevel")){
  
  
  # 1. Define and prep the recipe on the TRAINING data
  rec <- recipe(~ ., data = data[, OR_covariate_names]) %>%
    step_dummy(all_nominal_predictors()) %>%
    prep() # This "fits" the preprocessor
  
  # 2. Apply the *same* fitted recipe to both datasets
  design_matrix <- bake(rec, new_data = data[, OR_covariate_names])
  
  if(tunning){
    gamma_default <- dim(design_matrix)[2]
    svm_auto <- best.svm(x = design_matrix, y = data$logconcentration_plus_median, type='nu-regression', gamma = gamma_default*c(0.25,0.5,1,1.5,2,4,8), tunecontrol = tune.control(cross = 5))
  } else {
    svm_auto <- svm(x = design_matrix, y = data$logconcentration_plus_median, type='nu-regression')
  }
  # svm_auto <- svm(x = design_matrix, y = data$logconcentration_plus_median, type='nu-regression')
  # svm_auto <- best.svm(x = design_matrix, y = data$logconcentration_plus_median, type='nu-regression', gamma = (1/20)*c(0.5,1,1.5,2,4,8), tunecontrol = tune.control(cross = 5))
  # best.svm(x = design_matrix, y = data$logconcentration_plus_median, type='nu-regression')
  
  gpfit_RKHS = GpGp::fit_model(svm_auto$residuals,
                               locs = data[,c("longitude", "latitude")],
                               covfun_name = "matern_sphere",
                               convtol = 1e-03)
  
  krige_values_RKHS <- leave_one_out_kriging(locs = data[,c("longitude", "latitude")], 
                                             y_obs = svm_auto$residuals, 
                                             gp_model = gpfit_RKHS$covfun_name,
                                             gp_params = gpfit_RKHS$covparms,
                                             order = c("coordinate")) 
  
  pred <- predict(svm_auto, newdata = design_matrix)
  if(is.null(data_test)){
    krige_values_test_RKHS <- NULL
    pred_test <- NULL
  } else {
    design_matrix_test   <- bake(rec, new_data = data_test[, OR_covariate_names])
    pred_test <- predict(svm_auto, newdata = design_matrix_test)
    krige_values_test_RKHS <- GpGp::predictions(fit = gpfit_RKHS, 
                                                locs_pred = data_test[,c("longitude", "latitude")], 
                                                X_pred = rep(1,nrow(data_test)))
  }
  
  
  return(list(svm=svm_auto, gpfit=gpfit_RKHS,
              pred = pred, pred_test = pred_test, 
              krige_values = krige_values_RKHS, 
              krige_values_test = krige_values_test_RKHS,
              design_matrix = design_matrix,
              design_matrix_test = design_matrix_test))
}

split_nitrate_data <- function(data, p = 0.7, seed = 1998){
  set.seed(seed)
  data_full <- data
  # train_indices <- sample(nrow(data_full), nrow(data_full) * 0.7) # we would replace this with a stratified sample splitting
  # we partition to balance the distribution of nitrate
  train_indices <- createDataPartition(data$logconcentration_plus_median, list = FALSE, p = p) # from caret package
  data <- data_full[train_indices,]
  data_test <- data_full[-train_indices,]
  return(list(data=data, data_test=data_test))
}

# we also add log-transformation, and area in this function
load_nitrate_data <- function(file_path = "data/data_Nitrate_with_covar_median_well.csv", zero_inflated = FALSE){
  command <- paste("brctl download", shQuote(file_path))
  system(command)
  data <- read.csv(file_path)
  # Assuming 'data' is your data frame with a 'County' column.
  # Initialize the new 'area' column.
  data$area <- ""
  
  ## Northern Districts
  indices_northwest <- data$County %in% c("Douglas", "Bayfield", "Ashland", "Iron", 
                                          "Washburn", "Sawyer", "Burnett", "Polk", 
                                          "Barron", "Rusk", "St. Croix", "Dunn", 
                                          "Chippewa", "Pierce", "Eau Claire")
  
  data[indices_northwest, "area"] <- "Northwest"
  
  indices_northCentral <- data$County %in% c("Price", "Vilas", "Oneida", "Lincoln", 
                                             "Langlade", "Taylor", "Marathon", "Clark")
  data[indices_northCentral, "area"] <- "North Central"
  
  indices_northEast <- data$County %in% c("Florence", "Forest", "Marinette", "Oconto", 
                                          "Menominee", "Shawano", "Door", "Kewaunee")
  data[indices_northEast, "area"] <- "North East"
  
  ## Central Districts
  indices_westCentral <- data$County %in% c("Pepin", "Jackson", "Buffalo", "Trempealeau", 
                                            "La Crosse", "Monroe", "Eau Claire", "Dunn", "Pierce", "St. Croix")
  
  data[indices_westCentral, "area"] <- "West Central"
  
  indices_central <- data$County %in% c("Wood", "Portage", "Waupaca", "Juneau", 
                                        "Adams", "Waushara", "Marquette", "Green Lake")
  data[indices_central, "area"] <- "Central"
  
  indices_eastCentral <- data$County %in% c("Outagamie", "Winnebago", "Door", "Fond du Lac", 
                                            "Brown", "Calumet","Sheboygan","Manitowoc","Kewaunee")
  
  data[indices_eastCentral, "area"] <- "East Central"
  
  ## Southern Districts
  indices_southWest <- data$County %in% c("Vernon", "Crawford", "Grant", "Richland", 
                                          "Sauk", "Iowa", "Lafayette")
  data[indices_southWest, "area"] <- "South West"
  
  indices_southCentral <- data$County %in% c("Columbia", "Dodge", "Dane", "Jefferson", 
                                             "Green", "Rock")
  data[indices_southCentral, "area"] <- "South Central"
  
  indices_southEast <- data$County %in% c("Washington", "Waukesha", "Walworth", "Ozaukee", 
                                          "Milwaukee", "Racine", "Kenosha")
  data[indices_southEast, "area"] <- "South East"
  
  
  data$logWellDepth = log(data$WellDepth)
  data$logconcentration_plus_median = log(data$concentration_plus_median)
  if(!zero_inflated){
    data <- data[data$concentration_plus_median>0.5,]
  }
  return(data)
}

get_latest_record <- function(data){
  
  # by well ID
  latest_records <- data %>%
    # Make sure the date column is a true Date type for correct sorting
    mutate(SampleDate = as.Date(SampleDate)) %>%
    # Group the data by the well ID
    group_by(WellID) %>%
    # Within each group, keep only the row with the maximum date
    slice_max(order_by = SampleDate, n = 1) %>%
    # Optional: remove the grouping
    ungroup()
  
  
  # by well depth, longitude, latitude
  latest_records_depth <- latest_records %>%
    # Ensure the date column is a true Date type
    mutate(SampleDate = as.Date(SampleDate)) %>%
    # Group by all four identifying columns
    group_by(WellDepth, latitude, longitude) %>%
    # Keep the row with the latest date within each group
    slice_max(order_by = SampleDate, n = 1) %>%
    ungroup()
  
  return(latest_records_depth)
  
}

get_latest_and_hitorical_nitrate_record <- function(data){
  
  data_latest_oldest <- data %>% 
    mutate(SampleDate = as.Date(SampleDate)) %>%
    group_by(WellDepth, latitude, longitude, WellUse) %>%
    mutate(sampling_order = row_number(desc(SampleDate))) %>% 
    arrange(latitude, longitude, WellDepth, sampling_order) %>% 
    ungroup()
  
  data_latest <- data_latest_oldest[data_latest_oldest$sampling_order == 1,]
  data_historical <- data_latest_oldest[data_latest_oldest$sampling_order != 1,]
  return(list(latest_nitrate_record = data_latest, historical_nitrate_record = data_historical))
}


run_dc_algorithm <- function(params_initial,
                             k_matrix,
                             data,
                             T,
                             krige_adjust,
                             outcome_resid,
                             propensity_est,
                             lambda, 
                             smoothers,
                             cumint_smoothers,
                             trt_bounds,
                             threshold_val, 
                             kernel_bw, 
                             clip_epsilon,
                             clip_epsilon_bar,
                             tol = 0.001,
                             max_iter = 100){
  
  params_old <- params_initial
  
  for(k in 1:max_iter){
    
    cat(k,'-th iteration of DC algorithm','\n')
    
    sub_result <- optim(par = params_old, 
                        fn = compute_total_loss_smooth_dc_approx_RKHS, 
                        gr = d_compute_total_loss_smooth_dc_approx_RKHS, 
                        params_old = params_old,
                        K = k_matrix,
                        T = T,
                        krige_adjust = krige_adjust,
                        outcome_resid = outcome_resid,
                        propensity_est = propensity_est,
                        lambda = lambda, 
                        smoothers = smoothers,
                        cumint_smoothers = cumint_smoothers,
                        trt_bounds = trt_bounds,
                        threshold_val = threshold_val, 
                        kernel_bw= kernel_bw, 
                        clip_epsilon = clip_epsilon,
                        clip_epsilon_bar = clip_epsilon_bar,
                        surrogate_type = "Gaussian",
                        loss_type = "db",
                        method = "L-BFGS-B",
                        control = list(maxit =300, trace = 3))
    
    params <- sub_result$par
    
    # Check for convergence
    change <- sqrt(sum((params - params_old)^2))
    cat(sprintf("Iter: %d, Change: %f, Objective: %f\n", k, change, sub_result$value))
    
    if (change < tol) {
      cat("--- DC Algorithm Converged ---\n")
      converged <- TRUE
      final_iter <- k
      params_old <- params
      break
    }
    
    if(k==max_iter){
      converged <- FALSE
      final_iter <- k
      params_old <- params
      break
    }
    
    params_old <- params
  }
  return(list(par = params_old,
              optim_final = sub_result,
              iter = final_iter,
              converged = converged))  
}

compute_indirect_policy <- function(smoothers,trt_bounds,threshold_val,krige_adjust=NULL,spatial=FALSE){

  
  if(trt_bounds[1]>=trt_bounds[2]){
    stop("upper bound smaller or less than lower bound")
  }
  
  n <- length(smoothers)
  if(spatial){
    if(length(krige_adjust)!=length(smoothers)){
      stop("number of smoothers not equal numer of krige_adjust")
    }
    indirect_policy <- sapply(1:n, function(i) {
      if(smoothers[[i]](trt_bounds[1]) + krige_adjust[i] <= threshold_val) return(trt_bounds[1])
      if(smoothers[[i]](trt_bounds[2]) + krige_adjust[i]> threshold_val) return(trt_bounds[2])
      res <- uniroot(f = function(t) {smoothers[[i]](t) + krige_adjust[i] - threshold_val},
                     lower = trt_bounds[1],
                     upper = trt_bounds[2])
      return(res$root)
    })
  } else {
    indirect_policy <- sapply(1:n, function(i) {
      if(smoothers[[i]](trt_bounds[1]) <= threshold_val) return(trt_bounds[1])
      if(smoothers[[i]](trt_bounds[2]) > threshold_val) return(trt_bounds[2])
      res <- uniroot(f = function(t) {smoothers[[i]](t)- threshold_val},
                     lower = trt_bounds[1],
                     upper = trt_bounds[2])
      return(res$root)
    })
  }
  
}


get_direct_indirect_method_comparison <- function(data_split, indirect_policy, indirect_policy_test, plss_covariates_sf){
  
  plss_covariates <- plss_covariates_sf %>%
    mutate(
      longitude = st_coordinates(.)[,1],
      latitude = st_coordinates(.)[,2]
    ) %>%
    st_drop_geometry()
  
  indirect <- rbind(cbind(data_split$data[,c("Township", "Range", "Range.Direction", "Section","longitude", "latitude","concentration_plus_median","WellDepth")], indirect_policy = indirect_policy, test = rep(FALSE, length(indirect_policy))),
                    cbind(data_split$data_test[,c("Township", "Range", "Range.Direction", "Section","longitude", "latitude","concentration_plus_median","WellDepth")],  indirect_policy = indirect_policy_test, test = rep(TRUE, length(indirect_policy_test))))
  
  # indirect_sf <- st_as_sf(indirect,coords = c("longitude", "latitude"),crs = 4326)
  
  policy_comparison <- inner_join(indirect, plss_covariates, by = c("Township", "Range", "Range.Direction", "Section","longitude","latitude"))
  return(policy_comparison)
  # policy_comparison_sf <- st_as_sf(policy_comparison, coords = c("longitude", "latitude"),crs = 4326)
}

print_acc_mcc <- function(T, Y, policy, threshold_val){
  acc <- Metrics::accuracy(Y > threshold_val, policy  > T )
  mcc <- yardstick::mcc_vec(truth    = factor(Y      > threshold_val, levels = c(FALSE, TRUE)),
                             estimate = factor(policy > T,             levels = c(FALSE, TRUE)))
  cat('The accuracy is ',acc, '\n')
  cat('The mcc is ',mcc, '\n')
}


calculate_acc_mcc <- function(T, Y, policy, threshold_val){
  acc <- Metrics::accuracy(Y > threshold_val, policy  > T )
  mcc <- yardstick::mcc_vec(truth    = factor(Y      > threshold_val, levels = c(FALSE, TRUE)),
                             estimate = factor(policy > T,             levels = c(FALSE, TRUE)))
  return(list(acc = acc, mcc = mcc))
}

calculate_acc_mcc_two_sided_f1 <- function(T, Y, policy, threshold_val){
  acc <- Metrics::accuracy(Y > threshold_val, policy  > T )
  mcc <- yardstick::mcc_vec(truth = as.factor(Y > threshold_val), estimate = as.factor(policy  > T))
  
  
  score_pos <- Metrics::f1(Y > threshold_val, policy  > T)
  
  # 2. Calculate Inverse F1 (for class 0)
  # By doing (1 - variable), we turn 0s into 1s and 1s into 0s
  
  tp <- sum((Y > threshold_val) & (policy > T))
  fp <- sum((Y <= threshold_val) & (policy > T))
  tn <- sum((Y <= threshold_val) & (policy <= T))
  fn <- sum((Y > threshold_val) & (policy <= T))
  
  
  # 3. Average them for the Two-Sided (Macro) F1
  two_sided_f1 <- 0.5 * ((2* tp)/(2*tp + fp + fn) + (2* tn)/(2*tn + fn + fp))
  
  return(list(acc = acc, mcc = mcc, two_sided_f1 = two_sided_f1))
}

print_acc_mcc_two_sided_f1 <- function(T, Y, policy, threshold_val){
  acc <- Metrics::accuracy(Y > threshold_val, policy  > T )
  mcc <- yardstick::mcc_vec(truth = as.factor(Y > threshold_val), estimate = as.factor(policy  > T))
  
  
  score_pos <- Metrics::f1(Y > threshold_val, policy  > T)
  
  # 2. Calculate Inverse F1 (for class 0)
  # By doing (1 - variable), we turn 0s into 1s and 1s into 0s
  
  tp <- sum((Y > threshold_val) & (policy > T))
  fp <- sum((Y <= threshold_val) & (policy > T))
  tn <- sum((Y <= threshold_val) & (policy <= T))
  fn <- sum((Y > threshold_val) & (policy <= T))
  
  
  # 3. Average them for the Two-Sided (Macro) F1
  two_sided_f1 <- 0.5 * ((2* tp)/(2*tp + fp + fn) + (2* tn)/(2*tn + fn + fp))
  cat('acc: ',acc, '\n')
  cat('mcc: ',mcc, '\n')
  cat('two_sided_f1: ',two_sided_f1, '\n')
  return(0)
}






