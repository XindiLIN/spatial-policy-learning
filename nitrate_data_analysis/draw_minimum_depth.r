##### generate well depth for given threshold


## we use the data from 2020-2024 non-zero-inflated data to train the model
## we are considering using the data from 2025 as the test data
## But for now, we would stick to our original plan
## we first split the data into different areas, then combine them, it can benefit both computation due to the sample size and performance due to the
## heterogeneity of the data

## the covariate balancing has issue with whole central area data, and I don't know why
## the package is estimating the stabilized weight, which is not the propensity we need, so we try the regression-based propensity score

data <- load_nitrate_data(zero_inflated = FALSE)
data <- data[data$SampleYear>=2020,]
# data_2025_path <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/data_Nitrate_with_covar_median_well_2025.csv"
# data_test <- load_nitrate_data(file_path = data_2025_path, zero_inflated = FALSE)

threshold_lst <- c(log(5), log(10))
area_lst <- c("North Central","West Central","Northwest","Central","South West", "South Central","South East","East Central","North East")


for(threshold_val in threshold_lst){
  for(area in area_lst){
    
    
    # data_area <- data[data$area==area,]
    # data_test_area <- data_test[data_test$area==area,]
    # data_area_split <- list(data=data_area, data_test=data_test_area)
    
    
    data_area <- data[data$area==area,]
    data_area_split <- split_nitrate_data(data_area)
    
    
    outcome_regression_area <- outcome_regression_SVM(data = data_area_split$data, data_test = data_area_split$data_test,tunning = FALSE)
    
    
    
    cat('Training: ')
    Metrics::rmse(data_area_split$data$logconcentration_plus_median, outcome_regression_area$pred + outcome_regression_area$krige_values)
    Metrics::rmse(data_area_split$data$logconcentration_plus_median, outcome_regression_area$pred)
    cat('\n')
    cat('Test: ')
    Metrics::rmse(data_area_split$data_test$logconcentration_plus_median, outcome_regression_area$pred_test + outcome_regression_area$krige_values_test)
    Metrics::rmse(data_area_split$data_test$logconcentration_plus_median, outcome_regression_area$pred_test)
    
    
    ###### smoothers
    
    
    depth_range_area <- c(min(data_area_split$data$logWellDepth), max(data_area_split$data$logWellDepth))
    
    
    
    smoothers_area <- get_smoothers_RKHS(design_matrix = outcome_regression_area$design_matrix,
                                               design_matrix_test = outcome_regression_area$design_matrix_test,
                                               svm_auto = outcome_regression_area$svm,
                                               depth_range = depth_range_area)
    
    
    # The above code is very slow
    
    ###### indirect method
    
    
    
    indirect_policy_RKHS_area <- compute_indirect_policy(smoothers = smoothers_area$smoothers_RKHS, trt_bounds = depth_range_area, threshold_val = threshold_val, krige_adjust = outcome_regression_area$krige_values, spatial = TRUE)
    indirect_policy_test_RKHS_area <- compute_indirect_policy(smoothers = smoothers_area$smoothers_test_RKHS, trt_bounds = depth_range_area, threshold_val = threshold_val, krige_adjust = outcome_regression_area$krige_values_test, spatial = TRUE)
    
    
    cat('Indirect Training data:\n')
    print_acc_mcc(T = data_area_split$data$logWellDepth, 
                  Y = data_area_split$data$logconcentration_plus_median,
                  policy = indirect_policy_RKHS_area,
                  threshold_val = threshold_val)
    
    
    
    
    
    cat('Indirect Test data:\n')
    print_acc_mcc(T = data_area_split$data_test$logWellDepth, 
                  Y = data_area_split$data_test$logconcentration_plus_median,
                  policy = indirect_policy_test_RKHS_area,
                  threshold_val = threshold_val)
    
    # The mcc for threshold be log(10) is worse than that of log(5), but the accurary is better
    
    ###### propensity
    
    ####### covariate balance
    # They estimate the stabilized propensity but not the reciprocal of the conditional density, so they are not applicable

    # weights_cbps_area <- weightit(logWellDepth ~ StaticLevel+crop_type_combine+drainagecl+precipitation+cafolog,
    #                                     data = data_area_split$data, method = "cbps")

    # summary(weights_cbps_area$weights)
    
    # we would turn to estimate  conditional density
    
    propensity_regression_area <- propensity_SVM_regresion(data = data_area_split$data)
    
    
    ###### initialization
    
    initial_value_set <- c(depth_range_area[1] - 1, mean(depth_range_area), depth_range_area[2] + 1)
    kernel_bw_area <- 1.06 * sd(data_area_split$data$logWellDepth) *nrow(data_area_split$data)^(-1/5)
    
    
    direct_individual_multi_initial_area <- mcmapply(FUN = find_best_policy, 
                                                           i = 1:nrow(data_area_split$data), 
                                                           MoreArgs = list(initial_value_set = initial_value_set, 
                                                                           kernel_bw = kernel_bw_area,
                                                                           outcome_regression_object = outcome_regression_area, 
                                                                           data = data_area_split$data, 
                                                                           # weights = weights_cbps_area$weights,
                                                                           # weights = 1/gps_est,
                                                                           weights = 1/propensity_regression_area$gps_est,
                                                                           smoothers = smoothers_area, 
                                                                           depth_range = depth_range_area, 
                                                                           threshold_val = threshold_val, 
                                                                           clip_epsilon = 20), SIMPLIFY = TRUE)
    # saveRDS(direct_individual_multi_initial_area, "direct_individual_multi_initial_area.rds")
    # direct_individual_multi_initial_area <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/direct_individual_multi_initial_area.rds")
    
    print_acc_mcc(T = data_area_split$data$logWellDepth, 
                  Y = data_area_split$data$logconcentration_plus_median,
                  policy = direct_individual_multi_initial_area,
                  threshold_val = threshold_val)
    
    
    ###### select kernel bandwidth
    
    
    
    # there is no well use, and staticlevel 
    RKHS_covariate_names = c("StaticLevel", "crop_type_combine","drainagecl","precipitation","cafolog")
    rec <- recipe(~ ., data = data_area_split$data[, RKHS_covariate_names]) %>%
      step_dummy(all_nominal_predictors()) %>%
      prep() # This "fits" the preprocessor
    
    
    
    kernel_design_matrix <- bake(rec, new_data = data_area_split$data[, RKHS_covariate_names])
    kernel_design_matrix <- cbind(kernel_design_matrix, U = outcome_regression_area$krige_values)
    
    kernel_design_matrix_test <- bake(rec, new_data = data_area_split$data_test[, RKHS_covariate_names])
    kernel_design_matrix_test <- cbind(kernel_design_matrix_test, U = outcome_regression_area$krige_values_test)
    
    train_means <- colMeans(kernel_design_matrix)
    train_sds <- apply(kernel_design_matrix, 2, sd)
    
    kernel_design_matrix <- scale(kernel_design_matrix)
    kernel_design_matrix_test <- scale(kernel_design_matrix_test, center = train_means, scale = train_sds)
    
    
    
    for(m in c(-4,-3,-2,-1.5)){
      cat("m is: ",m,'\n')
      gamma <- 2^(m) * median(rdist(kernel_design_matrix))
      rbf <- rbfdot(sigma = 1/(gamma^2))
      k_matrix <- kernelMatrix(rbf, kernel_design_matrix)
      k_matrix <- cbind(rep(1,nrow(k_matrix)), k_matrix)
      
      k_matrix_test <- kernelMatrix(rbf, kernel_design_matrix_test, kernel_design_matrix)
      k_matrix_test <- cbind(rep(1,nrow(k_matrix_test)), k_matrix_test)
      
      initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_area, alpha = 0, lambda = 0)
      params_initial_glmnet <- coef(initial_glmnet)
      params_initial_glmnet <- as.numeric(params_initial_glmnet)
      cat('training: ')
      print_acc_mcc(T = data_area_split$data$logWellDepth, 
                    Y = data_area_split$data$logconcentration_plus_median, 
                    policy = k_matrix %*% params_initial_glmnet, 
                    threshold_val = threshold_val)
      
      cat('test: ')
      print_acc_mcc(T = data_area_split$data_test$logWellDepth, 
                    Y = data_area_split$data_test$logconcentration_plus_median, 
                    policy = k_matrix_test %*% params_initial_glmnet, 
                    threshold_val = threshold_val)
      
    }
    
    
    ###### non-DC algorithm
    
    # We will choose m is:  -3 
    
    
    # cat("lambda is: ",lambda,'\n')
    m <- -2
    lambda <- 0.25
    clip_epsilon <- 30
    gamma <- 2^(m) * median(fields::rdist(kernel_design_matrix))
    rbf <- rbfdot(sigma = 1/(gamma^2))
    k_matrix <- kernelMatrix(rbf, kernel_design_matrix)
    k_matrix <- cbind(rep(1,nrow(k_matrix)), k_matrix)
    
    k_matrix_test <- kernelMatrix(rbf, kernel_design_matrix_test, kernel_design_matrix)
    k_matrix_test <- cbind(rep(1,nrow(k_matrix_test)), k_matrix_test)
    
    initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_area, alpha = 0, lambda = 0)
    params_initial_glmnet <- coef(initial_glmnet)
    params_initial_glmnet <- as.numeric(params_initial_glmnet)
    
    kernel_optim_trt_DB_area <- optim(par = params_initial_glmnet, 
                                            fn = compute_total_loss_smooth_RKHS,
                                            gr = d_compute_total_loss_smooth_RKHS,
                                            K = k_matrix,
                                            T = data_area_split$data$logWellDepth,
                                            krige_adjust = outcome_regression_area$krige_values,
                                            outcome_resid = data_area_split$data$logconcentration_plus_median - outcome_regression_area$pred - outcome_regression_area$krige_values,
                                            # propensity_est = pmax(generalized_propensity_area$gps_est, 0.0001),
                                            # propensity_est = 1/weights_cbps_area$weights,
                                            propensity_est = gps_est, 
                                            lambda = lambda, 
                                            smoothers = smoothers_area$smoothers_RKHS,
                                            cumint_smoothers = smoothers_area$cumint_smoothers_RKHS,
                                            trt_bounds = depth_range_area,
                                            threshold_val = threshold_val, 
                                            kernel_bw = kernel_bw_area, 
                                            clip_epsilon = clip_epsilon,
                                            surrogate_type = "Gaussian",
                                            loss_type = "db",
                                            method = "L-BFGS-B",
                                            control = list(maxit = 1000, trace = 3)) 
    
    cat('Training: \n')
    print_acc_mcc(T = data_area_split$data$logWellDepth, 
                  Y = data_area_split$data$logconcentration_plus_median,
                  policy = k_matrix %*% kernel_optim_trt_DB_area$par,
                  threshold_val = threshold_val)
    
    # 11/18
    # The accuracy is  0.8986333 
    # The mcc is  0.6881898 
    
    cat('Test: \n')
    print_acc_mcc(T = data_area_split$data_test$logWellDepth, 
                  Y = data_area_split$data_test$logconcentration_plus_median, 
                  policy = k_matrix_test %*% kernel_optim_trt_DB_area$par, 
                  threshold_val = threshold_val)
    
    # 11.18
    # The accuracy is  0.8445982 
    # The mcc is  0.5060648 
    
    # Training: 
    # The accuracy is  0.7718797 
    # The mcc is  0.5568953
    # 
    # Test:
    # The accuracy is  0.7464252 
    # The mcc is  0.502661
    
    saveRDS(kernel_optim_trt_DB_area, 'kernel_optim_trt_DB_area.rds')
    
    
    ###### DC algorithm
    
    m <- -2
    lambda <- 0.25
    clip_epsilon <- 30
    clip_epsilon_bar <- 50
    tol <- 0.01
    gamma <- 2^(m) * median(fields::rdist(kernel_design_matrix))
    rbf <- rbfdot(sigma = 1/(gamma^2))
    k_matrix <- kernelMatrix(rbf, kernel_design_matrix)
    k_matrix <- cbind(rep(1,nrow(k_matrix)), k_matrix)
    
    k_matrix_test <- kernelMatrix(rbf, kernel_design_matrix_test, kernel_design_matrix)
    k_matrix_test <- cbind(rep(1,nrow(k_matrix_test)), k_matrix_test)
    
    initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_area, alpha = 0, lambda = 0)
    params_initial_glmnet <- coef(initial_glmnet)
    params_initial_glmnet <- as.numeric(params_initial_glmnet)
    
    max_iter <- 100
    
    params_old <- params_initial_glmnet
    
    dc_result <- run_dc_algorithm(params_initial = params_initial_glmnet,
                                 k_matrix = k_matrix,
                                 data,
                                 T = data_area_split$data$logWellDepth,
                                 krige_adjust = outcome_regression_area$krige_values,
                                 outcome_resid = data_area_split$data$logconcentration_plus_median - outcome_regression_area$pred - outcome_regression_area$krige_values,
                                 propensity_est = gps_est,
                                 lambda = lambda, 
                                 smoothers = smoothers_area$smoothers_RKHS,
                                 cumint_smoothers = smoothers_area$cumint_smoothers_RKHS,
                                 trt_bounds = depth_range_area,
                                 threshold_val = threshold_val, 
                                 kernel_bw= kernel_bw_area, 
                                 clip_epsilon = clip_epsilon,
                                 clip_epsilon_bar = clip_epsilon_bar,
                                 tol = 0.001,
                                 max_iter = 100)
    
    for(k in 1:max_iter){
      
      cat(k,'-th iteration of DC algorithm','\n')
      
      sub_result <- optim(par = params_old, 
                          fn = compute_total_loss_smooth_dc_approx_RKHS, 
                          gr = d_compute_total_loss_smooth_dc_approx_RKHS, 
                          params_old = params_old,
                          K = k_matrix,
                          T = data_area_split$data$logWellDepth,
                          krige_adjust = outcome_regression_area$krige_values,
                          outcome_resid = data_area_split$data$logconcentration_plus_median - outcome_regression_area$pred - outcome_regression_area$krige_values,
                          propensity_est = gps_est,
                          lambda = lambda, 
                          smoothers = smoothers_area$smoothers_RKHS,
                          cumint_smoothers = smoothers_area$cumint_smoothers_RKHS,
                          trt_bounds = depth_range_area,
                          threshold_val = threshold_val, 
                          kernel_bw= kernel_bw_area, 
                          clip_epsilon = clip_epsilon,
                          clip_epsilon_bar = clip_epsilon_bar,
                          surrogate_type = "Gaussian",
                          loss_type = "db",
                          method = "L-BFGS-B",
                          control = list(maxit =200, trace = 3))
      
      params <- sub_result$par
      
      # Check for convergence
      change <- sqrt(sum((params - params_old)^2))
      cat(sprintf("Iter: %d, Change: %f, Objective: %f\n", k, change, sub_result$value))
      
      if (change < tol) {
        cat("--- DC Algorithm Converged ---\n")
        params_old <- params
        break
      }
      params_old <- params
    }
    
    
    cat('Training: \n')
    print_acc_mcc(T = data_area_split$data$logWellDepth, 
                  Y = data_area_split$data$logconcentration_plus_median,
                  policy = k_matrix %*% dc_result$par,
                  threshold_val = threshold_val)
    
    # 11/18
    # The accuracy is  0.8986333 
    # The mcc is  0.6881898 
    
    cat('Test: \n')
    print_acc_mcc(T = data_area_split$data_test$logWellDepth, 
                  Y = data_area_split$data_test$logconcentration_plus_median, 
                  policy = k_matrix_test %*% dc_result$par,
                  threshold_val = threshold_val)
    
    
    # after the loop, we assign the sub_result as the kernel_optim_trt_DB_area to keep the variable name consistent
    kernel_optim_trt_DB_area <- dc_result$optim_final
    
    
    ###### visualize depth map
    
    
    file_path <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/plss_covariates.csv"
    command <- paste("brctl download", shQuote(file_path))
    system(command)
    plss_covariates <- read.csv(file_path)
    
    
    
    plss_covariates_area <- plss_covariates[plss_covariates$County %in% c("Vernon", "Crawford", "Grant", "Richland", 
                                                                                "Sauk", "Iowa", "Lafayette"),]
    
    
    
    plss_covariates_area$cafolog <- log(plss_covariates_area$cafo + 1)
    
    
    
    dim(plss_covariates_area)
    dim(na.omit(plss_covariates_area))
    
    # 5491   16
    # 3728   16
    
    plss_covariates_area <- na.omit(plss_covariates_area)
    
    
    # There are two plss section w with 'Vegtables' as their crop_type_combine, which is present in training data, we will just remove them
    
    
    plss_covariates_area <- plss_covariates_area[!plss_covariates_area$crop_type_combine %in% c("Vegtables"),]
    
    
    
    ####### Interpolate static water level
    
    
    gpfit_staticLevel_area <- GpGp::fit_model(data_area$StaticLevel,
                                                    locs = data_area[,c("longitude", "latitude")],
                                                    covfun_name = "matern_sphere")
    
    plss_covariates_area$StaticLevel <-  GpGp::predictions(fit = gpfit_staticLevel_area,
                                                                 locs_pred = plss_covariates_area[,c("longitude", "latitude")], 
                                                                 X_pred = rep(1,nrow(plss_covariates_area)))
    
    
    ####### get kriging values for plss
    
    
    plss_krige_values_area <- GpGp::predictions(fit = outcome_regression_area$gpfit, locs_pred = plss_covariates_area[,c("longitude", "latitude")], X_pred = rep(1,nrow(plss_covariates_area)))
    
    
    ####### build kernel matrices
    
    
    RKHS_covariate_names = c("StaticLevel", "crop_type_combine","drainagecl","precipitation","cafolog")
    rec <- recipe(~ ., data = data_area_split$data[, RKHS_covariate_names]) %>%
      step_dummy(all_nominal_predictors()) %>%
      prep() # This "fits" the preprocessor
    
    
    
    plss_kernel_design_matrix_area <- bake(rec, new_data = plss_covariates_area[,  RKHS_covariate_names])
    
    plss_kernel_design_matrix_area <- cbind(plss_kernel_design_matrix_area, U = plss_krige_values_area)
    plss_kernel_design_matrix_area <- scale(plss_kernel_design_matrix_area, center = train_means, scale = train_sds)
    
    
    
    plss_k_matrix_area <- kernelMatrix(rbf, plss_kernel_design_matrix_area, kernel_design_matrix)
    plss_k_matrix_area <- cbind(rep(1,nrow(plss_k_matrix_area)), plss_k_matrix_area)
    
    
    ####### generate map
    
    
    plss_covariates_area$policy <- plss_k_matrix_area %*% kernel_optim_trt_DB_area$par
    
    
    
    plss_covariates_area_sf <- st_as_sf(plss_covariates_area, coords = c('longitude','latitude'), crs=4326)
    
    counties <- counties(state = "WI", cb = TRUE, year = 2022) # 'cb = TRUE' for cartographic boundary (simplified) version
    counties = counties[,c("NAME","geometry")]
    colnames(counties) = c("County","geometry")
    
    ggplot()+
      geom_sf(data = plss_covariates_area_sf, aes(color = policy), size = 0.2) +
      geom_sf(data = plss_covariates_west_central_sf, aes(color = policy), size = 0.2) +
      geom_sf(data = plss_covariates_east_central_sf, aes(color = policy), size = 0.2) +
      geom_sf(data = plss_covariates_south_central_sf, aes(color = policy), size = 0.2) +
      geom_sf(data = plss_covariates_central_sf, aes(color = policy), size = 0.2) +
      geom_sf(data = counties,color = 'black') +
      scale_color_viridis_c(option = "plasma",
                            name = "Minimum Depth (Feet)",
                            labels = function(breaks) { round(exp(breaks), 0) })+
      labs(title = "5mg/L Threshold Well Depth", x = "Longitude", y = "Latitude") +
      theme_void()
    
    
    
    
    
    ggplot()+
      geom_sf(data = plss_covariates_area_sf, aes(color = U), size = 0.15) +
      geom_sf(data = counties,color = 'black')+
      scale_color_viridis_c(option = "plasma",
                            name = "Estimation of U",
                            labels = function(breaks) { round(exp(breaks), 0) })+
      labs(title = "5mg/L Threshold Well Depth", x = "Longitude", y = "Latitude") +
      theme_void()
    
    saveRDS(plss_covariates_area_sf, "plss_covariates_area_sf.rds")
    plss_covariates_area_sf$U <- plss_krige_values_area
    plss_covariates_area_sf <- readRDS("plss_covariates_area_sf.rds")
    
  }
}

