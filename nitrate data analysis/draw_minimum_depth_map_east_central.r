##### generate well depth at north east_central

###### outcome regression

data <- load_nitrate_data(zero_inflated = FALSE)
# data <- data[data$SampleYear>=2020, ]
data_east_central <- data[data$area=="East Central",]
data_east_central_split <- split_nitrate_data(data_east_central)


outcome_regression_east_central <- outcome_regression_SVM(data = data_east_central_split$data, data_test = data_east_central_split$data_test,tunning = FALSE)


cat('Training: ')
Metrics::rmse(data_east_central_split$data$logconcentration_plus_median, outcome_regression_east_central$pred + outcome_regression_east_central$krige_values)
Metrics::rmse(data_east_central_split$data$logconcentration_plus_median, outcome_regression_east_central$pred)
cat('\n')
cat('Test: ')
Metrics::rmse(data_east_central_split$data_test$logconcentration_plus_median, outcome_regression_east_central$pred_test + outcome_regression_east_central$krige_values_test)
Metrics::rmse(data_east_central_split$data_test$logconcentration_plus_median, outcome_regression_east_central$pred_test)



###### smoothers


depth_range_east_central <- c(min(data_east_central_split$data$logWellDepth), max(data_east_central_split$data$logWellDepth))



smoothers_east_central <- get_smoothers_RKHS(design_matrix = outcome_regression_east_central$design_matrix,
                                        design_matrix_test = outcome_regression_east_central$design_matrix_test,
                                        svm_auto = outcome_regression_east_central$svm,
                                        depth_range = depth_range_east_central)


# The above code is very slow, so we save the result


saveRDS(smoothers_east_central, "smoothers_east_central.rds")
smoothers_east_central <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_east_central.rds")


###### indirect method


threshold_val = log(5)


indirect_policy_RKHS_east_central <- compute_indirect_policy(smoothers = smoothers_east_central$smoothers_RKHS, trt_bounds = depth_range_east_central, threshold_val = threshold_val, krige_adjust = outcome_regression_east_central$krige_values, spatial = TRUE)

indirect_policy_test_RKHS_east_central <- compute_indirect_policy(smoothers = smoothers_east_central$smoothers_test_RKHS, trt_bounds = depth_range_east_central, threshold_val = threshold_val, krige_adjust = outcome_regression_east_central$krige_values_test, spatial = TRUE)



cat('Indirect Training data:\n')
print_acc_mcc(T = data_east_central_split$data$logWellDepth, 
              Y = data_east_central_split$data$logconcentration_plus_median,
              policy = indirect_policy_RKHS_east_central,
              threshold_val = log(5))



cat('Indirect Test data:\n')
print_acc_mcc(T = data_east_central_split$data_test$logWellDepth, 
              Y = data_east_central_split$data_test$logconcentration_plus_median,
              policy = indirect_policy_test_RKHS_east_central,
              threshold_val = log(5))



###### propensity

####### covariate balance

weights_cbps_east_central <- weightit(logWellDepth ~ StaticLevel+crop_type_combine+drainagecl+precipitation+cafolog, 
                                 data = data_east_central_split$data, method = "cbps")

summary(weights_cbps_east_central$weights)



###### initialization


initial_value_set <- c(depth_range_east_central[1] - 1, mean(depth_range_east_central), depth_range_east_central[2] + 1)
kernel_bw_east_central <- 1.06 * sd(data_east_central_split$data$logWellDepth) *nrow(data_east_central_split$data)^(-1/5)


direct_individual_multi_initial_east_central <- mcmapply(FUN = find_best_policy, 
                                                    i = 1:nrow(data_east_central_split$data), 
                                                    MoreArgs = list(initial_value_set = initial_value_set, 
                                                                    kernel_bw = kernel_bw_east_central,
                                                                    outcome_regression_object = outcome_regression_east_central, 
                                                                    data = data_east_central_split$data, 
                                                                    weights = weights_cbps_east_central$weights,
                                                                    smoothers = smoothers_east_central, 
                                                                    depth_range = depth_range_east_central, 
                                                                    threshold_val = log(5), 
                                                                    clip_epsilon = 20), SIMPLIFY = TRUE)

saveRDS(direct_individual_multi_initial_east_central, "direct_individual_multi_initial_east_central.rds")
direct_individual_multi_initial_east_central <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/direct_individual_multi_initial_east_central.rds")

print_acc_mcc(T = data_east_central_split$data$logWellDepth, 
              Y = data_east_central_split$data$logconcentration_plus_median,
              policy = direct_individual_multi_initial_east_central,
              threshold_val = log(5))



###### select kernel bandwidth


# there is no well use, and staticlevel 
RKHS_covariate_names = c("StaticLevel", "crop_type_combine","drainagecl","precipitation","cafolog")
rec <- recipe(~ ., data = data_east_central_split$data[, RKHS_covariate_names]) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep() # This "fits" the preprocessor



kernel_design_matrix <- bake(rec, new_data = data_east_central_split$data[, RKHS_covariate_names])
kernel_design_matrix <- cbind(kernel_design_matrix, U = outcome_regression_east_central$krige_values)

kernel_design_matrix_test <- bake(rec, new_data = data_east_central_split$data_test[, RKHS_covariate_names])
kernel_design_matrix_test <- cbind(kernel_design_matrix_test, U = outcome_regression_east_central$krige_values_test)

train_means <- colMeans(kernel_design_matrix)
train_sds <- apply(kernel_design_matrix, 2, sd)

kernel_design_matrix <- scale(kernel_design_matrix)
kernel_design_matrix_test <- scale(kernel_design_matrix_test, center = train_means, scale = train_sds)



for(m in c(-4,-3,-2,-1.5,-1)){
  cat("m is: ",m,'\n')
  gamma <- 2^(m) * median(rdist(kernel_design_matrix))
  rbf <- rbfdot(sigma = 1/(gamma^2))
  k_matrix <- kernelMatrix(rbf, kernel_design_matrix)
  k_matrix <- cbind(rep(1,nrow(k_matrix)), k_matrix)
  
  k_matrix_test <- kernelMatrix(rbf, kernel_design_matrix_test, kernel_design_matrix)
  k_matrix_test <- cbind(rep(1,nrow(k_matrix_test)), k_matrix_test)
  
  initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_east_central, alpha = 0, lambda = 0)
  params_initial_glmnet <- coef(initial_glmnet)
  params_initial_glmnet <- as.numeric(params_initial_glmnet)
  cat('training: ')
  print_acc_mcc(T = data_east_central_split$data$logWellDepth, 
                Y = data_east_central_split$data$logconcentration_plus_median, 
                policy = k_matrix %*% params_initial_glmnet, 
                threshold_val = log(5))
  
  cat('test: ')
  print_acc_mcc(T = data_east_central_split$data_test$logWellDepth, 
                Y = data_east_central_split$data_test$logconcentration_plus_median, 
                policy = k_matrix_test %*% params_initial_glmnet, 
                threshold_val = log(5))
  
}



###### non-DC algorithm

# We will choose m as -2 (no yet)


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

initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_east_central, alpha = 0, lambda = 0)
params_initial_glmnet <- coef(initial_glmnet)
params_initial_glmnet <- as.numeric(params_initial_glmnet)

kernel_optim_trt_DB_east_central <- optim(par = params_initial_glmnet, 
                                     fn = compute_total_loss_smooth_RKHS,
                                     gr = d_compute_total_loss_smooth_RKHS,
                                     K = k_matrix,
                                     T = data_east_central_split$data$logWellDepth,
                                     krige_adjust = outcome_regression_east_central$krige_values,
                                     outcome_resid = data_east_central_split$data$logconcentration_plus_median - outcome_regression_east_central$pred - outcome_regression_east_central$krige_values,
                                     # propensity_est = pmax(generalized_propensity_east_central$gps_est, 0.0001),
                                     propensity_est = 1/weights_cbps_east_central$weights,
                                     lambda = lambda, 
                                     smoothers = smoothers_east_central$smoothers_RKHS,
                                     cumint_smoothers = smoothers_east_central$cumint_smoothers_RKHS,
                                     trt_bounds = depth_range_east_central,
                                     threshold_val = log(5), 
                                     kernel_bw = kernel_bw_east_central, 
                                     clip_epsilon = clip_epsilon,
                                     surrogate_type = "Gaussian",
                                     loss_type = "db",
                                     method = "L-BFGS-B",
                                     control = list(maxit = 1000, trace = 3)) 




cat('Training: \n')
print_acc_mcc(T = data_east_central_split$data$logWellDepth, 
              Y = data_east_central_split$data$logconcentration_plus_median,
              policy = k_matrix %*% kernel_optim_trt_DB_east_central$par,
              threshold_val = log(5))


cat('Test: \n')
print_acc_mcc(T = data_east_central_split$data_test$logWellDepth, 
              Y = data_east_central_split$data_test$logconcentration_plus_median, 
              policy = k_matrix_test %*% kernel_optim_trt_DB_east_central$par, 
              threshold_val = log(5))



saveRDS(kernel_optim_trt_DB_east_central, 'kernel_optim_trt_DB_east_central.rds')

###### visualize depth map


file_path <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/plss_covariates.csv"
command <- paste("brctl download", shQuote(file_path))
system(command)
plss_covariates <- read.csv(file_path)

plss_covariates_east_central <- plss_covariates[plss_covariates$County %in% unique(data_east_central$County),]


plss_covariates_east_central$cafolog <- log(plss_covariates_east_central$cafo + 1)



dim(plss_covariates_east_central)
dim(na.omit(plss_covariates_east_central))


plss_covariates_east_central <- na.omit(plss_covariates_east_central)


plss_covariates_east_central <- plss_covariates_east_central[plss_covariates_east_central$crop_type_combine %in% unique(data_east_central$crop_type_combine),]

####### Interpolate static water level


gpfit_staticLevel_east_central <- GpGp::fit_model(data_east_central$StaticLevel,
                                             locs = data_east_central[,c("longitude", "latitude")],
                                             covfun_name = "matern_sphere")

plss_covariates_east_central$StaticLevel <-  GpGp::predictions(fit = gpfit_staticLevel_east_central,
                                                          locs_pred = plss_covariates_east_central[,c("longitude", "latitude")], 
                                                          X_pred = rep(1,nrow(plss_covariates_east_central)))


####### get kriging values for plss


plss_krige_values_east_central <- GpGp::predictions(fit = outcome_regression_east_central$gpfit, locs_pred = plss_covariates_east_central[,c("longitude", "latitude")], X_pred = rep(1,nrow(plss_covariates_east_central)))


####### build kernel matrices


RKHS_covariate_names = c("StaticLevel", "crop_type_combine","drainagecl","precipitation","cafolog")
rec <- recipe(~ ., data = data_east_central_split$data[, RKHS_covariate_names]) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep() # This "fits" the preprocessor



plss_kernel_design_matrix_east_central <- bake(rec, new_data = plss_covariates_east_central[,  RKHS_covariate_names])

plss_kernel_design_matrix_east_central <- cbind(plss_kernel_design_matrix_east_central, U = plss_krige_values_east_central)
plss_kernel_design_matrix_east_central <- scale(plss_kernel_design_matrix_east_central, center = train_means, scale = train_sds)


plss_k_matrix_east_central <- kernelMatrix(rbf, plss_kernel_design_matrix_east_central, kernel_design_matrix)
plss_k_matrix_east_central <- cbind(rep(1,nrow(plss_k_matrix_east_central)), plss_k_matrix_east_central)


####### generate map


plss_covariates_east_central$policy <- plss_k_matrix_east_central %*% kernel_optim_trt_DB_east_central$par
plss_covariates_east_central$U <-plss_krige_values_east_central

plss_covariates_east_central_sf <- st_as_sf(plss_covariates_east_central, coords = c('longitude','latitude'), crs=4326)

counties <- counties(state = "WI", cb = TRUE, year = 2022) # 'cb = TRUE' for cartographic boundary (simplified) version
counties = counties[,c("NAME","geometry")]
colnames(counties) = c("County","geometry")


ggplot()+
  geom_sf(data = plss_covariates_east_central_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_north_west_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_north_east_central_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_south_west_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_west_east_central_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_east_east_central_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_south_east_central_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_east_central_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = counties,color = 'black') +
  scale_color_viridis_c(option = "plasma",
                        name = "Minimum Depth (Feet)",
                        labels = function(breaks) { round(exp(breaks), 0) })+
  labs(title = "5mg/L Threshold Well Depth", x = "Longitude", y = "Latitude") +
  theme_void()



ggplot()+
  geom_sf(data = plss_covariates_east_central_sf, aes(color = U), size = 0.15) +
  geom_sf(data = counties,color = 'black')+
  scale_color_viridis_c(option = "plasma",
                        name = "Estimation of U",
                        labels = function(breaks) { round(exp(breaks), 0) })+
  labs(title = "5mg/L Threshold Well Depth", x = "Longitude", y = "Latitude") +
  theme_void()

plss_covariates_east_central_sf$U <- plss_krige_values_east_central
saveRDS(plss_covariates_east_central_sf, "plss_covariates_east_central_sf.rds")
plss_covariates_east_central_sf <- readRDS("plss_covariates_east_central_sf.rds")
