##### generate well depth at south west

###### outcome regression

data <- load_nitrate_data(zero_inflated = FALSE)
data_south_west <- data[data$area=="South West",]
data_south_west_split <- split_nitrate_data(data_south_west)



outcome_regression_south_west <- outcome_regression_SVM(data = data_south_west_split$data, data_test = data_south_west_split$data_test,tunning = FALSE)



cat('Training: ')
Metrics::rmse(data_south_west_split$data$logconcentration_plus_median, outcome_regression_south_west$pred + outcome_regression_south_west$krige_values)
Metrics::rmse(data_south_west_split$data$logconcentration_plus_median, outcome_regression_south_west$pred)
cat('\n')
cat('Test: ')
Metrics::rmse(data_south_west_split$data_test$logconcentration_plus_median, outcome_regression_south_west$pred_test + outcome_regression_south_west$krige_values_test)
Metrics::rmse(data_south_west_split$data_test$logconcentration_plus_median, outcome_regression_south_west$pred_test)


# Training: 
# 0.6429579
# 0.7089384
# Test:
# 0.6662577
# 0.7552599


###### smoothers


depth_range_south_west <- c(min(data_south_west_split$data$logWellDepth), max(data_south_west_split$data$logWellDepth))



smoothers_south_west <- get_smoothers_RKHS(design_matrix = outcome_regression_south_west$design_matrix,
                                              design_matrix_test = outcome_regression_south_west$design_matrix_test,
                                              svm_auto = outcome_regression_south_west$svm,
                                              depth_range = depth_range_south_west)


# The above code is very slow, so we save the result


saveRDS(smoothers_south_west, "smoothers_south_west.rds")
smoothers_south_west <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_south_west.rds")


###### indirect method


threshold_val = log(5)


indirect_policy_RKHS_south_west <- compute_indirect_policy(smoothers = smoothers_south_west$smoothers_RKHS, trt_bounds = depth_range_south_west, threshold_val = threshold_val, krige_adjust = outcome_regression_south_west$krige_values, spatial = TRUE)

indirect_policy_test_RKHS_south_west <- compute_indirect_policy(smoothers = smoothers_south_west$smoothers_test_RKHS, trt_bounds = depth_range_south_west, threshold_val = threshold_val, krige_adjust = outcome_regression_south_west$krige_values_test, spatial = TRUE)



cat('Indirect Training data:\n')
print_acc_mcc(T = data_south_west_split$data$logWellDepth, 
              Y = data_south_west_split$data$logconcentration_plus_median,
              policy = indirect_policy_RKHS_south_west,
              threshold_val = log(5))

# The accuracy is  0.8314351 
# The mcc is  0.4428726 

# Indirect Training data:
# The accuracy is  0.2314351 
# The mcc is  0.4428726 


cat('Indirect Test data:\n')
print_acc_mcc(T = data_south_west_split$data_test$logWellDepth, 
              Y = data_south_west_split$data_test$logconcentration_plus_median,
              policy = indirect_policy_test_RKHS_south_west,
              threshold_val = log(5))

# The accuracy is  0.8430016 
# The mcc is  0.4979457 

# Indirect Test data:
# The accuracy is  0.2430016 
# The mcc is  0.4979457 

###### propensity

####### covariate balance

weights_cbps_south_west <- weightit(logWellDepth ~ StaticLevel+crop_type_combine+drainagecl+precipitation+cafolog, 
                                       data = data_south_west_split$data, method = "cbps")

summary(weights_cbps_south_west$weights)


###### initialization


initial_value_set <- c(depth_range_south_west[1] - 1, mean(depth_range_south_west), depth_range_south_west[2] + 1)
kernel_bw_south_west <- 1.06 * sd(data_south_west_split$data$logWellDepth) *nrow(data_south_west_split$data)^(-1/5)


direct_individual_multi_initial_south_west <- mcmapply(FUN = find_best_policy, 
                                                          i = 1:nrow(data_south_west_split$data), 
                                                          MoreArgs = list(initial_value_set = initial_value_set, 
                                                                          kernel_bw = kernel_bw_south_west,
                                                                          outcome_regression_object = outcome_regression_south_west, 
                                                                          data = data_south_west_split$data, 
                                                                          weights = weights_cbps_south_west$weights,
                                                                          smoothers = smoothers_south_west, 
                                                                          depth_range = depth_range_south, 
                                                                          threshold_val = log(5), 
                                                                          clip_epsilon = 20), SIMPLIFY = TRUE)
saveRDS(direct_individual_multi_initial_south_west, "direct_individual_multi_initial_south_west.rds")
direct_individual_multi_initial_south_west <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/direct_individual_multi_initial_south_west.rds")

print_acc_mcc(T = data_south_west_split$data$logWellDepth, 
              Y = data_south_west_split$data$logconcentration_plus_median,
              policy = direct_individual_multi_initial_south_west,
              threshold_val = log(5))

# The accuracy is 0.9207289 
# The mcc is 0.7605969 

###### select kernel bandwidth



# there is no well use, and staticlevel 
RKHS_covariate_names = c("StaticLevel", "crop_type_combine","drainagecl","precipitation","cafolog")
rec <- recipe(~ ., data = data_south_west_split$data[, RKHS_covariate_names]) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep() # This "fits" the preprocessor



kernel_design_matrix <- bake(rec, new_data = data_south_west_split$data[, RKHS_covariate_names])
kernel_design_matrix <- cbind(kernel_design_matrix, U = outcome_regression_south_west$krige_values)

kernel_design_matrix_test <- bake(rec, new_data = data_south_west_split$data_test[, RKHS_covariate_names])
kernel_design_matrix_test <- cbind(kernel_design_matrix_test, U = outcome_regression_south_west$krige_values_test)

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
  
  initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_south_west, alpha = 0, lambda = 0)
  params_initial_glmnet <- coef(initial_glmnet)
  params_initial_glmnet <- as.numeric(params_initial_glmnet)
  cat('training: ')
  print_acc_mcc(T = data_south_west_split$data$logWellDepth, 
                Y = data_south_west_split$data$logconcentration_plus_median, 
                policy = k_matrix %*% params_initial_glmnet, 
                threshold_val = log(5))
  
  cat('test: ')
  print_acc_mcc(T = data_south_west_split$data_test$logWellDepth, 
                Y = data_south_west_split$data_test$logconcentration_plus_median, 
                policy = k_matrix_test %*% params_initial_glmnet, 
                threshold_val = log(5))
  
}

# m is:  -4 
# training: The accuracy is  0.9075171 
# The mcc is  0.7182911 
# test: The accuracy is  0.241405 
# The mcc is  0.4983094 
# m is:  -3 
# training: The accuracy is  0.2929385 
# The mcc is  0.6691536 
# test: The accuracy is  0.2451304 
# The mcc is  0.5120353 
# m is:  -2 
# training: The accuracy is  0.2603645 
# The mcc is  0.5539811 
# test: The accuracy is  0.2430016 
# The mcc is  0.5020877 
# m is:  -1.5 
# training: The accuracy is  0.2453303 
# The mcc is  0.4965617 
# test: The accuracy is  0.238744 
# The mcc is  0.4847051 



###### non-DC algorithm

# We will choose m is:  -3 


# cat("lambda is: ",lambda,'\n')
m <- -3
lambda <- 0.25
clip_epsilon <- 30
gamma <- 2^(m) * median(fields::rdist(kernel_design_matrix))
rbf <- rbfdot(sigma = 1/(gamma^2))
k_matrix <- kernelMatrix(rbf, kernel_design_matrix)
k_matrix <- cbind(rep(1,nrow(k_matrix)), k_matrix)

k_matrix_test <- kernelMatrix(rbf, kernel_design_matrix_test, kernel_design_matrix)
k_matrix_test <- cbind(rep(1,nrow(k_matrix_test)), k_matrix_test)

initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_south_west, alpha = 0, lambda = 0)
params_initial_glmnet <- coef(initial_glmnet)
params_initial_glmnet <- as.numeric(params_initial_glmnet)

kernel_optim_trt_DB_south_west <- optim(par = params_initial_glmnet, 
                                           fn = compute_total_loss_smooth_RKHS,
                                           gr = d_compute_total_loss_smooth_RKHS,
                                           K = k_matrix,
                                           T = data_south_west_split$data$logWellDepth,
                                           krige_adjust = outcome_regression_south_west$krige_values,
                                           outcome_resid = data_south_west_split$data$logconcentration_plus_median - outcome_regression_south_west$pred - outcome_regression_south_west$krige_values,
                                           # propensity_est = pmax(generalized_propensity_south_west$gps_est, 0.0001),
                                           propensity_est = 1/weights_cbps_south_west$weights,
                                           lambda = lambda, 
                                           smoothers = smoothers_south_west$smoothers_RKHS,
                                           cumint_smoothers = smoothers_south_west$cumint_smoothers_RKHS,
                                           trt_bounds = depth_range,
                                           threshold_val = log(5), 
                                           kernel_bw = kernel_bw_south_west, 
                                           clip_epsilon = clip_epsilon,
                                           surrogate_type = "Gaussian",
                                           loss_type = "db",
                                           method = "L-BFGS-B",
                                           control = list(maxit = 1000, trace = 3)) 




cat('Training: \n')
print_acc_mcc(T = data_south_west_split$data$logWellDepth, 
              Y = data_south_west_split$data$logconcentration_plus_median,
              policy = k_matrix %*% kernel_optim_trt_DB_south_west$par,
              threshold_val = log(5))
# 11/18
# The accuracy is  0.8986333 
# The mcc is  0.6881898 

cat('Test: \n')
print_acc_mcc(T = data_south_west_split$data_test$logWellDepth, 
              Y = data_south_west_split$data_test$logconcentration_plus_median, 
              policy = k_matrix_test %*% kernel_optim_trt_DB_south_west$par, 
              threshold_val = log(5))

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

saveRDS(kernel_optim_trt_DB_south_west, 'kernel_optim_trt_DB_south_west.rds')

###### visualize depth map


file_path <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/plss_covariates.csv"
command <- paste("brctl download", shQuote(file_path))
system(command)
plss_covariates <- read.csv(file_path)



plss_covariates_south_west <- plss_covariates[plss_covariates$County %in% c("Vernon", "Crawford", "Grant", "Richland", 
                                                                            "Sauk", "Iowa", "Lafayette"),]



plss_covariates_south_west$cafolog <- log(plss_covariates_south_west$cafo + 1)



dim(plss_covariates_south_west)
dim(na.omit(plss_covariates_south_west))

# 5491   16
# 3728   16

plss_covariates_south_west <- na.omit(plss_covariates_south_west)


# There are two plss section w with 'Vegtables' as their crop_type_combine, which is present in training data, we will just remove them


plss_covariates_south_west <- plss_covariates_south_west[!plss_covariates_south_west$crop_type_combine %in% c("Vegtables"),]



####### Interpolate static water level


gpfit_staticLevel_south_west <- GpGp::fit_model(data_south_west$StaticLevel,
                                                   locs = data_south_west[,c("longitude", "latitude")],
                                                   covfun_name = "matern_sphere")

plss_covariates_south_west$StaticLevel <-  GpGp::predictions(fit = gpfit_staticLevel_south_west,
                                                                locs_pred = plss_covariates_south_west[,c("longitude", "latitude")], 
                                                                X_pred = rep(1,nrow(plss_covariates_south_west)))


####### get kriging values for plss


plss_krige_values_south_west <- GpGp::predictions(fit = outcome_regression_south_west$gpfit, locs_pred = plss_covariates_south_west[,c("longitude", "latitude")], X_pred = rep(1,nrow(plss_covariates_south_west)))


####### build kernel matrices


RKHS_covariate_names = c("StaticLevel", "crop_type_combine","drainagecl","precipitation","cafolog")
rec <- recipe(~ ., data = data_south_west_split$data[, RKHS_covariate_names]) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep() # This "fits" the preprocessor



plss_kernel_design_matrix_south_west <- bake(rec, new_data = plss_covariates_south_west[,  RKHS_covariate_names])

plss_kernel_design_matrix_south_west <- cbind(plss_kernel_design_matrix_south_west, U = plss_krige_values_south_west)
plss_kernel_design_matrix_south_west <- scale(plss_kernel_design_matrix_south_west, center = train_means, scale = train_sds)



plss_k_matrix_south_west <- kernelMatrix(rbf, plss_kernel_design_matrix_south_west, kernel_design_matrix)
plss_k_matrix_south_west <- cbind(rep(1,nrow(plss_k_matrix_south_west)), plss_k_matrix_south_west)


####### generate map


plss_covariates_south_west$policy <- plss_k_matrix_south_west %*% kernel_optim_trt_DB_south_west$par



plss_covariates_south_west_sf <- st_as_sf(plss_covariates_south_west, coords = c('longitude','latitude'), crs=4326)

counties <- counties(state = "WI", cb = TRUE, year = 2022) # 'cb = TRUE' for cartographic boundary (simplified) version
counties = counties[,c("NAME","geometry")]
colnames(counties) = c("County","geometry")

ggplot()+
  geom_sf(data = plss_covariates_south_west_sf, aes(color = policy), size = 0.2) +
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
  geom_sf(data = plss_covariates_south_west_sf, aes(color = U), size = 0.15) +
  geom_sf(data = counties,color = 'black')+
  scale_color_viridis_c(option = "plasma",
                        name = "Estimation of U",
                        labels = function(breaks) { round(exp(breaks), 0) })+
  labs(title = "5mg/L Threshold Well Depth", x = "Longitude", y = "Latitude") +
  theme_void()

saveRDS(plss_covariates_south_west_sf, "plss_covariates_south_west_sf.rds")
plss_covariates_south_west_sf$U <- plss_krige_values_south_west
plss_covariates_south_west_sf <- readRDS("plss_covariates_south_west_sf.rds")
