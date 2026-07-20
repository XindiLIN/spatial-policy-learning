##### generate well depth at south center

###### outcome regression

data <- load_nitrate_data(zero_inflated = FALSE)
# data <- data[data$SampleYear>=2020, ]
data_south_central <- data[data$area=="South Central",]
data_south_central_split <- split_nitrate_data(data_south_central)



outcome_regression_south_central <- outcome_regression_SVM(data = data_south_central_split$data, data_test = data_south_central_split$data_test,tunning = FALSE)



cat('Training: ')
Metrics::rmse(data_south_central_split$data$logconcentration_plus_median, outcome_regression_south_central$pred + outcome_regression_south_central$krige_values)
Metrics::rmse(data_south_central_split$data$logconcentration_plus_median, outcome_regression_south_central$pred)
cat('\n')
cat('Test: ')
Metrics::rmse(data_south_central_split$data_test$logconcentration_plus_median, outcome_regression_south_central$pred_test + outcome_regression_south_central$krige_values_test)
Metrics::rmse(data_south_central_split$data_test$logconcentration_plus_median, outcome_regression_south_central$pred_test)

###### smoothers


depth_range_south_central <- c(min(data_south_central_split$data$logWellDepth), max(data_south_central_split$data$logWellDepth))



smoothers_south_central <- get_smoothers_RKHS(design_matrix = outcome_regression_south_central$design_matrix,
                                             design_matrix_test = outcome_regression_south_central$design_matrix_test,
                                             svm_auto = outcome_regression_south_central$svm,
                                             depth_range = depth_range_south_central)


# The above code is very slow, so we save the result


saveRDS(smoothers_south_central, "smoothers_south_central.rds")


###### indirect method


threshold_val = log(5)



indirect_policy_RKHS_south_central <- compute_indirect_policy(smoothers = smoothers_south_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, krige_adjust = outcome_regression_south_central$krige_values, spatial = TRUE)

indirect_policy_test_RKHS_south_central <- compute_indirect_policy(smoothers = smoothers_south_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, krige_adjust = outcome_regression_south_central$krige_values_test, spatial = TRUE)



cat('Indirect Training data:\n')
print_acc_mcc(T = data_south_central_split$data$logWellDepth, 
              Y = data_south_central_split$data$logconcentration_plus_median,
              policy = indirect_policy_RKHS_south_central,
              threshold_val = log(5))
# Indirect Training data:
# The accuracy is  0.7604464 
# The mcc is  0.5281033 


cat('Indirect Test data:\n')
print_acc_mcc(T = data_south_central_split$data_test$logWellDepth, 
              Y = data_south_central_split$data_test$logconcentration_plus_median,
              policy = indirect_policy_test_RKHS_south_central,
              threshold_val = log(5))

# The accuracy is  0.7530982 
# The mcc is  0.5118123 

###### propensity

####### covariate balance

weights_cbps_south_central <- weightit(logWellDepth ~ StaticLevel+crop_type_combine+drainagecl+precipitation+cafolog, 
                                      data = data_south_central_split$data, method = "cbps")

summary(weights_cbps_south_central$weights)


###### initialization


initial_value_set <- c(depth_range_south[1] - 1, mean(depth_range_south), depth_range_south[2] + 1)
kernel_bw_south_central <- 1.06 * sd(data_south_central_split$data$logWellDepth) *nrow(data_south_central_split$data)^(-1/5)


direct_individual_multi_initial_south_central <- mcmapply(FUN = find_best_policy, 
                                                         i = 1:nrow(data_south_central_split$data), 
                                                         MoreArgs = list(initial_value_set = initial_value_set, 
                                                                         kernel_bw = kernel_bw_south_central,
                                                                         outcome_regression_object = outcome_regression_south_central, 
                                                                         data = data_south_central_split$data, 
                                                                         weights = weights_cbps_south_central$weights,
                                                                         smoothers = smoothers_south_central, 
                                                                         depth_range = depth_range_south, 
                                                                         threshold_val = log(5), 
                                                                         clip_epsilon = 20), SIMPLIFY = TRUE)
saveRDS(direct_individual_multi_initial_south_central, "direct_individual_multi_initial_south_central.rds")


print_acc_mcc(T = data_south_central_split$data$logWellDepth, 
              Y = data_south_central_split$data$logconcentration_plus_median,
              policy = direct_individual_multi_initial_south_central,
              threshold_val = log(5))

# The accuracy is  0.9021369 
# The mcc is  0.8076363 

###### select kernel bandwidth



# there is no well use, and staticlevel 
RKHS_covariate_names = c("StaticLevel", "crop_type_combine","drainagecl","precipitation","cafolog")
rec <- recipe(~ ., data = data_south_central_split$data[, RKHS_covariate_names]) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep() # This "fits" the preprocessor



kernel_design_matrix <- bake(rec, new_data = data_south_central_split$data[, RKHS_covariate_names])
kernel_design_matrix <- cbind(kernel_design_matrix, U = outcome_regression_south_central$krige_values)

kernel_design_matrix_test <- bake(rec, new_data = data_south_central_split$data_test[, RKHS_covariate_names])
kernel_design_matrix_test <- cbind(kernel_design_matrix_test, U = outcome_regression_south_central$krige_values_test)

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
  
  initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_south_central, alpha = 0, lambda = 0)
  params_initial_glmnet <- coef(initial_glmnet)
  params_initial_glmnet <- as.numeric(params_initial_glmnet)
  cat('training: ')
  print_acc_mcc(T = data_south_central_split$data$logWellDepth, 
                Y = data_south_central_split$data$logconcentration_plus_median, 
                policy = k_matrix %*% params_initial_glmnet, 
                threshold_val = log(5))
  
  cat('test: ')
  print_acc_mcc(T = data_south_central_split$data_test$logWellDepth, 
                Y = data_south_central_split$data_test$logconcentration_plus_median, 
                policy = k_matrix_test %*% params_initial_glmnet, 
                threshold_val = log(5))
  
}

# m is:  -4 
# training: The accuracy is  0.8770927 
# The mcc is  0.7584321 
# test: The accuracy is  0.7219574 
# The mcc is  0.4451872 
# m is:  -3 
# training: The accuracy is  0.8335375 
# The mcc is  0.6729542 
# test: The accuracy is  0.7368923 
# The mcc is  0.4787644 
# m is:  -2 
# training: The accuracy is  0.7803185 
# The mcc is  0.5695388 
# test: The accuracy is  0.7356212 
# The mcc is  0.4775718 
# m is:  -1.5 
# training: The accuracy is  0.7535048 
# The mcc is  0.5195058 
# test: The accuracy is  0.7378456 
# The mcc is  0.483221 



###### non-DC algorithm


# cat("lambda is: ",lambda,'\n')
m <- -1.5
lambda <- 0.25
clip_epsilon <- 30
gamma <- 2^(m) * median(fields::rdist(kernel_design_matrix))
rbf <- rbfdot(sigma = 1/(gamma^2))
k_matrix <- kernelMatrix(rbf, kernel_design_matrix)
k_matrix <- cbind(rep(1,nrow(k_matrix)), k_matrix)

k_matrix_test <- kernelMatrix(rbf, kernel_design_matrix_test, kernel_design_matrix)
k_matrix_test <- cbind(rep(1,nrow(k_matrix_test)), k_matrix_test)

initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_south_central, alpha = 0, lambda = 0)
params_initial_glmnet <- coef(initial_glmnet)
params_initial_glmnet <- as.numeric(params_initial_glmnet)

kernel_optim_trt_DB_south_central <- optim(par = params_initial_glmnet, 
                                          fn = compute_total_loss_smooth_RKHS,
                                          gr = d_compute_total_loss_smooth_RKHS,
                                          K = k_matrix,
                                          T = data_south_central_split$data$logWellDepth,
                                          krige_adjust = outcome_regression_south_central$krige_values,
                                          outcome_resid = data_south_central_split$data$logconcentration_plus_median - outcome_regression_south_central$pred - outcome_regression_south_central$krige_values,
                                          # propensity_est = pmax(generalized_propensity_south_central$gps_est, 0.0001),
                                          propensity_est = 1/weights_cbps_south_central$weights,
                                          lambda = lambda, 
                                          smoothers = smoothers_south_central$smoothers_RKHS,
                                          cumint_smoothers = smoothers_south_central$cumint_smoothers_RKHS,
                                          trt_bounds = depth_range,
                                          threshold_val = log(5), 
                                          kernel_bw = kernel_bw_south_central, 
                                          clip_epsilon = clip_epsilon,
                                          surrogate_type = "Gaussian",
                                          loss_type = "db",
                                          method = "L-BFGS-B",
                                          control = list(maxit = 1000, trace = 3)) 




cat('Training: \n')
print_acc_mcc(T = data_south_central_split$data$logWellDepth, 
              Y = data_south_central_split$data$logconcentration_plus_median,
              policy = k_matrix %*% kernel_optim_trt_DB_south_central$par,
              threshold_val = log(5))
cat('Test: \n')
print_acc_mcc(T = data_south_central_split$data_test$logWellDepth, 
              Y = data_south_central_split$data_test$logconcentration_plus_median, 
              policy = k_matrix_test %*% kernel_optim_trt_DB_south_central$par, 
              threshold_val = log(5))

# Training: 
# The accuracy is  0.7718797 
# The mcc is  0.5568953
# 
# Test:
# The accuracy is  0.7464252 
# The mcc is  0.502661

###### visualize depth map


file_path <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/plss_covariates.csv"
command <- paste("brctl download", shQuote(file_path))
system(command)
plss_covariates <- read.csv(file_path)



plss_covariates_south_central <- plss_covariates[plss_covariates$County %in% c("Columbia", "Dodge", "Dane", "Jefferson", 
                                                                               "Green", "Rock"),]



plss_covariates_south_central$cafolog <- log(plss_covariates_south_central$cafo + 1)



dim(plss_covariates_south_central)
dim(na.omit(plss_covariates_south_central))



plss_covariates_south_central <- na.omit(plss_covariates_south_central)



####### Interpolate static water level


gpfit_staticLevel_south_central <- GpGp::fit_model(data_south_central$StaticLevel,
                                                  locs = data_south_central[,c("longitude", "latitude")],
                                                  covfun_name = "matern_sphere")

plss_covariates_south_central$StaticLevel <-  GpGp::predictions(fit = gpfit_staticLevel_south_central,
                                                               locs_pred = plss_covariates_south_central[,c("longitude", "latitude")], 
                                                               X_pred = rep(1,nrow(plss_covariates_south_central)))


####### get kriging values for plss


plss_krige_values_south_central <- GpGp::predictions(fit = outcome_regression_south_central$gpfit, locs_pred = plss_covariates_south_central[,c("longitude", "latitude")], X_pred = rep(1,nrow(plss_covariates_south_central)))


####### build kernel matrices


RKHS_covariate_names = c("StaticLevel", "crop_type_combine","drainagecl","precipitation","cafolog")
rec <- recipe(~ ., data = data_south_central_split$data[, RKHS_covariate_names]) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep() # This "fits" the preprocessor



plss_kernel_design_matrix_south_central <- bake(rec, new_data = plss_covariates_south_central[,  RKHS_covariate_names])

plss_kernel_design_matrix_south_central <- cbind(plss_kernel_design_matrix_south_central, U = plss_krige_values_south_central)
plss_kernel_design_matrix_south_central <- scale(plss_kernel_design_matrix_south_central, center = train_means, scale = train_sds)



plss_k_matrix_south_central <- kernelMatrix(rbf, plss_kernel_design_matrix_south_central, kernel_design_matrix)
plss_k_matrix_south_central <- cbind(rep(1,nrow(plss_k_matrix_south_central)), plss_k_matrix_south_central)


####### generate map


plss_covariates_south_central$policy <- plss_k_matrix_south_central %*% kernel_optim_trt_DB_south_central$par



plss_covariates_south_central_sf <- st_as_sf(plss_covariates_south_central, coords = c('longitude','latitude'), crs=4326)



ggplot()+
  geom_sf(data = plss_covariates_west_central_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_east_central_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_south_central_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_central_sf, aes(color = policy), size = 0.1) +
  geom_sf(data = counties,color = 'black') +
  scale_color_viridis_c(option = "plasma",
                        name = "Minimum Depth (Feet)",
                        labels = function(breaks) { round(exp(breaks), 0) })+
  labs(title = "5mg/L Threshold Well Depth", x = "Longitude", y = "Latitude") +
  theme_void()






ggplot()+
  geom_sf(data = plss_covariates_south_central_sf, aes(color = U), size = 0.15) +
  geom_sf(data = counties,color = 'black')+
  scale_color_viridis_c(option = "plasma",
                        name = "Estimation of U",
                        labels = function(breaks) { round(exp(breaks), 0) })+
  labs(title = "5mg/L Threshold Well Depth", x = "Longitude", y = "Latitude") +
  theme_void()


plss_covariates_south_central_sf$U = plss_krige_values_south_central

saveRDS(plss_covariates_south_central_sf, "plss_covariates_south_central_sf.rds")
plss_covariates_south_central_sf <- readRDS("plss_covariates_south_central_sf.rds")



