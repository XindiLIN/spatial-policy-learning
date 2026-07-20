##### generate well depth at north central

###### outcome regression

data <- load_nitrate_data(zero_inflated = FALSE)
data <- data[data$SampleYear>=2020, ]
data_south_east <- data[data$area=="South East",]
data_south_east_split <- split_nitrate_data(data_south_east)


outcome_regression_south_east <- outcome_regression_SVM(data = data_south_east_split$data, data_test = data_south_east_split$data_test,tunning = FALSE)


cat('Training: ')
Metrics::rmse(data_south_east_split$data$logconcentration_plus_median, outcome_regression_south_east$pred + outcome_regression_south_east$krige_values)
Metrics::rmse(data_south_east_split$data$logconcentration_plus_median, outcome_regression_south_east$pred)
cat('\n')
cat('Test: ')
Metrics::rmse(data_south_east_split$data_test$logconcentration_plus_median, outcome_regression_south_east$pred_test + outcome_regression_south_east$krige_values_test)
Metrics::rmse(data_south_east_split$data_test$logconcentration_plus_median, outcome_regression_south_east$pred_test)


# Training: 
# 0.6620767
# 0.6911073
# Test:
# 0.6821321
# 0.7264399


###### smoothers


depth_range_south_east <- c(min(data_south_east_split$data$logWellDepth), max(data_south_east_split$data$logWellDepth))



smoothers_south_east <- get_smoothers_RKHS(design_matrix = outcome_regression_south_east$design_matrix,
                                           design_matrix_test = outcome_regression_south_east$design_matrix_test,
                                           svm_auto = outcome_regression_south_east$svm,
                                           depth_range = depth_range_south_east)


# The above code is very slow, so we save the result


saveRDS(smoothers_south_east, "smoothers_south_east.rds")
smoothers_south_east <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_south_east.rds")


###### indirect method


threshold_val = log(5)


indirect_policy_RKHS_south_east <- compute_indirect_policy(smoothers = smoothers_south_east$smoothers_RKHS, trt_bounds = depth_range_south_east, threshold_val = threshold_val, krige_adjust = outcome_regression_south_east$krige_values, spatial = TRUE)

indirect_policy_test_RKHS_south_east <- compute_indirect_policy(smoothers = smoothers_south_east$smoothers_test_RKHS, trt_bounds = depth_range_south_east, threshold_val = threshold_val, krige_adjust = outcome_regression_south_east$krige_values_test, spatial = TRUE)



cat('Indirect Training data:\n')
print_acc_mcc(T = data_south_east_split$data$logWellDepth, 
              Y = data_south_east_split$data$logconcentration_plus_median,
              policy = indirect_policy_RKHS_south_east,
              threshold_val = threshold_val)

# The accuracy is  0.8704762 
# The mcc is  0.3837981 


cat('Indirect Test data:\n')
print_acc_mcc(T = data_south_east_split$data_test$logWellDepth, 
              Y = data_south_east_split$data_test$logconcentration_plus_median,
              policy = indirect_policy_test_RKHS_south_east,
              threshold_val = threshold_val)

# The accuracy is  0.8707665 
# The mcc is  0.3654044

###### propensity

####### covariate balance

weights_cbps_south_east <- weightit(logWellDepth ~ StaticLevel+crop_type_combine+drainagecl+precipitation+cafolog, 
                                    data = data_south_east_split$data, method = "cbps")

summary(weights_cbps_south_east$weights)

# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 3.471e-04 7.110e-01 8.113e-01 8.712e-01 9.088e-01 1.444e+01 

###### initialization

initial_value_set <- c(depth_range_south_east[1] - 1, mean(depth_range_south_east), depth_range_south_east[2] + 1)
kernel_bw_south_east <- 1.06 * sd(data_south_east_split$data$logWellDepth) *nrow(data_south_east_split$data)^(-1/5)


direct_individual_multi_initial_south_east <- mcmapply(FUN = find_best_policy, 
                                                       i = 1:nrow(data_south_east_split$data), 
                                                       MoreArgs = list(initial_value_set = initial_value_set, 
                                                                       kernel_bw = kernel_bw_south_east,
                                                                       outcome_regression_object = outcome_regression_south_east, 
                                                                       data = data_south_east_split$data, 
                                                                       weights = weights_cbps_south_east$weights,
                                                                       smoothers = smoothers_south_east, 
                                                                       depth_range = depth_range_south_east, 
                                                                       threshold_val = threshold_val, 
                                                                       clip_epsilon = 20), SIMPLIFY = TRUE)

saveRDS(direct_individual_multi_initial_south_east, "direct_individual_multi_initial_south_east.rds")
direct_individual_multi_initial_south_east <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/direct_individual_multi_initial_south_east.rds")

print_acc_mcc(T = data_south_east_split$data$logWellDepth, 
              Y = data_south_east_split$data$logconcentration_plus_median,
              policy = direct_individual_multi_initial_south_east,
              threshold_val = threshold_val)

# TThe accuracy is  0.9615238 
# The mcc is  0.8488292 

###### select kernel bandwidth


# there is no well use, and staticlevel 
RKHS_covariate_names = c("StaticLevel", "crop_type_combine","drainagecl","precipitation","cafolog")
rec <- recipe(~ ., data = data_south_east_split$data[, RKHS_covariate_names]) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep() # This "fits" the preprocessor



kernel_design_matrix <- bake(rec, new_data = data_south_east_split$data[, RKHS_covariate_names])
kernel_design_matrix <- cbind(kernel_design_matrix, U = outcome_regression_south_east$krige_values)

kernel_design_matrix_test <- bake(rec, new_data = data_south_east_split$data_test[, RKHS_covariate_names])
kernel_design_matrix_test <- cbind(kernel_design_matrix_test, U = outcome_regression_south_east$krige_values_test)

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
  
  initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_south_east, alpha = 0, lambda = 0)
  params_initial_glmnet <- coef(initial_glmnet)
  params_initial_glmnet <- as.numeric(params_initial_glmnet)
  cat('training: ')
  print_acc_mcc(T = data_south_east_split$data$logWellDepth, 
                Y = data_south_east_split$data$logconcentration_plus_median, 
                policy = k_matrix %*% params_initial_glmnet, 
                threshold_val = threshold_val)
  
  cat('test: ')
  print_acc_mcc(T = data_south_east_split$data_test$logWellDepth, 
                Y = data_south_east_split$data_test$logconcentration_plus_median, 
                policy = k_matrix_test %*% params_initial_glmnet, 
                threshold_val = threshold_val)
  
}

# m is:  -4 
# training: The accuracy is  0.9561905 
# The mcc is  0.826536 
# test: The accuracy is  0.8618538 
# The mcc is  0.3602784 
# m is:  -3 
# training: The accuracy is  0.9344762 
# The mcc is  0.7329126 
# test: The accuracy is  0.8743316 
# The mcc is  0.412891 
# m is:  -2 
# training: The accuracy is  0.904381 
# The mcc is  0.5899621 
# test: The accuracy is  0.868984 
# The mcc is  0.3923773 
# m is:  -1.5 
# training: The accuracy is  0.887619 
# The mcc is  0.5028002 
# test: The accuracy is  0.8752228 
# The mcc is  0.4241568 
# m is:  -1 
# training: The accuracy is  0.8838095 
# The mcc is  0.4766715 
# test: The accuracy is  0.8752228 
# The mcc is  0.4135167 



###### non-DC algorithm

# We will choose m as -1.5


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

initial_glmnet <- glmnet(x = k_matrix[,-1], y = direct_individual_multi_initial_south_east, alpha = 0, lambda = 0)
params_initial_glmnet <- coef(initial_glmnet)
params_initial_glmnet <- as.numeric(params_initial_glmnet)

kernel_optim_trt_DB_south_east <- optim(par = params_initial_glmnet, 
                                        fn = compute_total_loss_smooth_RKHS,
                                        gr = d_compute_total_loss_smooth_RKHS,
                                        K = k_matrix,
                                        T = data_south_east_split$data$logWellDepth,
                                        krige_adjust = outcome_regression_south_east$krige_values,
                                        outcome_resid = data_south_east_split$data$logconcentration_plus_median - outcome_regression_south_east$pred - outcome_regression_south_east$krige_values,
                                        # propensity_est = pmax(generalized_propensity_south_east$gps_est, 0.0001),
                                        propensity_est = 1/weights_cbps_south_east$weights,
                                        lambda = lambda, 
                                        smoothers = smoothers_south_east$smoothers_RKHS,
                                        cumint_smoothers = smoothers_south_east$cumint_smoothers_RKHS,
                                        trt_bounds = depth_range_south_east,
                                        threshold_val = threshold_val, 
                                        kernel_bw = kernel_bw_south_east, 
                                        clip_epsilon = clip_epsilon,
                                        surrogate_type = "Gaussian",
                                        loss_type = "db",
                                        method = "L-BFGS-B",
                                        control = list(maxit = 1000, trace = 3)) 




cat('Training: \n')
print_acc_mcc(T = data_south_east_split$data$logWellDepth, 
              Y = data_south_east_split$data$logconcentration_plus_median,
              policy = k_matrix %*% kernel_optim_trt_DB_south_east$par,
              threshold_val = threshold_val)
# 11/18
# The accuracy is  0.8887619 
# The mcc is  0.5121747 

cat('Test: \n')
print_acc_mcc(T = data_south_east_split$data_test$logWellDepth, 
              Y = data_south_east_split$data_test$logconcentration_plus_median, 
              policy = k_matrix_test %*% kernel_optim_trt_DB_south_east$par, 
              threshold_val = threshold_val)

# The accuracy is  0.8698752 
# The mcc is  0.4003615 


saveRDS(kernel_optim_trt_DB_south_east, 'kernel_optim_trt_DB_south_east.rds')

###### visualize depth map


file_path <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/plss_covariates.csv"
command <- paste("brctl download", shQuote(file_path))
system(command)
plss_covariates <- read.csv(file_path)



plss_covariates_south_east <- plss_covariates[plss_covariates$County %in% unique(data_south_east$County),]


plss_covariates_south_east$cafolog <- log(plss_covariates_south_east$cafo + 1)



dim(plss_covariates_south_east)
dim(na.omit(plss_covariates_south_east))

# 2672   16
# 1743   16

plss_covariates_south_east <- na.omit(plss_covariates_south_east)

# There are two plss section w with 'Vegtables' as their crop_type_combine, which is present in training data, we will just remove them

plss_covariates_south_east <- plss_covariates_south_east[plss_covariates_south_east$crop_type_combine %in% unique(data_south_east$crop_type_combine),]

####### Interpolate static water level


gpfit_staticLevel_south_east <- GpGp::fit_model(data_south_east$StaticLevel,
                                                locs = data_south_east[,c("longitude", "latitude")],
                                                covfun_name = "matern_sphere")

plss_covariates_south_east$StaticLevel <-  GpGp::predictions(fit = gpfit_staticLevel_south_east,
                                                             locs_pred = plss_covariates_south_east[,c("longitude", "latitude")], 
                                                             X_pred = rep(1,nrow(plss_covariates_south_east)))


####### get kriging values for plss


plss_krige_values_south_east <- GpGp::predictions(fit = outcome_regression_south_east$gpfit, locs_pred = plss_covariates_south_east[,c("longitude", "latitude")], X_pred = rep(1,nrow(plss_covariates_south_east)))


####### build kernel matrices


RKHS_covariate_names = c("StaticLevel", "crop_type_combine","drainagecl","precipitation","cafolog")
rec <- recipe(~ ., data = data_south_east_split$data[, RKHS_covariate_names]) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep() # This "fits" the preprocessor



plss_kernel_design_matrix_south_east <- bake(rec, new_data = plss_covariates_south_east[,  RKHS_covariate_names])

plss_kernel_design_matrix_south_east <- cbind(plss_kernel_design_matrix_south_east, U = plss_krige_values_south_east)
plss_kernel_design_matrix_south_east <- scale(plss_kernel_design_matrix_south_east, center = train_means, scale = train_sds)


plss_k_matrix_south_east <- kernelMatrix(rbf, plss_kernel_design_matrix_south_east, kernel_design_matrix)
plss_k_matrix_south_east <- cbind(rep(1,nrow(plss_k_matrix_south_east)), plss_k_matrix_south_east)


####### generate map


plss_covariates_south_east$policy <- plss_k_matrix_south_east %*% kernel_optim_trt_DB_south_east$par
plss_covariates_south_east$U <-plss_krige_values_south_east

plss_covariates_south_east_sf <- st_as_sf(plss_covariates_south_east, coords = c('longitude','latitude'), crs=4326)

counties <- counties(state = "WI", cb = TRUE, year = 2022) # 'cb = TRUE' for cartographic boundary (simplified) version
counties = counties[,c("NAME","geometry")]
colnames(counties) = c("County","geometry")


ggplot()+
  geom_sf(data = plss_covariates_south_east_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_north_east_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_north_west_sf, aes(color = policy), size = 0.2) +
  geom_sf(data = plss_covariates_north_central_sf, aes(color = policy), size = 0.2) +
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
  geom_sf(data = plss_covariates_south_east_sf, aes(color = U), size = 0.15) +
  geom_sf(data = counties,color = 'black')+
  scale_color_viridis_c(option = "plasma",
                        name = "Estimation of U",
                        labels = function(breaks) { round(exp(breaks), 0) })+
  labs(title = "5mg/L Threshold Well Depth", x = "Longitude", y = "Latitude") +
  theme_void()


saveRDS(plss_covariates_south_east_sf, "plss_covariates_south_east_sf.rds")
plss_covariates_south_east_sf <- readRDS("plss_covariates_south_east_sf.rds")
