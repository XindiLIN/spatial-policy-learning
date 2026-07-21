source('miscellaneous.r')


## Draw i.i.d. 

### Load Data

data <- load_nitrate_data(zero_inflated = FALSE)
data <- data[data$SampleYear>=2020, ]
# depth_range <- c(min(data$logWellDepth), max(data$logWellDepth))
depth_range <- c(min(data$logWellDepth), quantile(data$logWellDepth,0.999))

### Train SVM 

OR_covariate_names = c("logWellDepth","crop_type_combine","drainagecl","precipitation","cafolog","StaticLevel")

rec <- recipe(~ ., data = data[, OR_covariate_names]) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep() # This "fits" the preprocessor

design_matrix <- bake(rec, new_data = data[, OR_covariate_names])

svm_auto <- svm(x = design_matrix, y = data$logconcentration_plus_median, type='nu-regression')

### Make Prediction on PLSS

#### load PLSS

file_path <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/plss_covariates_static_no_na.csv"
command <- paste("brctl download", shQuote(file_path))
system(command)
plss_covariates <- read.csv(file_path)

#### Build Smoothers

plss_covariates$logWellDepth <- mean(depth_range)
plss_covariates_design_matrix <- bake(rec, new_data = plss_covariates[, OR_covariate_names])

smoothers_svm_whole <- get_smoothers_RKHS(design_matrix = plss_covariates_design_matrix,
                                          design_matrix_test = NULL,
                                          svm_auto = svm_auto,
                                          depth_range = depth_range,
                                          treatment_step = 0.1)

saveRDS(smoothers_svm_whole, "~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_svm_whole.rds")
smoothers_svm_whole <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_svm_whole.rds")
### Find Policy under i.i.d. Modeling

indirect_policy_RKHS_iid_whole <- compute_indirect_policy(smoothers = smoothers_svm_whole$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, krige_adjust = NULL, spatial = FALSE)

plss_covariates$indirect_policy_RKHS_iid_whole <- indirect_policy_RKHS_iid_whole

### Plot

plss_covariates_sf <- st_as_sf(plss_covariates, coords = c('longitude','latitude'), crs=4326)

counties <- counties(state = "WI", cb = TRUE, year = 2022) # 'cb = TRUE' for cartographic boundary (simplified) version
counties = counties[,c("NAME","geometry")]
colnames(counties) = c("County","geometry")


ggplot() +
  geom_sf(data = plss_covariates_sf, aes(color = indirect_policy_RKHS_iid_whole), size = 0.1) +
  geom_sf(data = counties,color = 'black') +
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        limits = c(log(5),log(2000)),
                        breaks = log(c(10, 100, 500, 1000, 2000)),
                        labels = function(breaks) { round(exp(breaks), 0)}) +
  labs(title = "5mg/L Threshold Well Depth: Non Spatial", x = "Longitude", y = "Latitude") +
  theme_void()





#### Indirect method at training and test data

area_lst <- c("North Central","West Central","Northwest","Central","South West", "South Central","South East","East Central","North East")
threshold_val = log(5)

data <- load_nitrate_data(zero_inflated = FALSE)
depth_range <- c(min(data$logWellDepth), max(data$logWellDepth))
# depth_range = c() # we would use a consistent depth_range from the original dataset

### North West
data <- load_nitrate_data(zero_inflated = FALSE)
data <- data[data$SampleYear>=2020, ]
data_north_west <- data[data$area=="Northwest",]
data_north_west_split <- split_nitrate_data(data_north_west)


### North Central
data <- load_nitrate_data(zero_inflated = FALSE)
data <- data[data$SampleYear>=2020, ]
data_north_central <- data[data$area=="North Central",]
data_north_central_split <- split_nitrate_data(data_north_central)

outcome_regression_north_central <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/outcome_regression_north_central.rds")
smoothers_north_central <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_north_central.rds")

indirect_policy_RKHS_north_central <- compute_indirect_policy(smoothers = smoothers_north_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, krige_adjust = outcome_regression_north_central$krige_values, spatial = TRUE)
indirect_policy_test_RKHS_north_central <- compute_indirect_policy(smoothers = smoothers_north_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, krige_adjust = outcome_regression_north_central$krige_values_test, spatial = TRUE)

data_north_central_split$data$indirect_policy_RKHS <- indirect_policy_RKHS_north_central
data_north_central_split$data_test$indirect_policy_test_RKHS_area <- indirect_policy_test_RKHS_north_central

### South West

data <- load_nitrate_data(zero_inflated = FALSE)
data_south_west <- data[data$area=="South West",]
data_south_west_split <- split_nitrate_data(data_south_west)

outcome_regression_south_west <- outcome_regression_SVM(data = data_south_west_split$data, data_test = data_south_west_split$data_test,tunning = FALSE)
depth_range_south_west <- c(min(data_south_west_split$data$logWellDepth), max(data_south_west_split$data$logWellDepth))
smoothers_south_west <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_south_west.rds")

indirect_policy_RKHS_south_west <- compute_indirect_policy(smoothers = smoothers_south_west$smoothers_RKHS, trt_bounds = depth_range_south_west, threshold_val = threshold_val, krige_adjust = outcome_regression_south_west$krige_values, spatial = TRUE)

indirect_policy_test_RKHS_south_west <- compute_indirect_policy(smoothers = smoothers_south_west$smoothers_test_RKHS, trt_bounds = depth_range_south_west, threshold_val = threshold_val, krige_adjust = outcome_regression_south_west$krige_values_test, spatial = TRUE)

data_south_west_split$data$indirect_policy_RKHS <- indirect_policy_RKHS_south_west
data_south_west_split$data_test$indirect_policy_test_RKHS_area <- indirect_policy_test_RKHS_south_west


file_path <- "/Users/xindilin/Desktop/2024 summer/groundwater_pesticide/data/plss_covariates.csv"
command <- paste("brctl download", shQuote(file_path))
system(command)
plss_covariates <- read.csv(file_path)

plss_covariates_south_west <- plss_covariates[plss_covariates$County %in% unique(data_south_west$County),]
plss_covariates_south_west <- na.omit(plss_covariates_south_west)




for(area in area_lst){
  data_area <- data[data$area==area,]
  data_area_split <- split_nitrate_data(data_area)
  
  
  outcome_regression_area <- outcome_regression_SVM(data = data_area_split$data, data_test = data_area_split$data_test,tunning = FALSE)
  
  depth_range_area <- c(min(data_area_split$data$logWellDepth), max(data_area_split$data$logWellDepth))
  smoothers_area
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
  
  
  
}
