### get indirect method

library(sf)
source("functions/miscellaneous.r")

threshold_val = log(5)

data <- load_nitrate_data(file_path = "data/data_Nitrate_with_covar.csv", zero_inflated = FALSE)

## load per-region PLSS covariates (east_central/west_central were saved with
## save() rather than saveRDS(), so they need load() instead of readRDS())
plss_covariates_south_west_sf    <- readRDS("nitrate_data_analysis/output/plss_covariates_south_west_sf.rds")
plss_covariates_south_central_sf <- readRDS("nitrate_data_analysis/output/plss_covariates_south_central_sf.rds")
plss_covariates_south_east_sf    <- readRDS("nitrate_data_analysis/output/plss_covariates_south_east_sf.rds")
load("nitrate_data_analysis/output/plss_covariates_west_central_sf.rds")
load("nitrate_data_analysis/output/plss_covariates_east_central_sf.rds")
plss_covariates_central_sf       <- readRDS("nitrate_data_analysis/output/plss_covariates_central_sf.rds")
plss_covariates_north_west_sf    <- readRDS("nitrate_data_analysis/output/plss_covariates_north_west_sf.rds")
plss_covariates_north_central_sf <- readRDS("nitrate_data_analysis/output/plss_covariates_north_central_sf.rds")
plss_covariates_north_east_sf    <- readRDS("nitrate_data_analysis/output/plss_covariates_north_east_sf.rds")

depth_range <- c(min(data$logWellDepth), max(data$logWellDepth))

## south west
data_south_west <- data[data$area=="South West",]
data_south_west_split <- split_nitrate_data(data_south_west)

smoothers_south_west <- readRDS("nitrate_data_analysis/output/smoothers_south_west.rds")
indirect_policy_RKHS_south_west <- compute_indirect_policy(smoothers = smoothers_south_west$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_south_west <- compute_indirect_policy(smoothers = smoothers_south_west$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_south_west <- get_direct_indirect_method_comparison(data_split = data_south_west_split, indirect_policy = indirect_policy_RKHS_south_west, 
                                                                      indirect_policy_test = indirect_policy_test_RKHS_south_west, plss_covariates = plss_covariates_south_west_sf)

## south central

data_south_central <- data[data$area=="South Central",]
data_south_central_split <- split_nitrate_data(data_south_central)

smoothers_south_central <- readRDS("nitrate_data_analysis/output/smoothers_south_central.rds")
indirect_policy_RKHS_south_central <- compute_indirect_policy(smoothers = smoothers_south_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_south_central <- compute_indirect_policy(smoothers = smoothers_south_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_south_central <- get_direct_indirect_method_comparison(data_split = data_south_central_split, indirect_policy = indirect_policy_RKHS_south_central, 
                                                                         indirect_policy_test = indirect_policy_test_RKHS_south_central, plss_covariates = plss_covariates_south_central_sf)

## south east

data_south_east <- data[data$area=="South East" & data$SampleYear>=2020,]
data_south_east_split <- split_nitrate_data(data_south_east)

smoothers_south_east <- readRDS("nitrate_data_analysis/output/smoothers_south_east.rds")
indirect_policy_RKHS_south_east <- compute_indirect_policy(smoothers = smoothers_south_east$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_south_east <- compute_indirect_policy(smoothers = smoothers_south_east$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_south_east <- get_direct_indirect_method_comparison(data_split = data_south_east_split, indirect_policy = indirect_policy_RKHS_south_east, 
                                                                      indirect_policy_test = indirect_policy_test_RKHS_south_east, plss_covariates = plss_covariates_south_east_sf)

## west central

data_west_central <- data[data$area=="West Central",]
data_west_central_split <- split_nitrate_data(data_west_central)

smoothers_west_central <- readRDS("nitrate_data_analysis/output/smoothers_west_central.rds")
indirect_policy_RKHS_west_central <- compute_indirect_policy(smoothers = smoothers_west_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_west_central <- compute_indirect_policy(smoothers = smoothers_west_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_west_central <- get_direct_indirect_method_comparison(data_split = data_west_central_split, indirect_policy = indirect_policy_RKHS_west_central, 
                                                                        indirect_policy_test = indirect_policy_test_RKHS_west_central, plss_covariates = plss_covariates_west_central_sf)

## central

data_central <- data[data$area=="Central" & data$SampleYear>=2020,]
data_central_split <- split_nitrate_data(data_central)

smoothers_central <- readRDS("nitrate_data_analysis/output/smoothers_central.rds")
indirect_policy_RKHS_central <- compute_indirect_policy(smoothers = smoothers_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_central <- compute_indirect_policy(smoothers = smoothers_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_central <- get_direct_indirect_method_comparison(data_split = data_central_split, indirect_policy = indirect_policy_RKHS_central, 
                                                                   indirect_policy_test = indirect_policy_test_RKHS_central, plss_covariates = plss_covariates_central_sf)
## east central

data_east_central <- data[data$area=="East Central",]
data_east_central_split <- split_nitrate_data(data_east_central)

smoothers_east_central <- readRDS("nitrate_data_analysis/output/smoothers_east_central.rds")
indirect_policy_RKHS_east_central <- compute_indirect_policy(smoothers = smoothers_east_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_east_central <- compute_indirect_policy(smoothers = smoothers_east_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_east_central <- get_direct_indirect_method_comparison(data_split = data_east_central_split, indirect_policy = indirect_policy_RKHS_east_central, 
                                                                        indirect_policy_test = indirect_policy_test_RKHS_east_central, plss_covariates = plss_covariates_east_central_sf)

## north west

data_north_west <- data[data$area=="Northwest" & data$SampleYear>=2020,]
data_north_west_split <- split_nitrate_data(data_north_west)

smoothers_north_west <- readRDS("nitrate_data_analysis/output/smoothers_north_west.rds")
indirect_policy_RKHS_north_west <- compute_indirect_policy(smoothers = smoothers_north_west$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_north_west <- compute_indirect_policy(smoothers = smoothers_north_west$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_north_west <- get_direct_indirect_method_comparison(data_split = data_north_west_split, indirect_policy = indirect_policy_RKHS_north_west, 
                                                                      indirect_policy_test = indirect_policy_test_RKHS_north_west, plss_covariates = plss_covariates_north_west_sf)

## north central

data_north_central <- data[data$area=="North Central" & data$SampleYear>=2020,]
data_north_central_split <- split_nitrate_data(data_north_central)

smoothers_north_central <- readRDS("nitrate_data_analysis/output/smoothers_north_central.rds")
indirect_policy_RKHS_north_central <- compute_indirect_policy(smoothers = smoothers_north_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_north_central <- compute_indirect_policy(smoothers = smoothers_north_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_north_central <- get_direct_indirect_method_comparison(data_split = data_north_central_split, indirect_policy = indirect_policy_RKHS_north_central, 
                                                                         indirect_policy_test = indirect_policy_test_RKHS_north_central, plss_covariates = plss_covariates_north_central_sf)

## north east

data_north_east <- data[data$area=="North East" & data$SampleYear>=2020,]
data_north_east_split <- split_nitrate_data(data_north_east)

smoothers_north_east <- readRDS("nitrate_data_analysis/output/smoothers_north_east.rds")
indirect_policy_RKHS_north_east <- compute_indirect_policy(smoothers = smoothers_north_east$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_north_east <- compute_indirect_policy(smoothers = smoothers_north_east$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_north_east <- get_direct_indirect_method_comparison(data_split = data_north_east_split, indirect_policy = indirect_policy_RKHS_north_east, 
                                                                      indirect_policy_test = indirect_policy_test_RKHS_north_east, plss_covariates = plss_covariates_north_east_sf)




# combine the policy in different regions
policy_comparison_combine <- rbind(policy_comparison_south_east, policy_comparison_north_east, policy_comparison_north_west,
                                   policy_comparison_north_central, policy_comparison_south_west, policy_comparison_west_central,
                                   policy_comparison_east_central, policy_comparison_south_central, policy_comparison_central)

policy_comparison_combine_sf <- st_as_sf(policy_comparison_combine,coords = c("longitude", "latitude"),crs = 4326)

policy_comparison_combine_sf_test <- policy_comparison_combine_sf[policy_comparison_combine_sf$test==1,]

global_min <- min(c(policy_comparison_combine_sf$indirect_policy, 
                    policy_comparison_combine_sf$policy), na.rm = TRUE)
global_max <- log(2000)
my_limits <- c(global_min, global_max)


