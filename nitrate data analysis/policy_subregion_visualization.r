# load('plss_covariates_central_sf.rds')
load('plss_covariates_east_central_sf.rds')
# load('plss_covariates_south_central_sf.rds')
load('plss_covariates_west_central_sf.rds')
# this file seems problematic
# load('plss_covariates_south_west_sf.rds')
## Reason: You used saveRDS to create the file, so you must use readRDS to open it.
## If you specifically want to use load() (so the variable simply appears in your environment without assignment), you must re-save the file using save()

plss_covariates_central_sf <- readRDS("plss_covariates_central_sf.rds")
plss_covariates_south_central_sf <- readRDS("plss_covariates_south_central_sf.rds")

plss_covariates_south_west_sf <- readRDS("plss_covariates_south_west_sf.rds")
plss_covariates_north_central_sf <- readRDS("plss_covariates_north_central_sf.rds")
plss_covariates_north_east_sf <- readRDS("plss_covariates_north_east_sf.rds")
plss_covariates_north_west_sf <- readRDS("plss_covariates_north_west_sf.rds")
plss_covariates_south_east_sf <- readRDS("plss_covariates_south_east_sf.rds")

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

plss_covariates_sf_combine <- rbind(plss_covariates_south_east_sf, plss_covariates_north_east_sf, plss_covariates_north_west_sf,
                                    plss_covariates_north_central_sf, plss_covariates_south_west_sf, plss_covariates_west_central_sf,
                                    plss_covariates_east_central_sf, plss_covariates_south_central_sf, plss_covariates_central_sf)


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

ggplot() +
  geom_sf(data = plss_covariates_sf_combine, aes(color = policy), size = 0.2) +
  geom_sf(data = counties,color = "black", size = 0.2)+
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        breaks = log(c(10, 100, 500, 1000, 2000)),
                        limits = my_limits,
                        labels = function(breaks) { round(exp(breaks), 0)}) +
  labs(title = "Safe, Cost Effective Policy for 5mg/L Threshold") +
  theme_void()


# Compare Policy Map with Covariates

# land use/crop type
ggplot() +
  geom_sf(data = plss_covariates_sf_combine[!(plss_covariates_sf_combine$crop_type_combine %in% c("Undefined","Vegtables","Cover Crop","Other Crop")),], aes(color=crop_type_combine ),size=0.2) +
  geom_sf(data = counties,color = 'black', size = 0.2) +  
  scale_color_manual(
    values = c(
      "Forest" = "darkgreen", 
      "Grass" = "lightgreen", 
      "Corn" = "red",  # Make "Corn" stand out with red
      "Developed" = "gray", 
      "Soybeans" = "blue",  # Make "Soybeans" stand out with blue
      "Other Crop" = "orange",  # Make "Other Crop" stand out with orange
      "Cover Crop" = "brown", 
      "Vegtables" = "purple"
    ),
    name = "Land Use"
  ) + 
  # scale_color_viridis_d(option = "plasma", name = "Crop Type") +  # Color scale
  # labs(title = "Log(Cafo)", x = "Longitude", y = "Latitude") +
  theme_void() +
  labs(title = "Land Use in Wisconsin", x = "Longitude", y = "Latitude") +
  guides(
    color = guide_legend(
      keywidth = unit(0.7, "cm"),  # Increase width of the legend squares
      keyheight = unit(0.7, "cm"),  # Increase height of the legend squares
      override.aes = list(shape = 15, size = 7)  # Set legend keys to squares (shape 15) and increase size
    )
  )

## CAFO



ggplot() +
  geom_sf(data = counties,color = 'black', size = 0.2) +  
  geom_sf(data = plss_covariates_sf_combine, aes(color= log(cafo+1)),size=0.2) +
  scale_color_viridis_c(option = "magma",name = "Ln(CAFO+1)") +  # Color scale
  labs(title = "CAFO") +
  theme_void()


## precipitation

ggplot() +
  geom_sf(data = counties,color = 'black', size = 0.2,fill = 'white') +  
  geom_sf(data = plss_covariates_sf_combine, aes(color= precipitation),size=0.2) +
  scale_color_gradient(low = 'white',  # Light blue
                       high = "skyblue", # Dark blue
                       name = "mm/day")+
  labs(title = "Precipitation", x = "Longitude", y = "Latitude") +
  theme_void()


## drainage level

plss_covariates_sf_combine$drainagecl = factor(plss_covariates_sf_combine$drainagecl,
                                               levels = c("Very poorly drained","Poorly drained","Somewhat poorly drained","Moderately well drained","Well drained","Somewhat excessively drained","Excessively drained"))


ggplot() +
  geom_sf(data = counties,color = 'black', size = 0.2,fill = 'lightgrey') +  
  geom_sf(data = plss_covariates_sf_combine, aes(color=drainagecl), size = 0.2) +
  scale_color_brewer(palette = "magma", name = "Soil Drainage Level",
                     labels = c(
                       "Very poorly drained" = "Very poor",  # Custom label
                       "Poorly drained" = "Poorly",                  # Custom label
                       "Somewhat poorly drained" = "Somewhat poor",              # Custom label
                       "Moderately well drained" = "Moderately well",         # Custom label
                       "Well drained" = "Well",
                       "Somewhat excessively drained" = "Somewhat excessive",
                       "Excessively drained" = "Excessive"
                     )) + 
  theme_void()+
  guides(
    color = guide_legend(
      keywidth = unit(0.7, "cm"),  # Increase width of the legend squares
      keyheight = unit(0.7, "cm"),  # Increase height of the legend squares
      override.aes = list(shape = 15, size = 7)  # Set legend keys to squares (shape 15) and increase size
    )
  )


## static water level

ggplot() +
  geom_sf(data = counties,color = 'black', size = 0.2,fill = 'white') +  
  geom_sf(data = plss_covariates_sf_combine, aes(color= StaticLevel),size=0.2) +
  scale_color_gradient(low = 'white',  # Light blue
                       high = "cadetblue4", # Dark blue
                       name = "Feet",
                       # breaks = log(c(10, 100, 200,400)),
                       # limits = c(log(1), log(500)),
                       # labels = function(breaks) { round(exp(breaks), 0)}
                       )+
  labs(title = "Static Water Level", x = "Longitude", y = "Latitude") +
  theme_void()

## current well depth



# Interpretation/visualization of the direct policy

## Regress Policy Over Covariates




plss_covariates_sf_combine$crop_type_combine <- factor(plss_covariates_sf_combine$crop_type_combine, levels = c(
  "Forest",
  "Undefined",
  "Corn",
  "Soybeans",
  "Grass",
  "Vegtables",
  "Developed"
))

policy_gam <- mgcv::gam(policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel + s(U) +  drainagecl + cafolog, data = plss_covariates_sf_combine)
# policy_gam <- mgcv::gam(policy ~ crop_type_combine * drainagecl + precipitation + StaticLevel   + s(U) + cafolog * drainagecl, data = plss_covariates_sf_combine)

# emm_results <- emmeans(policy_gam, ~ crop_type_combine + drainagecl)

## draw the land use plot

crop_coeffs <- tidy(policy_gam, parametric = TRUE) %>%
  filter(str_detect(term, "crop_type_combine") & !str_detect(term, "Undefined")) %>%
  mutate(
    # Calculate 95% Confidence Intervals
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,
    
    # Clean up term names: remove the prefix "crop_type_combine"
    term_clean = str_remove(term, "crop_type_combine")
  )

# 2. Create the Forest Plot
crop_plot <- ggplot(crop_coeffs, aes(x = estimate, y = term_clean)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "darkgreen") +
  
  # Add a vertical reference line at 0 (the value for "Forest")
  geom_vline(xintercept = 0, linetype = "dashed", color = "pink", size = 0.8) +
  
  labs(
    title = "Policy Relative to the 'Forest' base level",
    # subtitle = "Relative to the 'Forest' base level",
    x = "Well Depth (Log(Feet))",
    y = "Land Use"
  ) +
  theme(
    panel.grid.major = element_line(color = "gray10", size = 0.5, linetype = "dotted"),
    panel.grid.minor = element_line(color = "gray10", size = 0.25, linetype = "dotted"),
    
    # You are also setting the plot border/axis lines here:
    axis.line = element_line(color = "black", size = 0.5),
    # You might want to remove the border box drawn by theme_bw
    panel.border = element_blank()
  ) + theme_minimal()

print(crop_plot)
# 2. View the table (Predictions for every combination)
summary(emm_results)


## draw the drainage line
emm_drainage <- emmeans(policy_gam, ~ drainagecl)

# 2. Convert to data frame for ggplot2
emm_df <- as.data.frame(emm_drainage) 

# 3. Plot the ordered trend
drainage_plot <- ggplot(emm_df, aes(x = drainagecl, y = emmean)) +
  # Connect the means with a line to emphasize the trend 
  geom_line(aes(group = 1), color = "darkblue", size = 1) + 
  # Add points for the estimated means
  geom_point(size = 3, color = "darkblue") +
  # Add 95% confidence intervals
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "darkblue") +
  
  labs(
    title = "Predicted Policy Value Across Soil Drainage Classes (GAM)",
    subtitle = "Trend adjusted for smooth effects of continuous covariates",
    x = "Soil Drainage Class ",
    y = "Well Depth (log(Feet))"
  ) +
  theme_minimal() +
  # Rotate x-axis labels for readability
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(drainage_plot)


## plot effect of $U$

draw(policy_gam, select = "s(U)", rug = FALSE, ci_alpha = 0.2) + 
  geom_line(color = "#0479A8", linewidth = 1.0) +
  # coord_cartesian(xlim = c(15000, 66000)) +  # Set the visible range for the x-axis
  labs(
    title = "",
    x = "Kriging Value (Gaussian process)", # Set custom x-axis label
    y = "Depth Policy (Log Feet)"         # Set custom y-axis label
  )+
  # theme_minimal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
  )

## plot effect of precipitation
plot(ggpredict(policy_gam, terms = c("precipitation"))) + 
  labs(
    title = "Linear Effect of Precipitation",
    x = "precipitation (mm/day)",
    y = "Marginal Average of Policy (Log Feet)"
  ) +
  # scale_y_continuous(labels = exp)+
  theme_bw() # Change the whole theme

## plot effect of static water level
plot(ggpredict(policy_gam, terms = c("StaticLevel")))+ 
  labs(
    title = "Linear Effect of Static Level",
    x = "static level (feet)",
    y = "Marginal Average of Policy (Log Feet)",
  ) +
  theme_bw() # Change the whole theme

## plot effect of static CAFO

plot(ggpredict(policy_gam, terms = c("cafolog")))+ 
  labs(
    title = "Linear Effect of CAFO",
    x = "cafolog",
    y = "Marginal Average of Depth (Log Feet)",
  ) +
  theme_minimal() # Change the whole theme
### get indirect method

threshold_val = log(5)

data <- load_nitrate_data(zero_inflated = FALSE)

depth_range <- c(min(data$logWellDepth), max(data$logWellDepth))

## south west
data_south_west <- data[data$area=="South West",]
data_south_west_split <- split_nitrate_data(data_south_west)

smoothers_south_west <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_south_west.rds")
indirect_policy_RKHS_south_west <- compute_indirect_policy(smoothers = smoothers_south_west$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_south_west <- compute_indirect_policy(smoothers = smoothers_south_west$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_south_west <- get_direct_indirect_method_comparison(data_split = data_south_west_split, indirect_policy = indirect_policy_RKHS_south_west, 
                                      indirect_policy_test = indirect_policy_test_RKHS_south_west, plss_covariates = plss_covariates_south_west_sf)

## south central

data_south_central <- data[data$area=="South Central",]
data_south_central_split <- split_nitrate_data(data_south_central)

smoothers_south_central <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_south_central.rds")
indirect_policy_RKHS_south_central <- compute_indirect_policy(smoothers = smoothers_south_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_south_central <- compute_indirect_policy(smoothers = smoothers_south_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_south_central <- get_direct_indirect_method_comparison(data_split = data_south_central_split, indirect_policy = indirect_policy_RKHS_south_central, 
                                                                         indirect_policy_test = indirect_policy_test_RKHS_south_central, plss_covariates = plss_covariates_south_central_sf)

## south east

data_south_east <- data[data$area=="South East" & data$SampleYear>=2020,]
data_south_east_split <- split_nitrate_data(data_south_east)

smoothers_south_east <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_south_east.rds")
indirect_policy_RKHS_south_east <- compute_indirect_policy(smoothers = smoothers_south_east$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_south_east <- compute_indirect_policy(smoothers = smoothers_south_east$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_south_east <- get_direct_indirect_method_comparison(data_split = data_south_east_split, indirect_policy = indirect_policy_RKHS_south_east, 
                                                                      indirect_policy_test = indirect_policy_test_RKHS_south_east, plss_covariates = plss_covariates_south_east_sf)

## west central

data_west_central <- data[data$area=="West Central",]
data_west_central_split <- split_nitrate_data(data_west_central)

smoothers_west_central <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_west_central.rds")
indirect_policy_RKHS_west_central <- compute_indirect_policy(smoothers = smoothers_west_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_west_central <- compute_indirect_policy(smoothers = smoothers_west_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_west_central <- get_direct_indirect_method_comparison(data_split = data_west_central_split, indirect_policy = indirect_policy_RKHS_west_central, 
                                                                        indirect_policy_test = indirect_policy_test_RKHS_west_central, plss_covariates = plss_covariates_west_central_sf)

## central

data_central <- data[data$area=="Central",]
data_central_split <- split_nitrate_data(data_central)

smoothers_central <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_central.rds")
indirect_policy_RKHS_central <- compute_indirect_policy(smoothers = smoothers_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_central <- compute_indirect_policy(smoothers = smoothers_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_central <- get_direct_indirect_method_comparison(data_split = data_central_split, indirect_policy = indirect_policy_RKHS_central, 
                                                                   indirect_policy_test = indirect_policy_test_RKHS_central, plss_covariates = plss_covariates_central_sf)
## east central

data_east_central <- data[data$area=="East Central",]
data_east_central_split <- split_nitrate_data(data_east_central)

smoothers_east_central <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_east_central.rds")
indirect_policy_RKHS_east_central <- compute_indirect_policy(smoothers = smoothers_east_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_east_central <- compute_indirect_policy(smoothers = smoothers_east_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_east_central <- get_direct_indirect_method_comparison(data_split = data_east_central_split, indirect_policy = indirect_policy_RKHS_east_central, 
                                                                        indirect_policy_test = indirect_policy_test_RKHS_east_central, plss_covariates = plss_covariates_east_central_sf)

## north west

data_north_west <- data[data$area=="Northwest" & data$SampleYear>=2020,]
data_north_west_split <- split_nitrate_data(data_north_west)

smoothers_north_west <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_north_west.rds")
indirect_policy_RKHS_north_west <- compute_indirect_policy(smoothers = smoothers_north_west$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_north_west <- compute_indirect_policy(smoothers = smoothers_north_west$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_north_west <- get_direct_indirect_method_comparison(data_split = data_north_west_split, indirect_policy = indirect_policy_RKHS_north_west, 
                                                                      indirect_policy_test = indirect_policy_test_RKHS_north_west, plss_covariates = plss_covariates_north_west_sf)

## north central

data_north_central <- data[data$area=="North Central" & data$SampleYear>=2020,]
data_north_central_split <- split_nitrate_data(data_north_central)

smoothers_north_central <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_north_central.rds")
indirect_policy_RKHS_north_central <- compute_indirect_policy(smoothers = smoothers_north_central$smoothers_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)
indirect_policy_test_RKHS_north_central <- compute_indirect_policy(smoothers = smoothers_north_central$smoothers_test_RKHS, trt_bounds = depth_range, threshold_val = threshold_val, spatial = FALSE)

policy_comparison_north_central <- get_direct_indirect_method_comparison(data_split = data_north_central_split, indirect_policy = indirect_policy_RKHS_north_central, 
                                                                         indirect_policy_test = indirect_policy_test_RKHS_north_central, plss_covariates = plss_covariates_north_central_sf)

## north east

data_north_east <- data[data$area=="North East" & data$SampleYear>=2020,]
data_north_east_split <- split_nitrate_data(data_north_east)

smoothers_north_east <- readRDS("~/Desktop/policy learning/simulation/new_version_after_JSM/smoothers_north_east.rds")
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



# Nonsp Indirect method v.s. Our method

## analyze the indirect policy

# we relevel  the data

policy_comparison_combine_sf$drainagecl = factor(policy_comparison_combine_sf$drainagecl,
                                               levels = c("Very poorly drained","Poorly drained","Somewhat poorly drained","Moderately well drained","Well drained","Somewhat excessively drained","Excessively drained"))


policy_comparison_combine_sf$crop_type_combine <- factor(policy_comparison_combine_sf$crop_type_combine, levels = c(
  "Forest",
  "Undefined",
  "Corn",
  "Soybeans",
  "Grass",
  "Vegtables",
  "Developed"
))

policy_comparison_combine_sf$logWellDepth <- log(policy_comparison_combine_sf$WellDepth)

indirect_policy_gam <- mgcv::gam(indirect_policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel+drainagecl + cafolog, data = policy_comparison_combine_sf)
indirect_policy_gam <- mgcv::gam(policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel+drainagecl + cafolog, data = policy_comparison_combine_sf)


# 1. get the coefficients
crop_coeffs <- tidy(indirect_policy_gam, parametric = TRUE) %>%
  filter(str_detect(term, "crop_type_combine") & !str_detect(term, "Undefined")) %>%
  mutate(
    # Calculate 95% Confidence Intervals
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,
    
    # Clean up term names: remove the prefix "crop_type_combine"
    term_clean = str_remove(term, "crop_type_combine")
  )

#2. Create the Forest Plot
crop_plot <- ggplot(crop_coeffs[!(crop_coeffs$term_clean %in% c("Cover Crop")),], aes(x = estimate, y = term_clean)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "darkgreen") +
  
  # Add a vertical reference line at 0 (the value for "Forest")
  geom_vline(xintercept = 0, linetype = "dashed", color = "pink", size = 0.8) +
  
  labs(
    title = "Policy Relative to the 'Forest' base level",
    # subtitle = "Relative to the 'Forest' base level",
    x = "Well Depth (Log(Feet))",
    y = "Land Use"
  ) +
  theme(
    panel.grid.major = element_line(color = "gray10", size = 0.5, linetype = "dotted"),
    panel.grid.minor = element_line(color = "gray10", size = 0.25, linetype = "dotted"),
    
    # You are also setting the plot border/axis lines here:
    axis.line = element_line(color = "black", size = 0.5),
    # You might want to remove the border box drawn by theme_bw
    panel.border = element_blank()
  ) + theme_minimal()

print(crop_plot)
# 2. View the table (Predictions for every combination)
summary(emm_results)


## draw the drainage line
emm_drainage <- emmeans(indirect_policy_gam, ~ drainagecl)

# 2. Convert to data frame for ggplot2
emm_df <- as.data.frame(emm_drainage) 

# 3. Plot the ordered trend
drainage_plot <- ggplot(emm_df, aes(x = drainagecl, y = emmean)) +
  # Connect the means with a line to emphasize the trend 
  geom_line(aes(group = 1), color = "darkblue", size = 1) + 
  # Add points for the estimated means
  geom_point(size = 3, color = "darkblue") +
  # Add 95% confidence intervals
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "darkblue") +
  
  labs(
    title = "Predicted Policy Value Across Soil Drainage Classes (GAM)",
    subtitle = "Trend adjusted for smooth effects of continuous covariates",
    x = "Soil Drainage Class ",
    y = "Well Depth (log(Feet))"
  ) +
  theme_minimal() +
  # Rotate x-axis labels for readability
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(drainage_plot)

## compare the acc/mcc


calculate_acc_mcc_two_sided_f1(T = log(policy_comparison_combine_sf$WellDepth), 
                               Y = log(policy_comparison_combine_sf$concentration_plus_median), 
                               policy = policy_comparison_combine_sf$indirect_policy, 
                               threshold_val = log(5))

calculate_acc_mcc_two_sided_f1(T = log(policy_comparison_combine_sf_test$WellDepth), 
                               Y = log(policy_comparison_combine_sf_test$concentration_plus_median), 
                               policy = policy_comparison_combine_sf_test$indirect_policy, 
                               threshold_val = log(5))

calculate_acc_mcc_two_sided_f1(T = log(policy_comparison_combine_sf$WellDepth), 
                               Y = log(policy_comparison_combine_sf$concentration_plus_median), 
                               policy = policy_comparison_combine_sf$policy, 
                               threshold_val = log(5))

calculate_acc_mcc_two_sided_f1(T = log(policy_comparison_combine_sf_test$WellDepth), 
                  Y = log(policy_comparison_combine_sf_test$concentration_plus_median), 
                  policy = policy_comparison_combine_sf_test$policy, 
                  threshold_val = log(5))




## density plot

density_data_long <- st_drop_geometry(policy_comparison_combine_sf[policy_comparison_combine_sf$concentration_plus_median>=5,]) %>%
  dplyr::select(indirect_policy, policy, logWellDepth) %>%
  pivot_longer(
    cols = everything(), # Selects all columns 
    names_to = "Variable",
    values_to = "Value"
  )
  
# 3. Create the density plot
combined_density_plot <- density_data_long %>%
  ggplot(aes(x = Value, fill = Variable)) +
  geom_density(alpha = 0.5) + # alpha controls transparency for overlapping areas
  scale_x_continuous(
    name = "Depth (Feet)", # Change the X-axis title here
    breaks = log(c(20, 100, 500, 2000)) ,
    labels = function(x) {
      # Exponentiate the value and format to a clean number
      scales::number(exp(x), accuracy = 1) 
    }
  ) +
  labs(
    title = "Density Plot of  Well Depth Policy v.s. Observed Well Depth in Central Area",
    x = "Depth (Feet)",
    y = "Density",
    fill = "Variable"
  ) +
  theme_bw() 
# +
#   scale_fill_manual(values = c("logWellDepth" = "darkblue", "policy" = "darkgreen"),
#                     labels = c("logWellDepth" = "Well Depth", "policy" = "Policy"))

# Display the plot
print(combined_density_plot)


## visualization
ggplot() +
  geom_sf(data = policy_comparison_combine_sf_test, aes(color = pmin(indirect_policy,log(2000))), size = 0.5) +
  geom_sf(data = counties,color = "black", size = 0.5)+
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        breaks = log(c(10, 100, 500, 1000, 2000)),
                        limits = my_limits,
                        labels = function(breaks) { round(exp(breaks), 0)}) +
  labs(title = "Indirect Well Depth Policy") +
  theme_void()



ggplot() +
  geom_sf(data = policy_comparison_combine_sf_test, aes(color = pmin(policy,log(2000))), size = 0.5) +
  geom_sf(data = counties,color = "black", size = 0.5)+
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        breaks = log(c(10, 100, 500, 1000, 2000)),
                        limits = my_limits,
                        labels = function(breaks) { round(exp(breaks), 0)}) +
  labs(title = "Direct Well Depth Policy") +
  theme_void()

# 

