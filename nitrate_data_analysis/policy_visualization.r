library(tigris)
library(ggplot2)
library(tidyr)
library(tidymodels)
library(dplyr)
library(patchwork)
library(stringr)
library(emmeans)

source("functions/miscellaneous.r")


# Load PLSS covariates and the estimated policy at each cell 

wi_counties_path <- "nitrate_data_analysis/output/wi_counties.rds"
plss_covariates_combined_path <- "nitrate_data_analysis/output/plss_covariates_sf_combine.rds"


if (file.exists(plss_covariates_combined_path)) {
  plss_covariates_sf_combine <- readRDS(plss_covariates_combined_path)
} else {
  plss_covariates_central_sf       <- readRDS("nitrate_data_analysis/output/plss_covariates_central_sf.rds")
  load('nitrate_data_analysis/output/plss_covariates_east_central_sf.rds')
  load('nitrate_data_analysis/output/plss_covariates_west_central_sf.rds')
  plss_covariates_south_central_sf <- readRDS("nitrate_data_analysis/output/plss_covariates_south_central_sf.rds")
  plss_covariates_south_west_sf    <- readRDS("nitrate_data_analysis/output/plss_covariates_south_west_sf.rds")
  plss_covariates_north_central_sf <- readRDS("nitrate_data_analysis/output/plss_covariates_north_central_sf.rds")
  plss_covariates_north_east_sf    <- readRDS("nitrate_data_analysis/output/plss_covariates_north_east_sf.rds")
  plss_covariates_north_west_sf    <- readRDS("nitrate_data_analysis/output/plss_covariates_north_west_sf.rds")
  plss_covariates_south_east_sf    <- readRDS("nitrate_data_analysis/output/plss_covariates_south_east_sf.rds")
  
  plss_covariates_sf_combine <- rbind(plss_covariates_south_east_sf, plss_covariates_north_east_sf, plss_covariates_north_west_sf,
                                      plss_covariates_north_central_sf, plss_covariates_south_west_sf, plss_covariates_west_central_sf,
                                      plss_covariates_east_central_sf, plss_covariates_south_central_sf, plss_covariates_central_sf)
  
  saveRDS(plss_covariates_sf_combine, plss_covariates_combined_path)
}

if (file.exists(wi_counties_path)) {
  counties <- readRDS(wi_counties_path)
} else {
  counties <- tigris::counties(state = "WI", cb = TRUE, year = 2022)
  counties <- counties[, c("NAME", "geometry")]
  colnames(counties) <- c("County", "geometry")
  
  dir.create(dirname(wi_counties_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(counties, wi_counties_path)
}


# Policy visualization

## Estimated Required Depth

map_plot <- ggplot() +
  geom_sf(data = plss_covariates_sf_combine, aes(color = policy), size = 0.5) +
  geom_sf(data = counties,color = "black", size = 0.2)+
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        breaks = log(c(10, 100, 500, 1000, 2000)),
                        limits = log(c(5,2000)),
                        labels = function(breaks) { round(exp(breaks), 0)}) +
  theme_void() +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm")
  )

print(map_plot)



## Linear projection of estimated Required Depth over covariates (drainage level, land use)

plss_covariates_sf_combine$crop_type_combine <- factor(plss_covariates_sf_combine$crop_type_combine, levels = c(
  "Forest",
  "Undefined",
  "Corn",
  "Soybeans",
  "Grass",
  "Vegtables",
  "Developed"
))

plss_covariates_sf_combine$drainagecl <- factor(plss_covariates_sf_combine$drainagecl,
                                                 levels = c("Very poorly drained","Poorly drained","Somewhat poorly drained","Moderately well drained","Well drained","Somewhat excessively drained","Excessively drained"))

policy_gam <- mgcv::gam(policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel + s(U) +  drainagecl + cafolog, data = plss_covariates_sf_combine)

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
    x = "Well Depth (Log(Feet))",
    y = "Land Use"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

print(crop_plot)

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
    x = "Soil Drainage Class",
    y = "Well Depth (log(Feet))"
  ) +
  theme_minimal() +
  # Rotate x-axis labels for readability
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1)
  )

print(drainage_plot)

## combine the drainage plot (top-left), land-use plot (bottom-left), and
## spatial map (right) into one figure, matching the reference layout
fig_policy_combined <- (drainage_plot / crop_plot) | map_plot
fig_policy_combined <- fig_policy_combined + plot_layout(widths = c(1, 1.3))

print(fig_policy_combined)

ggsave("nitrate_data_analysis/figures/fig_policy_combined.png", plot = fig_policy_combined,
       width = 14, height = 8, dpi = 300, bg = "white")







policy_comparison_combine_sf <- st_as_sf(policy_comparison_combine,coords = c("longitude", "latitude"),crs = 4326)


policy_comparison_combine_sf_test <- policy_comparison_combine_sf[policy_comparison_combine_sf$test==1,]

global_min <- min(c(policy_comparison_combine_sf$indirect_policy, 
                    policy_comparison_combine_sf$policy), na.rm = TRUE)
global_max <- log(2000)
my_limits <- c(global_min, global_max)




# visualize the indirect policy


## interpolate the indirect method to PLSS grid

# this one use predicted at both train and test to interpolate
indirect_policy_Gp <- GpGp::fit_model(policy_comparison_combine_sf$indirect_policy, locs = st_coordinates(policy_comparison_combine_sf),covfun_name = "matern_sphere")
# this one use predicted at test to interpolate
indirect_policy_Gp_test <- GpGp::fit_model(policy_comparison_combine_sf_test$indirect_policy, locs = st_coordinates(policy_comparison_combine_sf_test),covfun_name = "matern_sphere")

indirect_policy_plss <- GpGp::predictions(fit = indirect_policy_Gp,locs_pred = st_coordinates(plss_covariates_sf_combine),X_pred = rep(1,nrow(plss_covariates_sf_combine)))


plss_covariates_sf_combine$indirect_policy <- indirect_policy_plss

## Estimated Required Depth

map_plot <- ggplot() +
  geom_sf(data = plss_covariates_sf_combine, aes(color = indirect_policy), size = 0.5) +
  geom_sf(data = counties,color = "black", size = 0.2)+
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        breaks = log(c(10, 100, 500, 1000, 2000)),
                        limits = log(c(5,2000)),
                        labels = function(breaks) { round(exp(breaks), 0)}) +
  theme_void() +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm")
  )

print(map_plot)



## Linear projection of estimated Required Depth (Indirect policy) over covariates (drainage level, land use)


indirect_policy_gam <- mgcv::gam(indirect_policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel +  drainagecl + cafolog, data = plss_covariates_sf_combine)

## draw the land use plot

crop_coeffs <- tidy(indirect_policy_gam, parametric = TRUE) %>%
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
    x = "Well Depth (Log(Feet))",
    y = "Land Use"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

print(crop_plot)

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
    x = "Soil Drainage Class",
    y = "Well Depth (log(Feet))"
  ) +
  theme_minimal() +
  # Rotate x-axis labels for readability
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1)
  )

print(drainage_plot)

## combine the drainage plot (top-left), land-use plot (bottom-left), and
## spatial map (right) into one figure, matching the reference layout
fig_policy_combined <- (drainage_plot / crop_plot) | map_plot
fig_policy_combined <- fig_policy_combined + plot_layout(widths = c(1, 1.3))

print(fig_policy_combined)








## we relevel the data


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

policy_comparison_combine_sf_test$drainagecl = factor(policy_comparison_combine_sf_test$drainagecl,
                                                 levels = c("Very poorly drained","Poorly drained","Somewhat poorly drained","Moderately well drained","Well drained","Somewhat excessively drained","Excessively drained"))




policy_comparison_combine_sf_test$crop_type_combine <- factor(policy_comparison_combine_sf_test$crop_type_combine, levels = c(
  "Forest",
  "Undefined",
  "Corn",
  "Soybeans",
  "Grass",
  "Vegtables",
  "Developed"
))




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



policy_comparison_combine_sf$logWellDepth <- log(policy_comparison_combine_sf$WellDepth)

indirect_policy_gam <- mgcv::gam(indirect_policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel+drainagecl + cafolog, data = policy_comparison_combine_sf)

indirect_policy_gam <- mgcv::gam(indirect_policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel+drainagecl + cafolog, data = policy_comparison_combine_sf_test)

indirect_policy_gam <- mgcv::gam(policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel+drainagecl + cafolog + s(U), data = policy_comparison_combine_sf)

indirect_policy_gam <- mgcv::gam(policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel+drainagecl + cafolog + s(U), data = policy_comparison_combine_sf_test)


# land use categories shared by both representations, in a common display order
crop_order <- c("Corn", "Soybeans", "Grass", "Vegtables", "Developed")

# 1. get the coefficients
crop_coeffs <- tidy(indirect_policy_gam, parametric = TRUE) %>%
  filter(str_detect(term, "crop_type_combine") & !str_detect(term, "Undefined")) %>%
  mutate(
    # Calculate 95% Confidence Intervals
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,

    # Clean up term names: remove the prefix "crop_type_combine"
    term_clean = str_remove(term, "crop_type_combine")
  ) %>%
  filter(term_clean %in% crop_order) %>%
  mutate(term_clean = factor(term_clean, levels = crop_order))

#2. Create the Forest Plot
crop_plot <- ggplot(crop_coeffs, aes(x = estimate, y = term_clean)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "darkgreen") +

  # Add a vertical reference line at 0 (the value for "Forest")
  geom_vline(xintercept = 0, linetype = "dashed", color = "pink", size = 0.8) +

  labs(
    title = "Coefficient (vs. Forest baseline)",
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

## draw the crop/land use marginal means (Forest included, as the baseline level)
## ggplot draws the last level of a discrete y-scale at the top, so Forest goes
## last here to appear as the first (top) row in the plot
crop_order_marginal <- c(crop_order, "Forest")

emm_crop <- emmeans(indirect_policy_gam, ~ crop_type_combine)

emm_crop_df <- as.data.frame(emm_crop) %>%
  filter(crop_type_combine %in% crop_order_marginal) %>%
  mutate(crop_type_combine = factor(crop_type_combine, levels = crop_order_marginal))

crop_marginal_plot <- ggplot(emm_crop_df, aes(x = emmean, y = crop_type_combine)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2, color = "darkgreen") +
  labs(
    title = "Marginal Mean",
    x = "Well Depth (log(Feet))",
    y = "Land Use"
  ) +
  theme_minimal()

print(crop_marginal_plot)

## same marginal means, with x/y axes swapped (land use on x-axis, well depth on y-axis)
crop_marginal_plot_vertical <- ggplot(emm_crop_df, aes(x = crop_type_combine, y = emmean)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "darkgreen") +
  labs(
    title = "Marginal Mean",
    x = "Land Use",
    y = "Well Depth (log(Feet))"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(crop_marginal_plot_vertical)

## side-by-side comparison of the two representations, sharing the same
## land-use order so differences in shape/ranking are easy to spot
crop_compare_plot <- crop_plot + crop_marginal_plot
print(crop_compare_plot)


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



# Nonsp Indirect method v.s. Our method (acc/mcc)


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


