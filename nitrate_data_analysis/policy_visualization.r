library(tigris)
library(ggplot2)
library(tidyr)
library(tidymodels)
library(dplyr)
library(patchwork)
library(stringr)
library(emmeans)

source("functions/miscellaneous.r")


# Load PLSS covariates and the estimated policy at each PLSS cell 

wi_counties_path <- "nitrate_data_analysis/output/wi_counties.rds"
plss_covariates_combined_path <- "nitrate_data_analysis/output/plss_covariates_sf_combine.rds"


if (file.exists(plss_covariates_combined_path)) {
  plss_covariates_sf_combine <- readRDS(plss_covariates_combined_path)
} else {
  # load and combine the direct policy
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
  
  # load the indirect policy estimation at observed locations
  policy_comparison_combine <- readRDS("nitrate_data_analysis/output/policy_comparison_combine.rds")
  
  policy_comparison_combine_sf <- st_as_sf(policy_comparison_combine,coords = c("longitude", "latitude"),crs = 4326)
  
  # policy_comparison_combine_sf_test <- policy_comparison_combine_sf[policy_comparison_combine_sf$test==1,]
  
  # interpolate the indirect policy
  ## because the indirect method can do grid search at observed locations, we interpolate the estimation at these locations to the PLSS grids.
  
  indirect_policy_Gp <- GpGp::fit_model(policy_comparison_combine_sf$indirect_policy, locs = st_coordinates(policy_comparison_combine_sf), covfun_name = "matern_sphere")
  # this one use predicted at test to interpolate
  # indirect_policy_Gp_test <- GpGp::fit_model(policy_comparison_combine_sf_test$indirect_policy, locs = st_coordinates(policy_comparison_combine_sf_test),covfun_name = "matern_sphere")
  
  indirect_policy_plss <- GpGp::predictions(fit = indirect_policy_Gp,locs_pred = st_coordinates(plss_covariates_sf_combine),X_pred = rep(1,nrow(plss_covariates_sf_combine)))
  
  
  plss_covariates_sf_combine$indirect_policy <- indirect_policy_plss
  
  ### re-level the categorical data for better figure illustration
  
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


# Estimated Required Depth

## direct method

direct_map_plot <- ggplot() +
  geom_sf(data = plss_covariates_sf_combine, aes(color = policy), size = 0.2) +
  geom_sf(data = counties,color = "black", size = 0.2)+
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        breaks = log(c(10, 100, 500, 1000, 2000)),
                        limits = log(c(5,2000)),
                        labels = function(breaks) { round(exp(breaks), 0)}) +
  labs(title = "Direct") +
  theme_void() +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

print(direct_map_plot)

## indirect method

indirect_map_plot <- ggplot() +
  geom_sf(data = plss_covariates_sf_combine, aes(color = indirect_policy), size = 0.2) +
  geom_sf(data = counties,color = "black", size = 0.2)+
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        breaks = log(c(10, 100, 500, 1000, 2000)),
                        limits = log(c(5,2000)),
                        labels = function(breaks) { round(exp(breaks), 0)}) +
  labs(title = "Non-spatial Indirect") +
  theme_void() +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

print(indirect_map_plot)

## combined map of direct vs. indirect, side by side sharing one legend
combined_map_plot <- direct_map_plot + indirect_map_plot + plot_layout(guides = "collect")

print(combined_map_plot)

## zoomed-in view over the high-required-depth hotspot in central Wisconsin (Central
## Sands region, identified as the sub-area with the highest average direct-method
## depth), to make the direct vs. indirect comparison easier to see in detail
zoom_xlim <- c(-90.2, -88.9)
zoom_ylim <- c(43.4, 44.6)

# geom_sf's `size` is a fixed rendering size, not a geographic one -- reusing
# direct_map_plot/indirect_map_plot as-is would draw points at the same physical size
# as in the whole-state view, which looks proportionally much sparser once zoomed into
# a small fraction of the state. Scale the point size up by the same linear factor the
# view was zoomed in by, so the points look the same *relative* size/density here as
# they do in the whole-state map.
zoom_scale <- sqrt(
  (diff(sf::st_bbox(plss_covariates_sf_combine)[c("xmin", "xmax")]) / diff(zoom_xlim)) *
  (diff(sf::st_bbox(plss_covariates_sf_combine)[c("ymin", "ymax")]) / diff(zoom_ylim))
)
zoom_point_size <- 0.2 * zoom_scale

direct_map_zoom <- ggplot() +
  geom_sf(data = plss_covariates_sf_combine, aes(color = policy), size = zoom_point_size) +
  geom_sf(data = counties, color = "black", size = 0.2) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim, expand = FALSE) +
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        breaks = log(c(10, 100, 500, 1000, 2000)),
                        limits = log(c(5,2000)),
                        labels = function(breaks) { round(exp(breaks), 0)}) +
  labs(title = "Direct") +
  theme_void() +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

indirect_map_zoom <- ggplot() +
  geom_sf(data = plss_covariates_sf_combine, aes(color = indirect_policy), size = zoom_point_size) +
  geom_sf(data = counties, color = "black", size = 0.2) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim, expand = FALSE) +
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        breaks = log(c(10, 100, 500, 1000, 2000)),
                        limits = log(c(5,2000)),
                        labels = function(breaks) { round(exp(breaks), 0)}) +
  labs(title = "Non-spatial Indirect") +
  theme_void() +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

fig1_zoom_combined <- direct_map_zoom + indirect_map_zoom + plot_layout(guides = "collect")

print(fig1_zoom_combined)

ggsave("nitrate_data_analysis/figures/fig1_zoom_central.png", plot = fig1_zoom_combined,
       width = 14, height = 6, dpi = 300, bg = "white")


# Regression of the estimated Required Depth over covariates 

policy_gam <- mgcv::gam(policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel + s(U) +  drainagecl + cafolog, data = plss_covariates_sf_combine)
indirect_policy_gam <- mgcv::gam(indirect_policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel + drainagecl + cafolog, data = plss_covariates_sf_combine)

## marginal mean of required depth of land use
crop_order <- c("Corn", "Soybeans", "Grass", "Vegtables", "Developed")
crop_order_marginal <- c(crop_order, "Forest")

emm_crop <- emmeans(policy_gam, ~ crop_type_combine)

## convert it to data frame for ggplot2
emm_crop_df <- as.data.frame(emm_crop) %>%
  filter(crop_type_combine %in% crop_order_marginal) %>%
  mutate(crop_type_combine = factor(crop_type_combine, levels = crop_order_marginal))


indirect_emm_crop <- emmeans(indirect_policy_gam, ~ crop_type_combine)

indirect_emm_crop_df <- as.data.frame(indirect_emm_crop) %>%
  filter(crop_type_combine %in% crop_order_marginal) %>%
  mutate(crop_type_combine = factor(crop_type_combine, levels = crop_order_marginal))


## marginal required depth of different drainage level

emm_drainage <- emmeans(policy_gam, ~ drainagecl)

indirect_emm_drainage <- emmeans(indirect_policy_gam, ~ drainagecl)

emm_drainage_df <- as.data.frame(emm_drainage)

indirect_emm_drainage_df <- as.data.frame(indirect_emm_drainage)


# plot marginal required depth of different land use

crop_marginal_combined_df <- bind_rows(
  emm_crop_df %>% mutate(method = "Direct"),
  indirect_emm_crop_df %>% mutate(method = "Non-spatial Indirect")
)

crop_marginal_plot <- ggplot(crop_marginal_combined_df, aes(x = crop_type_combine, y = emmean, color = crop_type_combine, shape = method)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), linewidth = 0.8, width = 0.2, position = position_dodge(width = 0.4)) +
  # positions stay on the model's log scale, but labels are converted back to feet
  scale_y_continuous(labels = function(x) round(exp(x), 0)) +
  # numeric 1-6 codes (matching crop_order_marginal's order) so the
  # ordering reads directly off the x-axis, alongside the color legend
  scale_x_discrete(labels = as.character(seq_along(crop_order_marginal))) +
  # same land-use color scheme as dataset_visualization.R's p_land_use; the
  # legend labels carry the same 1-6 numbering as the x-axis (matching
  # crop_order_marginal's order) so the correspondence lives in the legend,
  # not spelled out again in the axis title
  scale_color_manual(
    values = c(
      "Forest" = "darkgreen",
      "Grass" = "lightgreen",
      "Corn" = "red",
      "Developed" = "gray",
      "Soybeans" = "blue",
      "Other Crop" = "orange",
      "Cover Crop" = "brown",
      "Vegtables" = "purple"
    ),
    labels = c(
      "Corn" = "1 Corn",
      "Soybeans" = "2 Soybeans",
      "Grass" = "3 Grass",
      "Vegtables" = "4 Vegetables",
      "Developed" = "5 Developed",
      "Forest" = "6 Forest"
    ),
    name = "Land Use"
  ) +
  scale_shape_manual(values = c("Direct" = 16, "Non-spatial Indirect" = 17)) +
  labs(
    # title = "Marginal Mean",
    x = "Land Use",
    y = "Well Depth (Feet)",
    shape = "Method"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm")
    # legend.position = "inside",
    # legend.position.inside = c(0.85, 0.85),
    # legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )

print(crop_marginal_plot)



# marginal required depth of different drainage level


drainage_marginal_combined_df <- bind_rows(
  emm_drainage_df %>% mutate(method = "Direct"),
  indirect_emm_drainage_df %>% mutate(method = "Non-spatial Indirect")
)

drainage_plot <- ggplot(drainage_marginal_combined_df, aes(x = drainagecl, y = emmean, shape = method)) +
  # Connect the means with a line to emphasize the trend (grey, distinguished by
  # linetype rather than color, since color is reserved for drainage class below)
  geom_line(aes(group = method, linetype = method), color = "grey40", linewidth = 0.8, position = position_dodge(width = 0.3)) +
  # Add points for the estimated means
  geom_point(aes(color = drainagecl), size = 2.8, position = position_dodge(width = 0.3)) +
  # Add 95% confidence intervals
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = drainagecl), linewidth = 0.8, width = 0.2, position = position_dodge(width = 0.3)) +

  # positions stay on the model's log scale, but labels are converted back to feet
  scale_y_continuous(labels = function(x) round(exp(x), 0)) +
  # numeric 1-7 codes (matching the drainagecl factor's level order set above)
  # so the ordering reads directly off the x-axis, alongside the color legend
  scale_x_discrete(labels = as.character(1:7)) +
  # same drainage-level color scheme (sequential green) as dataset_visualization.R's
  # p_drainage; legend labels carry the same 1-7 numbering as the x-axis so the
  # correspondence lives in the legend, not spelled out again in the axis title
  scale_color_brewer(palette = "Greens", name = "Soil Drainage Level",
                     labels = c(
                       "Very poorly drained" = "1 Very poor",
                       "Poorly drained" = "2 Poorly",
                       "Somewhat poorly drained" = "3 Somewhat poor",
                       "Moderately well drained" = "4 Moderately well",
                       "Well drained" = "5 Well",
                       "Somewhat excessively drained" = "6 Somewhat excessive",
                       "Excessively drained" = "7 Excessive"
                     )) +
  scale_shape_manual(values = c("Direct" = 16, "Non-spatial Indirect" = 17)) +
  scale_linetype_manual(values = c("Direct" = "solid", "Non-spatial Indirect" = "dashed")) +
  # force the Method legend above the Soil Drainage Level legend; ggplot's
  # default stacking order isn't controlled by scale/aes() declaration order
  # here since shape+linetype merge into one "Method" guide
  guides(
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) +
  labs(
    y = "Well Depth (Feet)",
    x  = "Soil Drainage Level",
    shape = "Method",
    linetype = "Method"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm")
    # legend.position = "inside",
    # legend.position.inside = c(0.15, 0.85),
    # legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )

print(drainage_plot)




## Figure 1: direct map (left) and indirect map (right)
ggsave("nitrate_data_analysis/figures/fig1_map_combined.png", plot = combined_map_plot,
       width = 14, height = 6, dpi = 300, bg = "white")

## Figure 2: marginal mean of depth over land use
ggsave("nitrate_data_analysis/figures/fig2_marginal_land_use.png", plot = crop_marginal_plot,
       width = 10, height = 5, dpi = 300, bg = "white")

## Figure 3: marginal mean of depth over drainage level
ggsave("nitrate_data_analysis/figures/fig3_marginal_drainage.png", plot = drainage_plot,
       width = 10, height = 5, dpi = 300, bg = "white")





# 
# # visualize the indirect policy
# 
# 
# ## Estimated Required Depth
# 
# 
# 
# 
# ## Linear projection of estimated Required Depth (Indirect policy) over covariates (drainage level, land use)
# 
# 
# 
# 
# ## draw the land use plot
# 
# crop_coeffs <- tidy(indirect_policy_gam, parametric = TRUE) %>%
#   filter(str_detect(term, "crop_type_combine") & !str_detect(term, "Undefined")) %>%
#   mutate(
#     # Calculate 95% Confidence Intervals
#     lower = estimate - 1.96 * std.error,
#     upper = estimate + 1.96 * std.error,
#     
#     # Clean up term names: remove the prefix "crop_type_combine"
#     term_clean = str_remove(term, "crop_type_combine")
#   )
# 
# # 2. Create the Forest Plot
# crop_plot <- ggplot(crop_coeffs, aes(x = estimate, y = term_clean)) +
#   geom_point(size = 3, color = "darkgreen") +
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "darkgreen") +
#   
#   # Add a vertical reference line at 0 (the value for "Forest")
#   geom_vline(xintercept = 0, linetype = "dashed", color = "pink", size = 0.8) +
#   
#   labs(
#     x = "Well Depth (Log(Feet))",
#     y = "Land Use"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 16)
#   )
# 
# print(crop_plot)
# 
# ## draw the drainage line
# emm_drainage <- emmeans(indirect_policy_gam, ~ drainagecl)
# 
# # 2. Convert to data frame for ggplot2
# emm_df <- as.data.frame(emm_drainage)
# 
# # 3. Plot the ordered trend
# drainage_plot <- ggplot(emm_df, aes(x = drainagecl, y = emmean)) +
#   # Connect the means with a line to emphasize the trend
#   geom_line(aes(group = 1), color = "darkblue", size = 1) +
#   # Add points for the estimated means
#   geom_point(size = 3, color = "darkblue") +
#   # Add 95% confidence intervals
#   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "darkblue") +
#   
#   labs(
#     x = "Soil Drainage Class",
#     y = "Well Depth (log(Feet))"
#   ) +
#   theme_minimal() +
#   # Rotate x-axis labels for readability
#   theme(
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 16),
#     axis.text.x = element_text(size = 16, angle = 45, hjust = 1)
#   )
# 
# print(drainage_plot)
# 
# ## combine the drainage plot (top-left), land-use plot (bottom-left), and
# ## spatial map (right) into one figure, matching the reference layout
# fig_policy_combined <- (drainage_plot / crop_plot) | map_plot
# fig_policy_combined <- fig_policy_combined + plot_layout(widths = c(1, 1.3))
# 
# print(fig_policy_combined)
# 
# 
# 
# 
# 
# 
# 
# 
# ## we relevel the data
# 
# 
# policy_comparison_combine_sf$drainagecl = factor(policy_comparison_combine_sf$drainagecl,
#                                                levels = c("Very poorly drained","Poorly drained","Somewhat poorly drained","Moderately well drained","Well drained","Somewhat excessively drained","Excessively drained"))
# 
# 
# 
# 
# policy_comparison_combine_sf$crop_type_combine <- factor(policy_comparison_combine_sf$crop_type_combine, levels = c(
#   "Forest",
#   "Undefined",
#   "Corn",
#   "Soybeans",
#   "Grass",
#   "Vegtables",
#   "Developed"
# ))
# 
# policy_comparison_combine_sf_test$drainagecl = factor(policy_comparison_combine_sf_test$drainagecl,
#                                                  levels = c("Very poorly drained","Poorly drained","Somewhat poorly drained","Moderately well drained","Well drained","Somewhat excessively drained","Excessively drained"))
# 
# 
# 
# 
# policy_comparison_combine_sf_test$crop_type_combine <- factor(policy_comparison_combine_sf_test$crop_type_combine, levels = c(
#   "Forest",
#   "Undefined",
#   "Corn",
#   "Soybeans",
#   "Grass",
#   "Vegtables",
#   "Developed"
# ))
# 
# 
# 
# 
# ## visualization
# ggplot() +
#   geom_sf(data = policy_comparison_combine_sf_test, aes(color = pmin(indirect_policy,log(2000))), size = 0.5) +
#   geom_sf(data = counties,color = "black", size = 0.5)+
#   scale_color_viridis_c(option = "D",
#                         name = "Feet",
#                         breaks = log(c(10, 100, 500, 1000, 2000)),
#                         limits = my_limits,
#                         labels = function(breaks) { round(exp(breaks), 0)}) +
#   labs(title = "Indirect Well Depth Policy") +
#   theme_void()
# 
# 
# 
# ggplot() +
#   geom_sf(data = policy_comparison_combine_sf_test, aes(color = pmin(policy,log(2000))), size = 0.5) +
#   geom_sf(data = counties,color = "black", size = 0.5)+
#   scale_color_viridis_c(option = "D",
#                         name = "Feet",
#                         breaks = log(c(10, 100, 500, 1000, 2000)),
#                         limits = my_limits,
#                         labels = function(breaks) { round(exp(breaks), 0)}) +
#   labs(title = "Direct Well Depth Policy") +
#   theme_void()
# 
# 
# 
# # land use categories shared by both representations, in a common display order
# crop_order <- c("Corn", "Soybeans", "Grass", "Vegtables", "Developed")
# 
# # 1. get the coefficients
# crop_coeffs <- tidy(indirect_policy_gam, parametric = TRUE) %>%
#   filter(str_detect(term, "crop_type_combine") & !str_detect(term, "Undefined")) %>%
#   mutate(
#     # Calculate 95% Confidence Intervals
#     lower = estimate - 1.96 * std.error,
#     upper = estimate + 1.96 * std.error,
# 
#     # Clean up term names: remove the prefix "crop_type_combine"
#     term_clean = str_remove(term, "crop_type_combine")
#   ) %>%
#   filter(term_clean %in% crop_order) %>%
#   mutate(term_clean = factor(term_clean, levels = crop_order))
# 
# #2. Create the Forest Plot
# crop_plot <- ggplot(crop_coeffs, aes(x = estimate, y = term_clean)) +
#   geom_point(size = 3, color = "darkgreen") +
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "darkgreen") +
# 
#   # Add a vertical reference line at 0 (the value for "Forest")
#   geom_vline(xintercept = 0, linetype = "dashed", color = "pink", size = 0.8) +
# 
#   labs(
#     title = "Coefficient (vs. Forest baseline)",
#     x = "Well Depth (Log(Feet))",
#     y = "Land Use"
#   ) +
#   theme(
#     panel.grid.major = element_line(color = "gray10", size = 0.5, linetype = "dotted"),
#     panel.grid.minor = element_line(color = "gray10", size = 0.25, linetype = "dotted"),
# 
#     # You are also setting the plot border/axis lines here:
#     axis.line = element_line(color = "black", size = 0.5),
#     # You might want to remove the border box drawn by theme_bw
#     panel.border = element_blank()
#   ) + theme_minimal()
# 
# print(crop_plot)
# 
# 
# 
# 
# ## draw the drainage line
# emm_drainage <- emmeans(indirect_policy_gam, ~ drainagecl)
# 
# # 2. Convert to data frame for ggplot2
# emm_df <- as.data.frame(emm_drainage) 
# 
# # 3. Plot the ordered trend
# drainage_plot <- ggplot(emm_df, aes(x = drainagecl, y = emmean)) +
#   # Connect the means with a line to emphasize the trend 
#   geom_line(aes(group = 1), color = "darkblue", size = 1) + 
#   # Add points for the estimated means
#   geom_point(size = 3, color = "darkblue") +
#   # Add 95% confidence intervals
#   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "darkblue") +
#   
#   labs(
#     title = "Predicted Policy Value Across Soil Drainage Classes (GAM)",
#     subtitle = "Trend adjusted for smooth effects of continuous covariates",
#     x = "Soil Drainage Class ",
#     y = "Well Depth (log(Feet))"
#   ) +
#   theme_minimal() +
#   # Rotate x-axis labels for readability
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# print(drainage_plot)
# 
# 
# 
