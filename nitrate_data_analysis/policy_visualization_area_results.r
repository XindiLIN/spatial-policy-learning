library(tigris)
library(ggplot2)
library(tidyr)
library(tidymodels)
library(dplyr)
library(patchwork)
library(stringr)
library(emmeans)

source("functions/miscellaneous.r")

# ============================================================
# Draws the same three figures as policy_visualization.r, but sourced from
# run_all_regions.R's real output (nitrate_data_analysis/output/area_results/)
# instead of the old pre-refactor cache/per-region files.
# ============================================================

region_results_dir  <- "nitrate_data_analysis/output/area_results"
wi_counties_path     <- "nitrate_data_analysis/output/wi_counties.rds"

plss_direct   <- readRDS(file.path(region_results_dir, "plss_sf_combine_direct.rds"))
plss_indirect <- readRDS(file.path(region_results_dir, "plss_sf_combine_indirect_nonspatial.rds"))

# Both are built from the same prep$plss per region (see run_all_regions.R's
# run_region_pipeline()), so they should already be in 1:1 row correspondence --
# verify that rather than assuming it silently.
if (nrow(plss_direct) != nrow(plss_indirect)) {
  stop(sprintf(
    "plss_sf_combine_direct.rds (%d rows) and plss_sf_combine_indirect_nonspatial.rds (%d rows) don't match -- can't combine into one data frame.",
    nrow(plss_direct), nrow(plss_indirect)
  ))
}

# Both files use the column name "policy" for their own method's estimate --
# combine into one data frame with both, matching the shape the GAM-fitting
# code below expects (policy = direct, indirect_policy = indirect).
plss_covariates_sf_combine <- plss_direct
plss_covariates_sf_combine$indirect_policy <- plss_indirect$policy

# The direct method's PLSS scoring (region_pipeline_funcs.R's fit_direct_dc_policy())
# is an unclipped RKHS extrapolation -- unlike the indirect method, which is
# guaranteed to stay within [depth_range[1], depth_range[2]] by construction
# (compute_indirect_policy()'s uniroot() with boundary clamps). A small
# fraction (~0.9%) of PLSS grid points end up with implausible direct-method
# values (e.g. > 10,000 ft) as a result. But depth_range is per-region and can
# itself exceed 2000 ft (indirect's overall range across regions was 22 to
# 2825.5 ft), so the map's shared limits = log(c(5,2000)) already silently
# drops a few indirect points too -- truncate BOTH columns the same way here,
# so the marginal-effect GAMs (fig2/fig3) see the same effective range for
# both methods that the map already shows visually, rather than an
# inconsistent, method-specific treatment.
plss_covariates_sf_combine$policy          <- pmin(pmax(plss_covariates_sf_combine$policy, log(5)), log(2000))
plss_covariates_sf_combine$indirect_policy <- pmin(pmax(plss_covariates_sf_combine$indirect_policy, log(5)), log(2000))

# NEW pipeline output stores crop_type_combine/drainagecl as character, not
# factor -- relevel here the same way policy_visualization.r does, for
# consistent category ordering in the plots below.
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

if (file.exists(wi_counties_path)) {
  counties <- readRDS(wi_counties_path)
} else {
  counties <- tigris::counties(state = "WI", cb = TRUE, year = 2022)
  counties <- counties[, c("NAME", "geometry")]
  colnames(counties) <- c("County", "geometry")

  dir.create(dirname(wi_counties_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(counties, wi_counties_path)
}


# ============================================================
# Figure 1: estimated required depth map, direct vs. non-spatial indirect
# ============================================================

direct_map_plot <- ggplot() +
  geom_sf(data = plss_covariates_sf_combine, aes(color = policy), size = 0.2) +
  geom_sf(data = counties, color = "black", size = 0.2) +
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

indirect_map_plot <- ggplot() +
  geom_sf(data = plss_covariates_sf_combine, aes(color = indirect_policy), size = 0.2) +
  geom_sf(data = counties, color = "black", size = 0.2) +
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

## combined map of direct vs. indirect, side by side sharing one legend
combined_map_plot <- direct_map_plot + indirect_map_plot + plot_layout(guides = "collect")

print(combined_map_plot)


# ============================================================
# Regression of the estimated required depth over covariates, for the
# marginal-effect figures below
# ============================================================

policy_gam <- mgcv::gam(policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel + s(U) + drainagecl + cafolog, data = plss_covariates_sf_combine)
indirect_policy_gam <- mgcv::gam(indirect_policy ~ crop_type_combine + drainagecl + precipitation + StaticLevel + drainagecl + cafolog, data = plss_covariates_sf_combine)

## marginal mean of required depth by land use
crop_order <- c("Corn", "Soybeans", "Grass", "Vegtables", "Developed")
crop_order_marginal <- c(crop_order, "Forest")

emm_crop <- emmeans(policy_gam, ~ crop_type_combine)
emm_crop_df <- as.data.frame(emm_crop) %>%
  filter(crop_type_combine %in% crop_order_marginal) %>%
  mutate(crop_type_combine = factor(crop_type_combine, levels = crop_order_marginal))

indirect_emm_crop <- emmeans(indirect_policy_gam, ~ crop_type_combine)
indirect_emm_crop_df <- as.data.frame(indirect_emm_crop) %>%
  filter(crop_type_combine %in% crop_order_marginal) %>%
  mutate(crop_type_combine = factor(crop_type_combine, levels = crop_order_marginal))

## marginal mean of required depth by drainage level
emm_drainage <- emmeans(policy_gam, ~ drainagecl)
indirect_emm_drainage <- emmeans(indirect_policy_gam, ~ drainagecl)

emm_drainage_df <- as.data.frame(emm_drainage)
indirect_emm_drainage_df <- as.data.frame(indirect_emm_drainage)


# ============================================================
# Figure 2: marginal mean of required depth by land use
# ============================================================

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
  )

print(crop_marginal_plot)


# ============================================================
# Figure 3: marginal mean of required depth by drainage level
# ============================================================

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
  )

print(drainage_plot)


# ============================================================
# Save
# ============================================================

ggsave("nitrate_data_analysis/figures/fig1_map_combined.png", plot = combined_map_plot,
       width = 14, height = 6, dpi = 300, bg = "white")

ggsave("nitrate_data_analysis/figures/fig2_marginal_land_use.png", plot = crop_marginal_plot,
       width = 10, height = 5, dpi = 300, bg = "white")

ggsave("nitrate_data_analysis/figures/fig3_marginal_drainage.png", plot = drainage_plot,
       width = 10, height = 5, dpi = 300, bg = "white")
