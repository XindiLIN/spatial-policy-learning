library(tigris)
library(ggplot2)
library(tidyr)
library(tidymodels)
library(dplyr)
library(patchwork)
library(sf)

source("functions/miscellaneous.r")

nitrate_data <- load_nitrate_data(file_path = "data/data_Nitrate_with_covar.csv", zero_inflated = FALSE)
nitrate_data <- nitrate_data[nitrate_data$SampleYear>=2020,]

nitrate_data_sf <- st_as_sf(nitrate_data,coords = c("longitude","latitude"),crs=4326)

# Load data

wi_counties_path <- "nitrate_data_analysis/output/wi_counties.rds"
plss_covariates_combined_path <- "nitrate_data_analysis/output/plss_covariates_sf_combine.rds"

if (file.exists(plss_covariates_combined_path)) {
  plss_covariates_sf_combine <- readRDS(plss_covariates_combined_path)
} else {
  plss_covariates_central_sf       <- readRDS("nitrate_data_analysis/output/plss_covariates_central_sf.rds")
  plss_covariates_east_central_sf  <- readRDS("nitrate_data_analysis/output/plss_covariates_east_central_sf.rds")
  plss_covariates_west_central_sf  <- readRDS("nitrate_data_analysis/output/plss_covariates_west_central_sf.rds")
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

# Dataset Visualization

## shared panel styling: centered bold titles, and legend text/keys sized up to be
## comparable to the title/caption text rather than the (too-small) ggplot defaults
panel_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.6, "cm")
  )

## land use/crop type
p_land_use <- ggplot() +
  geom_sf(data = counties,color = 'black', size = 0.2) +
  geom_sf(data = plss_covariates_sf_combine[!(plss_covariates_sf_combine$crop_type_combine %in% c("Undefined","Vegtables","Cover Crop","Other Crop")),], aes(color=crop_type_combine ),size=0.05) +
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
  labs(title = "Land Use", x = "Longitude", y = "Latitude") +
  panel_theme +
  guides(
    color = guide_legend(
      keywidth = unit(0.7, "cm"),  # Increase width of the legend squares
      keyheight = unit(0.7, "cm"),  # Increase height of the legend squares
      override.aes = list(shape = 15, size = 7)  # Set legend keys to squares (shape 15) and increase size
    )
  )




## CAFO

p_cafo <- ggplot() +
  geom_sf(data = counties,color = 'black', size = 0.2) +
  geom_sf(data = plss_covariates_sf_combine, aes(color= log(cafo+1)),size=0.05) +
  scale_color_viridis_c(option = "magma",name = "Ln(CAFO+1)") +  # Color scale
  labs(title = "CAFO") +
  panel_theme


## precipitation

p_precip <- ggplot() +
  geom_sf(data = counties,color = 'black', size = 0.2,fill = 'white') +
  geom_sf(data = plss_covariates_sf_combine, aes(color= precipitation),size=0.05) +
  scale_color_gradient(low = 'white',  # Light blue
                       high = "skyblue", # Dark blue
                       name = "mm/day")+
  labs(title = "Precipitation", x = "Longitude", y = "Latitude") +
  panel_theme


## drainage level

plss_covariates_sf_combine$drainagecl = factor(plss_covariates_sf_combine$drainagecl,
                                               levels = c("Very poorly drained","Poorly drained","Somewhat poorly drained","Moderately well drained","Well drained","Somewhat excessively drained","Excessively drained"))


p_drainage <- ggplot() +
  geom_sf(data = counties,color = 'black', size = 0.2,fill = 'lightgrey') +
  geom_sf(data = plss_covariates_sf_combine, aes(color=drainagecl), size = 0.05) +
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
  labs(title = "Soil Drainage") +
  panel_theme +
  guides(
    color = guide_legend(
      keywidth = unit(0.7, "cm"),  # Increase width of the legend squares
      keyheight = unit(0.7, "cm"),  # Increase height of the legend squares
      override.aes = list(shape = 15, size = 7)  # Set legend keys to squares (shape 15) and increase size
    )
  )


## static water level

p_water_level <- ggplot() +
  geom_sf(data = counties,color = 'black', size = 0.2,fill = 'white') +
  geom_sf(data = plss_covariates_sf_combine, aes(color= StaticLevel),size=0.05) +
  scale_color_gradient(low = 'white',  # Light blue
                       high = "cadetblue4", # Dark blue
                       name = "Feet",
                       # breaks = log(c(10, 100, 200,400)),
                       # limits = c(log(1), log(500)),
                       # labels = function(breaks) { round(exp(breaks), 0)}
  )+
  labs(title = "Static Water Level", x = "Longitude", y = "Latitude") +
  panel_theme


## observed nitrate concentrations
p_nitrate <- ggplot() +
  geom_sf(data = counties,color = "black", size = 0.2)+
  geom_sf(data = nitrate_data_sf, aes(color = logconcentration_plus_median), size = 0.05) +
  scale_color_viridis_c(option = "plasma",
                        name = "mg/L",
                        breaks = log(c(1,3,7,20,50)),
                        limits = c(log(0.5),log(50)),
                        labels = function(breaks) { round(exp(breaks), 0)}) +  # Color scale
  labs(title = "Nitrate Levels") +
  panel_theme

## observed well depth
p_well_depth <- ggplot() +
  geom_sf(data = counties,color = "black", size = 0.2)+
  geom_sf(data = nitrate_data_sf, aes(color = WellDepth), size = 0.05) +
  scale_color_viridis_c(option = "D",
                        name = "Feet",
                        breaks = c(0,250,500,750,1000),
                        limits = c(0,1000)) +  # Color scale
  labs(title = "Well Depth") +
  panel_theme

## combine all 7 panels into one figure -- 3 columns, filled row by row, matching
## the reference layout (Nitrate Levels/Land Use/Soil Drainage, Well Depth/CAFO/
## Precipitation, Static Water Level alone in the last row)
fig_combined <- wrap_plots(
  list(p_nitrate, p_land_use, p_drainage,
       p_well_depth, p_cafo, p_precip,
       p_water_level),
  ncol = 3
)


# dir.create("nitrate_data_analysis/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("nitrate_data_analysis/figures/fig_combined.png", plot = fig_combined,
       width = 14, height = 12, dpi = 300, bg = "white")
