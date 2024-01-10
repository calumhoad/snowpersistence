# Compare snow cover with derived landsat trends from LandsatTS package
# Calum Hoad, 06/12/2023

# Import necessary libraries
library(terra)
library(dplyr)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)
library(cowplot)

# Read in the data ----

## Sentinel-2 ##

# All modelled NDVI values for Sentinel-2
s2.ndvi <- read.csv("../../data/sentinel-2/output/s2_all_models.csv") %>%
               st_as_sf(coords = c('X', 'Y'), crs = 32621)

# Snow data
s2.snow <- read.csv("../../data/uav/snow-metrics/blaesedalen-10m-auc-snowcover.csv") #%>%

# Join the snow and ndvi data
s2.data <- left_join(s2.ndvi, s2.snow, by = 'id') %>%
  select(-snow.persist, -X, -Y) %>%
  drop_na()


## Landsat ##

# Trend data from LandsatTS, with automatic screening of pixels
trends.auto <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_auto_7yr_trnds.csv')

# Trend data from LandsatTS, with manual screening of pixels
trends.manual <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_manual_7yr_trnds.csv')

# Growing season data from LandsatTS, with automatic screening of pixels
gs.auto <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_auto_7yr_gs_metric.csv')

# Growing season data from LandsatTS, with manual screening of pixels
gs.manual <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_manual_7yr_gs_metric.csv')

# Landsat (30m) snow cover data, derived from UAV imagery
ls.snow <- read.csv2('../../data/uav/snow-metrics/blaesedalen_30m_snowcover.csv') %>%
  rename(sample.id = id, 
         ndvi.max = 'ndvi_mx', 
         ndvi.max.doy = 'ndv_mx_dy')

# Checking match between Landsat and LandsatTS pixel centres

# Get data
lsatTS.pix.centres <- st_read( '../../data/lsatTS-output/pixel_centres.shp') %>%
  st_transform('epsg:32621') %>%
  filter(row_number() < 169)
lsat.pix.centres <- st_read('../../data/sentinel-2/output/ls_modelled_point_wide.shp')  

# Plot to check for match
ggplot() +
  geom_sf(data = lsatTS.pix.centres, aes(color = 'red', size = 2)) +#, aes(color = 'red')) +
  geom_sf(data = lsat.pix.centres, aes(color = 'blue', size = 1))#, aes(color = 'blue'))

# Create single dataset containing snow persistence, landsat sinlgle yr ndvi preds, LandsatTS preds ----

# Read in the output from LandsatTS
lsatTS <- st_read('../../data/lsatTS-output/blaesedalen/lsatTS_auto_7_gs.shp') %>%
  st_transform(crs = 32621) %>%
  filter(year == 2023)

# From landsat
lsat <- st_read('../../data/uav/snow-metrics/blaesedalen_30m_snowcover_andLS.shp')

# Check locations are consistent
ggplot() +
  geom_sf(data = lsat.all, aes(color = 'purple', size = 3)) +
  geom_sf(data = lsatTS, aes(color = 'red', size = 2)) +#, aes(color = 'red')) +
  geom_sf(data = lsat, aes(color = 'blue', size = 1))#, aes(color = 'blue'))

# Spatial join lsatTS to lsat
lsts.all <- st_join(lsat, lsatTS, left = TRUE) %>%
  rename(ls.ndvi.max = 'ndvi_mx.x', 
         ls.ndvi.max.doy = 'ndv_mx_dy', 
         lsts.ndvi.max = 'ndvi_mx.y', 
         lsts.ndvi.max.doy = 'ndv_mx_d',
         snow.persist = 'snw_prs') %>%
  dplyr::select(ls.ndvi.max, ls.ndvi.max.doy, 
         lsts.ndvi.max, lsts.ndvi.max.doy, 
         snow.persist, 
         geometry)

# Checking NaN values in dataset


lsat.na <- lsat.all %>% filter(is.na(snow.persist))

ggplot() + geom_sf(data = lsat.na) # Confirms NA values are outside the drone rast

################################################################################

###
# Plots of max.ndvi,doy, max.ndvi, and snow.persist ----
###

# Set some plot parameters here:
xdoy <- c(200, 270)
xndvi <- c(0.1, 0.7)
ysnow <- c(0, 1)
yndvi <- c(0.1, 0.7)

# SENTINEL-2, with Beck model applied
s2b.max.ndvi.doy.plot <- ggplot(s2.beck, aes(x = ndvi.max.doy, y = snow.persist)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max",
                        breaks = seq(0.3, 
                                     0.5, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0.1, 0.6),# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = '', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 max ndvi against snow persistence
s2b.max.ndvi.plot <- ggplot(s2.beck, aes(x = ndvi.max, y = snow.persist)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm') +
  xlim(xndvi) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "lsts.ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Sentinel-2 Beck\n\nSnow persistence', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2b.ndvi.metrics.plot <- ggplot(drop_na(s2.beck, snow.persist), aes(x = ndvi.max.doy, 
                                                                       y = ndvi.max)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.persist, size = snow.persist)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.persist",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0, 0.6),
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = '', 
       y = 'Maximum NDVI', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# SENTINEL-2, with modelling edited by Jakob
s2j.max.ndvi.doy.plot <- ggplot(s2.snow.ja, aes(x = ndvi.max.doy, y = snow.persist)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max",
                        breaks = seq(0.3, 
                                     0.5, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0.1, 0.6),# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = '', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 max ndvi against snow persistence
s2j.max.ndvi.plot <- ggplot(s2.snow.ja, aes(x = ndvi.max, y = snow.persist)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm') +
  xlim(xndvi) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "lsts.ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Sentinel-2 j\n\nSnow persistence', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2j.ndvi.metrics.plot <- ggplot(drop_na(s2.snow.ja, snow.persist), aes(x = ndvi.max.doy, 
                                                                      y = ndvi.max)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.persist, size = snow.persist)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.persist",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0, 0.6),
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = '', 
       y = 'Maximum NDVI', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

#SENTINEL-2 (ndvi-filtered)
s2f.max.ndvi.doy.plot <- ggplot(s2.snow.2, aes(x = ndvi.max.doy, y = snow.persist)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max",
                        breaks = seq(0.3, 
                                     0.5, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0.1, 0.6),# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = '', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 max ndvi against snow persistence
s2f.max.ndvi.plot <- ggplot(s2.snow.2, aes(x = ndvi.max, y = snow.persist)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm') +
  xlim(xndvi) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "lsts.ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Sentinel-2 filtered\n\nSnow persistence', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2f.ndvi.metrics.plot <- ggplot(drop_na(s2.snow.2, snow.persist), aes(x = ndvi.max.doy, 
                                                                   y = ndvi.max)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.persist, size = snow.persist)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.persist",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0, 0.6),
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = '', 
       y = 'Maximum NDVI', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 (non-ndvi filtered) max ndvi doy against snow persistence
s2.max.ndvi.doy.plot <- ggplot(s2.snow, aes(x = ndvi.max.doy, y = snow.persist)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max",
                        breaks = seq(0.3, 
                                     0.5, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0.1, 0.6),# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = '', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 max ndvi against snow persistence
s2.max.ndvi.plot <- ggplot(s2.snow, aes(x = ndvi.max, y = snow.persist)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm') +
  xlim(xndvi) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "lsts.ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Sentinel-2\n\nSnow persistence', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2.ndvi.metrics.plot <- ggplot(drop_na(s2.snow, snow.persist), aes(x = ndvi.max.doy, 
                                            y = ndvi.max)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.persist, size = snow.persist)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.persist",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0, 0.6),
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = '', 
       y = 'Maximum NDVI', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()


#LANDSAT
# max ndvi doy against snow persistence
ls.max.ndvi.doy.plot <- ggplot(ls.snow, aes(x = ndvi.max.doy, y = snow.persist)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max",
                        breaks = seq(0.3, 
                                     0.5, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0.1, 0.6),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = '', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# max ndvi against snow persistence
ls.max.ndvi.plot <- ggplot(ls.snow, aes(x = ndvi.max, y = snow.persist)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm') +
  xlim(xndvi) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "lsts.ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Landsat 8/9\n\nSnow persistence', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# Landsat ndvi metrics against eachother
ls.ndvi.metrics.plot <- ggplot(drop_na(ls.snow, snow.persist), aes(x = ndvi.max.doy, 
                                                                   y = ndvi.max)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.persist, size = snow.persist)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.persist",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0, 0.6),
                        na.value = 'yellow', # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) +  
  labs(x = '', 
       y = 'Maximum NDVI', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

###
# Plotting all three datasets
###

#LandsatTS
# max ndvi doy against snow persistence
lsts.max.ndvi.doy.plot <- ggplot(lsts.all, aes(x = lsts.ndvi.max.doy, y = snow.persist)) +
  geom_point(aes(color = lsts.ndvi.max)) +
  geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max",
                        breaks = seq(0.3, 
                                     0.5, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0.1, 0.6),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax',
                                                direction = 'vertical')) +
  labs(x = 'Maximum NDVI Day of Year', 
       y = '', 
       padding = 1) +
  theme_cowplot()

# max ndvi against snow persistence
lsts.max.ndvi.plot <- ggplot(lsts.all, aes(x = lsts.ndvi.max, y = snow.persist)) +
  geom_point(aes(color = lsts.ndvi.max.doy)) +
  geom_smooth(method = 'lm') +
  xlim(xndvi) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "lsts.ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = 'Maximum NDVI', 
       y = 'LandsatTS\n\nSnow persistence', 
       padding = 1) +
  theme(legend.position = 'bottom', 
        legend.justification = 'left', 
        legend.direction = 'horizontal') +
  theme_cowplot()

# LandsatTS NDVI metrics against eachother
lsts.ndvi.metrics.plot <- ggplot(drop_na(lsts.all, snow.persist), aes(x = lsts.ndvi.max.doy, 
                                                                   y = lsts.ndvi.max)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.persist, size = snow.persist)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.persist",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0, 0.6),
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = 'Maximum NDVI Day of Year', 
       y = 'Maximum NDVI', 
       padding = 1) +
  theme_cowplot()

# Arrange plots side by side
plots_combined <- s2.max.ndvi.plot + 
  s2.max.ndvi.doy.plot +
  s2.ndvi.metrics.plot +
  s2f.max.ndvi.plot +
  s2f.max.ndvi.doy.plot +
  s2f.ndvi.metrics.plot +
  s2j.max.ndvi.plot +
  s2j.max.ndvi.doy.plot +
  s2j.ndvi.metrics.plot +
  s2b.max.ndvi.plot +
  s2b.max.ndvi.doy.plot +
  s2b.ndvi.metrics.plot +
  #ls.max.ndvi.plot +
  #ls.max.ndvi.doy.plot +
  #ls.ndvi.metrics.plot +
  #lsts.max.ndvi.plot +
  #lsts.max.ndvi.doy.plot +
  #lsts.ndvi.metrics.plot +
  plot_layout(ncol = 3, nrow = 4)

plots_combined


# Map plots
# Plotting NDVI metrics as map
ggplot() +
  geom_sf(data = st_buffer(s2.data, dist = 5, endCapStyle = "SQUARE"),
          aes(fill = ndvi.max_b)) +
  scale_fill_viridis_c(limits = c(0.1, 0.6))

