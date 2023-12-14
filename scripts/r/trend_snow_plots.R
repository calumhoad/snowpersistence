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

# Sentinel-2 (30,) snow cover data, drived from UAV imagery
s2.snow <- read.csv2('../../data/uav/snow-metrics/blaesedalen_10m_snowcover.csv') %>%
  rename(ndvi.max = 'ndvi_mx', ndvi.max.doy = 'ndv_mx_') %>%
  mutate(ndvi.max.doy = yday(ndvi.max.doy))

# Join snow and NDVI metrics, if separate datasets ----
# Join data frames based on sample_id
joined.auto <- left_join(gs.auto, ls.snow, by = "sample.id")

# Filter gs datasets to only 2023 and max.ndvi.doy greater than 175
gs.auto <- gs.auto %>% filter(year == 2023 & ndvi.max.doy > 175)
gs.auto.joined <- left_join(gs.auto, snow, by = 'sample.id')

###
# Plots of max.ndvi,doy, max.ndvi, and snow.persist ----
###

#SENTINEL-2
# S2 max ndvi doy against snow persistence
s2.max.ndvi.doy.plot <- ggplot(s2.snow, aes(x = ndvi.max.doy, y = snow.persist)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm') +
  scale_color_viridis_c(name = "ndvi.max") +
  labs(x = 'Maximum NDVI Day of Year', 
       y = 'Snow persistence\n(unweighted average)', 
       padding = 1) +
  theme_cowplot()

# S2 max ndvi against snow persistence
s2.max.ndvi.plot <- ggplot(s2.snow, aes(x = ndvi.max, y = snow.persist)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm') +
  scale_color_viridis_c(name = "ndvi.max.doy") +
  labs(x = 'Maximum NDVI', 
       y = 'Sentinel-2\n\nSnow persistence\n(unweighted average)', 
       padding = 1) +
  theme_cowplot()
  
s2.ndvi.metrics.plot <- ggplot(drop_na(s2.snow, snow.persist), aes(x = ndvi.max.doy, 
                                            y = ndvi.max)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.persist, size = snow.persist)) +
  #geom_smooth(method = 'lm') +
  scale_color_viridis_c(name = "snow.persist") +
  labs(x = 'Maximum NDVI Day of Year', 
       y = 'Maximum NDVI', 
       padding = 1) +
  theme_cowplot()

#LANDSAT
# max ndvi doy against snow persistence
ls.max.ndvi.doy.plot <- ggplot(ls.snow, aes(x = ndvi.max.doy, y = snow.persist)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm') +
  scale_color_viridis_c(name = "ndvi.max") +
  labs(x = 'Maximum NDVI Day of Year', 
       y = 'Snow persistence\n(unweighted average)', 
       padding = 1) +
  theme_cowplot()

# max ndvi against snow persistence
ls.max.ndvi.plot <- ggplot(ls.snow, aes(x = ndvi.max, y = snow.persist)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm') +
  scale_color_viridis_c(name = "ndvi.max.doy") +
  labs(x = 'Maximum NDVI', 
       y = 'LANDSAT\n\nSnow persistence\n(unweighted average)', 
       padding = 1) +
  theme_cowplot()

ls.ndvi.metrics.plot <- ggplot(drop_na(ls.snow, snow.persist), aes(x = ndvi.max.doy, 
                                                                   y = ndvi.max)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.persist, size = snow.persist)) +
  #geom_smooth(method = 'lm') +
  scale_color_viridis_c(name = "snow.persist") +
  labs(x = 'Maximum NDVI Day of Year', 
       y = 'Maximum NDVI', 
       padding = 1) +
  theme_cowplot()


# Arrange plots side by side
plots_combined <- s2.max.ndvi.plot + 
  s2.max.ndvi.doy.plot +
  s2.ndvi.metrics.plot +
  ls.max.ndvi.plot +
  ls.max.ndvi.doy.plot +
  ls.ndvi.metrics.plot +
  plot_layout(ncol = 3, nrow = 2)

plots_combined

###
# Checking match between Landsat and LandsatTS pixel centres
###

# Get data
lsatTS.pix.centres <- st_read( '../../data/lsatTS-output/pixel_centres.shp') %>%
  st_transform('epsg:32621') %>%
  filter(row_number() < 169)
lsat.pix.centres <- st_read('../../data/sentinel-2/output/ls_modelled_point_wide.shp')  

# Plot to check for match
ggplot() +
  geom_sf(data = lsatTS.pix.centres, aes(color = 'red', size = 2)) +#, aes(color = 'red')) +
  geom_sf(data = lsat.pix.centres, aes(color = 'blue', size = 1))#, aes(color = 'blue'))
