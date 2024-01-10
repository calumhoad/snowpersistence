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
ysnow <- c(0, 25)
yndvi <- c(0.1, 0.7)

# SENTINEL-2, with Beck model applied
s2b.max.ndvi.doy.plot <- ggplot(s2.data, aes(x = ndvi.max.doy_b, y = snow.auc)) +
  geom_point(aes(color = ndvi.max_b)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(xdoy) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max_b",
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
s2b.max.ndvi.plot <- ggplot(s2.data, aes(x = ndvi.max_b, y = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy_b)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(xndvi) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max.doy_b",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Sentinel-2\nBeck (2006) double logistic\n\nSnow persistence', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2b.ndvi.metrics.plot <- ggplot(drop_na(s2.data, snow.auc), aes(x = ndvi.max.doy_b, 
                                                                       y = ndvi.max_b)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.auc, size = snow.auc)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.auc",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = ysnow,
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = '', 
       y = 'Maximum NDVI', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# SENTINEL-2, smoothed spline
s2s.max.ndvi.doy.plot <- ggplot(s2.data, aes(x = ndvi.max.doy_s, y = snow.auc)) +
  geom_point(aes(color = ndvi.max_s)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(xdoy) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max_s",
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
s2s.max.ndvi.plot <- ggplot(s2.data, aes(x = ndvi.max_s, y = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy_s)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(xndvi) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max.doy_s",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Sentinel-2\nsmoothed-spline\n\nSnow persistence', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2s.ndvi.metrics.plot <- ggplot(drop_na(s2.data, snow.auc), aes(x = ndvi.max.doy_s, 
                                                                      y = ndvi.max_s)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.auc, size = snow.auc)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.auc",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = ysnow,
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = '', 
       y = 'Maximum NDVI', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()


#SENTINEL-2, parabolic 2nd order polynomial

# S2 (non-ndvi filtered) max ndvi doy against snow persistence
s2p.max.ndvi.doy.plot <- ggplot(s2.data, aes(x = ndvi.max.doy_p, y = snow.auc)) +
  geom_point(aes(color = ndvi.max_p)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(xdoy) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max_p",
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
s2p.max.ndvi.plot <- ggplot(s2.data, aes(x = ndvi.max_p, y = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy_p)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(xndvi) +
  ylim(ysnow) +
  scale_color_viridis_c(name = "ndvi.max.doy_p",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Sentinel-2\n parabolic 2nd order poynomial\n\nSnow persistence', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2p.ndvi.metrics.plot <- ggplot(drop_na(s2.data, snow.auc), aes(x = ndvi.max.doy_p, 
                                            y = ndvi.max_p)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.auc, size = snow.auc)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.auc",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = ysnow,
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = '', 
       y = 'Maximum NDVI', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()



# Arrange plots side by side
plots_combined <- s2p.max.ndvi.plot + 
  s2p.max.ndvi.doy.plot +
  s2p.ndvi.metrics.plot +
  s2s.max.ndvi.plot +
  s2s.max.ndvi.doy.plot +
  s2s.ndvi.metrics.plot +
  s2b.max.ndvi.plot +
  s2b.max.ndvi.doy.plot +
  s2b.ndvi.metrics.plot +
  #s2b.max.ndvi.plot +
  #s2b.max.ndvi.doy.plot +
  #s2b.ndvi.metrics.plot +
  #ls.max.ndvi.plot +
  #ls.max.ndvi.doy.plot +
  #ls.ndvi.metrics.plot +
  #lsts.max.ndvi.plot +
  #lsts.max.ndvi.doy.plot +
  #lsts.ndvi.metrics.plot +
  plot_layout(ncol = 3, nrow = 3)

plots_combined


# Same plots again, but with the axis switched ----

# SENTINEL-2, with Beck model applied
s2b.max.ndvi.doy.plot <- ggplot(s2.data, aes(y = ndvi.max.doy_b, x = snow.auc)) +
  geom_point(aes(color = ndvi.max_b)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xdoy) +
  scale_color_viridis_c(name = "ndvi.max_b",
                        breaks = seq(0.3, 
                                     0.5, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0.1, 0.6),# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = 'Snow persistence', 
       y = '', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 max ndvi against snow persistence
s2b.max.ndvi.plot <- ggplot(s2.data, aes(y = ndvi.max_b, x = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy_b)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xndvi) +
  scale_color_viridis_c(name = "ndvi.max.doy_b",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = 'Snow persistence', 
       y = 'Sentinel-2\nBeck (2006) double logistic\n\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2b.ndvi.metrics.plot <- ggplot(drop_na(s2.data, snow.auc), aes(x = ndvi.max.doy_b, 
                                                                y = ndvi.max_b)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.auc, size = snow.auc)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.auc",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = ysnow,
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = 'Maximum NDVI DoY', 
       y = '', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# SENTINEL-2, smoothed spline
s2s.max.ndvi.doy.plot <- ggplot(s2.data, aes(y = ndvi.max.doy_s, x = snow.auc)) +
  geom_point(aes(color = ndvi.max_s)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xdoy) +
  scale_color_viridis_c(name = "ndvi.max_s",
                        breaks = seq(0.3, 
                                     0.5, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0.1, 0.6),# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = '\n\n———————————— Maximum NDVI DoY —————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()



# S2 max ndvi against snow persistence
s2s.max.ndvi.plot <- ggplot(s2.data, aes(y = ndvi.max_s, x = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy_s)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xndvi) +
  scale_color_viridis_c(name = "ndvi.max.doy_s",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Sentinel-2\nsmoothed-spline\n\n—————————————— Maximum NDVI ————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2s.ndvi.metrics.plot <- ggplot(drop_na(s2.data, snow.auc), aes(x = ndvi.max.doy_s, 
                                                                y = ndvi.max_s)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.auc, size = snow.auc)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.auc",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = ysnow,
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = '', 
       y = '\n—————————————— Maximum NDVI ————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()


#SENTINEL-2, parabolic 2nd order polynomial

# S2 (non-ndvi filtered) max ndvi doy against snow persistence
s2p.max.ndvi.doy.plot <- ggplot(s2.data, aes(y = ndvi.max.doy_p, x = snow.auc)) +
  geom_point(aes(color = ndvi.max_p)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xdoy) +
  scale_color_viridis_c(name = "ndvi.max_p",
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

s2p.max.ndvi.plot

# S2 max ndvi against snow persistence
s2p.max.ndvi.plot <- ggplot(s2.data, aes(y = ndvi.max_p, x = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy_p)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xndvi) +
  scale_color_viridis_c(name = "ndvi.max.doy_p",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Sentinel-2\n parabolic 2nd order poynomial\n\n\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2p.ndvi.metrics.plot <- ggplot(drop_na(s2.data, snow.auc), aes(x = ndvi.max.doy_p, 
                                                                y = ndvi.max_p)) +
  geom_point(position = 'jitter', alpha = 0.5, aes(color = snow.auc, size = snow.auc)) +
  #geom_smooth(method = 'lm') +
  xlim(xdoy) +
  ylim(yndvi) +
  scale_color_viridis_c(name = "snow.auc",
                        breaks = seq(0, 
                                     0.6, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = ysnow,
                        na.value = 'yellow',# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'Snow\nPersist',
                                                direction = 'vertical')) + 
  labs(x = 'Maximum NDVI DoY', 
       y = '', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()



# Arrange plots side by side
plots_combined <- s2p.max.ndvi.plot + 
  s2p.max.ndvi.doy.plot +
  s2p.ndvi.metrics.plot +
  s2s.max.ndvi.plot +
  s2s.max.ndvi.doy.plot +
  s2s.ndvi.metrics.plot +
  s2b.max.ndvi.plot +
  s2b.max.ndvi.doy.plot +
  s2b.ndvi.metrics.plot +
  #s2b.max.ndvi.plot +
  #s2b.max.ndvi.doy.plot +
  #s2b.ndvi.metrics.plot +
  #ls.max.ndvi.plot +
  #ls.max.ndvi.doy.plot +
  #ls.ndvi.metrics.plot +
  #lsts.max.ndvi.plot +
  #lsts.max.ndvi.doy.plot +
  #lsts.ndvi.metrics.plot +
  plot_layout(ncol = 3, nrow = 3)

plots_combined

# Map plots
# Plotting NDVI metrics as map
ggplot() +
  geom_sf(data = st_buffer(s2.data, dist = 5, endCapStyle = "SQUARE"),
          aes(fill = ndvi.max_b)) +
  scale_fill_viridis_c(limits = c(0.1, 0.6))

ndvi <- read.csv('../../data/sentinel-2/output/sentinel-2-ndvi-ts-pt-2023.csv') %>%
  st_as_sf(coords = c('X', 'Y'))

ndvi.join <- left_join(ndvi, s2.snow, by = 'id') %>%
  drop_na()

ndvi.sample <- ndvi.join %>% #filter(id == 300) %>%
  select(-X, -Y, -snow.av, -X2023.07.02, -X2023.07.12, -X2023.07.18, -X2023.07.26.y) %>%
  pivot_longer(!id & !geometry & !snow.auc, names_to = 'date', values_to = 'ndvi') %>%
  mutate(date = ifelse(date == 'X2023.06.26', yday('2023-06-26'),
                ifelse(date == 'X2023.07.08', yday('2023-07-08'), 
                       ifelse(date == 'X2023.07.26.x', yday('2023-07-26'), 
                              ifelse(date == 'X2023.07.29', yday('2023-07-29'), 
                                     ifelse(date == 'X2023.08.07', yday('2023-08-07'), 
                                            ifelse(date == 'X2023.08.08', yday('2023-08-08'), 
                                                   ifelse(date == 'X2023.08.17', yday('2023-08-17'), 
                                                          ifelse(date == 'X2023.09.23', yday('2023-09-23'), NA)))))))))

ggplot() +
  geom_line(data = ndvi.sample %>% filter(snow.auc < 2), aes(x = date, y = ndvi, group = id), color = '#9ebc9fff', size = 1) + #size = snow.auc)) +
  geom_line(data = ndvi.sample %>% filter(snow.auc > 15), aes(x = date, y = ndvi, group = id), color = '#d08c6dff', size = 2) +
  #scale_color_viridis_c() +
  labs( x = 'Day of Year', 
        y = 'NDVI Observation') +
  theme_cowplot()
# Un-used code ----
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


