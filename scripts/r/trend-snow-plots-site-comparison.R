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

# Blaesedalen ####

# All modelled NDVI values for Sentinel-2, Blaesedalen
s2.ndvi <- read.csv("../../data/sentinel-2/output/s2_all_models.csv") %>%
               st_as_sf(coords = c('X', 'Y'), crs = 32621)

# Snow data
s2.snow <- read.csv("../../data/uav/snow-metrics/blaesedalen-10m-auc-snowcover.csv") #%>%

# Join the snow and ndvi data
s2.data <- left_join(s2.ndvi, s2.snow, by = 'id') %>%
  select(-snow.persist, -X, -Y) %>%
  drop_na()

# Smoothed spline modelled data for S30
s30.ndvi <- read.csv('../../data/nasa-hls/blaesedalen/s30/output/s30_modelled_smoothed_spline_point_wide.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621)

# Snow data for S30
s30.snow <- read.csv("../../data/uav/snow-metrics/blaesedalen-30m-auc-snowcover.csv")

# Join the snow and ndvi data for S30
s30.data <- left_join(s30.ndvi, s30.snow, by = 'id') %>%
  select(-X, -Y) %>%
  drop_na()

# Kluane low ####
s2.ndvi.kl <- read.csv("../../data/sentinel-2/tidy-output/s2_kluane-low_smoothed_spline_point_wide.csv") %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)

# Snow, 10m
s2.snow.kl <- read.csv("../../data/uav/snow-metrics/kluane-low-10m-snow-cover.csv")

# Joined S2,kluane-low snow and NDVI
s2.kl.data <- left_join(s2.ndvi.kl, s2.snow.kl, by = 'id') %>%
  select(-X, -Y) %>%
  drop_na()

# S30 ndvi, kluane-low
s30.ndvi.kl <- read.csv("../../data/nasa-hls/output/s30_kluane-low-modelled_smoothed_spline_point_wide.csv") %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)

# Snow, 30m
s30.snow.kl <- read.csv("../../data/uav/snow-metrics/kluane-low-30m-snow-cover.csv")

# Joined S30, kluane-low snow and ndvi
s30.kl.data <- left_join(s30.ndvi.kl, s30.snow.kl, by = 'id') %>%
  select(-X, -Y) %>%
  drop_na()

# Kluane-high ####
s2.ndvi.kh <- read.csv("../../data/sentinel-2/tidy-output/s2_kluane-high_smoothed_spline_point_wide.csv") %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)

# Snow, 10m
s2.snow.kh <- read.csv("../../data/uav/snow-metrics/kluane-high-10m-snow-cover.csv")

# Joined S2, kluane-high snow and ndvi
s2.kh.data <- left_join(s2.ndvi.kh, s2.snow.kh, by = 'id') %>%
  select(-X, -Y) %>%
  drop_na()

# S30, kluane-high
s30.ndvi.kh <- read.csv("../../data/sentinel-2/tidy-output/s2_kluane-high_smoothed_spline_point_wide.csv") %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)

# snow, 30m
s30.snow.kh <- read.csv("../../data/uav/snow-metrics/kluane-high-30m-snow-cover.csv")

# Joined S30, kluane-high snow and ndvi
s30.kh.data <- left_join(s30.ndvi.kh, s30.snow.kh, by = 'id') %>%
  select(-X, -Y) %>%
  drop_na()

################################################################################

###
# Plots of max.ndvi,doy, max.ndvi, and snow.persist ----
###

# Set some plot parameters here:
xdoy <- c(200, 270)
xndvi <- c(0.1, 0.7)
ysnow <- c(0, 25)
yndvi <- c(0.1, 0.9)

# SENTINEL-2, with smoothed spline ----

# Blaesedalen
s2b.max.ndvi.doy.plot <- ggplot(s2.data, aes(y = ndvi.max.doy_s, x = snow.auc)) +
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
s2b.max.ndvi.plot <- ggplot(s2.data, aes(y = ndvi.max_s, x = snow.auc)) +
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
       y = 'Blaesedalen\nSentinel-2\nsmoothed-spline\n\n—————————————— Maximum NDVI ————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2b.ndvi.metrics.plot <- ggplot(drop_na(s2.data, snow.auc), aes(x = ndvi.max.doy_s, 
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

# Kluane-low ##
s2kl.max.ndvi.doy.plot <- ggplot(s2.kl.data, aes(y = ndvi.max.doy, x = snow.auc)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xdoy) +
  scale_color_viridis_c(name = "ndvi.max",
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
s2kl.max.ndvi.plot <- ggplot(s2.kl.data, aes(y = ndvi.max, x = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xndvi) +
  scale_color_viridis_c(name = "ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Kluane-low\nSentinel-2\nsmoothed-spline\n\n—————————————— Maximum NDVI ————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2kl.ndvi.metrics.plot <- ggplot(drop_na(s2.kl.data, snow.auc), aes(x = ndvi.max.doy, 
                                                                y = ndvi.max)) +
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


# Kluane-high ##
s2kh.max.ndvi.doy.plot <- ggplot(s2.kh.data, aes(y = ndvi.max.doy, x = snow.auc)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xdoy) +
  scale_color_viridis_c(name = "ndvi.max",
                        breaks = seq(0.3, 
                                     0.5, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0.1, 0.6),# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = 'Snow persistence', 
       y = '\n\n———————————— Maximum NDVI DoY —————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 max ndvi against snow persistence
s2kh.max.ndvi.plot <- ggplot(s2.kh.data, aes(y = ndvi.max, x = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xndvi) +
  scale_color_viridis_c(name = "ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = 'Snow persistence', 
       y = 'Kluane-high\nSentinel-2\nsmoothed-spline\n\n—————————————— Maximum NDVI ————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s2kh.ndvi.metrics.plot <- ggplot(drop_na(s2.kh.data, snow.auc), aes(x = ndvi.max.doy, 
                                                                    y = ndvi.max)) +
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
  labs(x = 'Maximum NDVI Day of Year', 
       y = '\n—————————————— Maximum NDVI ————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()


# Arrange plots side by side
plots_combined <-   s2b.max.ndvi.plot +
  s2b.max.ndvi.doy.plot +
  s2b.ndvi.metrics.plot +
  s2kl.max.ndvi.plot +
  s2kl.max.ndvi.doy.plot +
  s2kl.ndvi.metrics.plot +
  s2kh.max.ndvi.plot +
  s2kh.max.ndvi.doy.plot +
  s2kh.ndvi.metrics.plot +
  plot_layout(ncol = 3, nrow = 3)

plots_combined


# Plotting S30 data, smoothed spline ----
# Set some plot parameters here:
xdoy <- c(200, 270)
xndvi <- c(0.1, 0.95)
ysnow <- c(0, 25)
yndvi <- c(0.1, 0.9)

# Blaesedalen
s30b.max.ndvi.doy.plot <- ggplot(s30.data, aes(y = ndvi.max.doy, x = snow.auc)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xdoy) +
  scale_color_viridis_c(name = "ndvi.max",
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
s30b.max.ndvi.plot <- ggplot(s30.data, aes(y = ndvi.max, x = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xndvi) +
  scale_color_viridis_c(name = "ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Blaesedalen\nSentinel-2\nsmoothed-spline\n\n—————————————— Maximum NDVI ————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s30b.ndvi.metrics.plot <- ggplot(drop_na(s30.data, snow.auc), aes(x = ndvi.max.doy, 
                                                                y = ndvi.max)) +
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

# Kluane-low ##
s30kl.max.ndvi.doy.plot <- ggplot(s30.kl.data, aes(y = ndvi.max.doy, x = snow.auc)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xdoy) +
  scale_color_viridis_c(name = "ndvi.max",
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
s30kl.max.ndvi.plot <- ggplot(s30.kl.data, aes(y = ndvi.max, x = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xndvi) +
  scale_color_viridis_c(name = "ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = '', 
       y = 'Kluane-low\nSentinel-2\nsmoothed-spline\n\n—————————————— Maximum NDVI ————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()


# S2 ndvi metrics against each other  
s30kl.ndvi.metrics.plot <- ggplot(drop_na(s30.kl.data, snow.auc), aes(x = ndvi.max.doy, 
                                                                    y = ndvi.max)) +
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


# Kluane-high ##
s30kh.max.ndvi.doy.plot <- ggplot(s30.kh.data, aes(y = ndvi.max.doy, x = snow.auc)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xdoy) +
  scale_color_viridis_c(name = "ndvi.max",
                        breaks = seq(0.3, 
                                     0.5, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(0.1, 0.6),# Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = 'Snow persistence', 
       y = '\n\n———————————— Maximum NDVI DoY —————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 max ndvi against snow persistence
s30kh.max.ndvi.plot <- ggplot(s30.kh.data, aes(y = ndvi.max, x = snow.auc)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm', color = 'red') +
  xlim(ysnow) +
  ylim(xndvi) +
  scale_color_viridis_c(name = "ndvi.max.doy",
                        breaks = seq(220, 
                                     260, 
                                     length.out = 3),  # Adjust the number of breaks as needed
                        limits = c(210, 270),  # Set the limits to cover the entire range of lsts.ndvi.max.doy
                        guide = guide_colourbar(title = 'NDVI\nmax\nDoY',
                                                direction = 'vertical')) +
  labs(x = 'Snow persistence', 
       y = 'Kluane-high\nSentinel-2\nsmoothed-spline\n\n—————————————— Maximum NDVI ————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()

# S2 ndvi metrics against each other  
s30kh.ndvi.metrics.plot <- ggplot(drop_na(s30.kh.data, snow.auc), aes(x = ndvi.max.doy, 
                                                                    y = ndvi.max)) +
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
  labs(x = 'Maximum NDVI Day of Year', 
       y = '\n—————————————— Maximum NDVI ————————————\n', 
       padding = 1) +
  guides(color = 'none') +
  theme_cowplot()


# Arrange plots side by side
plots_combined <-   s30b.max.ndvi.plot +
  s30b.max.ndvi.doy.plot +
  s30b.ndvi.metrics.plot +
  s30kl.max.ndvi.plot +
  s30kl.max.ndvi.doy.plot +
  s30kl.ndvi.metrics.plot +
  s30kh.max.ndvi.plot +
  s30kh.max.ndvi.doy.plot +
  s30kh.ndvi.metrics.plot +
  plot_layout(ncol = 3, nrow = 3)

plots_combined



# Map plots ----
# Plotting NDVI metrics as map
ggplot() +
  geom_sf(data = st_buffer(s2.data, dist = 5, endCapStyle = "SQUARE"),
          aes(fill = ndvi.max_b)) +
  scale_fill_viridis_c(limits = c(0.1, 0.6))

# Read in raw NDVI data
ndvi <- read.csv('../../data/sentinel-2/output/sentinel-2-ndvi-ts-pt-2023.csv') %>%
  st_as_sf(coords = c('X', 'Y'))

# Join with snow data
ndvi.join <- left_join(ndvi, s2.snow, by = 'id') %>%
  drop_na()

# Format data for plotting
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

# Plot observations 
ggplot() +
  geom_line(data = ndvi.s %>% filter(snow.auc < 2), aes(x = doy, y = ndvi, group = id), color = '#9ebc9fff', alpha = 0.1, size = 1) + #size = snow.auc)) +
  geom_line(data = ndvi.s %>% filter(snow.auc > 15), aes(x = doy, y = ndvi, group = id), color = '#d08c6dff', alpha = 0.1, size = 2) +
  geom_point(data = s2.data %>% filter(snow.auc < 2), aes(x = ndvi.max.doy_s, y = ndvi.max_s), color = '#9ebc9fff', alpha = 0.2) +
  geom_point(data = s2.data %>% filter(snow.auc > 15), aes(x = ndvi.max.doy_s, y = ndvi.max_s), color = '#d08c6dff', alpha = 0.2) +
  #scale_color_viridis_c() +
  labs( x = 'Day of Year', 
        y = 'NDVI Observation') +
  theme_cowplot()

# Plot as map
ggplot() +
  geom_sf(data = st_buffer(ndvi.sample %>% filter(date == 229), dist = 5, endCapStyle = "SQUARE"), 
          aes(fill = ndvi)) +
  scale_fill_viridis_c(limits = c(0.1, 0.6))

# Bring in predicted NDVI from smoothed spline model
ndvi.s <- st_read('../../data/sentinel-2/output/s2_modelled_JA_point_long.shp') %>%
  select(id, doy_obs, ndv_prd)

# Rename variables
ndvi.s <- ndvi.s %>% rename(doy = 'doy_obs',
                            ndvi = 'ndv_prd')

# Join
ndvi.s <- left_join(ndvi.s, s2.snow, by = 'id') %>%
  select(id, doy, ndvi, snow.auc, geometry)

# Drop NAs
ndvi.s <- ndvi.s %>% drop_na()

# Join snow data 

# Plot observations 

# Create variables for filtering, indicative of lots of snow and little snow AUC
high.snow <- 10
low.snow <- 5

# Get the range of max doy from high and low snow
s2.data.low <- s2.data %>% filter(snow.auc < low.snow)
s2.data.high <- s2.data %>% filter(snow.auc > high.snow)

# Store min and max values
high.snow.max <- max(s2.data.high$ndvi.max.doy_s)
high.snow.min <- min(s2.data.high$ndvi.max.doy_s)
low.snow.max <- max(s2.data.low$ndvi.max.doy_s)
low.snow.min <- min(s2.data.low$ndvi.max.doy_s)

# Plot in the style of the conceptual plot from the beginning of this research proj.
ggplot() +
  geom_line(data = ndvi.s %>% filter(snow.auc < low.snow), aes(x = doy, y = ndvi, group = id), color = '#9ebc9fff', alpha = 0.1, size = 1) + #size = snow.auc)) +
  geom_line(data = ndvi.s %>% filter(snow.auc > high.snow), aes(x = doy, y = ndvi, group = id), color = '#d08c6dff', alpha = 0.1, size = 1) +
  geom_point(data = s2.data %>% filter(snow.auc < low.snow), aes(x = ndvi.max.doy_s, y = ndvi.max_s), color = '#9ebc9fff', alpha = 0.5) +
  geom_point(data = s2.data %>% filter(snow.auc > high.snow), aes(x = ndvi.max.doy_s, y = ndvi.max_s), color = '#d08c6dff', alpha = 0.5) +
  geom_segment(aes(x = high.snow.min, xend = high.snow.max, y = 0.07, yend = 0.07), color = '#d08c6dff', size = 1) +
  geom_segment(aes(x = low.snow.min, xend = low.snow.max, y = 0.08, yend = 0.08), color = '#9ebc9fff', size = 1) +
  #scale_color_viridis_c() +
  labs( x = 'Day of Year', 
        y = 'NDVI Observation') +
 # xlim(c(175, 260)) +
#  ylim(c(0, 0.6)) +
  coord_cartesian(xlim = c(175, 260), ylim = c(0, 0.6)) +
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


