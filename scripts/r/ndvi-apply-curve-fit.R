# Applying curve fit to NDVI data from all sites and sensors
# Calum Hoad, 25 Jan 2024

# Import necessary packages ----
library(terra)
library(dplyr)
library(rts)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(purrr)
library(broom)
library(viridis)
library(tidyverse)
library(bfast)
library(phenopix)
library(greenbrown)
library(greenbrown)
library(ggplot2)
library(cowplot)

# Source curve fitting functions from script
source('ndvi-curve-fitting-functions.R')


# Bring in the data ----

# SENTINEL-2 ## 
# Blaesedalen
s2.bl <- read.csv('../../data/ndvi/s2-blaesedalen-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621)
# Kluane-low
s2.kl <- read.csv('../../data/ndvi/s2-kluane-low-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)
# Kluane-high
s2.kh <- read.csv('../../data/ndvi/s2-kluane-high-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)

# NASA HLS S30 ###
s30.bl <- read.csv('../../data/ndvi/s30-blaesedalen-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621)

# Format the data ----

# Function for formatting the data, 
#   including setting all NDVI values < 0 to 0.
long_ndvi <- function(df) {
  df <- df %>%
    pivot_longer(!geometry & !id, names_to = 'doy', values_to = 'ndvi') %>%
    mutate(doy = sub('X', '', doy), 
           doy = sub('\\.', '-', doy), 
           doy = as_date(doy),
           doy = yday(doy)) %>%
    mutate(ndvi = ifelse(ndvi < 0, 0, ndvi)) %>%
    group_by(id)
}

# For S30, drop all data from Sept onwards due to quality error
s30.bl <- s30.bl %>% 
  select(-X2023.09.14, -X2023.09.22, -X2023.09.23, -X2023.10.03)

# Apply function
s2.bl <- long_ndvi(s2.bl)
s2.kl <- long_ndvi(s2.kl)
s2.kh <- long_ndvi(s2.kh)

s30.bl <- long_ndvi(s30.bl)

# Create a test dataset, where the outlying value on 2023-10-03 for Kluane is removed
s2.bl <- s2.bl %>% filter(doy != yday('2023-10-03'))

# Check the time raw time series data ----
ggplot() +
  geom_line(data = s2.bl, aes(x = doy, y = ndvi, group = id))
ggplot() +
  geom_line(data = s2.kl, aes(x = doy, y = ndvi, group = id))
ggplot() +
  geom_line(data = s2.kh, aes(x = doy, y = ndvi, group = id))
ggplot() +
  geom_line(data = s30.bl, aes(x = doy, y = ndvi, group = id))

# Fit smoothed spline model to the data ----
s2.bl.smooth <- s2.bl %>%
  group_modify(~model_ndvi_smoothedspline(.x))
s2.kl.smooth <- s2.kl %>%
  group_modify(~model_ndvi_smoothedspline(.x))
s2.kh.smooth <- s2.kh %>%
  group_modify(~model_ndvi_smoothedspline(.x))

s30.bl.smooth <- s30.bl %>%
  group_modify(~model_ndvi_smoothedspline(.x))

# Fit Beck double logistic to the data ----
s2.bl.beck <- s2.bl %>%
  group_modify(~model_ndvi_beck(.x))
s2.kl.beck <- s2.kl %>%
  group_modify(~model_ndvi_beck(.x))
s2.kh.beck <- s2.kh %>%
  group_modify(~model_ndvi_beck(.x))

# Format data 
s2.bl.smooth <- s2.bl.smooth %>%
  rename(ndvi.pred = 'ndvi.pred.doy.1')
s2.kl.smooth <- s2.kl.smooth %>%
  rename(ndvi.pred = 'ndvi.pred.doy.1')
s2.kh.smooth <- s2.kh.smooth %>%
  rename(ndvi.pred = 'ndvi.pred.doy.1')

s30.bl.smooth <- s30.bl.smooth %>%
  rename(ndvi.pred = 'ndvi.pred.doy.1')

# Join data to singular dataframes?

# Can join to singular sf dataframe, but Blaesedalen and Kluane have different 
#   projections (UTM08N and 21N). 

# Outputs ----

# Beck ##
# Blaesedalen
st_write(st_as_sf(s2.bl.beck),  '../../data/ndvi/s2-bl-beck-outlier-removed.csv', 
         layer_options = "GEOMETRY=AS_XY")
# Kluane-low
st_write(st_as_sf(s2.kl.beck),  '../../data/ndvi/s2-kl-beck.csv', 
         layer_options = "GEOMETRY=AS_XY")
# Kluane-high
st_write(st_as_sf(s2.kh.beck),  '../../data/ndvi/s2-kh-beck.csv', 
         layer_options = "GEOMETRY=AS_XY")
# Smoothed-spline ##
# Blaesedalen
st_write(st_as_sf(s2.bl.smooth),  '../../data/ndvi/s2-bl-smooth-outlier-removed.csv', 
         layer_options = "GEOMETRY=AS_XY")
# Kluane-low
st_write(st_as_sf(s2.kl.smooth),  '../../data/ndvi/s2-kl-smooth.csv', 
         layer_options = "GEOMETRY=AS_XY")
# Kluane-high
st_write(st_as_sf(s2.kh.smooth),  '../../data/ndvi/s2-kh-smooth.csv', 
         layer_options = "GEOMETRY=AS_XY")

# S30 Blaesedalen smooth
st_write(st_as_sf(s30.bl.smooth), '../../data/ndvi/s30-bl-smooth.csv',
         layer_options = 'GEOMETRY=AS_XY')


# Plot out the data for comparison of model fit ----
# 100 random pixels overview
plot.data <- s2.bl.smooth

ggplot(
  plot.data %>% filter(id %in% sample(unique(plot.data$id), 100)),
  aes(x = doy, colour = id, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred)) +
  theme_classic() +
  theme(legend.position = "none")

# 9 random pixels in detail
rand_id <- sample(plot.data %>% st_drop_geometry() %>% pull(id) %>% unique(), 9)

# For one model
plot.data <- s2.bl.smooth # Which model?

ggplot(
  plot.data %>% filter(id %in% rand_id),
  aes(x = doy, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred)) +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id) +
  theme_classic()

# For both models
ggplot() +
  # Beck
  geom_point(data = s2.bl.beck %>% filter(id %in% rand_id),
             aes(x = doy, y = ndvi, group = id)) +
  geom_line(data = s2.bl.beck %>% filter(id %in% rand_id),
            aes(x = doy,y = ndvi.pred, color = 'red')) +
  geom_point(data = s2.bl.beck %>% filter(id %in% rand_id),
             aes(x = ndvi.max.doy, y = ndvi.max), color = "blue") +
  # Smoothed spline
  geom_point(data = s2.bl.smooth %>% filter(id %in% rand_id),
             aes(x = doy, y = ndvi, group = id)) +
  geom_line(data = s2.bl.smooth %>% filter(id %in% rand_id),
            aes(x = doy, y = ndvi.pred, color = 'blue')) +
  geom_point(data = s2.bl.smooth %>% filter(id %in% rand_id),
             aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id) +
  theme_cowplot() +
  theme(legend.position = 'none')

# For both models
ggplot() +
  # Beck
  geom_point(data = s2.bl.beck %>% filter(id %in% rand_id),
             aes(x = doy, y = ndvi, group = id)) +
  geom_line(data = s2.bl.beck %>% filter(id %in% rand_id),
            aes(x = doy, y = ndvi.pred, color = 'red')) +
  geom_point(data = s2.bl.beck %>% filter(id %in% rand_id),
             aes(x = ndvi.max.doy, y = ndvi.max), color = "blue") +
  # Smoothed spline
  geom_point(data = s2.bl.smooth %>% filter(id %in% rand_id),
             aes(x = doy, y = ndvi, group = id)) +
  geom_line(data = s2.bl.smooth %>% filter(id %in% rand_id),
            aes(x = doy, y = ndvi.pred, color = 'blue')) +
  geom_point(data = s2.bl.smooth %>% filter(id %in% rand_id),
             aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id, labeller = label_parsed) +  # Parse label
  theme_cowplot() +
  theme(legend.position = "none")  # Remove legend

# Check effect of snow present in pixels on curve fitting ----

# Bring in the snow data
# Snow data
s2.bl.snow <- read.csv('../../data/snow/snow-cover-10m-blaesedalen.csv')

# Join data
s2.bl.joined <- left_join(s2.bl.smooth, s2.bl.snow, by = 'id') %>%
  drop_na(snow.auc) %>% # drop pixel ids for which we have no snow data
  # Create var for latest observed snow cover
  mutate(latest.snow = ifelse(X2023.07.26 != 0,'2023-07-26',
                              ifelse(X2023.07.18 != 0, '2023-07-18',
                                     ifelse(X2023.07.12 != 0, '2023-07-12',
                                            ifelse(X2023.07.02 != 0, '2023-07-02',
                                            '2023-06-01'))))) # arbitrary date


# For one model
plot.data <- s2.bl.joined %>% filter(latest.snow == '2023-06-01') # Which model?

# 9 random pixels in detail
rand_id <- sample(plot.data %>% st_drop_geometry() %>% pull(id) %>% unique(), 9)

ggplot(
  plot.data %>% filter(id %in% rand_id),
  aes(x = doy, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred)) +
  #geom_vline(xintercept = yday(plot.data$latest.snow), aes(color =  scale_color_viridis_c(snow.auc))) +
  # geom_segment(aes(x = yday(latest.snow) - 1, xend = yday(latest.snow) + 1, y = 0, yend = 0.5, color = snow.auc),
  #              size = 1) +
  scale_color_viridis_c() +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id) +
  theme_classic()



