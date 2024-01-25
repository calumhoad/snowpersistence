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
s2.bl <- read.csv('../../data/sentinel-2/tidy-output/s2-blaesedalen-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621)
# Kluane-low
s2.kl <- read.csv('../../data/sentinel-2/tidy-output/s2-kluane-low-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)
# Kluane-high
s2.kh <- read.csv('../../data/sentinel-2/tidy-output/s2-kluane-high-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)


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

# Apply function
s2.bl <- long_ndvi(s2.bl)
s2.kl <- long_ndvi(s2.kl)
s2.kh <- long_ndvi(s2.kh)


# Fit smoothed spline model to the data ----
s2.bl.smooth <- s2.bl %>%
  group_modify(~model_ndvi_smoothedspline(.x))
s2.kl.smooth <- s2.kl %>%
  group_modify(~model_ndvi_smoothedspline(.x))
s2.kh.smooth <- s2.kh %>%
  group_modify(~model_ndvi_smoothedspline(.x))

# Fit Beck double logistic to the data ----
s2.bl.beck <- s2.bl %>%
  group_modify(~model_ndvi_beck(.x))
s2.kl.beck <- s2.kl %>%
  group_modify(~model_ndvi_beck(.x))
s2.kh.beck <- s2.kh %>%
  group_modify(~model_ndvi_beck(.x))



# Quick quality control plots
# 100 random pixels overview
plot.data <- s2.smooth

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
ggplot(
  plot.data %>% filter(id %in% rand_id),
  aes(x = doy, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred)) +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id) +
  theme_classic()


# Outputs ---- 

# Wide format
s2.modelled.export.wide <- s2.modelled.ndvi %>%
  group_by(id) %>%
  filter(doy == 224) %>%
  dplyr::select(-doy, -ndvi.pred.doy.1, -ndvi.max.date, -ndvi, -ndvi.pred.doy)

st_write(st_as_sf(s2.modelled.export.wide),  '../../data/sentinel-2/tidy-output/s2_kluane-high_smoothed_spline_point_wide.csv', 
         layer_options = "GEOMETRY=AS_XY")

# Long format
s2.modelled.export.long <- s2.modelled.ndvi %>%
  rename(doy.obs = 'doy', 
         ndvi.obs = 'ndvi', 
         ndvi.pred = 'ndvi.pred.doy.1') %>%
  select(-ndvi.max.date)

st_write(st_as_sf(s2.modelled.export.long),  '../../data/sentinel-2/tidy-output/s2_kluane-high_modelled_smoothed_spline_point_long.csv', 
         layer_options = "GEOMETRY=AS_XY")

