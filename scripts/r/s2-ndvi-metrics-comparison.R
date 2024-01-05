# Script to compare NDVI max and NDVI max DoY between different models
# Calum Hoad, 5th Jan 2024 

library(tidyverse)
library(ggplot2)

# Bring in the data from the various model outputs ----

# Parabolic 2nd order polynomial
s2.para <- read.csv2('../../data/sentinel-2/output/s2_modelled_point_long.csv') %>%
    group_by(id) %>%
    filter(row_number() == 1) %>%
    select(id, geometry, ndvi.max, ndvi.max.doy) %>%
    rename(geometry.p = 'geometry', 
           ndvi.max.p = 'ndvi.max',
           ndvi.max.doy.p = 'ndvi.max.doy') %>%
    mutate(ndvi.max.doy.p = floor(ndvi.max.doy.p))

# Smoothed spline
s2.smooth <- read.csv2('../../data/sentinel-2/output/s2_modelled_JA_point_wide.csv') %>%
  select(id, geometry, ndvi.max, ndvi.max.doy) %>%
  rename(geometry.s = 'geometry', 
         ndvi.max.s = 'ndvi.max', 
         ndvi.max.doy.s = 'ndvi.max.doy')

# Beck double logistic
s2.beck <- read.csv('../../data/sentinel-2/output/s2_modelled_beck_point_wide.csv') %>%
  select(id, X, Y, ndvi.max, ndvi.max.doy) %>%
  rename(geometry.x.b = 'X',
         geometry.y.b = 'Y', 
         ndvi.max.b = 'ndvi.max', 
         ndvi.max.doy.b = 'ndvi.max.doy')


# Join the data from the various model outputs ----
s2.all <- left_join(s2.para, s2.smooth, by = 'id')
s2.all <- left_join(s2.all, s2.beck, by = 'id')

# Manual check of geometry values shows that all ids are consistent, drop geom 
# columns which are redundant
s2.all <- s2.all %>%
  select(-geometry.p, -geometry.s) %>%
  st_as_sf(coords = c("geometry.x.b", "geometry.y.b"), crs = 32621)


# Compare and plot the data ----
