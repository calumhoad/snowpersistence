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
           ndvi.max.doy.p = 'ndvi.max.doy')

# Smoothed spline
s2.smooth <- read.csv2('../../data/sentinel-2/output/s2_modelled_JA_point_wide.csv')

# Beck double logistic
s2.beck <- read.csv2('../../data/sentinel-2/output/s2_modelled_beck_point_wide.csv')

# Join the data from the various model outputs ----

# Compare and plot the data ----
