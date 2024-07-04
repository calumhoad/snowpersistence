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
source('a-ndvi-curve-fitting-functions.R')

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
#s2.bl <- s2.bl %>% filter(doy != yday('2023-10-03'))

# Check the raw time series data ----
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
s30.bl.beck <- s30.bl %>%
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
st_write(st_as_sf(s2.bl.beck),  '../../data/ndvi/s2-bl-beck.csv', 
         layer_options = "GEOMETRY=AS_XY")
# Kluane-low
st_write(st_as_sf(s2.kl.beck),  '../../data/ndvi/s2-kl-beck.csv', 
         layer_options = "GEOMETRY=AS_XY")
# Kluane-high
st_write(st_as_sf(s2.kh.beck),  '../../data/ndvi/s2-kh-beck.csv', 
         layer_options = "GEOMETRY=AS_XY")
# Smoothed-spline ##
# Blaesedalen
st_write(st_as_sf(s2.bl.smooth),  '../../data/ndvi/s2-bl-smooth.csv', 
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

