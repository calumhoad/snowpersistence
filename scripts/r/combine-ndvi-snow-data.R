# Combine snow and ndvi datasets to singlular datasets, read for analyses and plotting
# Calum Hoad, 1 Feb 2024

# Turn off scientific notation
options(scipen = 999)

# Import the necessary libraries ----
library(terra)
library(dplyr)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(tidyterra)
library(pbapply)
library(DescTools)


# Sentinel-2 ----

###
# Blaesedalen
###


# NDVI data
s2.bl.ndvi <- read.csv('../../data/ndvi/s2-bl-smooth.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621)

s2.bl.ndvi.no <- read.csv('../../data/ndvi/s2-bl-smooth-outlier-removed.csv') %>% # without outlier
  st_as_sf(coords = c('X', 'Y'), crs = 32621)

# Snow data
s2.bl.snow <- read.csv('../../data/snow/snow-cover-10m-blaesedalen.csv')

# Join data
s2.bl <- left_join(s2.bl.ndvi, s2.bl.snow, by = 'id') %>%
  drop_na(snow.auc) %>% # drop pixel ids for which we have no snow data
  select(-X, -Y, -ndvi.pred.doy, -X2023.07.02, -X2023.07.12, -X2023.07.18, -X2023.07.26)

s2.bl.no <- left_join(s2.bl.ndvi.no, s2.bl.snow, by = 'id') %>%
  drop_na(snow.auc) %>% # drop pixel ids for which we have no snow data
  select(-X, -Y, -ndvi.pred.doy, -X2023.07.02, -X2023.07.12, -X2023.07.18, -X2023.07.26)

# Output
st_write(st_as_sf(s2.bl), "../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv",
         layer_options = "GEOMETRY=AS_XY")

st_write(st_as_sf(s2.bl.no), "../../data/combined-ndvi-snow/s2-bl-smooth-joined-no-outliers.csv",
         layer_options = "GEOMETRY=AS_XY")

###
# Kluane Low
###

# NDVI data
s2.kl.ndvi <- read.csv('../../data/ndvi/s2-kl-smooth.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)

# Snow data
s2.kl.snow <- read.csv('../../data/snow/snow-cover-10m-kluane-low.csv')

# Join data
s2.kl <- left_join(s2.kl.ndvi, s2.kl.snow, by = 'id') %>%
  drop_na(snow.auc) %>% # drop pixel ids for which we have no snow data
  select(-X, -Y, -ndvi.pred.doy, -X2022.06.29, -X2022.07.05, -X2022.07.18, -X2022.08.01, -X2022.08.14)

# Output
st_write(st_as_sf(s2.kl), "../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv",
         layer_options = "GEOMETRY=AS_XY")

###
# Kluane High
###

# NDVI data
s2.kh.ndvi <- read.csv('../../data/ndvi/s2-kh-smooth.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)

# Snow data
s2.kh.snow <- read.csv("../../data/snow/snow-cover-10m-kluane-high.csv")

# Join data
s2.kh <- left_join(s2.kh.ndvi, s2.kh.snow, by = 'id') %>%
  drop_na(snow.auc) %>%
  select(-X, -Y, -ndvi.pred.doy, -X2022.07.09, -X2022.07.19, -X2022.07.29, -X2022.08.04, -X2022.08.13)

# Output
st_write(st_as_sf(s2.kh), "../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv", 
         layer_options = "GEOMETRY=AS_XY")


# NASA HLS S30 ----

###
# Blaesedalen
###

# NDVI
s30.bl.ndvi <- read.csv('../../data/ndvi/s30-bl-smooth.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621)

# Snow data
s30.bl.snow <- read.csv('../../data/snow/snow-cover-30m-blaesedalen.csv')

# Join data
s30.bl <- left_join(s30.bl.ndvi, s30.bl.snow, by = 'id') %>%
  drop_na(snow.auc) %>%
  select(-X, -Y, -ndvi.pred.doy, -X2023.07.02, -X2023.07.12, -X2023.07.18, -X2023.07.26)

# Output
st_write(st_as_sf(s30.bl), '../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv', 
         layer_options = "GEOMETRY=AS_XY")
