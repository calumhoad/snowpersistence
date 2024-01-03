# Regression analyses for snow cover, maxNDVI, maxNDVI DoY
# Calum Hoad, 3rd Jan 2024 

library(ggplot2)    # Graphics library
library(sf)         # Spatial data types and handling
library(mapview)    # Visualize spatial data
library(spdep)      # Diagnosing spatial dependence
library(spatialreg) # Spatial lag and spatial error model
library(tidyverse)

# Import the data ----

# Sentinel-2 (30,) snow cover data, derived from UAV imagery
s2.snow <- read.csv2('../../data/uav/snow-metrics/blaesedalen_10m_snowcover.csv') %>%
  rename(ndvi.max = 'ndvi_mx', ndvi.max.doy = 'ndv_mx_') %>%
  mutate(ndvi.max.doy = yday(ndvi.max.doy))

s2.snow$geometry <- lapply(strsplit(s2.snow$geometry, "\\), \\("), function(x) {
  # Remove 'list(' and ')' from the coordinates and split by ','
  coords <- gsub("list\\(|\\)", "", x)
  coords <- strsplit(coords, ", ")
  # Convert coordinates to numeric and combine as matrix
  matrix(as.numeric(unlist(coords)), ncol = 2, byrow = TRUE)
})
  
# Convert the data frame to sf object
data_sf <- st_sf(data, coords = c("V1", "V2"))
