# Sentinel-2 NDVI Metrics
# Calum Hoad, 05/12/2023

# Import necessary libraries
library(terra)
library(dplyr)
library(rts)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)

# Path to directory containing S2 images
date.dir <- list.files('../../data/sentinel-2/2023-08-17/', full.names = TRUE)
s2.stacked <- rast(date.dir) %>%
  project('epsg:32621')

plot(s2.stacked)

s2.stacked
# OR

# S2 paths
date <- rast('path')

# S2 paths to list
dates <- c(date1, date2, date3, date4...)

# Crop S2 data to same area as UAV imagery
crop_to <- function(x) {
  terra::crop(x, uavImage)
}

s2.cropped <- lapply(dates, crop_to)


# Apply function to calculate NDVI ----
s2_ndvi <- function(x) {
  ndvi <- (x$green-x$nir3)/(x$green+x$nir3)
  names(rast1ndsi) <- "ndvi"
}

s2.ndvi <- lapply(s2.cropped, s2_ndvi)