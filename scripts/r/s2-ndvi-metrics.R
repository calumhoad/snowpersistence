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

# Create a list of the S2 R10m files for each S2 scene
d20230626 <- list.files('../../data/sentinel-2/20230626/S2B_MSIL2A_20230626T153819_N0509_R011_T21WXT_20230626T183843.SAFE/GRANULE/L2A_T21WXT_A032927_20230626T153935/IMG_DATA/R10m/', full.names = TRUE)
d20230708 <- list.files('../../data/sentinel-2/20230708/S2A_MSIL2A_20230708T152811_N0509_R111_T21WXT_20230708T214952.SAFE/GRANULE/L2A_T21WXT_A042007_20230708T152945/IMG_DATA/R10m/', full.names = TRUE)
d20230726 <- list.files('../../data/sentinel-2/20230726/S2B_MSIL2A_20230726T153819_N0509_R011_T21WXT_20230726T184037.SAFE/GRANULE/L2A_T21WXT_A033356_20230726T153828/IMG_DATA/R10m/', full.names = TRUE)
d20230729 <- list.files('../../data/sentinel-2/20230729/S2B_MSIL2A_20230729T154819_N0509_R054_T21WXT_20230729T181552.SAFE/GRANULE/L2A_T21WXT_A033399_20230729T154940/IMG_DATA/R10m/', full.names = TRUE)
d20230807 <- list.files('../../data/sentinel-2/20230807/S2A_MSIL2A_20230807T152811_N0509_R111_T21WXS_20230807T212701.SAFE/GRANULE/L2A_T21WXS_A042436_20230807T153242/IMG_DATA/R10m/', full.names = TRUE)
d20230808 <- list.files('../../data/sentinel-2/20230808/S2B_MSIL2A_20230808T154819_N0509_R054_T21WXS_20230914T102422.SAFE/GRANULE/L2A_T21WXS_A033542_20230808T154852/IMG_DATA/R10m/', full.names = TRUE)
d20230817 <- list.files('../../data/sentinel-2/20230817/S2A_MSIL2A_20230817T152941_N0509_R111_T21WXT_20230817T214159.SAFE/GRANULE/L2A_T21WXT_A042579_20230817T153311/IMG_DATA/R10m/', full.names = TRUE)
d20230923 <- list.files('../../data/sentinel-2/20230923/S2A_MSIL2A_20230922T154951_N0509_R054_T21WXT_20230922T234400.SAFE/GRANULE/L2A_T21WXT_A043094_20230922T155005/IMG_DATA/R10m/', full.names = TRUE)

# Get the Sentinel-2 10m bands as raster stack, project + crop to extent of UAV imagery
# As function?

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s2.data <- list(d20230626, d20230708)


# Get uAV imagery over plot to use for cropping
uav <- project(rast('../../data/uav/M3M-exports/5cm/20230702-clipped-5cm-div128.tif'), 'epsg:32621')

d20230726

# Create function for stacking rasters from lists of filepaths, then cropping to extent of UAV imagery
import_s2 <- function(x) {
    x <- rast(x) %>%
      #project('epsg:32621') %>%
      crop(uav)
}

# Try applying the function 
s2.data.import <- lapply(s2.data, import_s2)

plot(s2.data.import[[2]])

# Attempt outside of function
s2.stacked <- rast(d20230726) %>%
  #project('epsg:32621') %>%
  crop(project(rast('../../data/uav/M3M-exports/5cm/20230702-clipped-5cm-div128.tif'), 'epsg:32621'))

plot(s2.stacked)
# Rename bands to make reading logical
names(s2.stacked) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')


# Check the data
plot(s2.stacked)


# Apply function to calculate NDVI ----
s2_ndvi <- function(x) {
  x <- (x$green-x$nir)/(x$green+x$nir)
  # names(x) <- c('ndvi')
}

s2.ndvi <- s2_ndvi(s2.stacked)

# Rename band to 'ndvi'
names(s2.ndvi) <- c('ndvi')

plot(s2.ndvi)

# Convert raster to points
s2.ndvi.points <- st_as_sf(as.points(s2.ndvi, values = TRUE))

# Use lapply to iterate through the dataframe by row (pixel.id) and

