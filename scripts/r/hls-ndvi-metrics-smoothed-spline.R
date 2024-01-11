# NASA HLS 30m NDVI metrics
# Calum Hoad, 11 Jan 2024

# Calculate NDVI maximum and NDVI maximum DoY by fitting a smoothed spline
# to every pixel in a HLS S30 time series, then extracting
# the vertex.

# Help writing model functions for smoothed spline from Jakob Assmann. 

# Import necessary libraries
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
library(pbapply)
#library(janitor)

# Part 1: Reads in Sentinel-2 data as a terra rast and clips to aoi extent,
#   calculates NDVI per pixel, then converts raster to point geom.
# Part 2: Reads in Sentinel-2 NDVI timeseries as sf point, then fits curves

###
# PART 1
###

# Create a list of the S2 R10m files for each S2 scene
d20230406 <- list.files('../../data/nasa-hls/s30/2023-04-06/', full.names = TRUE)
d20230501 <- list.files('../../data/nasa-hls/s30/2023-05-01/', full.names = TRUE)
d20230516 <- list.files('../../data/nasa-hls/s30/2023-05-16/', full.names = TRUE)
d20230522 <- list.files('../../data/nasa-hls/s30/2023-05-22/', full.names = TRUE)
d20230608 <- list.files('../../data/nasa-hls/s30/2023-06-08/', full.names = TRUE)
d20230626 <- list.files('../../data/nasa-hls/s30/2023-06-26/', full.names = TRUE)
d20230708 <- list.files('../../data/nasa-hls/s30/2023-07-08/', full.names = TRUE)
d20230726 <- list.files('../../data/nasa-hls/s30/2023-07-26/', full.names = TRUE)
d20230729 <- list.files('../../data/nasa-hls/s30/2023-07-29/', full.names = TRUE)
d20230807 <- list.files('../../data/nasa-hls/s30/2023-08-07/', full.names = TRUE)
d20230808 <- list.files('../../data/nasa-hls/s30/2023-08-08/', full.names = TRUE)
d20230914 <- list.files('../../data/nasa-hls/s30/2023-09-14/', full.names = TRUE)
d20230922 <- list.files('../../data/nasa-hls/s30/2023-09-22/', full.names = TRUE)
d20231003 <- list.files('../../data/nasa-hls/s30/2023-10-03/', full.names = TRUE)

d20230406[1:13]
# List of imagery dates, for later use
dates <- c('2023-04-06',
           '2023-05-01', 
           '2023-05-16',
           '2023-05-22', 
           '2023-06-08', 
           '2023-06-26', 
           '2023-07-08', 
           '2023-07-26', 
           '2023-07-29', 
           '2023-08-07', 
           '2023-08-08', 
           '2023-09-14', 
           '2023-09-22', 
           '2023-10-03')

# Get the relevant HLS bands as raster stack, project + crop to extent of UAV imagery

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s30.data <- list(d20230406,
                d20230501,
                d20230516,
                d20230522,
                d20230608,
                d20230626,
                d20230708,
                d20230726,
                d20230729,
                d20230807,
                d20230808,
                d20230914,
                d20230922,
                d20231003)

# Get uAV imagery over plot to use for cropping - RE-EXPORT UAV IMAGERY SO RE-PROJECT IS AVOIDED
uav <- project(rast('../../data/uav/M3M-exports/5cm/20230702-clipped-5cm-div128.tif'), 'epsg:32621')

# Create function for stacking rasters from lists of filepaths, then cropping to extent of UAV imagery
import_s30 <- function(x) {
  x <- rast(x[1:13]) %>%
    crop(uav)
  #x <- set.names(object = x, nm = c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp'))
}

# Import the S30 data using function
s30.data.import <- pblapply(s30.data, import_s30)

# Function to calculate NDVI (B8A = narrow NIR, B4 = red)
s30_ndvi <- function(x) {
  x <- (x$NIR_Narrow-x$Red)/(x$NIR_Narrow+x$Red)
}

# Apply function to stacked S30 raster list
s30.ndvi <- pblapply(s30.data.import, s30_ndvi)

# Stack NDVI rasters into single spatRast
s30.ndvi <- rast(s30.ndvi)

# Name raster layers with dates
names(s30.ndvi) <- dates

# Extract raster time series to points
s30.ndvi.points <- st_as_sf(as.points(s30.ndvi, values = TRUE)) %>%
  mutate(id = row_number())

# Write out the extracted points
st_write(s30.ndvi.points, '../../data/nasa-hls/s30/output/s30-ndvi-ts-pt-2023.csv', 
         layer_options = "GEOMETRY=AS_XY")
