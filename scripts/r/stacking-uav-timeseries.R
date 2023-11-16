# Script to create raster stacks of UAV orthomosaics in Terra
# Calum Hoad, 16/11/2023

# Installs
install.packages('rts')

# Library imports
library(terra)
library(dplyr)
library(rts)

# Paths to the data
one <- ('../../data/uav/MAIA-exports/20220629/20220629-div32768_clipped.tif')
two <- ('../../data/uav/MAIA-exports/20220705/20220705-div32768_clipped.tif') 
three <- ('../../data/uav/MAIA-exports/20220718/20220718-div32768_clipped.tif')
four <- ('../../data/uav/MAIA-exports/20220814/20220814-div32768_clipped.tif')

# Data to raster
rast1 <- rast(one)
rast2 <- rast(two)
rast3 <- rast(three)
rast4 <- rast(four)

# Resample the raster data, to the same spatial resolution
rast1 <- resample(rast1, rast3, method = "bilinear")
rast2 <- resample(rast2, rast3, method = "bilinear")
rast4 <- resample(rast4, rast3, method = "bilinear")

# For each raster, calculate the spectral indices
rast1ndvi <- (rast1$Band_8-rast1$Band_4)/((rast1$Band_8+rast1$Band_4)+.0001)
names(rast1ndvi) <- "ndvi"

rast2ndvi <- (rast2$Band_8-rast2$Band_4)/((rast2$Band_8+rast2$Band_4)+.0001)
names(rast2ndvi) <- "ndvi"

rast3ndvi <- (rast3$Band_8-rast3$Band_4)/((rast3$Band_8+rast3$Band_4)+.0001)
names(rast3ndvi) <- "ndvi"

rast4ndvi <- (rast4$Band_8-rast4$Band_4)/((rast4$Band_8+rast4$Band_4)+.0001)
names(rast4ndvi) <- "ndvi"

rast.list <- list(rast1, rast2, rast3, rast4)



plot(rast1$ndvi)

# Function to calculate NDVI for MAIA/Sentinel-2 data (B8 = NIR, B4 = Red)
#   THIS DOESN'T WORK, IMPROVE UNDERSTANDING OF FUNCS IN R 
terra.ndvi <- function(raster) {
  ndvi_rast <- (raster$Band_8 - raster$Band_4)/((raster$Band_8 + raster$Band_4)+.0001)
  names(ndvi_rast) <- "ndvi"
  raster <- c(raster, ndvi_rast)
  raster
}

lapply(rast.list, terra.ndvi)

# Create list of all rasters
all_rasters <- rast(list(rast1, rast2, rast3, rast4))

# raster to raster sstack
rast_stack <- c(all_rasters)

rast_stack

# Name the rasters in the stack
names(rast_stack) <- c("20220629", "20220705", "20220718", "20220814")

length(rast_stack)

# Try plotting to check the raster stack is working
plotRGB(rast4, r = 4, b = 3, g = 2, stretch = 'lin')

# Attempting to use rts package ----
ts.stack <- c(rast1ndvi, rast2ndvi, rast3ndvi, rast4ndvi)   # Create stack
d <- c("2022-06-29","2022-07-05","2022-07-18","2022-08-14") # List image dates
d <- as.Date(d)                                             # Make dates date format

rt <- rts(ts.stack, d)                                      # Make raster time series

# Plot raster time series
plot(rt)

# Try to find a cell close to the middle of the plot
cellFromRowCol(rt, 3000, 3000)

# Extract time series for a given cell location
t <- rt[20069309]
head(t)
plot(t)

# Try this again
t <- rt[20079309]
head(t)
plot(t)

