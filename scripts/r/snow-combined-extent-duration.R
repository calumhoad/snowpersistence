# Script for calculating area-under-curve snow cover
# Calum Hoad, 09 Jan 2024

### Purpose of script
# There are two principle metrics for snow cover:
# 1) Snow cover duration (SCD), defined as the number of days for which a given
#   area remains above a given threshold for snow cover.
# 2) Snow cover extent (SCE), the percentage of a given area covered by snow at
#   a given point in time.
#
# For the analysis of snow's relationship with key NDVI metrics, we need a snow
# metric which encapsulates both the extent of snow within EO pixels and its
# evolution across time (i.e. a pixel with a 10 x 10 cm patch of snow which melts
# out completely on July 26th is not the same as a pixel with a 200 x 200 cm patch
# which melts out on the dame date)
#
# In order to calculate a more meaningful snow metric for this study, this script
# will calculate snow persistence as the area-under-the-curve of snow cover extent
# (y) across time (x).
###

# Import the necessary libraries ----
library(terra)
library(dplyr)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(tidyterra)
library(pbapply)


# Bring in and format the data ----

# Paths to Blaesedalen data, read in as SpatRast
t1 <- project(rast('../../data/uav/M3M-exports/5cm/20230702-clipped-5cm-div128.tif'), 'epsg:32621')
t2 <- project(rast('../../data/uav/M3M-exports/5cm/20230712-clipped-5cm-div128.tif'), 'epsg:32621') 
t3 <- project(rast('../../data/uav/M3M-exports/5cm/20230718-clipped-5cm-div128.tif'), 'epsg:32621')
t4 <- project(rast('../../data/uav/M3M-exports/5cm/20230726-clipped-5cm-div128.tif'), 'epsg:32621')

# Assign band names
s2.bands <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3') # MAIA-S2
m3.bands <- c('green', 'nir', 're', 'red', 'RGB-R', 'RGB-G', 'RGB-B') # Mavic M3M bands 

names(t1) <- m3.bands
names(t2) <- m3.bands
names(t3) <- m3.bands
names(t4) <- m3.bands

# Create a list of the UAV rasters which can be passed to functions
bl <- list(t1, t2, t3, t4)

# Plot out the first raster in RGB as logic check
ggplot() +
  geom_spatraster_rgb(data = t4, r = 5, g = 6, b = 7, max_col_value = 0.6)

# Calculate snow cover at each time step ----

# Function to calculate altered NDSI
calc_ndsi <- function(x) {
  x <- (x$green - x$nir)
}

# Apply calc_ndsi over raster list
bl.ndsi <- pblapply(bl, calc_ndsi)

# Plot out the NDSI rasters
ggplot() +
  geom_spatraster(data = bl.ndsi[[2]], 
                  na.rm = TRUE) +
  scale_fill_viridis_c(limits = c(-0.3, 1))

# Create function to classify snow free (NDSI < 0) and snow covered (> 0) pixels
class_snow <- function(x) {
  x <- terra::classify(x, rbind(c(-2, 0, 0), c(0, 2, 1)))
}

# Apply class_snow to NDSI rasters
bl.snow <- pblapply(bl.ndsi, class_snow)

# Stack rasters
bl.snow <- rast(bl.snow)

plot(bl.snow)

# Using only raw red-band reflectance values

# Get only the red band
t1r <- t1['red']
t2r <- t2['red']
t3r <- t3['red']
t4r <- t4['red']

# Red band only rasters to list, for use with function
bl.red <- list(t1r, t2r, t3r, t4r)

# Plot to check values
plot(rast(bl.red))

# Function to classify snow free as 0, where red reflectance < 0.4, 
#   and snow covered as 1, where red reflectance > 0.4. 
class_snow_red <- function(x) {
  x <- terra::classify(x, rbind(c(-1, 0.4, 0), c(0.4, 2, 1)))
}

# Run function across list of red-band only rasters
bl.r.snow <- pblapply(bl.red, class_snow_red)

# Plot to check output
plot(rast(bl.r.snow))

# Rename snow class rasters
names(bl.snow[[1]]) <- ('snow.t1')
names(bl.snow[[2]]) <- ('snow.t2')
names(bl.snow[[3]]) <- ('snow.t3')
names(bl.snow[[4]]) <- ('snow.t4')

# Create raster where every pix = 1, to facilitate later pixel count via sum
num.pixels <- terra::classify(bl.snow[[2]], rbind(c(-2, 2, 1)))
# Assign name to this
names(num.pixels) <- c('tot.pixels')


# Get grid matching spatial resolution of EO data ----

