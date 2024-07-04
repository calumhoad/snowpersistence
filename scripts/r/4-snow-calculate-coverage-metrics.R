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


# Bring in and format the data ----

# Blaesedalen ###
blt1 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-02-5cm-clipped.tif')
blt2 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-12-5cm-clipped.tif')
blt3 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-18-5cm-clipped.tif')
blt4 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-26-5cm-clipped.tif')

# Dates of imagery
bld1 <- '2023-07-02'
bld2 <- '2023-07-12'
bld3 <- '2023-07-18'
bld4 <- '2023-07-26'

# Kluane Low ###
klt1 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-06-29-5cm-clipped.tif')
klt2 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-07-05-5cm-clipped.tif')
klt3 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-07-18-5cm-clipped.tif')
klt4 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-08-01-5cm-clipped.tif')
klt5 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-08-14-5cm-clipped.tif')

# Dates
kld1 <- '2022-06-29'
kld2 <- '2022-07-05'
kld3 <- '2022-07-18'
kld4 <- '2022-08-01'
kld5 <- '2022-08-14'

# Kluane High ###
kht1 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-09-5cm-clipped.tif')
kht2 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-19-5cm-clipped.tif')
kht3 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-29-5cm-clipped.tif')
kht4 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-08-04-5cm-clipped.tif')
kht5 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-08-13-5cm-clipped.tif')

# Dates
khd1 <- '2022-07-09'
khd2 <- '2022-07-19'
khd3 <- '2022-07-29'
khd4 <- '2022-08-04'
khd5 <- '2022-08-13'

# Assign band names
s2.bands <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3') # MAIA-S2
m3.bands <- c('green', 'nir', 're', 'red', 'RGB-R', 'RGB-G', 'RGB-B') # Mavic M3M bands 


# Assign band names to spatRaster layers

# Blaesedalen
names(blt1) <- m3.bands
names(blt2) <- m3.bands
names(blt3) <- m3.bands
names(blt4) <- m3.bands


# Kluane low
names(klt1) <- s2.bands
names(klt2) <- s2.bands
names(klt3) <- s2.bands
names(klt4) <- s2.bands
names(klt5) <- s2.bands

# Kluane high
names(kht1) <- s2.bands
names(kht2) <- s2.bands
names(kht3) <- s2.bands
names(kht4) <- s2.bands
names(kht5) <- s2.bands

# Plot out the first raster in RGB as logic check
ggplot() +
  geom_spatraster_rgb(data = blt3, r = 4, g = 3, b = 2, max_col_value = 0.6)

## Using only raw red-band reflectance values ##

# Get only the red band
band.filter <- 'red'

# Blaesedalen
blt1 <- blt1[band.filter]
blt2 <- blt2[band.filter]
blt3 <- blt3[band.filter]
blt4 <- blt4[band.filter]

# Kluane low
klt1 <- klt1[band.filter]
klt2 <- klt2[band.filter]
klt3 <- klt3[band.filter]
klt4 <- klt4[band.filter]
klt5 <- klt5[band.filter]

# Kluane High
kht1 <- kht1[band.filter]
kht2 <- kht2[band.filter]
kht3 <- kht3[band.filter]
kht4 <- kht4[band.filter]
kht5 <- kht5[band.filter]

# Create a list of the UAV rasters which can be passed to functions
bl <- list(blt1, blt2, blt3, blt4)
kl <- list(klt1, klt2, klt3, klt4, klt5)
kh <- list(kht1, kht2, kht3, kht4, kht5)


# Calculate snow cover at each time step ----

## Using NDSI ##

# NOTE:
#   Using NDSI does not work well for either of the UAV data sets in this analysis. 
#   As neither M3M nor MAIA has SWIR bands, the SWIR band is substituted for the 
#   NIR band in NDSI calculation. The result is no-where near as clean an index
#   with confusion around vegetation and other landscape features. Due to this,
#   the classification of choice uses the Red band, code is further below.

# # Function to calculate altered NDSI
# calc_ndsi <- function(x) {
#   #x <- (x$green - x$nir) / (x$green + x$nir) # For M3M
#   x <- (x$green - x$nir3) / (x$green + x$nir3) # For S2
# }
# 
# # Apply calc_ndsi over raster list
# bl.ndsi <- pblapply(bl, calc_ndsi)
# 
# # Plot out the NDSI rasters
# ggplot() +
#   geom_spatraster(data = bl.ndsi[[2]], 
#                   na.rm = TRUE) +
#   scale_fill_viridis_c(limits = c(-0.3, 1))
# 
# plot(rast(bl.ndsi), breaks = c(0, 1))
# 
# # Create function to classify snow free (NDSI < 0) and snow covered (> 0) pixels
# class_snow <- function(x) {
#   x <- terra::classify(x, rbind(c(-2, -0.21, 0), c(-0.21, 2, 1)))
# }
# 
# # Apply class_snow to NDSI rasters
# bl.snow <- pblapply(bl.ndsi, class_snow)
# 
# # Stack rasters
# bl.snow <- rast(bl.snow)
# 
# plot(bl.snow)
# 
# # Rename snow class rasters
# names(bl.snow[[1]]) <- d1
# names(bl.snow[[2]]) <- d2
# names(bl.snow[[3]]) <- d3
# names(bl.snow[[4]]) <- d4
# names(bl.snow[[5]]) <- d5

# Blaesedalen ----
# Plot to check values
plot(rast(bl), breaks = c(0, 0.3, 2))
# Function to classify snow free as 0, where red reflectance < 0.4, 
#   and snow covered as 1, where red reflectance > 0.4. 
bl_class_snow_red <- function(x) {
  x <- terra::classify(x, rbind(c(-1, 0.3, 0), c(0.3, 2, 1)))
}

# Run function across list of red-band only rasters
bl.snow <- pblapply(bl, bl_class_snow_red)

# Name the raster layers with the date they contain data from
names(bl.snow[[1]]) <- bld1
names(bl.snow[[2]]) <- bld2
names(bl.snow[[3]]) <- bld3
names(bl.snow[[4]]) <- bld4

# Kluane ----
# Plot to check values
plot(rast(kl), breaks = c(0, 0.4, 2))

# Function to classify snow free as 0, where red reflectance < 0.4, 
#   and snow covered as 1, where red reflectance > 0.4. 
kl_class_snow_red <- function(x) {
  x <- terra::classify(x, rbind(c(-1, 0.4, 0), c(0.4, 2, 1)))
}

# Run function across list of red-band only rasters at Kluane low
kl.snow <- pblapply(kl, kl_class_snow_red)

# Name the raster layers with the date they contain data from
names(kl.snow[[1]]) <- kld1
names(kl.snow[[2]]) <- kld2
names(kl.snow[[3]]) <- kld3
names(kl.snow[[4]]) <- kld4
names(kl.snow[[5]]) <- kld5

# Run function across list of red-band only rasters at Kluane low
kh.snow <- pblapply(kh, kl_class_snow_red)

# Name the raster layers with the date they contain data from
names(kh.snow[[1]]) <- khd1
names(kh.snow[[2]]) <- khd2
names(kh.snow[[3]]) <- khd3
names(kh.snow[[4]]) <- khd4
names(kh.snow[[5]]) <- khd5


# Use manual gap-filling and masking to tidy the classification at Kluane plots ----

# Convert to rast
#bl.r.snow <- rast(bl.r.snow)
#plot(bl.r.snow)
# Export classification rasters
#writeRaster(bl.r.snow, '../../data/uav/snow-classification/blaesedalen-red-0-3.tif')


# OPEN CLASSIFICATIONS IN GIS SOFTWARE AND CREATE: 
# 1) snow-fill.shp, with polygons where there are false negative snow classifications,
# 2) snow-mask.shp, with polygons masking out areas where there are false positive classifications.

# Import shape files of manual snow fill (filling false negative) and manual 
# snow mask (masking false positive) areas 
snow.fill <- st_read('../../data/uav/snow-classification/kluane-snow-polygons.shp')
snow.mask <- st_read('../../data/uav/snow-classification/kluane-snow-free-polygons.shp')

# Function which uses manually created polygons to reclassify and reduce error
manual_snow_class <- function(x) {
  class.rast <- x
  rast.date <- names(x)# get date from name of rast
  # Create rasters from polygon data  
  rast.snow.fill <- rasterize(snow.fill %>% filter(date == rast.date), x, background = 0)
  rast.snow.mask <- rasterize(snow.mask %>% filter(date == rast.date), x, background = 0)
  # reclassify as 1 if pixel is in snow.fill polygons
  x <- ifel(class.rast == 1 | rast.snow.fill == 1, 1, 0)
  # reclassify as 0 if pixel is in snow.mask polygons, else previous value
  x <- ifel(rast.snow.mask == 1, 0, x)
  # reinstate NA pixels from original raster
  x <- ifel(is.na(class.rast), NA, x)
}

# Run function over classified red raster 
kl.snow <- pblapply(kl.snow, manual_snow_class)
kh.snow <- pblapply(kh.snow, manual_snow_class)

# Stack list of snow rasters into single raster
bl.snow <- rast(bl.snow)
kl.snow <- rast(kl.snow)
kh.snow <- rast(kh.snow)

# Checking outputs
plot(kl.snow)
plotRGB(kl.snow, r = kld1, g = kld2, b = kld3, scale = 1)
plot(kh.snow)
plotRGB(kh.snow, r = khd1, g = khd2, b = khd3, scale = 1)

# Create raster for each plot where every pix = 1, to facilitate later pixel count via sum
bl.num.pixels <- terra::classify(bl.snow[[2]], rbind(c(-2, 2, 1)))
kl.num.pixels <- terra::classify(kl.snow[[2]], rbind(c(-2, 2, 1)))
kh.num.pixels <- terra::classify(kh.snow[[3]], rbind(c(-2, 2, 1)))

# Assign name to these count rasters
names(bl.num.pixels) <- c('tot.pixels')
names(kl.num.pixels) <- c('tot.pixels')
names(kh.num.pixels) <- c('tot.pixels')

# Get grid matching spatial resolution of EO data ----

# HLS S30 Bring in the list of pixel centres from HLSS30
s30.data <- read_csv('../../data/ndvi/s30-blaesedalen-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  select(id, geometry)

# buffer square to polygon
s30.poly <- st_buffer(s30.data, 15, endCapStyle = "SQUARE")

# plot to check
ggplot() + geom_sf(data = s30.poly) + geom_sf(data = s30.data)


# Sentinel-2, bring in the list of pixel centres from S2 script
# Blaesedalen
bl.s2.centres <- read_csv('../../data/ndvi/s2-blaesedalen-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  select(id, geometry)
# Kluane-high
kh.s2.centres <- read_csv('../../data/ndvi/s2-kluane-high-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  select(id, geometry)
# Kluane-low
kl.s2.centres <- read_csv('../../data/ndvi/s2-kluane-low-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  select(id, geometry)

# buffer square to polygon 
bl.s2.poly <- st_buffer(bl.s2.centres, 5, endCapStyle = "SQUARE")
kl.s2.poly <- st_buffer(kl.s2.centres, 5, endCapStyle = "SQUARE")
kh.s2.poly <- st_buffer(kh.s2.centres, 5, endCapStyle = "SQUARE")

# Use pixel polygons to extract the sum of the reclassed snow pixels ----

# SENTINEL-2
# Extract number of pixels which are snow covered from classified raster stack
bl.s2.snow.cover <- terra::extract(bl.snow, bl.s2.poly, fun = 'sum', ID = TRUE, bind = TRUE) 

# Extract number of pixels per polygon, using num.pixels raster (all pix = 1)
bl.s2.snow.cover <- terra::extract(bl.num.pixels, bl.s2.snow.cover, fun = 'sum', ID = TRUE, 
                                bind = TRUE)

# And the same for the Kluane plots...
kl.s2.snow.cover <- terra::extract(kl.snow, kl.s2.poly, fun = 'sum', ID = TRUE, bind = TRUE)
kl.s2.snow.cover <- terra::extract(kl.num.pixels, kl.s2.snow.cover, fun = 'sum', ID = TRUE, bind = TRUE)

kh.s2.snow.cover <- terra::extract(kh.snow, kh.s2.poly, fun = 'sum', ID = TRUE, bind = TRUE)
kh.s2.snow.cover <- terra::extract(kh.num.pixels, kh.s2.snow.cover, fun = 'sum', ID = TRUE, bind = TRUE)

# And for HLSS30 data at Blaesedalen...
# Extract number of pixels which are snow covered from classified raster stack
bl.s30.snow.cover <- terra::extract(bl.snow, s30.poly, fun = 'sum', ID = TRUE, bind = TRUE)

# Extract number of pixels per polygon, using num.pixels raster (all pix = 1)
bl.s30.snow.cover <- terra::extract(bl.num.pixels, bl.s30.snow.cover, fun = 'sum', ID = TRUE, bind = TRUE)


# Calculate snow cover as percentage ----
check <- st_as_sf(bl.s2.snow.cover)

# Functions to convert to sf, calculate percentages
bl.convert.snow <- function(data){
  st_as_sf(data) %>%
    select(id, geometry, tot.pixels, !!d1, !!d2, !!d3, !!d4) %>% #, !!d5) %>%
    # Calculate percentage cover per S2 pixel
    mutate(!!d1 := .data[[d1]]/tot.pixels, 
           !!d2 := .data[[d2]]/tot.pixels, 
           !!d3 := .data[[d3]]/tot.pixels, 
           !!d4 := .data[[d4]]/tot.pixels) %>%#, 
    #!!d5 := .data[[d5]]/tot.pixels) %>%
    # Prevent impossible increases in snowcover between time steps
    mutate(!!d2 := ifelse(.data[[d2]] > .data[[d1]], .data[[d1]], .data[[d2]])) %>%
    mutate(!!d3 := ifelse(.data[[d3]] > .data[[d2]], .data[[d2]], .data[[d3]])) %>% # If snow greater in next step,
    mutate(!!d4 := ifelse(.data[[d4]] > .data[[d3]], .data[[d3]], .data[[d4]])) %>% # take value from last step.
    #mutate(!!d5 := ifelse(.data[[d5]] > .data[[d4]], .data[[d4]], .data[[d5]])) %>%
    # Calculate average cover over the season
    mutate(snow.av = (0.25*.data[[d1]]) + (0.25*.data[[d2]]) + (0.25*.data[[d3]]) + (0.25*.data[[d4]]))# + (0.2*.data[[d5]])) #%>%
  # Rename time points with actual dates
  # rename('2023-07-02' = snow.t1,
  #        '2023-07-12' = snow.t2,
  #        '2023-07-18' = snow.t3,
  #        '2023-07-26' = snow.t4)
}

kl.convert.snow <- function(data){
  st_as_sf(data) %>%
    select(id, geometry, tot.pixels, !!d1, !!d2, !!d3, !!d4, !!d5) %>%
    # Calculate percentage cover per S2 pixel
    mutate(!!d1 := .data[[d1]]/tot.pixels, 
           !!d2 := .data[[d2]]/tot.pixels, 
           !!d3 := .data[[d3]]/tot.pixels, 
           !!d4 := .data[[d4]]/tot.pixels, 
           !!d5 := .data[[d5]]/tot.pixels) %>%
    # Prevent impossible increases in snowcover between time steps
    mutate(!!d2 := ifelse(.data[[d2]] > .data[[d1]], .data[[d1]], .data[[d2]])) %>%
    mutate(!!d3 := ifelse(.data[[d3]] > .data[[d2]], .data[[d2]], .data[[d3]])) %>% # If snow greater in next step,
    mutate(!!d4 := ifelse(.data[[d4]] > .data[[d3]], .data[[d3]], .data[[d4]])) %>% # take value from last step.
    mutate(!!d5 := ifelse(.data[[d5]] > .data[[d4]], .data[[d4]], .data[[d5]])) %>%
    # Calculate average cover over the season
    mutate(snow.av = (0.25*.data[[d1]]) + (0.25*.data[[d2]]) + (0.25*.data[[d3]]) + (0.25*.data[[d4]]) + (0.2*.data[[d5]])) #%>%
  # Rename time points with actual dates
  # rename('2023-07-02' = snow.t1,
  #        '2023-07-12' = snow.t2,
  #        '2023-07-18' = snow.t3,
  #        '2023-07-26' = snow.t4)
}

# Apply function to all datasets
bl.s2.snow.ts <- bl.convert.snow(bl.s2.snow.cover)
kl.s2.snow.ts <- kl.convert.snow(kl.s2.snow.cover)
kh.s2.snow.ts <- kl.convert.snow(kh.s2.snow.cover)
bl.s30.snow.ts <- bl.convert.snow(bl.s30.snow.cover)


# Interpolate between observations and calculate the area under the curve ----

# Format the data
format.snow.data <- function(extracted.data) { #filter(extracted.data, id == 189) %>%
  extracted.data %>%
  select(-tot.pixels) %>%
  pivot_longer(!id & !geometry & !snow.av, names_to = 'date', values_to = "snow.cover") %>%
  drop_na() %>%
  group_by(id)
}

# Apply data formatting
bl.s2.snow.ts <- format.snow.data(bl.s2.snow.ts)
kl.s2.snow.ts <- format.snow.data(kl.s2.snow.ts)
kh.s2.snow.ts <- format.snow.data(kh.s2.snow.ts)
bl.s30.snow.ts <- format.snow.data(bl.s30.snow.ts)

# Plot the data to check
ggplot() +
  geom_point(data = snow.data, aes(x = date, y = snow.cover)) +
  geom_line(data = snow.data, aes(x = date, y = snow.cover, group = id))

# Function for AUC calculation
calc_auc <- function(data) {
  # Calculate area under curve
  auc <- AUC(x = lubridate::yday(data$date), 
             y = data$snow.cover, 
             from = start,
             to = end, 
             method = 'trapezoid')
  # If all snow.cover values are 0 (at no point is there any snow in pixel),
  #   auc will return as NA. Replace all NA values with 0.
  auc <- ifelse(is.na(auc), 0, auc)
  # Write AUC value back to df
  data <- data %>%
    mutate(snow.auc = auc)
}

## BLAESEDALEN, S2
# Set the time period for which the area under the curve will be calculated
start <- lubridate::yday(bld1)
end <- lubridate::yday(bld4)

# Iterate auc calc over the bl dataframe
bl.s2.snow.auc <- bl.s2.snow.ts %>%
  group_modify(~ calc_auc(.x))

## BLAESEDALEN, S30 
# Set the time period for which the area under the curve will be calculated
start <- lubridate::yday(bld1)
end <- lubridate::yday(bld4)

# Iterate auc calc over the bl dataframe
bl.s30.snow.auc <- bl.s30.snow.ts %>%
  group_modify(~ calc_auc(.x))

## KLUANE LOW, S2
# Set the time period for which the area under the curve will be calculated
start <- lubridate::yday(kld1)
end <- lubridate::yday(kld5)

# Iterate auc calc over the bl dataframe
kl.s2.snow.auc <- kl.s2.snow.ts %>%
  group_modify(~ calc_auc(.x))

## KLUANE HIGH, S2
# Set the time period for which the area under the curve will be calculated
start <- lubridate::yday(khd1)
end <- lubridate::yday(khd5)

# Iterate auc calc over the bl dataframe
kh.s2.snow.auc <- kh.s2.snow.ts %>%
  group_modify(~ calc_auc(.x))

# Plot the output ----

# Histogram
ggplot() +
  geom_histogram(data = snow.auc, aes(x = snow.auc), binwidth = 0.5)

# Map of snow.auc values
ggplot() +
  geom_sf(data = st_as_sf(snow.auc), aes(fill = snow.auc)) +
  scale_fill_viridis_c()

# Plot the snow.av against the snow.auc
ggplot() +
  geom_point(data = snow.auc, aes(x = snow.auc, y = snow.av))

# Write the output ----

# Format to wide data for output
format.output.data <- function(data) {
  pivot_wider(data, id_cols = c(id, geometry, snow.auc, snow.av), 
              names_from = date, 
              values_from = c(snow.cover)) %>%
    st_as_sf() %>%
    st_centroid()
}

# Apply to all datasets
bl.s2.snow.auc <- format.output.data(bl.s2.snow.auc)
bl.s30.snow.auc <- format.output.data(bl.s30.snow.auc)
kl.s2.snow.auc <- format.output.data(kl.s2.snow.auc)
kh.s2.snow.auc <- format.output.data(kh.s2.snow.auc)


# Write all output files
st_write(bl.s2.snow.auc, "../../data/snow/snow-cover-10m-blaesedalen.csv",
         layer_options = "GEOMETRY=AS_XY")
st_write(bl.s30.snow.auc, "../../data/snow/snow-cover-30m-blaesedalen.csv",
         layer_options = "GEOMETRY=AS_XY")
st_write(kl.s2.snow.auc, "../../data/snow/snow-cover-10m-kluane-low.csv",
         layer_options = "GEOMETRY=AS_XY")
st_write(kh.s2.snow.auc, "../../data/snow/snow-cover-10m-kluane-high.csv",
         layer_options = "GEOMETRY=AS_XY")












