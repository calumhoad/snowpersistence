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
t1 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-02-5cm-clipped.tif')
t2 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-12-5cm-clipped.tif')
t3 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-18-5cm-clipped.tif')
t4 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-26-5cm-clipped.tif')

# Dates of imagery
d1 <- '2023-07-02'
d2 <- '2023-07-12'
d3 <- '2023-07-18'
d4 <- '2023-07-26'

# Kluane Low ###
t1 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-06-29-5cm-clipped.tif')
t2 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-07-05-5cm-clipped.tif')
t3 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-07-18-5cm-clipped.tif')
t4 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-08-01-5cm-clipped.tif')
t5 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-08-14-5cm-clipped.tif')

# Dates
d1 <- '2022-06-29'
d2 <- '2022-07-05'
d3 <- '2022-07-18'
d4 <- '2022-08-01'
d5 <- '2022-08-14'

# Kluane High ###
t1 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-09-5cm-clipped.tif')
t2 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-19-5cm-clipped.tif')
t3 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-29-5cm-clipped.tif')
t4 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-08-04-5cm-clipped.tif')
t5 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-08-13-5cm-clipped.tif')

# Dates
d1 <- '2022-07-09'
d2 <- '2022-07-19'
d3 <- '2022-07-29'
d4 <- '2022-08-04'
d5 <- '2022-08-13'

# Assign band names
s2.bands <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3') # MAIA-S2
m3.bands <- c('green', 'nir', 're', 'red', 'RGB-R', 'RGB-G', 'RGB-B') # Mavic M3M bands 

# Choose which sensor
bands <- s2.bands

# Assign band names to spatRaster layers
names(t1) <- bands
names(t2) <- bands
names(t3) <- bands
names(t4) <- bands
names(t5) <- bands

# Create a list of the UAV rasters which can be passed to functions
bl <- list(t1, t2, t3, t4, t5)

# Plot out the first raster in RGB as logic check
ggplot() +
  geom_spatraster_rgb(data = t3, r = 4, g = 3, b = 2, max_col_value = 0.6)


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


## Using only raw red-band reflectance values ##

# Get only the red band
band.filter <- 'red'

t1r <- t1[band.filter]
t2r <- t2[band.filter]
t3r <- t3[band.filter]
t4r <- t4[band.filter]
t5r <- t5[band.filter]

# Red band only rasters to list, for use with function
bl.red <- list(t1r, t2r, t3r, t4r, t5r)

# Plot to check values
plot(rast(bl.red), breaks = c(0, 0.4, 2))

# Function to classify snow free as 0, where red reflectance < 0.4, 
#   and snow covered as 1, where red reflectance > 0.4. 
class_snow_red <- function(x) {
  x <- terra::classify(x, rbind(c(-1, 0.3, 0), c(0.3, 2, 1)))
}

# Run function across list of red-band only rasters
bl.r.snow <- pblapply(bl.red, class_snow_red)

# Name the raster layers with the date they contain data from
names(bl.r.snow[[1]]) <- d1
names(bl.r.snow[[2]]) <- d2
names(bl.r.snow[[3]]) <- d3
names(bl.r.snow[[4]]) <- d4
names(bl.r.snow[[5]]) <- d5


# Use manual gap-filling and masking to tidy the classification ----

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
snow <- pblapply(bl.r.snow, manual_snow_class)

# Stack list of snow rasters into single raster
snow <- rast(snow)

# Checking outputs
plot(snow)
plotRGB(bl.r.snow, r = d1, g = d2, b = d3, scale = 1)

# Create raster where every pix = 1, to facilitate later pixel count via sum
num.pixels <- terra::classify(bl.r.snow[[2]], rbind(c(-2, 2, 1)))
# Assign name to this
names(num.pixels) <- c('tot.pixels')


# Get grid matching spatial resolution of EO data ----

# Note:
#   Need to edit this in order to get Kluane data.

# HLS S30 Bring in the list of pixel centres from LandsatTS
s30.data <- read_csv('../../data/nasa-hls/s30/output/s30_modelled_smoothed_spline_point_wide.csv') %>%
                   st_as_sf(coords = c('X', 'Y'), crs = 32621)
# Kluane-low
s30.data <- read_csv('../../data/nasa-hls/output/s30-kluane-low-ndvi-ts-pt-2023.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)
# Kluane-high
s30.data <- read_csv('../../data/nasa-hls/output/s30-kluane-high-ndvi-ts-pt-2023.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)

# buffer square to polygon
s30.poly <- st_buffer(s30.data, 15, endCapStyle = "SQUARE")

# plot to check
ggplot() + geom_sf(data = s30.poly) + geom_sf(data = s30.data)


# Sentinel-2, bring in the list of pixel centres from S2 script
# Blaesedalen
s2.centres <- read_csv('../../data/ndvi/s2-blaesedalen-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  select(id, geometry)
# Kluane-high
s2.centres <- read_csv('../../data/sentinel-2/tidy-output/s2-kluane-high-ndvi-ts-pt-2023.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  select(id, geometry)
# Kluane-low
s2.centres <- read_csv('../../data/ndvi/s2-kluane-low-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  select(id, geometry)

# buffer square to polygon 
s2.poly <- st_buffer(s2.centres, 5, endCapStyle = "SQUARE")

# plot to check
ggplot() + geom_sf(data = s2.poly) + geom_sf(data = s2.centres)  



# Use pixel polygons to extract the sum of the reclassed snow pixels ----

# SENTINEL-2
# Extract number of pixels which are snow covered from classified raster stack
s2.r.snow.cover <- terra::extract(snow, s2.poly, fun = 'sum', ID = TRUE, bind = TRUE) # red

# Extract number of pixels per polygon, using num.pixels raster (all pix = 1)
s2.r.snow.cover <- terra::extract(num.pixels, s2.r.snow.cover, fun = 'sum', ID = TRUE, 
                                bind = TRUE)
# NASA HLS S30
# Extract number of pixels which are snow covered from classified raster stack
s30.snow.cover <- terra::extract(snow, s30.poly, fun = 'sum', ID = TRUE, bind = TRUE)

# Extract number of pixels per polygon, using num.pixels raster (all pix = 1)
s30.snow.cover <- terra::extract(num.pixels, s30.snow.cover, fun = 'sum', ID = TRUE, bind = TRUE)


# Calculate snow cover as percentage ----
check <- st_as_sf(s2.r.snow.cover)
# Note: Should we just drop t5, as there's never any snow in it and it makes the
#   script much more complex. 

# Convert the spatvector to an sf
extracted.data <- st_as_sf(s2.r.snow.cover) %>%
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
  mutate(snow.av = (0.2*.data[[d1]]) + (0.2*.data[[d2]]) + (0.2*.data[[d3]]) + (0.2*.data[[d4]]) + (0.2*.data[[d5]])) #%>%
  # Rename time points with actual dates
  # rename('2023-07-02' = snow.t1,
  #        '2023-07-12' = snow.t2,
  #        '2023-07-18' = snow.t3,
  #        '2023-07-26' = snow.t4)

# Interpolate between observations and calculate the area under the curve ----

# Format the data (S2)
snow.data <- extracted.data %>% #filter(extracted.data, id == 189) %>%
  select(-tot.pixels) %>%
  pivot_longer(!id & !geometry & !snow.av, names_to = 'date', values_to = "snow.cover") %>%
  drop_na() %>%
  group_by(id)

# Format the data (S30)
snow.data <- extracted.data %>% #filter(extracted.data, id == 189) %>%
  select(-tot.pixels) %>%
  pivot_longer(!id & !geometry & !snow.av, names_to = 'date', values_to = "snow.cover") %>%
  drop_na() %>%
  group_by(id)

# Plot the data
ggplot() +
  geom_point(data = snow.data, aes(x = date, y = snow.cover)) +
  geom_line(data = snow.data, aes(x = date, y = snow.cover, group = id))

# Create predicted values for July 1st and 31st?
# Draw a line between t1 and t2,
# predict t0 (July 1st)
# Draw a line between t3 and t4, 
# predict t5 (July 31st)
# Store points in df?

# Set the time period for which the area under the curve will be calculated
start <- lubridate::yday(d1)
end <- lubridate::yday(d4)

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

# Iterate this over the dataframe
snow.auc <- snow.data %>%
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
output.data <- pivot_wider(snow.auc, id_cols = c(id, geometry, snow.auc, snow.av), 
                           names_from = date, 
                           values_from = c(snow.cover)) %>%
  st_as_sf() %>%
  st_centroid()

# Write the file 
st_write(output.data, "../../data/snow/snow-cover-10m-kluane-low.csv",
         layer_options = "GEOMETRY=AS_XY")


# Troubleshooting ----

# There is something going wrong, where some of the 02nd July pixels have NA 
# values for snow cover. Investigating this below:

# Get the pixels which have NA values
test.2 <- output.data %>%
  filter(is.na(`2023-07-02`))

# Plot the pixels missing values as red
ggplot() +
  geom_sf(data = st_as_sf(snow.auc), fill = 'blue') +
  geom_sf(data = st_as_sf(test.2), fill = 'red') +
  geom_sf(data = output.data)
  
# All the pixels which are missing values for 2nd July are in the NE of the plot.
# The extents should be exactly the same for all of the data, why are some missing?

extracted.data












