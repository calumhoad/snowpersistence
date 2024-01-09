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

# Rename snow class rasters
names(bl.snow[[1]]) <- ('snow.t1')
names(bl.snow[[2]]) <- ('snow.t2')
names(bl.snow[[3]]) <- ('snow.t3')
names(bl.snow[[4]]) <- ('snow.t4')

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

# Rename snow class rasters
names(bl.r.snow[[1]]) <- ('snow.t1')
names(bl.r.snow[[2]]) <- ('snow.t2')
names(bl.r.snow[[3]]) <- ('snow.t3')
names(bl.r.snow[[4]]) <- ('snow.t4')

# Stack classed rasters
bl.r.snow <- rast(bl.r.snow)

# Plot to check output
plot(bl.r.snow)

# Create raster where every pix = 1, to facilitate later pixel count via sum
num.pixels <- terra::classify(bl.snow[[2]], rbind(c(-2, 2, 1)))
# Assign name to this
names(num.pixels) <- c('tot.pixels')


# Get grid matching spatial resolution of EO data ----


# Landsat, Bring in the list of pixel centres from LandsatTS
landsat_centres <- st_read('../../data/lsatTS-output/pixel_centres.shp') %>%
  mutate(site = str_split(sample_id, "_") %>% 
           sapply(`[`, 1)) %>%
  filter(site == 'blaesedalen') %>%
  st_transform(crs = 32621)

# buffer square to polygon
pixel_poly <- st_buffer(landsat_centres, 15, endCapStyle = "SQUARE")

# plot to check
ggplot() + geom_sf(data = landsat_centres) + geom_sf(data = pixel_poly)


# Sentinel-2, bring in the list of pixel centres from S2 script
s2.centres <- st_read('../../data/sentinel-2/output/s2_modelled_point.shp')

# buffer square to polygon 
s2.poly <- st_buffer(s2.centres, 5, endCapStyle = "SQUARE")

# plot to check
ggplot() + geom_sf(data = s2.poly) + geom_sf(data = s2.centres)  



# Use pixel polygons to extract the sum of the reclassed snow pixels ----

# SENTINEL-2
# Extract number of pixels which are snow covered from classified raster stack
s2.r.snow.cover <- terra::extract(bl.r.snow, s2.poly, fun = 'sum', ID = TRUE, bind = TRUE) # red

# Extract number of pixels per polygon, using num.pixels raster (all pix = 1)
s2.r.snow.cover <- terra::extract(num.pixels, s2.r.snow.cover, fun = 'sum', ID = TRUE, 
                                bind = TRUE)
# LANDSAT
# Extract number of pixels which are snow covered from classified raster stack
ls.snow.cover <- terra::extract(bl.snow, ls.poly, fun = 'sum', ID = TRUE, bind = TRUE)

# Extract number of pixels per polygon, using num.pixels raster (all pix = 1)
ls.snow.cover <- terra::extract(num.pixels, ls.snow.cover, fun = 'sum', ID = TRUE, bind = TRUE)


# Calculate snow cover as percentage ----

# Convert the spatvector to an sf
extracted.data <- st_as_sf(s2.r.snow.cover) %>%
  # Calculate percentage cover per S2 pixel
  mutate(snow.t1 = snow.t1/tot.pixels, 
         snow.t2 = snow.t2/tot.pixels, 
         snow.t3 = snow.t3/tot.pixels, 
         snow.t4 = snow.t4/tot.pixels) %>%
  # Rename time points with actual dates
  rename('2023-07-02' = snow.t1,
         '2023-07-12' = snow.t2,
         '2023-07-18' = snow.t3,
         '2023-07-26' = snow.t4)

# Interpolate between observations and calculate the area under the curve ----

# Get a sample dataset to work with (one sample id)
snow.data <- extracted.data %>% #filter(extracted.data, id == 189) %>%
  select(-ndvi_mx, -ndv_mx_, -tot.pixels) %>%
  pivot_longer(!id & !geometry, names_to = 'date', values_to = "snow.cover") %>%
  drop_na() %>%
  group_by(id)

# Plot the sample data
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
start <- lubridate::yday('2023-07-02')
end <- lubridate::yday('2023-07-26')

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

# Write the output ----

output.data <- pivot_wider(snow.auc, id_cols = c(id, geometry, snow.auc), 
                           names_from = date, 
                           values_from = c(snow.cover))

# There is something going wrong, where some of the 02nd July pixels have NA 
# values for snow cover. Investigating this below:

# Get the pixels which have NA values
test.2 <- output.data %>%
  filter(is.na(`2023-07-02`))

# Plot the pixels missing values as red
ggplot() +
  geom_sf(data = st_as_sf(snow.auc), fill = 'blue') +
  geom_sf(data = st_as_sf(test.2), fill = 'red') 
  
# All the pixels which are missing values for 2nd July are in the NE of the plot.
# The extents should be exactly the same for all of the data, why are some missing?

extracted.data












