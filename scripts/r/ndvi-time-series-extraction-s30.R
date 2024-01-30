# NASA HLS 30m NDVI metrics
# Calum Hoad, 29 Jan 2024

# Extract NDVI time series from NASA HLS S30 data, over Kluane
# High, Kluane Low, and Blaesedalen plot locations.

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
library(tidyterra)
#library(janitor)

###
# NASA HLS S30
###

# Blaesedalen ###

# Create a list of the S30 files for each S30 scene
d20230406 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-04-06/', full.names = TRUE)
d20230501 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-05-01/', full.names = TRUE)
d20230516 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-05-16/', full.names = TRUE)
d20230522 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-05-22/', full.names = TRUE)
d20230608 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-06-08/', full.names = TRUE)
d20230626 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-06-26/', full.names = TRUE)
d20230708 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-07-08/', full.names = TRUE)
d20230726 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-07-26/', full.names = TRUE)
d20230729 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-07-29/', full.names = TRUE)
d20230807 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-08-07/', full.names = TRUE)
d20230808 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-08-08/', full.names = TRUE)
d20230914 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-09-14/', full.names = TRUE)
d20230922 <- list.files('../../data/nasa-hls/blaesedalen/s30/X2023-09-22/', full.names = TRUE)
d20230923 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-09-23/', full.names = TRUE)
d20231003 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-10-03/', full.names = TRUE)

# List of imagery dates, for later use
bl.dates <- c('2023-04-06',
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
           #'2023-09-22',
           '2023-09-23', 
           '2023-10-03')

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s30.bl.data <- list(d20230406,
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
                 #d20230922,
                 d20230923,
                 d20231003)

# Kluane ###

# Create a list of the S30 files for each S30 scene
d20220404 <- list.files('../../data/nasa-hls/kluane/s30/2022-04-04/', full.names = TRUE)
d20220512 <- list.files('../../data/nasa-hls/kluane/s30/2022-05-12/', full.names = TRUE)
d20220529 <- list.files('../../data/nasa-hls/kluane/s30/2022-05-29/', full.names = TRUE)
d20220608 <- list.files('../../data/nasa-hls/kluane/s30/2022-06-08/', full.names = TRUE)
d20220708 <- list.files('../../data/nasa-hls/kluane/s30/2022-07-08/', full.names = TRUE)
d20220721 <- list.files('../../data/nasa-hls/kluane/s30/2022-07-21/', full.names = TRUE)
d20220726 <- list.files('../../data/nasa-hls/kluane/s30/2022-07-26/', full.names = TRUE)
d20220812 <- list.files('../../data/nasa-hls/kluane/s30/2022-08-12/', full.names = TRUE)
d20220916 <- list.files('../../data/nasa-hls/kluane/s30/2022-09-16/', full.names = TRUE)
d20220924 <- list.files('../../data/nasa-hls/kluane/s30/2022-09-24/', full.names = TRUE)
d20220929 <- list.files('../../data/nasa-hls/kluane/s30/2022-09-29/', full.names = TRUE)
d20221006 <- list.files('../../data/nasa-hls/kluane/s30/2022-10-06/', full.names = TRUE)


# List of imagery dates, for later use
k.dates <- c('2022-04-04',
           '2022-05-12',
           '2022-05-29',
           '2022-06-08',
           '2022-07-08',
           '2022-07-21',
           '2022-07-26',
           '2022-08-12',
           '2022-09-16',
           '2022-09-24', 
           '2022-09-29', 
           '2022-10-06')

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s30.k.data <- list(d20220404,
                d20220512,
                d20220529,
                d20220608,
                d20220708,
                d20220721,
                d20220726,
                d20220812,
                d20220916,
                d20220924, 
                d20220929, 
                d20221006)


# Import UAV data for each site, to give area for NDVI time series extraction

# Blaesedalen
bl.uav <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-02-5cm-clipped.tif')

# Kluane, low, get uav imagery
kl.uav <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-07-05-5cm-clipped.tif')

# Kluane, high, get UAV imagery
kh.uav <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-09-5cm-clipped.tif')


# Create function for stacking rasters from lists of filepaths, then cropping to extent of UAV imagery
import_s30 <- function(x, plot) {
  x <- rast(x[1:13]) %>%
    crop(plot)
}

# Import the data
s30.bl.import <- pblapply(s30.bl.data, import_s30, plot = bl.uav) # Blaesedalen
s30.kl.import <- pblapply(s30.k.data, import_s30, plot = kl.uav) # Kluane-low
s30.kh.import <- pblapply(s30.k.data, import_s30, plot = kh.uav) # Kluane-high

# Check import function works by plotting rasters from list
ggplot() +
  geom_spatraster_rgb(data = s30.bl.import[[13]], r = 4, g = 3, b = 2, 
                      max_col_value = 0.2, 
                      interpolate = FALSE)

ggplot() +
  geom_spatraster(data = s30.bl.import[[13]], r = 4, g = 3, b = 2, 
                      max_col_value = 0.2, 
                      interpolate = FALSE) +
  facet_wrap(~lyr)

# What is going on with weird NDVI values?
crs(s30.bl.import[[14]])

d2023
test <- crop(rast(c(d20230922[[4]], d20230922[[13]])), bl.uav)
test$NDVI <- (test$NIR_Narrow-test$Red)/(test$NIR_Narrow+test$Red)
ggplot() +
  geom_spatraster(data = test) +
  scale_color_viridis_c() +
  facet_wrap(~lyr)
# Function to calculate NDVI (B8A = narrow NIR, B4 = red)
s30_ndvi <- function(x) {
  x <- (x$NIR_Narrow-x$Red)/(x$NIR_Narrow+x$Red)
}

# Apply function to calculate S2 NDVI
s30.bl.ndvi <- pblapply(s30.bl.import, s30_ndvi)
s30.kl.ndvi <- pblapply(s30.kl.import, s30_ndvi)
s30.kh.ndvi <- pblapply(s30.kh.import, s30_ndvi)

# Stack NDVI rasters into single spatRast
s30.bl.ndvi <- rast(s30.bl.ndvi)
s30.kl.ndvi <- rast(s30.kl.ndvi)
s30.kh.ndvi <- rast(s30.kh.ndvi)

# Name spatRaster layers with dates
names(s30.bl.ndvi) <- bl.dates
names(s30.kl.ndvi) <- k.dates
names(s30.kh.ndvi) <- k.dates

# Plot to check
ggplot() +
  geom_spatraster(data = s30.bl.ndvi) +
  scale_fill_viridis_c() +
  facet_wrap(~lyr)

plot(s30.bl.ndvi)

# Extract raster time series to points
s30.bl.ndvi.points <- st_as_sf(as.points(s30.bl.ndvi, values = TRUE)) %>%
  mutate(id = row_number())

s30.kl.ndvi.points <- st_as_sf(as.points(s30.kl.ndvi, values = TRUE)) %>%
  mutate(id = row_number())

s30.kh.ndvi.points <- st_as_sf(as.points(s30.kh.ndvi, values = TRUE)) %>%
  mutate(id = row_number())


# Synthesise 0 values for days late in year, due to assymetry of the datasets ----

# Using the same dates as for the S2 data

# Blaesedalen ###
origin <- '2023-01-01'
d1 <- as.character(as_date(290, origin = origin)) # S30 data not available, but S2 data available, and 0 value predominant
d2 <- as.character(as_date(329, origin = origin)) # as per S2 synthetic imagery dates, and d3 and d4 the same
d3 <- as.character(as_date(331, origin = origin))
d4 <- as.character(as_date(344, origin = origin))

# Assign 0 NDVI value to the synthetic imagery dates generated above
s30.bl.ndvi.points <- s30.bl.ndvi.points %>%
  mutate(!!d1 := 0, 
         !!d2 := 0, 
         !!d3 := 0,
         !!d4 := 0)

# Kluane ###
origin <- '2022-01-01'
d1 <- as.character(as_date(315, origin = origin)) 
d2 <- as.character(as_date(324, origin = origin)) 
d3 <- as.character(as_date(346, origin = origin))

# Assign 0 NDVI value to the synthetic imagery dates generated above
s30.kl.ndvi.points <- s30.kl.ndvi.points %>%
  mutate(!!d1 := 0,
         !!d2 := 0, 
         !!d3 := 0)
s30.kh.ndvi.points <- s30.kl.ndvi.points %>%
  mutate(!!d1 := 0,
         !!d2 := 0, 
         !!d3 := 0)

# Checking the data is logical ----
test <- s30.bl.ndvi.points %>% 
  pivot_longer(!id & !geometry, names_to = 'date', values_to = 'ndvi') %>%
  group_by(id)

ggplot() +
  geom_line(data = test, aes(x = as_date(date), y = ndvi, group = id))

ggplot() +
  geom_spatraster_rgb(data = rast(s30.bl.data[[14]]), r = 4, g = 3, b = 2, max_col_value = .3 ) #+
geom_sf(data = s30.bl.ndvi.points, aes(color = 'red', size = 4))

# Outputs ---- 

# Wide format
s30.modelled.export.wide <- s30.modelled.ndvi %>%
  group_by(id) %>%
  filter(doy == 220) %>%
  dplyr::select(-doy, -ndvi.pred.doy.1, -ndvi.max.date, -ndvi, -ndvi.pred.doy)

st_write(st_as_sf(s30.modelled.export.wide),  '../../data/nasa-hls/output/s30_kluane-high-modelled_smoothed_spline_point_wide.csv', 
         layer_options = "GEOMETRY=AS_XY")

# Long format
s30.modelled.export.long <- s30.modelled.ndvi %>%
  rename(doy.obs = 'doy', 
         ndvi.obs = 'ndvi', 
         ndvi.pred = 'ndvi.pred.doy.1') %>%
  select(-ndvi.max.date)

st_write(st_as_sf(s30.modelled.export.long),  '../../data/nasa-hls/output/s30_kluane-high-modelled_smoothed_spline_point_long.csv', 
         layer_options = "GEOMETRY=AS_XY")
