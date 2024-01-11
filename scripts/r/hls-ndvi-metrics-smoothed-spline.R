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


###
# PART 2
###

# Apply smoothed spline to every pixel time series ----

# If Part 1 of this script has not been run, read in the data
s30.ndvi.points <- read.csv('../../data/nasa-hls/s30/output/s30-ndvi-ts-pt-2023.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621)

# Get a dataframe of points from the raster
s30.ndvi.long <- s30.ndvi.points %>%
  pivot_longer(!geometry & !id, names_to = 'doy', values_to = 'ndvi') %>%
  mutate(doy = sub('X', '', doy), 
         doy = sub('\\.', '-', doy), 
         doy = as_date(doy),
         doy = yday(doy)) %>%
  group_by(id)

# Filter out obs with ndvi < 0.1, and groups where there are less than 5 obs
s30.ndvi.long <- s30.ndvi.long %>% 
  filter(ndvi >= 0.1) %>%
  filter(n_distinct(doy) >= 5)


# Function for fitting parabolic 2nd order polynomial model
model_fit <- function(df) {
  # Using a spline smoother
  smooth.spline(x = df$doy, y = df$ndvi, spar = 0.5)
}

# Function for calculation of vertex
# https://quantifyinghealth.com/plot-a-quadratic-function-in-r/
find_vertex = function(model_fit) {
  # Get model coefficients
  a = model_fit$coefficients[3]
  b = model_fit$coefficients[2]
  c = model_fit$coefficients[1]
  # Determine x coordinate of vertex
  vertex_x = -b / (2 * a)
  # Calculate y coordinate of vertex
  vertex_y = a * vertex_x**2 + b * vertex_x + c
  # Strip attributes and return as data.frame
  return(data.frame(
    x = as.numeric(vertex_x),
    y = as.numeric(vertex_y)
  ))
}

# Define function to model, find vertex, and precict values
model_ndvi <- function(data) {
  
  ### This doesn't work, have filtered above instead
  # Filter the dataset, to retain only records where NDVI >= 0.1
  # data <- data %>%
  #  filter(ndvi >= 0.1) %>%
  #  filter(n_distinct(doy) >= 5)
  
  # Use function to fit model
  model <- model_fit(data)
  
  # Generate predictions for curve plotting (for time-period doy 130-300)
  pred <- predict(model, data.frame(doy = 130:300))
  # pred <- unlist(pred$y)
  
  # use function to find vertex (linear model)
  # vertex <- find_vertex(model)
  # find vertex based on predictions (spline smoother)
  vertex <- data.frame(
    x = pred$x[pred$y == max(pred$y)], 
    y = pred$y[pred$y == max(pred$y)]
  )
  
  # Write necessary values back to df
  data <- suppressMessages(full_join(data, data.frame(
    doy = 130:300,
    ndvi.max = vertex$y,
    ndvi.max.doy = vertex$x,
    ndvi.pred = pred
  ))) %>% 
    # add missing geometries
    mutate(
      geometry = st_geometry(data[1, ])
    )
  return(data)
}

# Following tutorial here: https://data-se.netlify.app/2018/12/10/new-split-apply-combine-variant-in-dplyr-group-split/

# Apply model_ndvi to data using group_modify
s30.modelled.ndvi <- s30.ndvi.long %>%
  group_modify(~ model_ndvi(.x)) %>%
  mutate(ndvi.max.date = as_date(ndvi.max.doy),
         # Calculate fractional day of year, by adding decimal to date
         ndvi.max.doy = yday(ndvi.max.date) + (ndvi.max.doy - floor(ndvi.max.doy)))

# Quick quality control plots
# 100 random pixels overview
ggplot(
  s30.modelled.ndvi %>% filter(id %in% sample(unique(s30.modelled.ndvi$id), 100)),
  aes(x = doy, colour = id, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred.doy.1)) +
  theme_classic() +
  theme(legend.position = "none")

# 9 random pixels in detail
rand_id <- sample(s30.modelled.ndvi %>% st_drop_geometry() %>% pull(id) %>% unique(), 9)
#rand_id <- c(172, 343, 493, 839, 884, 1026, 1097, 1139, 1165)
ggplot(
  s30.modelled.ndvi %>% filter(id %in% rand_id),
  aes(x = doy, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred.doy.1)) +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id) +
  theme_classic()

# Outputs ---- 

# Wide format
s30.modelled.export.wide <- s30.modelled.ndvi %>%
  group_by(id) %>%
  filter(doy == 220) %>%
  dplyr::select(-doy, -ndvi.pred.doy.1, -ndvi.max.date, -ndvi, -ndvi.pred.doy)

st_write(st_as_sf(s30.modelled.export.wide),  '../../data/nasa-hls/s30/output/s30_modelled_smoothed_spline_point_wide.csv', 
         layer_options = "GEOMETRY=AS_XY")

# Long format
s30.modelled.export.long <- s30.modelled.ndvi %>%
  rename(doy.obs = 'doy', 
         ndvi.obs = 'ndvi', 
         ndvi.pred = 'ndvi.pred.doy.1') %>%
  select(-ndvi.max.date)

st_write(st_as_sf(s30.modelled.export.long),  '../../data/nasa-hls/s30/output/s30_modelled_smoothed_spline_point_long.csv', 
         layer_options = "GEOMETRY=AS_XY")
