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
library(purrr)
library(broom)
library(viridis)
#library(janitor)

# Part 1: Reads in Sentinel-2 data as a terra rast and clips to aoi extent,
#   calculates NDVI per pixel, then converts raster to point geom.
# Part 2: Reads in Sentinel-2 NDVI timeseries as sf point, then fits curves

#### @Jakob - skip to part 2 ####

###
# PART 1
###

# Create a list of the S2 R10m files for each S2 scene
d20230626 <- list.files('../../data/sentinel-2/20230626/S2B_MSIL2A_20230626T153819_N0509_R011_T21WXT_20230626T183843.SAFE/GRANULE/L2A_T21WXT_A032927_20230626T153935/IMG_DATA/R10m/', full.names = TRUE)
d20230708 <- list.files('../../data/sentinel-2/20230708/S2A_MSIL2A_20230708T152811_N0509_R111_T21WXT_20230708T214952.SAFE/GRANULE/L2A_T21WXT_A042007_20230708T152945/IMG_DATA/R10m/', full.names = TRUE)
d20230726 <- list.files('../../data/sentinel-2/20230726/S2B_MSIL2A_20230726T153819_N0509_R011_T21WXT_20230726T184037.SAFE/GRANULE/L2A_T21WXT_A033356_20230726T153828/IMG_DATA/R10m/', full.names = TRUE)
d20230729 <- list.files('../../data/sentinel-2/20230729/S2B_MSIL2A_20230729T154819_N0509_R054_T21WXT_20230729T181552.SAFE/GRANULE/L2A_T21WXT_A033399_20230729T154940/IMG_DATA/R10m/', full.names = TRUE)
d20230807 <- list.files('../../data/sentinel-2/20230807/S2A_MSIL2A_20230807T152811_N0509_R111_T21WXS_20230807T212701.SAFE/GRANULE/L2A_T21WXS_A042436_20230807T153242/IMG_DATA/R10m/', full.names = TRUE)
d20230808 <- list.files('../../data/sentinel-2/20230808/S2B_MSIL2A_20230808T154819_N0509_R054_T21WXS_20230914T102422.SAFE/GRANULE/L2A_T21WXS_A033542_20230808T154852/IMG_DATA/R10m/', full.names = TRUE)
d20230817 <- list.files('../../data/sentinel-2/20230817/S2A_MSIL2A_20230817T152941_N0509_R111_T21WXT_20230817T214159.SAFE/GRANULE/L2A_T21WXT_A042579_20230817T153311/IMG_DATA/R10m/', full.names = TRUE)
d20230923 <- list.files('../../data/sentinel-2/20230923/S2A_MSIL2A_20230922T154951_N0509_R054_T21WXT_20230922T234400.SAFE/GRANULE/L2A_T21WXT_A043094_20230922T155005/IMG_DATA/R10m/', full.names = TRUE)

# List of imagery dates, for later use
dates <- c('2023-06-26', 
           '2023-07-08',
           '2023-07-26', 
           '2023-07-29', 
           '2023-08-07', 
           '2023-08-08', 
           '2023-08-17', 
           '2023-09-23')

# Get the Sentinel-2 10m bands as raster stack, project + crop to extent of UAV imagery
# As function?

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s2.data <- list(d20230626, 
                d20230708, 
                d20230726, 
                d20230729, 
                d20230807, 
                d20230808, 
                d20230817, 
                d20230923)

# Get uAV imagery over plot to use for cropping - RE-EXPORT UAV IMAGERY SO RE-PROJECT IS AVOIDED
uav <- project(rast('../../data/uav/M3M-exports/5cm/20230702-clipped-5cm-div128.tif'), 'epsg:32621')

# Create function for stacking rasters from lists of filepaths, then cropping to extent of UAV imagery
import_s2 <- function(x) {
    x <- rast(x) %>%
      crop(uav)
    #x <- set.names(object = x, nm = c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp'))
}

# Import the data
s2.data.import <- lapply(s2.data, import_s2)

# Check import function works by plotting rasters from list
plot(s2.data.import[[6]])

# Rename bands to make reading logical
names(s2.data.import[[1]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[2]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[3]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[4]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[5]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[6]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[7]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[8]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')

# Apply function to calculate NDVI ----
s2_ndvi <- function(x) {
  x <- (x$nir-x$red)/(x$nir+x$red)
  # names(x) <- c('ndvi')
}

s2.ndvi <- lapply(s2.data.import, s2_ndvi)

# Stack NDVI rasters into single spatRast
s2.ndvi <- rast(s2.ndvi)

# Name spatRaster layers with dates
names(s2.ndvi) <- dates

# Extract raster time series to points
s2.ndvi.points <- st_as_sf(as.points(s2.ndvi, values = TRUE)) %>%
  mutate(id = row_number())

# Write out the extracted points
st_write(s2.ndvi.points, '../../data/sentinel-2/output/sentinel-2-ndvi-ts-pt-2023.csv', 
         layer_options = "GEOMETRY=AS_XY")


###
# PART 2
###

# Apply parabolic 2nd order polynomial to every pixel in the df ----

# If Part 1 of this script has not been run, read in the data
s2.ndvi.points <- read.csv('../../data/sentinel-2/output/sentinel-2-ndvi-ts-pt-2023.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621)

# Get a dataframe of points from the raster
s2.ndvi.long <- s2.ndvi.points %>%
  pivot_longer(!geometry & !id, names_to = 'doy', values_to = 'ndvi') %>%
  mutate(doy = sub('X', '', doy), 
         doy = sub('\\.', '-', doy), 
         doy = as_date(doy),
         doy = yday(doy)) %>%
  group_by(id)

# Filter out obs with ndvi < 0.1, and groups where there are less than 5 obs
s2.ndvi.long <- s2.ndvi.long %>% 
  # filter(ndvi >= 0.1) %>%
  filter(n_distinct(doy) >= 5)


# Function for fitting parabolic 2nd order polynomial model
model_fit <- function(df) {
  # Using a linear model polynom order 2
  lm(data = df, ndvi ~ poly(doy, 2, raw = T))
  # Using a spline smoother
  #smooth.spline(x = df$doy, y = df$ndvi, spar = 0.5)
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
  vertex <- find_vertex(model)
  # find vertex based on predictions (spline smoother)
  #vertex <- data.frame(
  #  x = pred$x[pred$y == max(pred$y)], 
  #  y = pred$y[pred$y == max(pred$y)]
  #)
  
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
s2.modelled.ndvi <- s2.ndvi.long %>%
  group_modify(~ model_ndvi(.x)) %>%
  mutate(ndvi.max.date = as_date(ndvi.max.doy),
         # Calculate fractional day of year, by adding decimal to date
         ndvi.max.doy = yday(ndvi.max.date) + (ndvi.max.doy - floor(ndvi.max.doy)))

# Quick quality control plots
# 100 random pixels overview
ggplot(
  s2.modelled.ndvi %>% filter(id %in% sample(unique(s2.modelled.ndvi$id), 100)),
  aes(x = doy, colour = id, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred)) +
  theme_classic() +
  theme(legend.position = "none")

# 9 random pixels in detail
rand_id <- sample(s2.modelled.ndvi %>% st_drop_geometry() %>% pull(id) %>% unique(), 9)
ggplot(
  s2.modelled.ndvi %>% filter(id %in% rand_id),
  aes(x = doy, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred)) +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id) +
  theme_classic()

# Outputs ---- 

# Wide format
s2.modelled.export.wide <- s2.modelled.ndvi %>%
  group_by(id) %>%
  filter(row_number() == 1) %>%
  select(!doy & !ndvi & !ndvi.pred)

st_write(st_as_sf(s2.modelled.export.wide),  '../../data/sentinel-2/output/s2_modelled_point_wide_filterndvi.shp')
write.csv2(s2.modelled.export.wide,  '../../data/sentinel-2/output/s2_modelled_point_wide_filterndvi.csv')

# Long format
s2.modelled.export.long <- s2.modelled.ndvi %>%
  rename(doy.obs = 'doy', 
         ndvi.obs = 'ndvi')

st_write(st_as_sf(s2.modelled.export.long),  '../../data/sentinel-2/output/s2_modelled_point_long.shp')
write.csv2(s2.modelled.export.long,  '../../data/sentinel-2/output/s2_modelled_point_long_filteredndvi.csv')
