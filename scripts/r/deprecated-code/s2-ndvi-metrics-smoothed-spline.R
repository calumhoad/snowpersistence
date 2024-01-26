# Sentinel-2 NDVI Metrics
# Calum Hoad, 05/12/2023

# Calculate NDVI maximum and NDVI maximum DoY by fitting a smoothed spline
# to every pixel in a Sentinel-2 time series, then extracting
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
library(tidyterra)
#library(janitor)

# Part 1: Reads in Sentinel-2 data as a terra rast and clips to aoi extent,
#   calculates NDVI per pixel, then converts raster to point geom.
# Part 2: Reads in Sentinel-2 NDVI timeseries as sf point, then fits curves

###
# PART 1
###

# Import ***Blaesedalen*** Imagery

# Create a list of the S2 R10m files for each S2 scene
d20230406 <- list.files('../../data/sentinel-2/imagery/20230406/S2A_MSIL2A_20230406T151941_N0509_R068_T21WXS_20230406T214452.SAFE/GRANULE/L2A_T21WXS_A040677_20230406T152237/IMG_DATA/R10m/', full.names = TRUE)
d20230501 <- list.files('../../data/sentinel-2/imagery/20230501/S2B_MSIL2A_20230501T151809_N0509_R068_T21WXS_20230501T173400.SAFE/GRANULE/L2A_T21WXS_A032126_20230501T152117/IMG_DATA/R10m/', full.names = TRUE)
d20230516 <- list.files('../../data/sentinel-2/imagery/20230516/S2A_MSIL2A_20230516T151801_N0509_R068_T21WXS_20230516T215553.SAFE/GRANULE/L2A_T21WXS_A041249_20230516T152031/IMG_DATA/R10m/', full.names = TRUE)
d20230522 <- list.files('../../data/sentinel-2/imagery/20230522/S2A_MSIL2A_20230522T153811_N0509_R011_T21WXS_20230522T231356.SAFE/GRANULE/L2A_T21WXS_A041335_20230522T154223/IMG_DATA/R10m/', full.names = TRUE)
d20230525 <- list.files('../../data/sentinel-2/imagery/20230525/S2A_MSIL2A_20230525T154941_N0509_R054_T22WDB_20230525T233156.SAFE/GRANULE/L2A_T22WDB_A041378_20230525T155122/IMG_DATA/R10m/', full.names = TRUE)
d20230608 <- list.files('../../data/sentinel-2/imagery/20230608/S2A_MSIL2A_20230608T152811_N0509_R111_T21WXS_20230608T215653.SAFE/GRANULE/L2A_T21WXS_A041578_20230608T152919/IMG_DATA/R10m/', full.names = TRUE)
d20230626 <- list.files('../../data/sentinel-2/imagery/20230626/S2B_MSIL2A_20230626T153819_N0509_R011_T21WXT_20230626T183843.SAFE/GRANULE/L2A_T21WXT_A032927_20230626T153935/IMG_DATA/R10m/', full.names = TRUE)
d20230708 <- list.files('../../data/sentinel-2/imagery/20230708/S2A_MSIL2A_20230708T152811_N0509_R111_T21WXT_20230708T214952.SAFE/GRANULE/L2A_T21WXT_A042007_20230708T152945/IMG_DATA/R10m/', full.names = TRUE)
d20230726 <- list.files('../../data/sentinel-2/imagery/20230726/S2B_MSIL2A_20230726T153819_N0509_R011_T21WXT_20230726T184037.SAFE/GRANULE/L2A_T21WXT_A033356_20230726T153828/IMG_DATA/R10m/', full.names = TRUE)
d20230729 <- list.files('../../data/sentinel-2/imagery/20230729/S2B_MSIL2A_20230729T154819_N0509_R054_T21WXT_20230729T181552.SAFE/GRANULE/L2A_T21WXT_A033399_20230729T154940/IMG_DATA/R10m/', full.names = TRUE)
d20230807 <- list.files('../../data/sentinel-2/imagery/20230807/S2A_MSIL2A_20230807T152811_N0509_R111_T21WXS_20230807T212701.SAFE/GRANULE/L2A_T21WXS_A042436_20230807T153242/IMG_DATA/R10m/', full.names = TRUE)
d20230808 <- list.files('../../data/sentinel-2/imagery/20230808/S2B_MSIL2A_20230808T154819_N0509_R054_T21WXS_20230914T102422.SAFE/GRANULE/L2A_T21WXS_A033542_20230808T154852/IMG_DATA/R10m/', full.names = TRUE)
d20230817 <- list.files('../../data/sentinel-2/imagery/20230817/S2A_MSIL2A_20230817T152941_N0509_R111_T21WXT_20230817T214159.SAFE/GRANULE/L2A_T21WXT_A042579_20230817T153311/IMG_DATA/R10m/', full.names = TRUE)
d20230923 <- list.files('../../data/sentinel-2/imagery/20230923/S2A_MSIL2A_20230922T154951_N0509_R054_T21WXT_20230922T234400.SAFE/GRANULE/L2A_T21WXT_A043094_20230922T155005/IMG_DATA/R10m/', full.names = TRUE)
d20231003 <- list.files('../../data/sentinel-2/imagery/20231003/S2A_MSIL2A_20231003T152051_N0509_R068_T21WXT_20231003T202504.SAFE/GRANULE/L2A_T21WXT_A043251_20231003T152052/IMG_DATA/R10m/', full.names = TRUE)
d20231011 <- list.files('../../data/sentinel-2/imagery/20231011/S2B_MSIL2A_20231011T153059_N0509_R111_T21WXS_20231011T190635.SAFE/GRANULE/L2A_T21WXS_A034457_20231011T153057/IMG_DATA/R10m/', full.names = TRUE)
d20231017 <- list.files('../../data/sentinel-2/imagery/20231017/S2B_MSIL2A_20231017T155149_N0509_R054_T21WXT_20231017T201439.SAFE/GRANULE/L2A_T21WXT_A034543_20231017T155218/IMG_DATA/R10m/', full.names = TRUE)

# List of imagery dates, for later use
dates <- c('2023-04-06', 
           '2023-05-01', 
           '2023-05-16', 
           '2023-05-22', 
           #'2023-05-25', 
           '2023-06-08', 
           #'2023-06-25', 
           '2023-06-26', 
           '2023-07-08',
           '2023-07-26', 
           '2023-07-29', 
           '2023-08-07', 
           '2023-08-08', 
           '2023-08-17', 
           '2023-09-23',
           '2023-10-03',
           '2023-10-11', 
           '2023-10-17')

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s2.data <- list(d20230406,
                d20230501,
                d20230516,
                d20230522,
                #d20230525, UTM22N
                d20230608,
                #d20230625, UTM22N
                d20230626, 
                d20230708, 
                d20230726, 
                d20230729, 
                d20230807, 
                d20230808, 
                d20230817, 
                d20230923,
                d20231003,
                d20231011, 
                d20231017)

# Import ***Kluane*** Imagery
d20220404 <- list.files('../../data/sentinel-2/kluane/20220404/S2A_MSIL2A_20220404T203021_N0400_R114_T08VLN_20220404T231548.SAFE/GRANULE/L2A_T08VLN_A035432_20220404T203754/IMG_DATA/R10m/', full.names = TRUE)
d20220512 <- list.files('../../data/sentinel-2/kluane/20220512/S2B_MSIL2A_20220512T204019_N0400_R014_T08VLN_20220513T130008.SAFE/GRANULE/L2A_T08VLN_A027067_20220512T204030/IMG_DATA/R10m/', full.names = TRUE)
d20220529 <- list.files('../../data/sentinel-2/kluane/20220529/S2B_MSIL2A_20220529T202849_N0400_R114_T08VLN_20220529T225514.SAFE/GRANULE/L2A_T08VLN_A027310_20220529T203220/IMG_DATA/R10m/', full.names = TRUE)
d20220608 <- list.files('../../data/sentinel-2/kluane/20220608/S2B_MSIL2A_20220608T202849_N0400_R114_T08VLN_20220608T224737.SAFE/GRANULE/L2A_T08VLN_A027453_20220608T203026/IMG_DATA/R10m/', full.names = TRUE)
d20220708 <- list.files('../../data/sentinel-2/kluane/20220708/S2B_MSIL2A_20220708T202849_N0400_R114_T08VLN_20220708T225301.SAFE/GRANULE/L2A_T08VLN_A027882_20220708T202849/IMG_DATA/R10m/', full.names = TRUE)
#d20220718 <- list.files('../../data/sentinel-2/kluane/20220718/S2B_MSIL2A_20220718T202849_N0400_R114_T08VLN_20220718T225258.SAFE/GRANULE/L2A_T08VLN_A028025_20220718T202849/IMG_DATA/R10m/', full.names = TRUE)
d20220721 <- list.files('../../data/sentinel-2/kluane/20220721/S2B_MSIL2A_20220721T204029_N0400_R014_T08VLN_20220721T231153.SAFE/GRANULE/L2A_T08VLN_A028068_20220721T204620/IMG_DATA/R10m/', full.names = TRUE)
d20220726 <- list.files('../../data/sentinel-2/kluane/20220726/S2A_MSIL2A_20220726T204031_N0400_R014_T08VLN_20220727T004510.SAFE/GRANULE/L2A_T08VLN_A037048_20220726T204216/IMG_DATA/R10m/', full.names = TRUE)
d20220812 <- list.files('../../data/sentinel-2/kluane/20220812/S2A_MSIL2A_20220812T202901_N0400_R114_T08VLN_20220813T031902.SAFE/GRANULE/L2A_T08VLN_A037291_20220812T202856/IMG_DATA/R10m/', full.names = TRUE)
d20220916 <- list.files('../../data/sentinel-2/kluane/20220916/S2B_MSIL2A_20220916T203109_N0400_R114_T08VLN_20220916T231624.SAFE/GRANULE/L2A_T08VLN_A028883_20220916T203411/IMG_DATA/R10m/', full.names = TRUE)
d20220924 <- list.files('../../data/sentinel-2/kluane/20220924/S2A_MSIL2A_20220924T204211_N0400_R014_T08VLN_20220925T004557.SAFE/GRANULE/L2A_T08VLN_A037906_20220924T204210/IMG_DATA/R10m/', full.names = TRUE)
d20220929 <- list.files('../../data/sentinel-2/kluane/20220929/S2B_MSIL2A_20220929T204239_N0400_R014_T08VLN_20220929T232306.SAFE/GRANULE/L2A_T08VLN_A029069_20220929T204451/IMG_DATA/R10m/', full.names = TRUE)
d20221006 <- list.files('../../data/sentinel-2/kluane/20221006/S2B_MSIL2A_20221006T203319_N0400_R114_T08VLN_20221007T000659.SAFE/GRANULE/L2A_T08VLN_A029169_20221006T203321/IMG_DATA/R10m/', full.names = TRUE)
#d20221108 <- list.files('../../data/sentinel-2/kluane/20221108/S2B_MSIL2A_20221108T204649_N0400_R014_T08VLN_20221108T212814.SAFE/GRANULE/L2A_T08VLN_A029641_20221108T204652/IMG_DATA/R10m/', full.names = TRUE)

# List of imagery dates, for later use
dates <- c('2022-04-04', 
           '2022-05-12', 
           '2022-05-29', 
           '2022-06-08', 
           '2022-07-08', 
           #'2022-07-18', 
           '2022-07-21', 
           '2022-07-26', 
           '2022-08-12', 
           '2022-09-16', 
           '2022-09-24', 
           '2022-09-29', 
           '2022-10-06') 
           #'2022-11-08'

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s2.data <- list(d20220404, 
                d20220512, 
                d20220529, 
                d20220608, 
                d20220708,
                #d20220718,
                d20220721,
                d20220726,
                d20220812, 
                d20220916,
                d20220924,
                d20220929,
                d20221006)
                #d20221108) # Unuseable? Clouds

# Blaesedalen, Get uAV imagery over plot to use for cropping - RE-EXPORT UAV IMAGERY SO RE-PROJECT IS AVOIDED
uav <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-02-5cm-clipped.tif')

# Kluane, low, get uav imagery
uav <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-07-05-5cm-clipped.tif')

# Kluane, high, get UAV imagery
uav <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-09-5cm-clipped.tif')

# Create function for stacking rasters from lists of filepaths, then cropping to extent of UAV imagery
import_s2 <- function(x) {
    x <- rast(x) %>%
      crop(uav)
    #x <- set.names(object = x, nm = c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp'))
}

plot(uav)

# Import the data
s2.data.import <- pblapply(s2.data, import_s2)

# Check import function works by plotting rasters from list
plotRGB(s2.data.import[[7]], r = 4*0.7, g = 3*0.7, b = 2*0.7)

# Set list of band names for Sentinel-2
s2.bands <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')

# Function to rename bands (gives warnings, but seems to work)
rename_bands <- function(x) {
  names(x) <- s2.bands
  update(x, names = TRUE)
}

# Apply function
s2.data.import <- pblapply(s2.data.import, rename_bands)


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

# Plot to check
ggplot() +
  geom_spatraster(data = s2.ndvi) +
  scale_fill_viridis_c() +
  facet_wrap(~lyr)

# Extract raster time series to points
s2.ndvi.points <- st_as_sf(as.points(s2.ndvi, values = TRUE)) %>%
  mutate(id = row_number())

# Synthesise 0 values for days late in year
# Which dates to synthesise data for?
s2.ndvi.long <- s2.ndvi.points %>%
  pivot_longer(!geometry & !id, names_to = 'doy', values_to = 'ndvi') %>%
  mutate(doy = sub('X', '', doy), 
         doy = sub('\\.', '-', doy), 
         doy = as_date(doy),
         doy = yday(doy)) %>%
  group_by(doy)

summarise_at(s2.ndvi.long, 
             vars(ndvi),
             list(mean.ndvi = mean))

diff <- (220-96)

end <- 220 + diff

possible.dates <- seq(290, 344, by = 3)

possible.dates <- sample(x = possible.dates, size = 3)

possible.dates

origin <- '2023-01-01'

d1 <- as.character(as_date(possible.dates[[1]], origin = origin))
d2 <- as.character(as_date(possible.dates[[2]], origin = origin))
d3 <- as.character(as_date(possible.dates[[3]], origin = origin)) 

s2.ndvi.synthetic <- s2.ndvi.points %>%
  mutate(!!d1 := 0, 
         !!d2 := 0, 
         !!d3 := 0)

s2.ndvi.synth.long <- s2.ndvi.synthetic %>%
  pivot_longer(!geometry & !id, names_to = 'doy', values_to = 'ndvi') %>%
  mutate(#doy = sub('X', '', doy), 
         #doy = sub('\\.', '-', doy), 
         doy = as_date(doy),
         doy = yday(doy)) %>%
  group_by(id)

ggplot() +
  geom_point(data = s2.ndvi.synth.long, aes(x = doy, y = ndvi))
  
# Write out the extracted points
st_write(s2.ndvi.points, '../../data/sentinel-2/tidy-output/s2-kluane-high-ndvi-ts-pt-2023.csv', 
         layer_options = "GEOMETRY=AS_XY")


###
# PART 2
###

# Apply smoothed spline to every pixel time series ----

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
  geom_line(aes(y = ndvi.pred.doy.1)) +
  theme_classic() +
  theme(legend.position = "none")

# 9 random pixels in detail
rand_id <- sample(s2.modelled.ndvi %>% st_drop_geometry() %>% pull(id) %>% unique(), 9)
#rand_id <- c(172, 343, 493, 839, 884, 1026, 1097, 1139, 1165)
ggplot(
  s2.modelled.ndvi %>% filter(id %in% rand_id),
  aes(x = doy, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred.doy.1)) +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id) +
  theme_classic()

# Outputs ---- 

# Wide format
s2.modelled.export.wide <- s2.modelled.ndvi %>%
  group_by(id) %>%
  filter(doy == 224) %>%
  dplyr::select(-doy, -ndvi.pred.doy.1, -ndvi.max.date, -ndvi, -ndvi.pred.doy)

st_write(st_as_sf(s2.modelled.export.wide),  '../../data/sentinel-2/tidy-output/s2_kluane-high_smoothed_spline_point_wide.csv', 
         layer_options = "GEOMETRY=AS_XY")

# Long format
s2.modelled.export.long <- s2.modelled.ndvi %>%
  rename(doy.obs = 'doy', 
         ndvi.obs = 'ndvi', 
         ndvi.pred = 'ndvi.pred.doy.1') %>%
  select(-ndvi.max.date)

st_write(st_as_sf(s2.modelled.export.long),  '../../data/sentinel-2/tidy-output/s2_kluane-high_modelled_smoothed_spline_point_long.csv', 
         layer_options = "GEOMETRY=AS_XY")
