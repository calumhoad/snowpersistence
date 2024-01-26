# Check Sentinel-2 data over Kluane
# Calum Hoad, 12 Jan 2024

# Found that imagery from July 18th and November 8th is the issue. July imagery
# contains cloud, November issue has sun angle so low that mountains are casting
# shadows into the plot

# Import necessary libraries
library(terra)
library(dplyr)
library(rts)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(tidyterra)
library(purrr)
library(broom)
library(viridis)
library(pbapply)
#library(janitor)

# Import ***Kluane*** Imagery
d20220404 <- list.files('../../data/sentinel-2/kluane/20220404/S2A_MSIL2A_20220404T203021_N0400_R114_T08VLN_20220404T231548.SAFE/GRANULE/L2A_T08VLN_A035432_20220404T203754/IMG_DATA/R10m/', full.names = TRUE)
d20220512 <- list.files('../../data/sentinel-2/kluane/20220512/S2B_MSIL2A_20220512T204019_N0400_R014_T08VLN_20220513T130008.SAFE/GRANULE/L2A_T08VLN_A027067_20220512T204030/IMG_DATA/R10m/', full.names = TRUE)
d20220529 <- list.files('../../data/sentinel-2/kluane/20220529/S2B_MSIL2A_20220529T202849_N0400_R114_T08VLN_20220529T225514.SAFE/GRANULE/L2A_T08VLN_A027310_20220529T203220/IMG_DATA/R10m/', full.names = TRUE)
d20220608 <- list.files('../../data/sentinel-2/kluane/20220608/S2B_MSIL2A_20220608T202849_N0400_R114_T08VLN_20220608T224737.SAFE/GRANULE/L2A_T08VLN_A027453_20220608T203026/IMG_DATA/R10m/', full.names = TRUE)
d20220708 <- list.files('../../data/sentinel-2/kluane/20220708/S2B_MSIL2A_20220708T202849_N0400_R114_T08VLN_20220708T225301.SAFE/GRANULE/L2A_T08VLN_A027882_20220708T202849/IMG_DATA/R10m/', full.names = TRUE)
d20220718 <- list.files('../../data/sentinel-2/kluane/20220718/S2B_MSIL2A_20220718T202849_N0400_R114_T08VLN_20220718T225258.SAFE/GRANULE/L2A_T08VLN_A028025_20220718T202849/IMG_DATA/R10m/', full.names = TRUE)
d20220721 <- list.files('../../data/sentinel-2/kluane/20220721/S2B_MSIL2A_20220721T204029_N0400_R014_T08VLN_20220721T231153.SAFE/GRANULE/L2A_T08VLN_A028068_20220721T204620/IMG_DATA/R10m/', full.names = TRUE)
d20220726 <- list.files('../../data/sentinel-2/kluane/20220726/S2A_MSIL2A_20220726T204031_N0400_R014_T08VLN_20220727T004510.SAFE/GRANULE/L2A_T08VLN_A037048_20220726T204216/IMG_DATA/R10m/', full.names = TRUE)
d20220812 <- list.files('../../data/sentinel-2/kluane/20220812/S2A_MSIL2A_20220812T202901_N0400_R114_T08VLN_20220813T031902.SAFE/GRANULE/L2A_T08VLN_A037291_20220812T202856/IMG_DATA/R10m/', full.names = TRUE)
d20220916 <- list.files('../../data/sentinel-2/kluane/20220916/S2B_MSIL2A_20220916T203109_N0400_R114_T08VLN_20220916T231624.SAFE/GRANULE/L2A_T08VLN_A028883_20220916T203411/IMG_DATA/R10m/', full.names = TRUE)
d20220924 <- list.files('../../data/sentinel-2/kluane/20220924/S2A_MSIL2A_20220924T204211_N0400_R014_T08VLN_20220925T004557.SAFE/GRANULE/L2A_T08VLN_A037906_20220924T204210/IMG_DATA/R10m/', full.names = TRUE)
d20220929 <- list.files('../../data/sentinel-2/kluane/20220929/S2B_MSIL2A_20220929T204239_N0400_R014_T08VLN_20220929T232306.SAFE/GRANULE/L2A_T08VLN_A029069_20220929T204451/IMG_DATA/R10m/', full.names = TRUE)
d20221006 <- list.files('../../data/sentinel-2/kluane/20221006/S2B_MSIL2A_20221006T203319_N0400_R114_T08VLN_20221007T000659.SAFE/GRANULE/L2A_T08VLN_A029169_20221006T203321/IMG_DATA/R10m/', full.names = TRUE)
d20221108 <- list.files('../../data/sentinel-2/kluane/20221108/S2B_MSIL2A_20221108T204649_N0400_R014_T08VLN_20221108T212814.SAFE/GRANULE/L2A_T08VLN_A029641_20221108T204652/IMG_DATA/R10m/', full.names = TRUE)

# List of imagery dates, for later use
dates <- c('2022-04-04', 
           '2022-05-12', 
           '2022-05-29', 
           '2022-06-08', 
           '2022-07-08', 
           '2022-07-18', 
           '2022-07-21', 
           '2022-07-26', 
           '2022-08-12', 
           '2022-09-16', 
           '2022-09-24', 
           '2022-09-29', 
           '2022-10-06', 
           '2022-11-08')

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s2.data <- list(d20220404, 
                d20220512, 
                d20220529, 
                d20220608, 
                d20220708,
                #d20220718, # Cloudy?
                d20220721, 
                d20220726,
                d20220812, 
                d20220916,
                d20220924,
                d20220929,
                d20221006,
                d20221108) # Unuseable? Clouds + shadow from mountains


# Kluane, low, get uav imagery
uav.low <- project(rast('../../data/uav/MAIA-exports/20220705/20220705-div32768_clipped.tif'), crs = crs('epsg:32608'))

# Kluane, high, get UAV imagery
uav.high <- project(rast('../../data/uav/MAIA-exports/20220729/20220729-div32768.tif'), crs('epsg:32608'))

# Get spatExtent and convert to sf
kl.low <- st_as_sf(vect(ext(uav.low)), crs = crs('epsg:32608'))
kl.high <- st_as_sf(vect(ext(uav.high)), crs = crs('epsg:32608'))

# Draw a polygon around the area encompassing both sites
kl.both <- st_union(kl.high, kl.low) %>%
  st_as_sf(st_bbox(), crs = 'epsg:32608') %>%
  st_buffer(dist = 500, endCapStyle = "SQUARE")


# Create function for stacking rasters from lists of filepaths, then cropping to extent of UAV imagery
import_s2 <- function(x) {
  x <- rast(x) %>%
    crop(kl.both)
  #x <- set.names(object = x, nm = c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp'))
}

# Use function to import the data
s2.data.import <- pblapply(s2.data, import_s2)

# Assign CRS to polygons
kl.high <- st_set_crs(kl.high, 'epsg:32608') 
kl.low <- st_set_crs(kl.low, 'epsg:32608')

# Plot data with plot extents overlain
ggplot() +
  geom_spatraster_rgb(data = s2.data.import[[13]], r = 4, g = 3, b = 2, max_col_value = 15000) +
  geom_sf(data = kl.high, color = 'red', fill = NA) +
  geom_sf(data = kl.low, color = 'blue', fill = NA)

# Indentified cloud covered plots
ggplot() +
  geom_spatraster_rgb(data = s2.data.import[[6]], r = 4, g = 3, b = 2, max_col_value = 6000) +
  geom_sf(data = kl.high, color = 'red', fill = NA) +
  geom_sf(data = kl.low, color = 'blue', fill = NA)

