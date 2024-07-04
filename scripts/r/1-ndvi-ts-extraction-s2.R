# Sentinel-2 NDVI Metrics
# Calum Hoad, 25/01/2024

# Extract NDVI time series from Sentinel 2 data, over Kluane
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


###
# SENTINEL-2
###

# Import ***Blaesedalen*** Imagery

# Create a list of the S2 R10m files for each S2 scene
d20230406 <- list.files('../../data/sentinel-2/imagery/20230406/S2A_MSIL2A_20230406T151941_N0509_R068_T21WXS_20230406T214452.SAFE/GRANULE/L2A_T21WXS_A040677_20230406T152237/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230501 <- list.files('../../data/sentinel-2/imagery/20230501/S2B_MSIL2A_20230501T151809_N0509_R068_T21WXS_20230501T173400.SAFE/GRANULE/L2A_T21WXS_A032126_20230501T152117/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230516 <- list.files('../../data/sentinel-2/imagery/20230516/S2A_MSIL2A_20230516T151801_N0509_R068_T21WXS_20230516T215553.SAFE/GRANULE/L2A_T21WXS_A041249_20230516T152031/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230522 <- list.files('../../data/sentinel-2/imagery/20230522/S2A_MSIL2A_20230522T153811_N0509_R011_T21WXS_20230522T231356.SAFE/GRANULE/L2A_T21WXS_A041335_20230522T154223/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230608 <- list.files('../../data/sentinel-2/imagery/20230608/S2A_MSIL2A_20230608T152811_N0509_R111_T21WXS_20230608T215653.SAFE/GRANULE/L2A_T21WXS_A041578_20230608T152919/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230626 <- list.files('../../data/sentinel-2/imagery/20230626/S2B_MSIL2A_20230626T153819_N0509_R011_T21WXT_20230626T183843.SAFE/GRANULE/L2A_T21WXT_A032927_20230626T153935/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230708 <- list.files('../../data/sentinel-2/imagery/20230708/S2A_MSIL2A_20230708T152811_N0509_R111_T21WXT_20230708T214952.SAFE/GRANULE/L2A_T21WXT_A042007_20230708T152945/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230726 <- list.files('../../data/sentinel-2/imagery/20230726/S2B_MSIL2A_20230726T153819_N0509_R011_T21WXT_20230726T184037.SAFE/GRANULE/L2A_T21WXT_A033356_20230726T153828/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230729 <- list.files('../../data/sentinel-2/imagery/20230729/S2B_MSIL2A_20230729T154819_N0509_R054_T21WXT_20230729T181552.SAFE/GRANULE/L2A_T21WXT_A033399_20230729T154940/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230807 <- list.files('../../data/sentinel-2/imagery/20230807/S2A_MSIL2A_20230807T152811_N0509_R111_T21WXS_20230807T212701.SAFE/GRANULE/L2A_T21WXS_A042436_20230807T153242/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230808 <- list.files('../../data/sentinel-2/imagery/20230808/S2B_MSIL2A_20230808T154819_N0509_R054_T21WXS_20230914T102422.SAFE/GRANULE/L2A_T21WXS_A033542_20230808T154852/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230817 <- list.files('../../data/sentinel-2/imagery/20230817/S2A_MSIL2A_20230817T152941_N0509_R111_T21WXT_20230817T214159.SAFE/GRANULE/L2A_T21WXT_A042579_20230817T153311/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20230923 <- list.files('../../data/sentinel-2/imagery/20230923/S2A_MSIL2A_20230922T154951_N0509_R054_T21WXT_20230922T234400.SAFE/GRANULE/L2A_T21WXT_A043094_20230922T155005/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20231003 <- list.files('../../data/sentinel-2/imagery/20231003/S2A_MSIL2A_20231003T152051_N0509_R068_T21WXT_20231003T202504.SAFE/GRANULE/L2A_T21WXT_A043251_20231003T152052/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20231011 <- list.files('../../data/sentinel-2/imagery/20231011/S2B_MSIL2A_20231011T153059_N0509_R111_T21WXS_20231011T190635.SAFE/GRANULE/L2A_T21WXS_A034457_20231011T153057/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20231017 <- list.files('../../data/sentinel-2/imagery/20231017/S2B_MSIL2A_20231017T155149_N0509_R054_T21WXT_20231017T201439.SAFE/GRANULE/L2A_T21WXT_A034543_20231017T155218/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)

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
           '2023-08-17', 
           '2023-09-23',
           '2023-10-03',
           '2023-10-11', 
           '2023-10-17')

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s2.bl.data <- list(d20230406,
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
                d20230817, 
                d20230923,
                d20231003,
                d20231011, 
                d20231017)


# Import ***Kluane*** Imagery

# List the files
d20220404 <- list.files('../../data/sentinel-2/kluane/20220404/S2A_MSIL2A_20220404T203021_N0400_R114_T08VLN_20220404T231548.SAFE/GRANULE/L2A_T08VLN_A035432_20220404T203754/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20220512 <- list.files('../../data/sentinel-2/kluane/20220512/S2B_MSIL2A_20220512T204019_N0400_R014_T08VLN_20220513T130008.SAFE/GRANULE/L2A_T08VLN_A027067_20220512T204030/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20220529 <- list.files('../../data/sentinel-2/kluane/20220529/S2B_MSIL2A_20220529T202849_N0400_R114_T08VLN_20220529T225514.SAFE/GRANULE/L2A_T08VLN_A027310_20220529T203220/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20220608 <- list.files('../../data/sentinel-2/kluane/20220608/S2B_MSIL2A_20220608T202849_N0400_R114_T08VLN_20220608T224737.SAFE/GRANULE/L2A_T08VLN_A027453_20220608T203026/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20220708 <- list.files('../../data/sentinel-2/kluane/20220708/S2B_MSIL2A_20220708T202849_N0400_R114_T08VLN_20220708T225301.SAFE/GRANULE/L2A_T08VLN_A027882_20220708T202849/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20220721 <- list.files('../../data/sentinel-2/kluane/20220721/S2B_MSIL2A_20220721T204029_N0400_R014_T08VLN_20220721T231153.SAFE/GRANULE/L2A_T08VLN_A028068_20220721T204620/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20220726 <- list.files('../../data/sentinel-2/kluane/20220726/S2A_MSIL2A_20220726T204031_N0400_R014_T08VLN_20220727T004510.SAFE/GRANULE/L2A_T08VLN_A037048_20220726T204216/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20220812 <- list.files('../../data/sentinel-2/kluane/20220812/S2A_MSIL2A_20220812T202901_N0400_R114_T08VLN_20220813T031902.SAFE/GRANULE/L2A_T08VLN_A037291_20220812T202856/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20220916 <- list.files('../../data/sentinel-2/kluane/20220916/S2B_MSIL2A_20220916T203109_N0400_R114_T08VLN_20220916T231624.SAFE/GRANULE/L2A_T08VLN_A028883_20220916T203411/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20220924 <- list.files('../../data/sentinel-2/kluane/20220924/S2A_MSIL2A_20220924T204211_N0400_R014_T08VLN_20220925T004557.SAFE/GRANULE/L2A_T08VLN_A037906_20220924T204210/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20220929 <- list.files('../../data/sentinel-2/kluane/20220929/S2B_MSIL2A_20220929T204239_N0400_R014_T08VLN_20220929T232306.SAFE/GRANULE/L2A_T08VLN_A029069_20220929T204451/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)
d20221006 <- list.files('../../data/sentinel-2/kluane/20221006/S2B_MSIL2A_20221006T203319_N0400_R114_T08VLN_20221007T000659.SAFE/GRANULE/L2A_T08VLN_A029169_20221006T203321/IMG_DATA/R10m/', pattern = "\\.jp2$", full.names = TRUE)

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
s2.k.data <- list(d20220404, 
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
import_s2 <- function(x, plot) {
    x <- rast(x) %>%
      crop(plot)
}


# Import the data
s2.bl.import <- pblapply(s2.bl.data, import_s2, plot = bl.uav) # Blaesedalen
s2.kl.import <- pblapply(s2.k.data, import_s2, plot = kl.uav) # Kluane-low
s2.kh.import <- pblapply(s2.k.data, import_s2, plot = kh.uav) # Kluane-high

# Check import function works by plotting rasters from list
plotRGB(s2.bl.import[[7]], r = 4*0.7, g = 3*0.7, b = 2*0.7)

# Set list of band names for Sentinel-2
s2.bands <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')

# Function to rename bands (gives warnings, but seems to work)
rename_bands <- function(x, bands) {
  names(x) <- bands
  update(x, names = TRUE)
}

# Apply function to rename bands
s2.bl.import <- pblapply(s2.bl.import, rename_bands, bands = s2.bands)
s2.kl.import <- pblapply(s2.kl.import, rename_bands, bands = s2.bands)
s2.kh.import <- pblapply(s2.kh.import, rename_bands, bands = s2.bands)


# Create function to calculate S2 NDVI ----
s2_ndvi <- function(x) {
  x <- (x$nir-x$red)/(x$nir+x$red)
  # names(x) <- c('ndvi')
}

# Apply function to calculate S2 NDVI
s2.bl.ndvi <- pblapply(s2.bl.import, s2_ndvi)
s2.kl.ndvi <- pblapply(s2.kl.import, s2_ndvi)
s2.kh.ndvi <- pblapply(s2.kh.import, s2_ndvi)


# Stack NDVI rasters into single spatRast
s2.bl.ndvi <- rast(s2.bl.ndvi)
s2.kl.ndvi <- rast(s2.kl.ndvi)
s2.kh.ndvi <- rast(s2.kh.ndvi)

# Name spatRaster layers with dates
names(s2.bl.ndvi) <- bl.dates
names(s2.kl.ndvi) <- k.dates
names(s2.kh.ndvi) <- k.dates

# Plot to check
ggplot() +
  geom_spatraster(data = s2.bl.ndvi) +
  scale_fill_viridis_c() +
  facet_wrap(~lyr)

# Extract raster time series to points
s2.bl.ndvi.points <- st_as_sf(as.points(s2.bl.ndvi, values = TRUE)) %>%
  mutate(id = row_number())

s2.kl.ndvi.points <- st_as_sf(as.points(s2.kl.ndvi, values = TRUE)) %>%
  mutate(id = row_number())

s2.kh.ndvi.points <- st_as_sf(as.points(s2.kh.ndvi, values = TRUE)) %>%
  mutate(id = row_number())

# Plotting to check consisency of ts
test <- s2.bl.ndvi.points %>% 
  pivot_longer(!id & !geometry, names_to = 'date', values_to = 'ndvi') %>%
  group_by(id)

ggplot() +
  geom_line(data = test, aes(x = as_date(date), y = ndvi, group = id))

# Synthesise 0 values for days late in year, due to assymetry of the datasets ----

# Which dates to synthesise data for?
mid.season = 220 # assume that doy 220 is in approx the middle of the season

# Blaesedalen ###
origin = '2023-01-01'
# Get first and last imagery dates
first.bl <- yday(first(bl.dates))
last.bl <- yday(last(bl.dates))

# Difference between mid-season and first image, in days
first.diff <- mid.season - first.bl
last.diff <- last.bl - mid.season

# Work out what a symmetrical last image day would be, if we assume centre is mid.season
synth.last.bl <- mid.season + first.diff # DOY = 344

# Take the 'symmetrical' last image date as our synthetic last image day
d3 <- as.character(as_date(synth.last.bl, origin = origin)) 

# Get a range of possible S2 imagery dates, between the last image and the 
# symmetrical last day of the season
possible.dates <- seq(last.bl, synth.last.bl, by = 3)

# Sample a further two possible dates for synthetic imagery from this range
# possible.dates <- sample(x = possible.dates, size = 2)
possible.dates <- c(329, 331)

# Convert these doy numbers to real dates
origin <- '2023-01-01'
d1 <- as.character(as_date(possible.dates[[1]], origin = origin))
d2 <- as.character(as_date(possible.dates[[2]], origin = origin))

# Assign 0 NDVI value to the synthetic imagery dates generated above
s2.bl.ndvi.points <- s2.bl.ndvi.points %>%
  mutate(!!d1 := 0, 
         !!d2 := 0, 
         !!d3 := 0)

# Kluane ###
origin = '2022-01-01'
# First and last imagery dates
first.k <- yday(first(k.dates))
last.k <- yday(last(k.dates))  

# Difference between mid-season and first image
k.first.diff <- mid.season - first.k
k.last.diff <- last.k - mid.season

# Work out what a symmetrical last image day would be
synth.last.k <- mid.season + k.first.diff
synth.last.k
# Take the 'symmetrical' last image date as our synthetic last image day
origin <- '2022-01-01'
d3 <- as.character(as_date(synth.last.k, origin = origin)) 

# Get a range of possible S2 imagery dates, between the last image and the 
# symmetrical last day of the season
possible.dates <- seq(last.k, synth.last.k, by = 3)

# Sample a further two possible dates for synthetic imagery from this range
#k.possible.dates <- sample(x = possible.dates, size = 2)
k.possible.dates <- c(315, 324)

# Convert these doy numbers to real dates
d1 <- as.character(as_date(possible.dates[[1]], origin = origin))
d2 <- as.character(as_date(possible.dates[[2]], origin = origin))

# Assign 0 NDVI value to the synthetic imagery dates generated above
s2.kl.ndvi.points <- s2.kl.ndvi.points %>%
  mutate(!!d1 := 0, 
         !!d2 := 0, 
         !!d3 := 0)
s2.kh.ndvi.points <- s2.kh.ndvi.points %>%
  mutate(!!d1 := 0, 
         !!d2 := 0, 
         !!d3 := 0)


# Write out the NDVI time series ----
st_write(s2.bl.ndvi.points, '../../data/ndvi/s2-blaesedalen-ndvi-ts-pt.csv', 
         layer_options = "GEOMETRY=AS_XY")
st_write(s2.kl.ndvi.points, '../../data/ndvi/s2-kluane-low-ndvi-ts-pt.csv', 
         layer_options = "GEOMETRY=AS_XY")
st_write(s2.kh.ndvi.points, '../../data/ndvi/s2-kluane-high-ndvi-ts-pt.csv', 
         layer_options = "GEOMETRY=AS_XY")
