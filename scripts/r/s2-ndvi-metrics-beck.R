# Sentinel-2 NDVI Metrics
# Calum Hoad, 05/12/2023

# Turn off scientific notation
options(scipen = 999)

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
library(tidyverse)
library(bfast)
library(phenopix)
library(greenbrown)
library(greenbrown)
library(ggplot2)
library(cowplot)

# Part 1: Reads in Sentinel-2 data as a terra rast and clips to aoi extent,
#   calculates NDVI per pixel, then converts raster to point geom.
# Part 2: Reads in Sentinel-2 NDVI timeseries as sf point, then fits curves

#### @Jakob - skip to part 2 ####

###
# PART 1
###

# Create a list of the S2 R10m files for each S2 scene
d20230406 <- list.files('../../data/sentinel-2/imagery/20230406/S2A_MSIL2A_20230406T151941_N0509_R068_T21WXS_20230406T214452.SAFE/GRANULE/L2A_T21WXS_A040677_20230406T152237/IMG_DATA/R10m/', full.names = TRUE)
d20230501 <- list.files('../../data/sentinel-2/imagery/20230501/S2B_MSIL2A_20230501T151809_N0509_R068_T21WXS_20230501T173400.SAFE/GRANULE/L2A_T21WXS_A032126_20230501T152117/IMG_DATA/R10m/', full.names = TRUE)
d20230516 <- list.files('../../data/sentinel-2/imagery/20230516/S2A_MSIL2A_20230516T151801_N0509_R068_T21WXS_20230516T215553.SAFE/GRANULE/L2A_T21WXS_A041249_20230516T152031/IMG_DATA/R10m/', full.names = TRUE)
d20230522 <- list.files('../../data/sentinel-2/imagery/20230522/S2A_MSIL2A_20230522T153811_N0509_R011_T21WXS_20230522T231356.SAFE/GRANULE/L2A_T21WXS_A041335_20230522T154223/IMG_DATA/R10m/', full.names = TRUE)
d20230525 <- list.files('../../data/sentinel-2/imagery/20230525/S2A_MSIL2A_20230525T154941_N0509_R054_T22WDB_20230525T233156.SAFE/GRANULE/L2A_T22WDB_A041378_20230525T155122/IMG_DATA/R10m/', full.names = TRUE)
d20230608 <- list.files('../../data/sentinel-2/imagery/20230608/S2A_MSIL2A_20230608T152811_N0509_R111_T21WXS_20230608T215653.SAFE/GRANULE/L2A_T21WXS_A041578_20230608T152919/IMG_DATA/R10m/', full.names = TRUE)
d20230625 <- list.files('../../data/sentinel-2/imagery/20230626/S2B_MSIL2A_20230626T153819_N0509_R011_T21WXT_20230626T183843.SAFE/GRANULE/L2A_T21WXT_A032927_20230626T153935/IMG_DATA/R10m/', full.names = TRUE)
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

# Get the Sentinel-2 10m bands as raster stack, project + crop to extent of UAV imagery
# As function?

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
names(s2.data.import[[1]])  <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[2]])  <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[3]])  <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[4]])  <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[5]])  <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[6]])  <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[7]])  <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[8]])  <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[9]])  <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[10]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[11]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[12]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[13]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[14]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[15]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[16]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')

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
st_write(s2.ndvi.points, '../../data/sentinel-2/output/sentinel-2-ndvi-ext-ts-pt-2023.csv', 
         layer_options = "GEOMETRY=AS_XY")


###
# PART 2
###

# Apply parabolic 2nd order polynomial to every pixel in the df ----

# If Part 1 of this script has not been run, read in the data
s2.ndvi.points <- read.csv('../../data/sentinel-2/output/sentinel-2-ndvi-ext-ts-pt-2023.csv') %>%
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
  mutate(ndvi = ifelse(ndvi < 0, 0, ndvi)) %>% # reassign values < 0 as 0, per Beck 
  # filter(ndvi >= 0.1) %>%
  filter(n_distinct(doy) >= 5)

# Shorter test dataset
s2.ndvi.test <- s2.ndvi.long %>% filter(id %in% rand_id)
s2.ndvi.test

#### QUESTION ####
#
# Should have code here to re-assign all negative NDVI values in the time series
# to 0, as per the Beck paper?
#
##################

# In order to get Beck to work, need to get more data days from either end of the 
# growing season. Then need to build the function into the group_map function
# so that it iterates over each pixel id

# Apply Beck, using script provided by Jakob ###
# Fit curve using greenbrown package function using Beck et al. 2006 Eq 3
# and including weighting of "overestimated NDVI values" (see paper for details)
# Note: greenbrown also allows for fitting with Elmore et al. 2011 Eq 4
#       but I tried (!) and with thin samples as we have, the gentle sloping
#       summer phenology introduced by the addtional term can not reliably
#       estimated using the low frequency (gap-y) data that we have. For that
#       use "FitDoubleLogElmore()" with the same syntax. 


# Function for fitting Beck ----
fit_beck  <- function(df) {
  double_log_model <- FitDoubleLogBeck(
    x = df$ndvi,
    t = df$doy,
    weighting = TRUE, # Does this default to TRUE? Re-run model with this explicit.
    #tout = seq(1, 12, length = 365), # Time steps of output, test this?
    plot = TRUE)
}

# Function to generate NDVI predictions using Beck
predict.dbl_log_model <- function(model_object, doys_to_predict) {
  eval(
    model_object$formula,
    c(
      list(t = doys_to_predict),
      split(model_object$params, names(model_object$params))
    )
  )
}

# Predict data for whole year and extract peak ndvi and associated doy
year_in_doys <- 1:365

model_ndvi_beck <- function(data) {
  
  # Fit Beck model to data for pixel
  model <- fit_beck(data)
  
  # Get geometry of id
  geom <- st_geometry(data[1, ])
  
  # Predict values and write back to original dataframe
  data <- suppressMessages(full_join(data, data.frame(
    doy = year_in_doys,
    ndvi_pred = predict.dbl_log_model(model, year_in_doys)
  ))) %>%
    arrange(doy) %>%
    mutate(ndvi_max = max(ndvi_pred)) %>%
    mutate(ndvi_max_doy = doy[which(ndvi_pred == ndvi_max[1])][1]) %>%
    mutate(geometry = geom)
  return(data)
}

# Apply Beck functions
s2.modelled.ndvi.beck <- s2.ndvi.test %>%
  group_modify(~ model_ndvi_beck(.x))

# Quick quality control plots
# 100 random pixels overview
ggplot(
  s2.modelled.ndvi.beck %>% filter(id %in% sample(unique(s2.modelled.ndvi.beck$id), 100)),
  aes(x = doy, colour = id, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi_pred)) +
  theme_classic() +
  theme(legend.position = "none")

# 9 random pixels in detail
rand_id <- sample(s2.modelled.ndvi.beck %>% st_drop_geometry() %>% pull(id) %>% unique(), 9)
ggplot(
  s2.modelled.ndvi.beck %>% filter(id %in% rand_id),
  aes(x = doy, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi_pred)) +
  geom_point(aes(x = ndvi_max_doy, y = ndvi_max), color = "red") +
  facet_wrap(~id) +
  theme_classic()

# Outputs ---- 

# Wide format
s2.modelled.export.wide <- s2.modelled.ndvi.beck %>%
  group_by(id) %>%
  filter(row_number() == 1) %>%
  dplyr::select(!doy & !ndvi & !ndvi_pred) %>%
  rename(ndvi.max = 'ndvi_max', 
         ndvi.max.doy = 'ndvi_max_doy')

st_write(st_as_sf(s2.modelled.export.wide),  
         '../../data/sentinel-2/output/s2_modelled_beck_point_wide.csv',
         layer_options = "GEOMETRY=AS_XY")

write.csv2(s2.modelled.export.wide,  '../../data/sentinel-2/output/s2_modelled_JA_point_wide.csv')

# Long format
s2.modelled.export.long <- s2.modelled.ndvi.beck %>%
  rename(ndvi.pred = 'ndvi_pred',
         ndvi.max = 'ndvi_max', 
         ndvi.max.doy = 'ndvi_max_doy')

st_write(st_as_sf(s2.modelled.export.long),  
         '../../data/sentinel-2/output/s2_modelled_beck_point_long.csv',
         layer_options = "GEOMETRY=AS_XY")

st_write(s2.ndvi.points, '../../data/sentinel-2/output/sentinel-2-ndvi-ext-ts-pt-2023.csv', 
         layer_options = "GEOMETRY=AS_XY")


write.csv2(s2.modelled.export.long,  '../../data/sentinel-2/output/s2_modelled_JA_point_long.csv')
