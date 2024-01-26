# Script for calculating average snow cover per pixel
# Calum Hoad, 04/12/2023

# Import necessary libraries ----
library(terra)
library(dplyr)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(pbapply)


# Paths to the data, read in as SpatRast ----
one <- rast('../../data/uav/MAIA-exports/20220629/20220629-div32768_clipped.tif')
two <- rast('../../data/uav/MAIA-exports/20220705/20220705-div32768_clipped.tif') 
three <- rast('../../data/uav/MAIA-exports/20220718/20220718-div32768_clipped.tif')
four <- rast('../../data/uav/MAIA-exports/20220814/20220814-div32768_clipped.tif')

# Paths to blaesedalen data
one <- project(rast('../../data/uav/M3M-exports/5cm/20230702-clipped-5cm-div128.tif'), 'epsg:32621')
two <- project(rast('../../data/uav/M3M-exports/5cm/20230712-clipped-5cm-div128.tif'), 'epsg:32621') 
three <- project(rast('../../data/uav/M3M-exports/5cm/20230718-clipped-5cm-div128.tif'), 'epsg:32621')
four <- project(rast('../../data/uav/M3M-exports/5cm/20230726-clipped-5cm-div128.tif'), 'epsg:32621')

# Resample if not already same res and extent
one <- resample(one, three, method = 'bilinear')
two <- resample(two, three, method = 'bilinear')
four <- resample(four, three, method = 'bilinear')


# Rename bands in raster to specifications of sensor (MAIA Sentinel-2) ----

# MAIA Sentinel-2
names(one) <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3')
names(two) <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3')
names(three) <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3')
names(four) <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3')

# DJI M3M
names(one) <- c('green', 'nir', 're', 'red', 'RGB-R', 'RGB-G', 'RGB-B')
names(two) <- c('green', 'nir', 're', 'red', 'RGB-R', 'RGB-G', 'RGB-B')
names(three) <- c('green', 'nir', 're', 'red', 'RGB-R', 'RGB-G', 'RGB-B')
names(four) <- c('green', 'nir', 're', 'red', 'RGB-R', 'RGB-G', 'RGB-B')

# Create a list of the rasters which can be passed to functions
bl <- list(one, two, three, four)


# Calculate snow cover ----

# Function to calculate NDSI
calc_ndsi <- function(x) {
  x <- (x$green-x$nir)/(x$green+x$nir)
}

# Apply function to raster list using lapply or pblapply
bl <- pblapply(bl, calc_ndsi)

# Create function to classify snow free (< 0) and snow covered (> 0) pixels
class_snow <- function(x) {
  x <- terra::classify(x, rbind(c(-2, 0, 0), c(0, 2, 1)))
}

# Apply function 
bl.snow <- pblapply(bl, class_snow)

# Stack rasters
bl.snow <- rast(bl.snow)

# Rename all rasters in list
names(bl.snow[[1]]) <- ('snow.1')
names(bl.snow[[2]]) <- ('snow.2')
names(bl.snow[[3]]) <- ('snow.3')
names(bl.snow[[4]]) <- ('snow.4')

# Create raster where every pix = 1, to facilitate later pixel count via sum
num.pixels <- terra::classify(bl.snow[[2]], rbind(c(-2, 2, 1)))

# Assign name  
names(num.pixels) <- c('tot.pixels')


# Get grid matching spatial resolution of EO data ----

# Landsat, Bring in the list of pixel centres from LandsatTS
landsat_centres <- st_read('../../data/lsatTS-output/pixel_centres.shp') %>%
  mutate(site = str_split(sample_id, "_") %>% 
           sapply(`[`, 1)) %>%
  filter(site == 'blaesedalen') %>%
  st_transform(crs = 32621)

pixel_poly <- st_buffer(pixel_centres, 15, endCapStyle = "SQUARE")

ggplot() + geom_sf(data = pixel_centres) + geom_sf(data = pixel_poly)

# Sentinel-2, bring in the list of pixel centres from S2 script

# without NDVI filter 
s2.centres <- st_read('../../data/sentinel-2/output/s2_modelled_point.shp')
# with NDVI filter
s2.centres <- st_read('../../data/sentinel-2/output/s2_modelled_point_wide_filterndvi.shp')

s2.poly <- st_buffer(s2.centres, 5, endCapStyle = "SQUARE")

ggplot() + geom_sf(data = s2.poly) + geom_sf(data = s2.centres)  


# Bring in alternative landsat centres from the ls metrics script
ls.centres <- st_read('../../data/sentinel-2/output/ls_modelled_point_wide.shp')

ls.poly <- st_buffer(ls.centres, 15, endCapStyle = "SQUARE")

plot(ls.poly)


# Use pixel polygons to sum the value of the reclassed snow pixels ----

# SENTINEL-2
# Extract number of pixels which are snow covered from classified raster stack
s2.snow.cover <- terra::extract(bl.snow, s2.poly, fun = 'sum', ID = TRUE, bind = TRUE)

# Extract number of pixels per polygon, using num.pixels raster (all pix = 1)
s2.snow.cover <- terra::extract(num.pixels, s2.snow.cover, fun = 'sum', ID = TRUE, 
                              bind = TRUE)
# LANDSAT
# Extract number of pixels which are snow covered from classified raster stack
ls.snow.cover <- terra::extract(bl.snow, ls.poly, fun = 'sum', ID = TRUE, bind = TRUE)

# Extract number of pixels per polygon, using num.pixels raster (all pix = 1)
ls.snow.cover <- terra::extract(num.pixels, ls.snow.cover, fun = 'sum', ID = TRUE, bind = TRUE)


# Calculation of average snow statistic ----
# Calculate average snow coverage per pixel over whole period of obs
extracted_data <- st_as_sf(s2.snow.cover) %>%
  #rename(tot.pixels = pixels) %>%
  mutate(snow.persist = (0.25 * (snow.1/tot.pixels)) +
                        (0.25 * (snow.2/tot.pixels)) +
                        (0.25 * (snow.3/tot.pixels)) +
                        (0.25 * (snow.4/tot.pixels))) %>%
  mutate(snow.1 = (snow.1/tot.pixels)) %>% # And percent cover at each time point
  mutate(snow.2 = (snow.2/tot.pixels)) %>%
  mutate(snow.3 = (snow.3/tot.pixels)) %>%
  mutate(snow.4 = (snow.4/tot.pixels)) %>%
  rename('2023-07-02' = snow.1,
         '2023-07-12' = snow.2,
         '2023-07-18' = snow.3,
         '2023-07-26' = snow.4)

# Output the data to csv
# SENTINEL-2
write.csv2(extracted_data, '../../data/uav/snow-metrics/blaesedalen_10m_snowcover_filterndvi.csv')
# LANDSAT
write.csv2(extracted_data, '../../data/uav/snow-metrics/blaesedalen_30m_snowcover.csv')
write.csv2(extracted_data, '../../data/uav/snow-metrics/blaesedalen_30m_snowcover_andLS.csv')
st_write(extracted_data, '../../data/uav/snow-metrics/blaesedalen_30m_snowcover_andLS.shp')
# Format data for plotting
extracted_data_long <- pivot_longer(extracted_data, 
                                    !sample_id & !site & !tot.pixels, 
                                    names_to = 'date', 
                                    values_to = 'percent.cover') %>%
  mutate(date = ifelse(date == 'ndsi.1', as_date('2023-07-02'),
                            ifelse(date == 'ndsi.2', as_date('2023-07-12'), 
                                   ifelse(date == 'ndsi.3', as_date('2023-07-18'), 
                                          ifelse(date == 'ndsi.4', as_date('2023-07-26'), NaN)))))

ggplot(extracted_data_long, aes(x = as_date(date), y = percent.cover, color = sample_id, group_by(sample_id))) +
  geom_line(show.legend = FALSE)


# Compare Landsat imagery with pixel polygons ----
landsat <- ('../../data/landsat-imagery/kluane_test.tif')

landsat <- rast(landsat)

landsat <- crop(landsat, one)

plotRGB(landsat, 4, 3, 2, scale = 1.1, stretch = 'lin', smooth = FALSE)
points(pixel_centres)
lines(pixel_poly)


single_point <- st_read('../../data/lsatTS-output/pixel_centres.shp') %>%
  mutate(site = str_split(sample_id, "_") %>% 
           sapply(`[`, 1)) %>%
  filter(site == 'blaesedalen')# %>%
 # filter(sample_id == 'blaesedalen_1')

single_poly <- mutate(single_point, geometry =  buffer_square(single_point, length = 30))
plot(single_poly)

options(scipen = 999)


