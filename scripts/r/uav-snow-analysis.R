# function for calculating altered NDSI
# Calum Hoad, 04/12/2023

install.packages('rgdal')

# Import necessary libraries
library(terra)
library(dplyr)
library(rts)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)

# Paths to the data, read in as SpatRast
one <- rast('../../data/uav/MAIA-exports/20220629/20220629-div32768_clipped.tif')
two <- rast('../../data/uav/MAIA-exports/20220705/20220705-div32768_clipped.tif') 
three <- rast('../../data/uav/MAIA-exports/20220718/20220718-div32768_clipped.tif')
four <- rast('../../data/uav/MAIA-exports/20220814/20220814-div32768_clipped.tif')

# Paths to blaesedalen data
one <- rast('../../data/uav/M3M-exports/5cm/20230702-clipped-5cm-div128.tif')
two <- rast('../../data/uav/M3M-exports/5cm/20230712-clipped-5cm-div128.tif') 
three <- rast('../../data/uav/M3M-exports/5cm/20230718-clipped-5cm-div128.tif')
four <- rast('../../data/uav/M3M-exports/5cm/20230726-clipped-5cm-div128.tif')

# Resample if not already same res and extent
one <- resample(one, three, method = 'bilinear')
two <- resample(two, three, method = 'bilinear')
four <- resample(four, three, method = 'bilinear')

# Rename bands in raster to specifications of sensor (MAIA Sentinel-2)

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


# Calculate the NDSI
rast1ndsi <- (one$green-one$nir)/(one$green+one$nir)
names(rast1ndsi) <- "ndsi.1"
rast2ndsi <- (two$green-two$nir)/(two$green+two$nir)
names(rast2ndsi) <- "ndsi.2"
rast3ndsi <- (three$green-three$nir)/(three$green+three$nir)
names(rast3ndsi) <- "ndsi.3"
rast4ndsi <- (four$green-four$nir)/(four$green+four$nir)
names(rast4ndsi) <- "ndsi.4"

# Reclassify every pixel as 1, to facilitate pixel count via sum
num_pixels <- terra::classify(rast1ndsi, rbind(c(-1, 1, 1)))
names(num_pixels) <- c('pixels')

# Classify, where < 0.1 not snow, where > 0.1 snow
classNDSIhigher <- terra::classify(rast1ndsi, rbind(c(-1, 0.15, 0), c(0.15, 1, 1)))
classNDSIlower <- terra::classify(rast1ndsi, rbind(c(-1, 0.1, 0), c(0.1, 1, 1)))
classNDSIlowest <- terra::classify(rast1ndsi, rbind(c(-1, 0, 0), c(0, 1, 1)))
plot(classNDSIhigher)
plot(classNDSIlower)
plot(classNDSIlowerst)

class_one <- terra::classify(rast1ndsi, rbind(c(-1, 0, 0), c(0, 1, 1)))
class_two <- terra::classify(rast2ndsi, rbind(c(-1, 0, 0), c(0, 1, 1)))
class_three <- terra::classify(rast3ndsi, rbind(c(-1, 0, 0), c(0, 1, 1)))
class_four <- terra::classify(rast4ndsi, rbind(c(-1, 0, 0), c(0, 1, 1)))

# Plot RGB to check snow classification
plotRGB(one, 4, 3, 2, scale = 0.6, stretch = 'lin')
plot(class_one)

# Get grid matching spatial resolution of EO data ----

# Bring in the list of pixel centres 
pixel_centres <- st_read('../../data/lsatTS-output/pixel_centres.shp') %>%
  mutate(site = str_split(sample_id, "_") %>% 
           sapply(`[`, 1)) %>%
  filter(site == 'blaesedalen') %>%
  st_transform(crs = 32622)

pixel_poly <- st_buffer(pixel_centres, 15, endCapStyle = "SQUARE")

ggplot() + geom_sf(data = pixel_centres) + geom_sf(data = pixel_poly)

 
# Use pixel polygons to sum the value of the reclassed snow pixels ----
extract_sum1 <- terra::extract(class_one, pixel_poly, fun = 'sum', ID = TRUE, 
                              bind = TRUE)
extract_sum2 <- terra::extract(class_two, extract_sum1, fun = 'sum', ID = TRUE, 
                                bind = TRUE)
extract_sum3 <- terra::extract(class_three, extract_sum2, fun = 'sum', ID = TRUE, 
                               bind = TRUE)
extract_sum4 <- terra::extract(class_four, extract_sum3, fun = 'sum', ID = TRUE, 
                               bind = TRUE)

# Extract sum from raster where every pixel = 1, in order to get tot.pix 
extract_sum5 <- terra::extract(num_pixels, extract_sum4, fun = 'sum', ID = TRUE, 
                              bind = TRUE)


# Calculation of average snow statistic ----
# Calculate average snow coverage per pixel over whole period of obs
extracted_data <- as.data.frame(extract_sum5) %>%
  rename(tot.pixels = pixels) %>%
  mutate(snow.persist = (0.25 * (ndsi.1/tot.pixels)) +
                        (0.25 * (ndsi.2/tot.pixels)) +
                        (0.25 * (ndsi.3/tot.pixels)) +
                        (0.25 * (ndsi.4/tot.pixels))) %>%
  mutate(ndsi.1 = (ndsi.1/tot.pixels)) %>% # And percent cover at each time point
  mutate(ndsi.2 = (ndsi.2/tot.pixels)) %>%
  mutate(ndsi.3 = (ndsi.3/tot.pixels)) %>%
  mutate(ndsi.4 = (ndsi.4/tot.pixels))

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

# Question for Jakob - how do we get consistent pixel polys for landsat data?

single_point <- st_read('../../data/lsatTS-output/pixel_centres.shp') %>%
  mutate(site = str_split(sample_id, "_") %>% 
           sapply(`[`, 1)) %>%
  filter(site == 'blaesedalen')# %>%
 # filter(sample_id == 'blaesedalen_1')

single_poly <- mutate(single_point, geometry =  buffer_square(single_point, length = 30))
plot(single_poly)
