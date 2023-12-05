# function for calculating altered NDSI
# Calum Hoad, 04/12/2023

# Install
devtools::install_github("tlhenvironment/buffeRs")

R.Version()

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

one <- resample(one, three, method = 'bilinear')
two <- resample(two, three, method = 'bilinear')
four <- resample(four, three, method = 'bilinear')

# Rename bands in raster to specifications of sensor (MAIA Sentinel-2)
names(one) <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3')
names(two) <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3')
names(three) <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3')
names(four) <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3')

# Calculate the NDSI
rast1ndsi <- (one$green-one$nir3)/(one$green+one$nir3)
names(rast1ndsi) <- "ndsi.1"
rast2ndsi <- (two$green-two$nir3)/(two$green+two$nir3)
names(rast2ndsi) <- "ndsi.2"
rast3ndsi <- (three$green-three$nir3)/(three$green+three$nir3)
names(rast3ndsi) <- "ndsi.3"
rast4ndsi <- (four$green-four$nir3)/(four$green+four$nir3)
names(rast4ndsi) <- "ndsi.4"

# Reclassify every pixel as 1, to facilitate pixel count via sum
num_pixels <- terra::classify(rast1ndsi, rbind(c(-1, 1, 1)))
names(num_pixels) <- 'pixels'
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

# Get grid matching spatial resolution of EO data ----

# Bring in the list of pixel centres 
pixel_centres <- st_read('../../data/lsatTS-output/pixel_centres.shp') %>%
  mutate(site = str_split(sample_id, "_") %>% 
           sapply(`[`, 1)) %>%
  filter(site == 'kluanelow')

# Buffer the pixel centres to a landsat pixel size using square buffer ----
st_buffer(pixel_centres, dist = 15, endCapStyle = 'SQUARE') # needs projection

# Convert to terra vect
pixel_centres <- vect(pixel_centres)

# buffer using terra
pixel_poly <- buffer(pixel_centres, width = 15, quadsegs = 1)

 lines(pixel_poly)
 points(pixel_centres)
 
# Use pixel polygons to sum the value of the reclassed snow pixels ----
extract_sum1 <- terra::extract(class_one, pixel_poly, fun = 'sum', ID = TRUE, 
                              bind = TRUE)
extract_sum2 <- terra::extract(class_two, extract_sum1, fun = 'sum', ID = TRUE, 
                                bind = TRUE)
extract_sum3 <- terra::extract(class_three, extract_sum2, fun = 'sum', ID = TRUE, 
                               bind = TRUE)
extract_sum4 <- terra::extract(class_four, extract_sum3, fun = 'sum', ID = TRUE, 
                               bind = TRUE)



# Repeat the above step for each time point, and add a new variable for each time

# Need to get a sum of the pixels contained within each 
extract_sum5 <- terra::extract(num_pixels, extract_sum4, fun = 'sum', ID = TRUE, 
                              bind = TRUE)


# Calculation of average snow statistic ----
# Do this as spatvector or as a raster? A dataframe. 
extracted_data <- as.data.frame(extract_sum5) %>%
  rename(tot.pixels = green) %>%
  mutate(snow.persist = (0.25 * (ndsi.1/tot.pixels)) +
                        (0.25 * (ndsi.2/tot.pixels)) +
                        (0.25 * (ndsi.3/tot.pixels)) +
                        (0.25 * (ndsi.4/tot.pixels)))
  
  # Plot snow persistence
ggplot(extracted_data, aes(x = sample_id, y = snow.persist)) +
  geom_bar(stat = 'identity')

kl120 <- st_read('../../data/lsatTS-output/pixel_centres.shp') %>%
  mutate(site = str_split(sample_id, "_") %>% 
           sapply(`[`, 1)) %>%
  filter(sample_id == 'kluanelow_120')

plot(class_one)
points(kl120)
lines(pixel_poly)

# Compare Landsat imagery with pixel polygons ----
landsat <- ('../../data/landsat-imagery/kluane_test.tif')

landsat <- rast(landsat)

landsat <- crop(landsat, one)

plotRGB(landsat, 4, 3, 2, scale = 1.1, stretch = 'lin', smooth = FALSE)
points(pixel_centres)
lines(pixel_poly)

# Question for Jakob - how do we get consistent pixel polys for landsat data?
