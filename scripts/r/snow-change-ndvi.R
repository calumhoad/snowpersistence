# Script for calculating area-under-curve snow cover
# Calum Hoad, 09 Jan 2024

# Turn off scientific notation
options(scipen = 999)

# Import the necessary libraries ----
library(terra)
library(dplyr)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(tidyterra)
library(pbapply)
library(DescTools)
library(lubridate)
library(cowplot)

# Bring in and format the data ----

# Load the snow points shp
snow.points <- vect('../../data/snow/points/snow-points.shp')

# Now do the NDVI side

# Blaesedalen ###
blt1 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-02-5cm-clipped.tif')
blt2 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-12-5cm-clipped.tif')
blt3 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-18-5cm-clipped.tif')
blt4 <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-26-5cm-clipped.tif')

# Dates of imagery
bld1 <- '2023-07-02'
bld2 <- '2023-07-12'
bld3 <- '2023-07-18'
bld4 <- '2023-07-26'

# Blaesedalen
names(blt1) <- m3.bands
names(blt2) <- m3.bands
names(blt3) <- m3.bands
names(blt4) <- m3.bands

# Create a list of the UAV rasters which can be passed to functions
bl <- list(blt1, blt2, blt3, blt4)

# Function to calculate the NDVI
calc_ndvi_m3m <- function(x){
  x <- (x$nir - x$red) / (x$nir + x$red) # For M3M
}

bl.ndvi <- pblapply(bl, calc_ndvi_m3m)

bl.ndvi.rast <- rast(bl.ndvi)

names(bl.ndvi.rast) <- c(bld1, bld2, bld3, bld4)


# extract from rast to points
snow.ndvi <- terra::extract(bl.ndvi.rast, snow.points)

# pivot
snow.ndvi <- snow.ndvi %>% pivot_longer(cols = !ID, names_to = 'date', values_to = 'ndvi')

ggplot() +
  geom_spatraster(data = bl.ndvi.rast[[4]]) +
  geom_sf(data = snow.points)

# Plot
snow.plot <- ggplot() +
  geom_line(data = snow.ndvi, aes(x = lubridate::as_date(date), y = ndvi, group = ID), alpha = 0.2) +
  scale_x_continuous(breaks = c(as_date(bld1), as_date(bld2), as_date(bld3), as_date(bld4)),
                     labels = c('July 2nd', 'July 12th', 'July 18th', 'July 26th')) + 
  ylab('NDVI') +
  xlab('Date') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))

cowplot::save_plot('../../plots/supplementary/snow-ndvi.png', snow.plot, 
                   base_height = 70, base_width = 140, bg = 'white', units = 'mm')



























































