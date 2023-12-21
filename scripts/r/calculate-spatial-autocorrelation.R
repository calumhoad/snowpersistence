# Calculate spatial autocorrelation of the dataset
# Calum Hoad, 2023/12/20

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
library(spdep)
#library(janitor)


# Read in data ----

# Sentinel-2 (30,) snow cover data, drived from UAV imagery
s2.snow <- st_read('../../data/sentinel-2/output/s2_modelled_point_wide.shp')

s2.snow.polygons <- st_buffer(s2.snow, dist = 15, endCapStyle = "SQUARE")
plot(s2.snow.polygons)

ggplot() + geom_sf(data = s2.snow.polygons)


# Pre-requisites for calculating Moran's Index ----

# Neighbour list
nb <- knn2nb(knearneigh(st_coordinates(s2.snow), k = 5))

# Compute weightings from the neighbours list
weights <- nb2listw(nb, style = 'W')


# Compute Moran's Index ----

###
# Following tutorial here (https://rspatial.org/analysis/3-spauto.html)
###

# sf to spatVector
s2.snow.polygons <- vect(s2.snow.polygons)

# Matrix of adjacent pixels, with method = Queen (as in chess)
ww <- terra::adjacent(s2.snow.polygons, 'queen', pairs = FALSE)

# Autocor function from the terra package
ac <- autocor(s2.snow.polygons$ndvi_mx, ww, 'moran')

# Create 99 more Moran's tests, using random data
m <- sapply(1:99, function(i) {
  autocor(sample(s2.snow.polygons$ndvi_mx), ww, 'moran')
})

# Pull out a p-value, representative of the number of times the Moran's Index
# on the random data exceeds or is equal to the Moran's test on the standard data
pval <- sum(m >= ac)/100

# Moran scatter plot
n <- length(s2.snow.polygons)

ms <- cbind(id = rep(1:n, each = n), y = rep(s2.snow.polygons$ndvi_mx, each = n), value = as.vector(wm * s2.snow.polygons$ndvi_mx))
