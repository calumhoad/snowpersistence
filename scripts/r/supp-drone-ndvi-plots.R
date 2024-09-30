# Supplementary plots of drone derived NDVI
# Calum Hoad, 26th September 2024, calum.hoad@ed.ac.uk

# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(mapview)
library(spdep)
library(spatialreg)
library(stargazer)
library(cowplot)
library(terra)
library(magick)
library(ggnewscale)

# Bring in the data
# Turn off scientific notation
options(scipen = 999)

# Bring in and format the data ----

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

# Kluane Low ###
klt1 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-06-29-5cm-clipped.tif')
klt2 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-07-05-5cm-clipped.tif')
klt3 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-07-18-5cm-clipped.tif')
klt4 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-08-01-5cm-clipped.tif')
klt5 <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-08-14-5cm-clipped.tif')

# Dates
kld1 <- '2022-06-29'
kld2 <- '2022-07-05'
kld3 <- '2022-07-18'
kld4 <- '2022-08-01'
kld5 <- '2022-08-14'

# Kluane High ###
kht1 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-09-5cm-clipped.tif')
kht2 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-19-5cm-clipped.tif')
kht3 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-29-5cm-clipped.tif')
kht4 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-08-04-5cm-clipped.tif')
kht5 <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-08-13-5cm-clipped.tif')

# Dates
khd1 <- '2022-07-09'
khd2 <- '2022-07-19'
khd3 <- '2022-07-29'
khd4 <- '2022-08-04'
khd5 <- '2022-08-13'

# Assign band names
s2.bands <- c('violet', 'blue', 'green', 'red', 're1', 're2', 'nir1', 'nir2', 'nir3') # MAIA-S2
m3.bands <- c('green', 'nir', 're', 'red', 'RGB-R', 'RGB-G', 'RGB-B') # Mavic M3M bands 


# Assign band names to spatRaster layers

# Blaesedalen
names(blt1) <- m3.bands
names(blt2) <- m3.bands
names(blt3) <- m3.bands
names(blt4) <- m3.bands


# Kluane low
names(klt1) <- s2.bands
names(klt2) <- s2.bands
names(klt3) <- s2.bands
names(klt4) <- s2.bands
names(klt5) <- s2.bands

# Kluane high
names(kht1) <- s2.bands
names(kht2) <- s2.bands
names(kht3) <- s2.bands
names(kht4) <- s2.bands
names(kht5) <- s2.bands

# Create a list of the UAV rasters which can be passed to functions
bl <- list(blt1, blt2, blt3, blt4)
kl <- list(klt1, klt2, klt3, klt4, klt5)
kh <- list(kht1, kht2, kht3, kht4, kht5)

# Function to calculate the NDVI
calc_ndvi_m3m <- function(x){
  x <- (x$nir - x$red) / (x$nir + x$red) # For M3M
}

calc_ndvi_maia <- function(x) {
  x <- (x$nir2 - x$red) / (x$nir2+ x$red) 
}

# Apply the NDVI function over the lists
kl.ndvi <- pblapply(kl, calc_ndvi_maia)
kh.ndvi <- pblapply(kh, calc_ndvi_maia)
bl.ndvi <- pblapply(bl, calc_ndvi_m3m)

# Convert back to rasters per site
kl.ndvi.rast <- rast(kl.ndvi)
kh.ndvi.rast <- rast(kh.ndvi)
bl.ndvi.rast <- rast(bl.ndvi)

# Name layers with dates
names(kl.ndvi.rast) <- c(kld1, kld2, kld3, kld4, kld5)
names(kh.ndvi.rast) <- c(khd1, khd2, khd3, khd4, khd5)
names(bl.ndvi.rast) <- c(bld1, bld2, bld3, bld4)

# Plotting function
plot_ndvi <- function(x, rows, cols){
  ndvi.plot <- ggplot() +
  geom_spatraster(data = x) +
  scale_fill_viridis_c(limits = c(0, 1),
                       oob = scales::oob_squish) +
  facet_wrap(~lyr, nrow = rows, ncol = cols) +
  theme_cowplot() +
  theme(
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.title.y = element_blank(),  # Remove y-axis label
    axis.text.x = element_blank(),   # Remove x-axis text
    axis.text.y = element_blank()    # Remove y-axis text
  )
  return(ndvi.plot)
}

# Creating plots
kl.plots <- plot_ndvi(kl.ndvi.rast, cols = 2, rows = 3)
kh.plots <- plot_ndvi(kh.ndvi.rast, cols = 2, rows = 3)
bl.plots <- plot_ndvi(bl.ndvi.rast, cols = 2, rows = 2)

cowplot::save_plot('../../plots/figures/uav-ndvi-plots-kl.png', kl.plots, base_height = 240, base_width = 180, units = 'mm', bg = 'white')
cowplot::save_plot('../../plots/figures/uav-ndvi-plots-kh.png', kh.plots, base_height = 240, base_width = 180, units = 'mm', bg = 'white')
cowplot::save_plot('../../plots/figures/uav-ndvi-plots-bl.png', bl.plots, base_height = 200, base_width = 180, units = 'mm', bg = 'white')
