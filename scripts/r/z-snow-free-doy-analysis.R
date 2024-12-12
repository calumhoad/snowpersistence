# Get snow free DOY from MODIS MOD10A1 data
# Calum Hoad, 23/11/2024

# Libraries
library(ncdf4)
library(terra)
library(stars)
library(ncmeta)
library(tidyverse)
library(raster)
library(lubridate)
library(tidyterra)
library(lubridate)
library(cowplot)

# Get coords of the sites
# Set the coords to extract data from the cdf files
coords <- data.frame(row.names = c('KL', 'BL'),
                     lon = c(-138.40488817583486, -53.45637381046565),
                     #lat = c(60.97933154339103, 69.29345972224118)) # Moved point northward, to avoid cell boundary.
                     lat = c(61.1, 69.29345972224118))    

# Load the MODIS MOD10A1 data
bl.mod.files <- list.files('../../data/mod-snow/MOD10A1_61-20241123_135820/', full.names = T)
kl.mod.files <- list.files('../../data/mod-snow/MOD10A1_61-20241123_190751/', full.names = T)

# Load the Arclim data
arclim.sse <- rast('../../data/arclim/arclim_SSE.nc')
arclim.mean <- rast('../../data/arclim/arclim_means.tif')

# Rename layers in Arclim mean
names(arclim.mean) <- c('GSL', 'GDD', 'FGS', 'FDD', 'ROS', 'WWE', 'WWI',
                        'HWMI', 'VPDI', 'SWI', 'SSL', 'SSO', 'SSE', 'HWE', 
                        'TAVG', 'PRA', 'SFA', 'WSA')

# Arclim extraction ----
# Set point for Blaesedalen
bl.point <- vect(coords[2,], crs = 'EPSG:4326') 
kl.point <- vect(coords[1,], crs = 'EPSG:4326') 

# Check for location
plot(crop(arclim.mean$SSE, buffer(bl.point, 100000)))
plot(bl.point, add = T)

plot(crop(arclim.mean$SSE, buffer(kl.point, 100000)))
plot(kl.point, add = T)

# Extract annual SSE from the Arclim annual SSE dataset
bl.arclim <- terra::extract(arclim.sse, bl.point) %>%
  tibble() %>%
  pivot_longer(cols = !ID, names_to = 'year', values_to = 'SSE')

kl.arclim <- terra::extract(arclim.sse, kl.point) %>%
  tibble() %>%
  pivot_longer(cols = !ID, names_to = 'year', values_to = 'SSE')

ggplot() +
  geom_line(data = bl.arclim, aes(x = year, y = SSE, group = ID), colour = 'blue') +
  geom_line(data = kl.arclim, aes(x = year, y = SSE, group = ID), colour = 'red') +
  theme_cowplot()

# Extract mean SSW 1991-2021 from the Arclim dataset
bl.mean.arclim <- terra::extract(arclim.mean, bl.point)
kl.mean.arclim <- terra::extract(arclim.mean, kl.point)


# MODIS extraction ----
# Set point for Blaesedalen
bl.point <- vect(coords[2,], crs = 'EPSG:4326') %>% project(rast(bl.mod.files[[1]]))
kl.point <- vect(coords[1,], crs = 'EPSG:4326') %>% project(rast(kl.mod.files[[1]]))

# Check point is within the first rast
ggplot() +
  geom_spatraster(data = rast(bl.mod.files[[1]])$NDSI_Snow_Cover %>% project('EPSG:4326')) +
  geom_spatvector(data = bl.point)

ggplot() +
  geom_spatraster(data = rast(kl.mod.files[[4]])$NDSI_Snow_Cover %>% project('EPSG:4326')) +
  geom_spatvector(data = kl.point)

# Dataframe to hold the data
extracted.data <- data.frame()
# Iterate through the data objects, extract the data, bind
for (file in bl.mod.files) {
  r <- rast(file)
  doy <- sub(".*A\\d{4}(\\d{3}).*", "\\1", file)
  snow.data <- data.frame(terra::extract(r, bl.point))
  snow.data$doy <- doy
  extracted.data <- bind_rows(extracted.data, snow.data)
}

# Dataframe to hold the data
kl.extracted.data <- data.frame()
# Iterate through the data objects, extract the data, bind
for (file in kl.mod.files) {
  r <- rast(file)
  doy <- sub(".*A\\d{4}(\\d{3}).*", "\\1", file)
  snow.data <- data.frame(terra::extract(r, kl.point))
  snow.data$doy <- as.numeric(doy)
  kl.extracted.data <- bind_rows(kl.extracted.data, snow.data)
}

# Filter to values < 100 and convert doy to date object
extracted.data <- extracted.data %>% 
  mutate(doy = as.numeric(doy)) %>%
  mutate(date = ymd(paste0("2023-01-01")) + days(as.numeric(doy) - 1)) %>%
  mutate(doy.corrected = yday(date)) %>%
           filter(NDSI_Snow_Cover <= 100)

# Filter to values < 100 and convert doy to date object
kl.extracted.data <- kl.extracted.data %>% 
  mutate(doy = as.numeric(doy)) %>%
  mutate(date = ymd(paste0("2022-01-01")) + days(as.numeric(doy) - 1)) %>%
  mutate(doy.corrected = yday(date)) %>%
  filter(NDSI_Snow_Cover <= 100)

# Plot out the time series
snow.plot <- ggplot() +
  geom_line(data = extracted.data, aes(x = doy, y = NDSI_Snow_Cover, group = ID), 
            colour = 'blue') +
  geom_line(data = kl.extracted.data, aes(x = doy, y = NDSI_Snow_Cover, group = ID),
            colour = 'red') +
  #geom_hline(yintercept = 40, linewidth = 0.5, colour = 'grey') +
  xlab('Day of Year') +
  ylab('MOD10A1 NDSI Snow Cover (%)') +
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, 365, by = 10)) +  # Adjust for appropriate range and interval
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save out the plot
cowplot::save_plot('../../plots/figures/modis-mod10A1-ndsi-snow-cover-both-sites.png',
                   snow.plot, 
                   base_height = 100,
                   base_width = 180,
                   units = 'mm',
                   bg = 'white')
