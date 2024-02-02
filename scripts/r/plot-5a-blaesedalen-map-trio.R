# Plot 5a for paper: Mapped values at Blaesedalen
# Calum Hoad, 1 Feb 2024

library(ggplot2)
library(dplyr)
library(sf)
library(mapview)
library(spdep)
library(spatialreg)
library(stargazer)
library(cowplot)
library(tidyr)

# Data ----
# Read in the data and reduce to unique id
# Blaesedalen
s2.bl <- read.csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  st_buffer(dist = 5, endCapStyle = "SQUARE")

# Maps ----

ndvi.max <- ggplot() +
  geom_sf(data = s2.bl, aes(fill = ndvi.max)) +
  scale_fill_viridis_c()

snow.auc <- ggplot() +
  geom_sf(data = s2.bl, aes(fill = snow.auc)) +
  scale_fill_viridis_c()

ndvi.max.doy <- ggplot() +
  geom_sf(data = s2.bl, aes(fill = ndvi.max.doy)) +
  scale_fill_viridis_c()

combined <- plot_grid(snow.auc, 
                      ndvi.max,
                      ndvi.max.doy,
                      ncol = 3, 
                      align = 'h')

combined
