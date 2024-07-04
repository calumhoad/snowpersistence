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

# Blaesedalen S30
s30.bl <- read.csv('../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs= 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  st_buffer(dist = 15, endCapStyle = "SQUARE")

# Maps ----

ndvi.max <- ggplot() +
  geom_sf(data = s2.bl, aes(fill = ndvi.max)) +
  scale_fill_viridis_c()

snow.auc <- ggplot() +
  geom_sf(data = s2.bl, aes(fill = snow.auc)) +
  scale_fill_viridis_c() +
  theme_void()

ndvi.max.doy <- ggplot() +
  geom_sf(data = s2.bl, aes(fill = ndvi.max.doy)) +
  scale_fill_viridis_c() +
  theme_void()

s30.ndvi.max.doy <- ggplot() +
  geom_sf(data = s30.bl, aes(fill = ndvi.max.doy)) +
  scale_fill_viridis_c() +
  theme_void()

combined <- plot_grid(ndvi.max.doy,
                      snow.auc, 
                      s30.ndvi.max.doy,
                      ncol = 3, 
                      align = 'h')

combined

cowplot::save_plot('../../plots/figures/figure-5av1.png', combined, base_height = 140, base_width = 180, units = 'mm')

s30.ndvi.max.doy
