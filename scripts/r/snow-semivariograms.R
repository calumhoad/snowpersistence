# Semivariograms for snow
# Calum Hoad, 7 March 2024

library(ggplot2)
library(dplyr)
library(sf)
library(mapview)
library(spdep)
library(spatialreg)
library(stargazer)
library(cowplot)
library(broom)
library(tidyr)
library(tidyverse)
library(tidyverse)
library(sp)
library(gstat)
library(spdep)
library(spatialreg)
library(ggplot2)
library(cowplot)
library(colorspace)
library(pbapply)
library(broom)
library(kableExtra)

install.packages('gstat')
# Turn off scientific notation
options(scipen = 999)


# Data ----
# Read in the data and reduce to unique id
# Blaesedalen
s2.bl <- read.csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Kluane low
s2.kl <- read.csv('../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Kluane high
s2.kh <- read.csv('../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Linear models

# Variograms ----

# Fit variogram snow
vario_snow <- as_Spatial(s2.kl) %>% as("SpatialPointsDataFrame") %>%
  variogram(snow.auc ~ 1, data = ., cutoff = 130, width = 10)
vario_snow.max_fit <- fit.variogram(vario_snow, model = vgm(model = "Sph")) 


# Visualise results
plot <- ggplot(data = vario_snow) +
    geom_point(aes(x = dist, y = gamma)) +
    geom_line(aes(x = dist, y = gamma), 
              data = variogramLine(vario_snow.max_fit, 
                                   dist_vector = seq(10,130,10))) +
    geom_vline(xintercept = vario_snow.max_fit$range) +
    annotate("text", x = vario_snow.max_fit$range, y  = Inf, 
             label = paste0(" range = ", round(vario_snow.max_fit$range, 1), " m"),
             hjust = 0, vjust = 1.5) +
    scale_x_continuous(limits = c(10, 130), breaks = seq(10,130,10)) +
    labs(x = "lag distance (m)", y = "semivariance (gamma)") +
    theme_cowplot()
plot


cowplot::save_plot('../../plots/snow-semivariograms/kluane-low.png', plot, base_height = 140, base_width = 260, units = 'mm', bg = 'white')
