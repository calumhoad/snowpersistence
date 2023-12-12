# Compare snow cover with derived landsat trends from LandsatTS package
# Calum Hoad, 06/12/2023

# Import necessary libraries
library(terra)
library(dplyr)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)

# Read in the data
trends.auto <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_auto_7yr_trnds.csv')
trends.manual <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_manual_7yr_trnds.csv')
gs.auto <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_auto_7yr_gs_metric.csv')
gs.manual <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_manual_7yr_gs_metric.csv')

view(gs.auto)

ls.snow <- read.csv2('../../data/uav/snow-metrics/blaesedalen_30m_snowcover.csv') %>%
  rename(sample.id = sample_id)

s2.snow <- read.csv2('../../data/uav/snow-metrics/blaesedalen_10m_snowcover.csv') %>%
  rename(ndvi.max = 'ndvi_mx', ndvi.max.doy = 'ndv_mx_') %>%
  mutate(ndvi.max.doy = yday(ndvi.max.doy))
# Join data frames based on sample_id
joined.auto <- left_join(gs.auto, ls.snow, by = "sample.id")

# Filter gs datasets to only 2023
gs.auto <- gs.auto %>% filter(year == 2023 & ndvi.max.doy > 175)
gs.auto.joined <- left_join(gs.auto, snow, by = 'sample.id')



# S2 max ndvi doy against snow persistence
s2.max.ndvi.doy.plot <- ggplot(s2.snow, aes(x = ndvi.max.doy, y = snow.persist)) +
  geom_point(aes(color = ndvi.max)) +
  geom_smooth(method = 'lm') +
  scale_color_viridis_c(name = "ndvi.max")

# S2 max ndvi against snow persistence
s2.max.ndvi.plot <- ggplot(s2.snow, aes(x = ndvi.max, y = snow.persist)) +
  geom_point(aes(color = ndvi.max.doy)) +
  geom_smooth(method = 'lm') +
  scale_color_viridis_c(name = "ndvi.max.doy") 

# Arrange plots side by side
plots_combined <- s2.max.ndvi.plot + s2.max.ndvi.doy.plot +
  plot_layout(ncol = 2)

plots_combined




