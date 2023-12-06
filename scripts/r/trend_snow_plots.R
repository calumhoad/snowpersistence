# Compare snow cover with derived landsat trends from LandsatTS package
# Calum Hoad, 06/12/2023

# Import necessary libraries
library(terra)
library(dplyr)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)

# Read in the data
trends.auto <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_auto_7yr_trnds.csv')
trends.manual <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_manual_7yr_trnds.csv')
gs.auto <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_auto_7yr_gs_metric.csv')
gs.manual <- read.csv2('../../data/lsatTS-output/blaesedalen/blaesedalen_manual_7yr_gs_metric.csv')

view(gs.auto)

snow <- read.csv2('../../data/uav/snow-metrics/blaesedalen_30m_snowcover.csv') %>%
  rename(sample.id = sample_id)

# Join data frames based on sample_id
joined.manual <- left_join(trends.manual, snow, by = "sample.id")

# Filter gs datasets to only 2023
gs.auto <- gs.auto %>% filter(year == 2023 & ndvi.max.doy > 175)
gs.auto.joined <- left_join(gs.auto, snow, by = 'sample.id')

