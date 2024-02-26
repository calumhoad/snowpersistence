# Script to bring in LandsatTS data and plot against snow

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

# Load in the data ----

# Trend data
lsat.trnds <- read_csv("../../data/landsat-ts/trends/bl-auto-7-trends.csv") %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
  st_transform(crs = 32621) %>%
  st_buffer(dist = 15, endCapStyle = "SQUARE")

# Check the data
ggplot() +
  geom_sf(data = lsat.trnds)

# Snow data
lsat.snow <- read_csv("../../data/landsat-ts/snow/snow-cover-30m-blaesedalen-landsatts.csv")


# Join the data ----
lsat.trnds <- left_join(lsat.trnds, lsat.snow, by = 'sample.id')


# Plot the data
ggplot() +
  geom_point(data = lsat.trnds, aes(x = snow.auc, y = slope, color = trend.cat))

ggplot(data = lsat.trnds) +
  geom_boxplot(data = lsat.trnds, aes(snow.auc, trend.cat, fill = trend.cat)) +
  scale_fill_manual(values = c('#DB9B3B', '#F1EDE7', '#518663'))

# Statistical testing ----
browning <- filter(lsat.trnds, trend.cat == 'browning') %>%
  drop_na()
notrend <- filter(lsat.trnds, trend.cat == 'no_trend')

wilcox.test(browning$snow.auc, notrend$snow.auc)

ggplot() +
  geom_histogram(data = notrend, aes(snow.auc))
