# Work out the difference between typical peak-NDVI timing and
# magnitude for high and low snow persistence pixels
# Calum Hoad, 08/05/2024

# Dependencies ----
library(tidyverse)
library(sf)
library(ggplot2)
library(cowplot)
library(colorspace)
library(pbapply)
library(broom)
library(INLA)
library(rlang)
library(gstat)

## Load data ----
# Blasedalen
s2.bl <- read_csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Kluane low
s2.kl <- read_csv('../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Kluane high
s2.kh <- read_csv('../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Blaesedalen, S30
s30.bl <- read_csv('../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv',
                   show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Filter the data to obtain the upper and lower quartile of snow.auc
upper <- function(data) {
  lower.quartile <- quantile(s2.bl$snow.auc, probs = 0.25)
  upper.quartile <- quantile(s2.bl$snow.auc, probs = 0.75)
  data <- data %>% filter(snow.auc >= upper.quartile)
}

lower <- function(data) {
  lower.quartile <- quantile(s2.bl$snow.auc, probs = 0.25)
  upper.quartile <- quantile(s2.bl$snow.auc, probs = 0.75)
  data <- data %>% filter(snow.auc <= lower.quartile)
}

bl.max.diff <- (mean(lower(s2.bl)$ndvi.max) - mean(upper(s2.bl)$ndvi.max))
bl.max.diff
kl.max.diff <- (mean(lower(s2.kl)$ndvi.max) - mean(upper(s2.kl)$ndvi.max))
kl.max.diff
kh.max.diff <- (mean(lower(s2.kh)$ndvi.max) - mean(upper(s2.kh)$ndvi.max))
kh.max.diff
all.max.mean <- mean(c(bl.max.diff, kl.max.diff, kh.max.diff))
all.max.mean

bl.doy.diff <- (mean(lower(s2.bl)$ndvi.max.doy) - mean(upper(s2.bl)$ndvi.max.doy))
bl.doy.diff
kl.doy.diff <- (mean(lower(s2.kl)$ndvi.max.doy) - mean(upper(s2.kl)$ndvi.max.doy))
kl.doy.diff
kh.doy.diff <- (mean(lower(s2.kh)$ndvi.max.doy) - mean(upper(s2.kh)$ndvi.max.doy))
kh.doy.diff
all.doy.mean <- mean(c(bl.doy.diff, kl.doy.diff, kh.doy.diff))
all.doy.mean
