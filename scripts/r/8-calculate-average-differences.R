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
  data.no.zero <- data %>% filter(snow.auc != 0)
  upper.quartile <- quantile(data.no.zero$snow.auc, probs = 0.9)
  data <- data %>% filter(snow.auc >= upper.quartile)
  return(data)
}

lower <- function(data) {
  lower.quartile <- quantile(data$snow.auc, probs = 0.1)
  #data <- data %>% filter(snow.auc <= lower.quartile)
  data <- data %>% filter(snow.auc == 0)
  return(data)
}

# Calculate the difference between the average for the upper and lower quartile
bl.max.diff <- (mean(lower(s2.bl)$ndvi.max) - mean(upper(s2.bl)$ndvi.max))
bl.max.diff
kl.max.diff <- (mean(lower(s2.kl)$ndvi.max) - mean(upper(s2.kl)$ndvi.max))
kl.max.diff
kh.max.diff <- (mean(lower(s2.kh)$ndvi.max) - mean(upper(s2.kh)$ndvi.max))
kh.max.diff

# Average the difference between all plots
all.max.mean <- mean(c(bl.max.diff, kl.max.diff, kh.max.diff))
all.max.mean

# Calculate the difference between the average for the upper and lower quartile
bl.doy.diff <- (mean(lower(s2.bl)$ndvi.max.doy) - mean(upper(s2.bl)$ndvi.max.doy))
bl.doy.diff
kl.doy.diff <- (mean(lower(s2.kl)$ndvi.max.doy) - mean(upper(s2.kl)$ndvi.max.doy))
kl.doy.diff
kh.doy.diff <- (mean(lower(s2.kh)$ndvi.max.doy) - mean(upper(s2.kh)$ndvi.max.doy))
kh.doy.diff

# Average the difference between all plots
all.doy.mean <- mean(c(bl.doy.diff, kl.doy.diff, kh.doy.diff))
all.doy.mean

# Create a dataframe from the results
max.diff <- c(bl = bl.max.diff, kl = kl.max.diff, kh = kh.max.diff, av.diff = all.max.mean)
doy.diff <- c(bl = bl.doy.diff, kl = kl.doy.diff, kh = kh.doy.diff, av.diff = all.doy.mean)

diff.df <- data.frame(max.diff, doy.diff) 

# Comparison of S2 and HLSS30 model outputs
s2.peak.ols <- -0.005
s2.peak.inla <- -0.002
s2.doy.ols <- 0.38
s2.doy.inla <- 0.24

s30.peak.ols <- -0.01
s30.peak.inla <- -0.01
s30.doy.ols <- 1.399
s30.doy.inla <- 1.385

peak.ols.diff <- s30.peak.ols / s2.peak.ols
peak.inla.diff <- s30.peak.inla / s2.peak.inla
doy.ols.diff <- s30.doy.ols / s2.doy.ols
doy.inla.diff <- s30.doy.inla / s2.doy.ols

doy.ols.diff

# As percentage
(peak.ols.diff/s2.peak.ols)*2
(peak.inla.diff/s2.peak.inla)*2
(doy.ols.diff/s2.doy.ols)*100
(doy.inla.diff/s2.doy.inla)*100
