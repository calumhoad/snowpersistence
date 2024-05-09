# Fit all models, output tables
# Calum Hoad, calum.hoad@ed.ac.uk, 09/05/2024

# Dependencies ----
library(tidyverse)
library(sf)
library(ggplot2)
library(cowplot)
library(pbapply)
library(INLA)
library(gt)
library(stargazer)
library(gtsummary)
library(broom)

# Data prep ----

## Load data
# Blasedalen
s2.bl <- read_csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "S2 BL", range = round(51.6/10), colour = '#4984BF')

# Kluane low
s2.kl <- read_csv('../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "S2 KL", range = round(61.9/10), colour = '#F5A40C')

# Kluane high
s2.kh <- read_csv('../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "S2 KH", range = round(38.5/10), colour = '#F23835')

# Blaesedalen, S30
s30.bl <- read_csv("../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv",
                   show_col_types = FALSE
) %>%
  st_as_sf(coords = c("X", "Y"), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "S30 BL", range = round(110/30), colour = '#4984BF')

# Combine into one data object (keeping order of Calum's plots)
data_list <- list(s2.kl, s2.kh, s2.bl, s30.bl)


# Spearmans rho ----

# Simple correlation tables (from Jakob's script)
bind_rows(
  map(data_list, function(site_data) {
    tibble(
      site_name = unique(site_data$site),
      var_x = "snow.auc",
      var_y = "ndvi.max",
      spearmans_rho = round(cor(site_data$snow.auc, site_data$ndvi.max, method = "spearman"), 3),
      p_value = round(cor.test(site_data$snow.auc, site_data$ndvi.max, method = "spearman")$p.value, 3)
    )
  }),
  map(data_list, function(site_data) {
    tibble(
      site_name = unique(site_data$site),
      var_x = "snow.auc",
      var_y = "ndvi.max.doy",
      spearmans_rho = round(cor(site_data$snow.auc, site_data$ndvi.max.doy, method = "spearman"), 3),
      p_value = round(cor.test(site_data$snow.auc, site_data$ndvi.max.doy, method = "spearman")$p.value, 3)
    )
  })
) %>% gt() %>%
  gtsave("../../plots/stats-tables/spearmans_table.png")


# Linear models, OLS regression ----

# peak-NDVI magnitude
lm.bl.max <- lm(data = s2.bl, formula = ndvi.max ~ snow.auc)
lm.kl.max <- lm(data = s2.kl, formula = ndvi.max ~ snow.auc)
lm.kh.max <- lm(data = s2.kh, formula = ndvi.max ~ snow.auc)

stargazer(lm.kl.max, lm.kh.max, lm.bl.max, type = 'html',
          column.labels = c("Kluane Low", "Kluane High", "Blaesedalen"), 
          out = '../../plots/stats-tables/linear-models-max.html')

# peak-NDVI magnitude, KH (y = ln(x + 1))
ln.kh.max <- lm(data = s2.kh, formula = ndvi.max ~ log(snow.auc + 1))

stargazer(ln.kh.max, type = 'html',
          column.labels = "Kluane High", 
          out = '../../plots/stats-tables/linear-models-kh-ln-max.html')

# peak-NDVI timing
lm.bl.doy <- lm(data = s2.bl, formula = ndvi.max.doy ~ snow.auc)
lm.kl.doy <- lm(data = s2.kl, formula = ndvi.max.doy ~ snow.auc)
lm.kh.doy <- lm(data = s2.kh, formula = ndvi.max.doy ~ snow.auc)

stargazer(lm.kl.doy, lm.kh.doy, lm.bl.doy, type = 'html',
          column.labels = c("Kluane Low", "Kluane High", "Blaesedalen"), 
          out = '../../plots/stats-tables/linear-models-doy.html')

# peak-NDVI timing, KH (y = ln(x + 1))
ln.kh.doy <- lm(data = s2.kh, formula = ndvi.max.doy ~ log(snow.auc + 1))

stargazer(ln.kh.doy, type = 'html',
          column.labels = "Kluane High", 
          out = '../../plots/stats-tables/linear-models-kh-ln-doy.html')


# INLA Matern 2D models ----
