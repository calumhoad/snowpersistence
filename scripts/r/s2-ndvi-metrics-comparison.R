# Script to compare NDVI max and NDVI max DoY between different models
# Calum Hoad, 5th Jan 2024 

library(tidyverse)
library(ggplot2)
library(cowplot)

# Bring in the data from the various model outputs ----

# Parabolic 2nd order polynomial
s2.para <- read.csv2('../../data/sentinel-2/output/s2_modelled_point_long.csv') %>%
    group_by(id) %>%
    filter(row_number() == 1) %>%
    select(id, geometry, ndvi.max, ndvi.max.doy) %>%
    rename(geometry.p = 'geometry', 
           ndvi.max_p = 'ndvi.max',
           ndvi.max.doy_p = 'ndvi.max.doy') %>%
    mutate(ndvi.max.doy_p = floor(ndvi.max.doy_p))

# Smoothed spline
s2.smooth <- read.csv2('../../data/sentinel-2/output/s2_modelled_JA_point_wide.csv') %>%
  select(id, geometry, ndvi.max, ndvi.max.doy) %>%
  rename(geometry.s = 'geometry', 
         ndvi.max_s = 'ndvi.max', 
         ndvi.max.doy_s = 'ndvi.max.doy')

# Beck double logistic
s2.beck <- read.csv('../../data/sentinel-2/output/s2_modelled_beck_point_wide.csv') %>%
  select(id, X, Y, ndvi.max, ndvi.max.doy) %>%
  rename(geometry.x.b = 'X',
         geometry.y.b = 'Y', 
         ndvi.max_b = 'ndvi.max', 
         ndvi.max.doy_b = 'ndvi.max.doy')


# Join the data from the various model outputs ----
s2.all <- left_join(s2.para, s2.smooth, by = 'id')
s2.all <- left_join(s2.all, s2.beck, by = 'id')

# Manual check of geometry values shows that all ids are consistent, drop geom 
# columns which are redundant
s2.all <- s2.all %>%
  select(-geometry.p, -geometry.s) %>%
  st_as_sf(coords = c("geometry.x.b", "geometry.y.b"), crs = 32621)


# Compare and plot the data ----
s2.all.long <- s2.all %>%
  pivot_longer(
      cols = starts_with("ndvi.max"),
      names_to = c(".value", "model"),
      names_sep = "_",
      values_to = "ndvi.max"
    )

# Plot each pixel's ndvi.max and ndvi.max.doy for all three models, 
# and the point for each model with a line
ggplot(s2.all.long %>% filter(id %in% sample(unique(s2.all.long$id), 20)) %>% group_by(id)) +
  geom_line(aes(x = ndvi.max.doy, y = ndvi.max, group = id)) +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max, color = model)) +
  xlim(c(190, 260)) +
  ylim(c(0, 0.7)) +
  theme(cowplot())



ggplot(s2.all %>% filter(id %in% sample(unique(s2.all$id), 100))) +
  geom_point(aes(x = ndvi.max.doy.p, y = ndvi.max.p, colour = 'red')) +
  geom_point(aes(x = ndvi.max.doy.s, y = ndvi.max.s, colour = 'blue')) +
  geom_point(aes(x = ndvi.max.doy.b, y = ndvi.max.b, colour = 'violet')) +
  xlim(c(180, 260)) +
  ylim(c(0, 0.7))








