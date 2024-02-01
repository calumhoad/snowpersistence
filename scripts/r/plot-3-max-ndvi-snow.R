# Plot 2 for paper - NDVI max is influenced by snow persistence
# Calum Hoad, 1 Feb 2024

library(ggplot2)
library(dplyr)
library(sf)
library(mapview)
library(spdep)
library(spatialreg)
library(stargazer)
library(cowplot)

# Data ----
# Read in the data and reduce to unique id
# Blaesedalen
s2.bl <- read.csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Kluane low
s2.kl <- read.csv('../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Kluane high
s2.kh <- read.csv('../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() ==1) %>%
  ungroup()

# Models ----

###
# Blaesedalen Spatial Error Models
###

#s2.bl.nb <- dnearneigh(s2.bl, d1 = 0, d2 = 10)
s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s2.bl.lw <- nb2listw(s2.bl.nb, style = 'W', zero.policy = FALSE)

sem.s2.bl.max <- errorsarlm(ndvi.max ~ snow.auc,
                            data = s2.bl, 
                            listw = s2.bl.lw,
                            zero.policy = FALSE)

sem.s2.bl.doy <- errorsarlm(ndvi.max.doy ~ snow.auc,
                            data = s2.bl,
                            listw = s2.bl.lw, 
                            zero.policy = FALSE)


###
# Kluane low Spatial Error Models
###

s2.kl.nb <- poly2nb(st_buffer(s2.kl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
s2.kl.lw <- nb2listw(s2.kl.nb, style = 'W', zero.policy = TRUE)

sem.s2.kl.max <- errorsarlm(ndvi.max ~ snow.auc, 
                            data = s2.kl, 
                            listw = s2.kl.lw, 
                            zero.policy = TRUE)

sem.s2.kl.doy <- errorsarlm(ndvi.max.doy ~ snow.auc, 
                            data = s2.kl, 
                            listw = s2.kl.lw, 
                            zero.policy = TRUE)

###
# Kluane high
###

s2.kh.nb <- poly2nb(st_buffer(s2.kh, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
s2.kh.lw <- nb2listw(s2.kh.nb, style = 'W', zero.policy = TRUE)

sem.s2.kh.max <- errorsarlm(ndvi.max ~ snow.auc, 
                            data = s2.kh, 
                            listw = s2.kh.lw,
                            zero.policy = TRUE)

sem.s2.kh.doy <- errorsarlm(ndvi.max.doy ~ snow.auc,
                            data = s2.kh, 
                            listw = s2.kh.lw, 
                            zero.policy = TRUE)


# Predictions ----
# Generate predicted values based on model

# Blaesedalen
bl.pred <- data.frame(snow.auc = seq(min(s2.bl$snow.auc), 
                                     max(s2.bl$snow.auc), 
                                     0.5))

bl.pred <- bl.pred %>% 
  mutate(ndvi = predict(sem.s2.bl.max, newdata = bl.pred))# %>%
  

# kluane low
kl.pred <- data.frame(snow.auc = seq(min(s2.kl$snow.auc),
                                     max(s2.kl$snow.auc),
                                     0.5))

kl.pred <- kl.pred %>%
  mutate(ndvi = predict(sem.s2.kl.max, newdata = kl.pred))

# Kluane high
kh.pred <- data.frame(snow.auc = seq(min(s2.kh$snow.auc),
                                     max(s2.kh$snow.auc),
                                     0.5))
kh.pred <- kh.pred %>%
  mutate(ndvi = predict(sem.s2.kh.max, newdata = kh.pred))


# Plot ----

bl <- ggplot() +
  geom_point(data = s2.bl, aes(x = snow.auc, y = ndvi.max)) +
  geom_line(data = bl.pred, aes(x = snow.auc, y = ndvi), color = 'red') +
  xlab('Snow persistence (area under curve)') +
  ylab('Maximum NDVI') +
  theme_cowplot()

kl <- ggplot() +
  geom_point(data = s2.kl, aes(x = snow.auc, y = ndvi.max)) +
  geom_line(data = kl.pred, aes(x = snow.auc, y = ndvi), color = 'red') +
  xlab('Snow persistence (area under curve)') +
  ylab('Maximum NDVI') +
  theme_cowplot()

kh <- ggplot() +
  geom_point(data = s2.kh, aes(x = snow.auc, y = ndvi.max)) +
  geom_line(data = kh.pred, aes(x = snow.auc, y = ndvi), color = 'red') +
  xlab('Snow persistence (area under curve)') +
  ylab('Maximum NDVI') +
  theme_cowplot()

combined.plots <- plot_grid(bl,
                            kl,
                            kh,
                            ncol = 1, 
                            align = 'v')

