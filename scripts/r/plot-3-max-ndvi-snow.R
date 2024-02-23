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

# Models ----

###
# Blaesedalen Spatial Error Models
###

s2.bl.nb <- dnearneigh(s2.bl, d1 = 0, d2 = 60)
#s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
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

stargazer(sem.s2.bl.max, sem.s2.bl.doy, type = 'html', 
          out = '../../data/statistical-output/sem-blaesedalen-nb-60.html')

###
# Kluane low Spatial Error Models
###

s2.kl.nb <- dnearneigh(s2.kl, d1 = 0, d2 = 60)
#s2.kl.nb <- poly2nb(st_buffer(s2.kl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
s2.kl.lw <- nb2listw(s2.kl.nb, style = 'W', zero.policy = TRUE)

sem.s2.kl.max <- errorsarlm(ndvi.max ~ snow.auc, 
                            data = s2.kl, 
                            listw = s2.kl.lw, 
                            zero.policy = TRUE)

sem.s2.kl.doy <- errorsarlm(ndvi.max.doy ~ snow.auc, 
                            data = s2.kl, 
                            listw = s2.kl.lw, 
                            zero.policy = TRUE)

stargazer(sem.s2.kl.max, sem.s2.kl.doy, type = 'html',
          out = '../../data/statistical-output/sem-kluane-low-nb-60.html')

###
# Kluane high
###
s2.kh.nb <- dnearneigh(s2.kh, d1 = 0, d2 = 60)
#s2.kh.nb <- poly2nb(st_buffer(s2.kh, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
s2.kh.lw <- nb2listw(s2.kh.nb, style = 'W', zero.policy = TRUE)

sem.s2.kh.max <- errorsarlm(ndvi.max ~ snow.auc, 
                            data = s2.kh, 
                            listw = s2.kh.lw,
                            zero.policy = TRUE)

sem.s2.kh.doy <- errorsarlm(ndvi.max.doy ~ snow.auc,
                            data = s2.kh, 
                            listw = s2.kh.lw, 
                            zero.policy = TRUE)

stargazer(sem.s2.kh.max, sem.s2.kh.doy, type = 'html', 
          out = '../../data/statistical-output/sem-kluane-high-nb-60.html')


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

# Blaesedalen
bl.coefficients <- sem.s2.bl.max$coefficients
bl.equation <- sprintf("y = %.3f + %.3f * x", bl.coefficients[1], bl.coefficients[2])

bl <- ggplot() +
  geom_point(data = s2.bl, aes(x = snow.auc, y = ndvi.max), alpha = 0.5, color = '#4984BF') +
  geom_line(data = bl.pred, aes(x = snow.auc, y = ndvi), color = 'black', linewidth = 1) +
  annotate("text", x = max(s2.bl$snow.auc), y = 0.6, label = bl.equation, hjust = 1, vjust = 1) +
  xlab('') +
  ylab('') +
  theme_cowplot()

# Kluane low
kl.coefficients <- sem.s2.kl.max$coefficients
kl.equation <- sprintf("y = %.3f + %.3f * x", kl.coefficients[1], kl.coefficients[2])

kl <- ggplot() +
  geom_point(data = s2.kl, aes(x = snow.auc, y = ndvi.max), alpha = 0.5, color = '#F5A40C') +
  geom_line(data = kl.pred, aes(x = snow.auc, y = ndvi), color = 'black', linewidth = 1) +
  annotate("text", x = max(s2.kl$snow.auc), y = 0.6, label = kl.equation, hjust = 1, vjust = 1) + 
   xlab('') +
  ylab('') +
  theme_cowplot()

# Kluane high
kh.coefficients <- sem.s2.kh.max$coefficients
kh.equation <- sprintf("y = %.3f + %.3f * x", kh.coefficients[1], kh.coefficients[2])

kh <- ggplot() +
  geom_point(data = s2.kh, aes(x = snow.auc, y = ndvi.max), alpha = 0.5, color = '#F23835') +
  geom_line(data = kh.pred, aes(x = snow.auc, y = ndvi), color = 'black', linewidth = 1) +
  annotate("text", x = max(s2.kh$snow.auc), y = 0.6, label = kh.equation, hjust = 1, vjust = 1) +  
  xlab('') +
  ylab('') +
  theme_cowplot()

# Combine Kluane plots to single row
kluane.plot <- plot_grid(kl, kh, 
                         ncol = 2, 
                         align = 'h')

                         labels = c('(b) ', '(c) '))

# Add Blaesedalen plot as extra row
combined.plots <- plot_grid(bl,
                            kluane.plot,
                            nrow = 2, 
                            align = 'h') 
                            labels = c('(a) ', '', ''))

# Show plots
combined.plots

# Save plots
cowplot::save_plot('../../plots/figures/figure-3v5.png', combined.plots, base_height = 140, base_width = 180, units = 'mm')
