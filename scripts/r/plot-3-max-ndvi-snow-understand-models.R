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
library(broom)

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

dist2 <- 200

s2.bl.nb <- dnearneigh(s2.bl, d1 = 0, d2 = dist2)
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

s2.kl.nb <- dnearneigh(s2.kl, d1 = 0, dist2)
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

stargazer(sem.s2.kl.max, sem.s2.kl.doy, type = 'text')#,
          out = '../../data/statistical-output/sem-kluane-low-nb-60.html')

###
# Kluane high
###
s2.kh.nb <- dnearneigh(s2.kh, d1 = 0, dist2)
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

# Extract confidence intervals
bl.ci <- confint(sem.s2.bl.max, level = 0.95)
bl.ci.low.int <- bl.ci[2,1]
bl.ci.low.slo <- bl.ci[3,1]
bl.ci.high.int <- bl.ci[2,2]
bl.ci.high.slo <- bl.ci[3,2]

# Use slope and intercept of high and low CI to predict values
bl.y.low <- (bl.ci.low.slo * bl.pred$snow.auc) + bl.ci.low.int  # Low CI line
bl.y.high <- (bl.ci.high.slo * bl.pred$snow.auc) + bl.ci.high.int

bl.model <- tidy(sem.s2.bl.max, conf.int = TRUE, conf.level = 0.95)

str(sem.s2.bl.max)

# Remove snow.auc which occur in both the original data and predict
common.ids <- intersect(bl.pred$snow.auc, s2.bl$snow.auc)
bl.pred <- bl.pred %>% filter(!snow.auc %in% common.ids)

# Predict
bl.pred <- bl.pred %>% 
  mutate(ndvi = predict(sem.s2.bl.max, newdata = bl.pred))#,
         #ci.low = bl.m#, interval = "confidence", level = 0.9))


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
  geom_ribbon(aes(x = bl.pred$snow.auc, ymin = y.low, ymax = y.high), fill = '#4984BF', alpha = 0.15, color = NA) +
  annotate("text", x = max(s2.bl$snow.auc), y = 0.6, label = bl.equation, hjust = 1, vjust = 1) +
  xlab('') +
  ylab('') +
  theme_cowplot()

bl
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

                         #labels = c('(b) ', '(c) '))

# Add Blaesedalen plot as extra row
combined.plots <- plot_grid(bl,
                            kluane.plot,
                            nrow = 2, 
                            align = 'h') 
                            #labels = c('(a) ', '', ''))

# Show plots
combined.plots

# Save plots
cowplot::save_plot('../../plots/figures/figure-3v5.png', combined.plots, base_height = 140, base_width = 180, units = 'mm')


# Checking what is going on with models ----

# Plot predictor and response variables and model parameters to understand data better (Jakob code)
plot_grid(
  ggplot() +
    geom_sf(data = s2.kl, aes(colour = ndvi.max, size = snow.auc)) +
    scale_colour_viridis_c(),
  ggplot() +
    geom_sf(data = s2.kl %>% mutate(
      fitted = fitted(sem.s2.kl.max)),
      mapping = aes(colour = fitted, size = snow.auc)) +
    scale_colour_viridis_c()
)

# Same plot, but for Blaesedalen data
plot_grid(
  ggplot() +
    geom_sf(data = s2.kl, aes(colour = ndvi.max, size = snow.auc)) +
    scale_colour_viridis_c(),
  ggplot() +
    geom_sf(data = s2.kl %>% mutate(
      fitted = fitted(sem.s2.kl.max)),
      mapping = aes(colour = fitted, size = snow.auc)) +
    scale_colour_viridis_c(), 
  ggplot() +
    geom_sf(data = s2.bl, aes(colour = ndvi.max, size = snow.auc)) +
    scale_colour_viridis_c(),
  ggplot() +
    geom_sf(data = s2.bl %>% mutate(
      fitted = fitted(sem.s2.bl.max)),
      mapping = aes(colour = fitted, size = snow.auc)) +
    scale_colour_viridis_c()
)

# Scatter plots of actual and fitted values
plot_grid(
ggplot() +
  geom_point(data = s2.kl, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
  geom_point(data = s2.kl %>% mutate(fitted = fitted(sem.s2.kl.max)), aes(x = snow.auc, y = fitted), colour = 'blue'),
ggplot() +
  geom_point(data = s2.kh, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
  geom_point(data = s2.kh %>% mutate(fitted = fitted(sem.s2.kh.max)), aes(x = snow.auc, y = fitted), colour = 'blue'),
ggplot() +
  geom_point(data = s2.bl, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
  geom_point(data = s2.bl %>% mutate(fitted = fitted(sem.s2.bl.max)), aes(x = snow.auc, y = fitted), colour = 'blue'), 
labels = c('KL', 'KH', 'BL')
)
