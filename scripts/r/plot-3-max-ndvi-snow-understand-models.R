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


# Create linear models for later testing with Moran's I ----

# Maximum NDVI and NDVI DoY is a function of snow persistence
# Blaesedalen
lm.s2.bl.max <- lm(s2.bl$ndvi.max ~ s2.bl$snow.auc)
lm.s2.bl.doy <- lm(s2.bl$ndvi.max.doy ~ s2.bl$snow.auc)
# kluane low
lm.s2.kl.max <- lm(s2.kl$ndvi.max ~ s2.kl$snow.auc)
lm.s2.kl.doy <- lm(s2.kl$ndvi.max.doy ~ s2.kl$snow.auc)
# Kluane high
lm.s2.kh.max <- lm(s2.kh$ndvi.max ~ s2.kh$snow.auc)
lm.s2.kh.doy <- lm(s2.kh$ndvi.max.doy ~ s2.kh$snow.auc)


# Sensitivity plots for determining d2 ----
sensit.df <- data.frame()
# Function for calculating moran's I with increasing neighbourhood dist, plotting result
sensit_check <- function(site.name, model, data) {
  # Empty vector for storage of Moran's I values
  moran.I <- c()
  for (d in seq(10, 400, 10)) {
    s2.sens.nb <- dnearneigh(data, d1 = 0, d2 = d)
    s2.sens.lw <- nb2listw(s2.sens.nb, style = 'W', zero.policy = TRUE)
    moran <- moran.mc(model$residuals, s2.sens.lw, 
                      nsim = 999, zero.policy = TRUE)
    moran.I <- c(moran.I, moran$statistic)
  }
  
  moran.I <- data.frame(moran = moran.I,
                        distance = seq(10, 400, 10))
  
  ggplot(data = moran.I, aes(x = distance, y = moran)) +
    geom_point(aes(y = moran)) +
    geom_line(aes(y = moran))
  
  sensit.df <<- moran.I %>%
    mutate(site = site.name)
  
}

sensit.df
# Run function on each model
# Blaesedalen
sensit_check('blaesedalen', lm.s2.bl.max, s2.bl)
sensit.bl <- sensit.df
# Kluane low
sensit_check('kluane-low', lm.s2.kl.max, s2.kl)
sensit.kl <- sensit.df
# Kluane high
sensit_check('kluane-high', lm.s2.kh.max, s2.kh)
sensit.kh <- sensit.df

# Bind the sensitivity dataframe for each site together
sensit.all <- rbind(sensit.bl, sensit.kl, sensit.kh)
sensit.all
# Output this df to prevent running sensitivity checks again (slow)
write_csv(sensit.all, '../../data/statistical-output/morans-i-sensitivity-neighbourhood.csv')

# Combined plot of the sensitivity check for Moran's I
ggplot() +
  geom_line(data = sensit.all, aes(x = distance, y = moran, group = site, colour = site)) #+
  #geom_vline(aes(x = dist2))

# read in the moran's i data with increasing distance of neighbourhood, sensitivity check
sensit.data <- read_csv('../../data/statistical-output/morans-i-sensitivity-neighbourhood.csv')

# Spatial Error Models ----

# Set maximum distance of neighbourhood
dist2 <- 200

# For loop to iterate through increasing neighbourhood dist and produce plots
for (dist2 in seq(10, 400, 10)) {

###
# Blaesedalen Spatial Error Models
###
s2.bl.nb <- dnearneigh(s2.bl, d1 = 0, d2 = dist2)
#s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s2.bl.lw <- nb2listw(s2.bl.nb, style = 'W', zero.policy = FALSE)

sem.s2.bl.max <- errorsarlm(ndvi.max ~ snow.auc,
                            data = s2.bl, 
                            listw = s2.bl.lw,
                            zero.policy = FALSE)

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


# Examine model summaries ----
stargazer(sem.s2.bl.max, sem.s2.kl.max, sem.s2.kh.max, type = 'text')


# Checking what is going on with models ----

# Scatter plots of actual and fitted values
out.plot <- plot_grid(
ggplot() +
  geom_point(data = s2.kl, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
  geom_point(data = s2.kl %>% mutate(fitted = fitted(sem.s2.kl.max)), aes(x = snow.auc, y = fitted), colour = 'blue') +
  ylim(0, 0.8),
ggplot() +
  geom_point(data = s2.kh, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
  geom_point(data = s2.kh %>% mutate(fitted = fitted(sem.s2.kh.max)), aes(x = snow.auc, y = fitted), colour = 'blue') +
  ylim(0, 0.6),
ggplot() +
  geom_point(data = s2.bl, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
  geom_point(data = s2.bl %>% mutate(fitted = fitted(sem.s2.bl.max)), aes(x = snow.auc, y = fitted), colour = 'blue') +
  ylim(0, 0.55),
ggplot() +
  geom_line(data = sensit.all, aes(x = distance, y = moran, group = site, colour = site)) +
  geom_vline(xintercept = dist2),
labels = c('KL', 'KH', 'BL', paste0('nb =', dist2, 'm'))
)

out.plot

#cowplot::save_plot(paste0('../../plots/moran-sensitivity/neighbourhood-', dist2, '-metres.png'), out.plot, base_height = 140, base_width = 260, units = 'mm')


# Maps of observed and fitted data
# Kluane Low
kl.maps <- plot_grid(
  ggplot() +
    geom_sf(data = s2.kl, aes(colour = ndvi.max, size = snow.auc)) +
    scale_colour_viridis_c(),
  #scale_colour_viridis_c(limits = c(0.1, 0.6)),
  ggplot() +
    geom_sf(data = s2.kl %>% mutate(
      fitted = fitted(sem.s2.kl.max)),
      mapping = aes(colour = fitted, size = snow.auc)) +
    scale_colour_viridis_c(), 
  labels = c('KL', dist2)
)

cowplot::save_plot(paste0('../../plots/moran-sensitivity/kl-maps/neighbourhood-', dist2, '-metres.png'), kl.maps, base_height = 140, base_width = 260, units = 'mm')

# Kluane high
kh.maps <- plot_grid(
  ggplot() +
    geom_sf(data = s2.kh, aes(colour = ndvi.max, size = snow.auc)) +
    scale_colour_viridis_c(),
  #scale_colour_viridis_c(limits = c(0.1, 0.6)),
  ggplot() +
    geom_sf(data = s2.kh %>% mutate(
      fitted = fitted(sem.s2.kh.max)),
      mapping = aes(colour = fitted, size = snow.auc)) +
    scale_colour_viridis_c(),
  labels = c('KH', dist2)
)

cowplot::save_plot(paste0('../../plots/moran-sensitivity/kh-maps/neighbourhood-', dist2, '-metres.png'), kh.maps, base_height = 140, base_width = 260, units = 'mm')

# Blaesedalen
bl.maps <- plot_grid(
  ggplot() +
    geom_sf(data = s2.bl, aes(colour = ndvi.max, size = snow.auc)) +
    scale_colour_viridis_c(),
  #scale_colour_viridis_c(limits = c(0.1, 0.6)),
  ggplot() +
    geom_sf(data = s2.bl %>% mutate(
      fitted = fitted(sem.s2.bl.max)),
      mapping = aes(colour = fitted, size = snow.auc)) +
    scale_colour_viridis_c(), 
  labels = c('BL', dist2)
)

cowplot::save_plot(paste0('../../plots/moran-sensitivity/bl-maps/neighbourhood-', dist2, '-metres.png'), bl.maps, base_height = 140, base_width = 260, units = 'mm')

}

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
    #scale_colour_viridis_c(limits = c(0.1, 0.6)),
  ggplot() +
    geom_sf(data = s2.kl %>% mutate(
      fitted = fitted(sem.s2.kl.max)),
      mapping = aes(colour = fitted, size = snow.auc)) +
      scale_colour_viridis_c(),
    #scale_colour_viridis_c(limits = c(0.1, 0.6)), 
  ggplot() +
    geom_sf(data = s2.bl, aes(colour = ndvi.max, size = snow.auc)) +
    scale_colour_viridis_c(),
    #scale_colour_viridis_c(limits = c(0.1, 0.6)),
  ggplot() +
    geom_sf(data = s2.bl %>% mutate(
      fitted = fitted(sem.s2.bl.max)),
      mapping = aes(colour = fitted, size = snow.auc)) +
      scale_colour_viridis_c()
    #scale_colour_viridis_c(limits = c(0.1, 0.6))
)

# Sensitivity plot on its own
sens.plot <- ggplot() +
  geom_line(data = sensit.all, aes(x = distance, y = moran, group = site, colour = site)) #+
  #geom_vline(xintercept = dist2),


cowplot::save_plot(paste0('../../plots/moran-sensitivity/sensitivity-plot.png'), sens.plot, base_height = 140, base_width = 260, units = 'mm')
