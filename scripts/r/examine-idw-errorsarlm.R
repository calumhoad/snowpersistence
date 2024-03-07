# Check effect of implementing inverse distance weighting on spatial error models
# Calum Hoad, 7 March 2024

library(ggplot2)
library(dplyr)
library(sf)
library(mapview)
library(spdep)
library(spatialreg)
library(stargazer)
library(cowplot)
library(broom)
library(tidyr)
library(tidyverse)

# Turn off scientific notation
options(scipen = 999)


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

# read in the moran's i data with increasing distance of neighbourhood, sensitivity check
lm.moran.data <- read_csv('../../data/statistical-output/morans-i-sensitivity-neighbourhood.csv')

# Create an empty dataframe to store model summary values
var.names <- c('site', 'weights', 'dist', 'aic', 'loglike', 'sig2', 'walds', 'waldp', 'lograts', 'logratp')
aic.df <- data.frame(matrix(ncol = length(var.names), nrow = 0))
colnames(aic.df) <- var.names

# Create linear models for later use
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


# Set dist2 manually, if running only part of loop
dist2 <- 60

# For loop to iterate through increasing neighbourhood dist and produce plots
for (dist2 in seq(10, 200, 10)) {
  
  ###
  # Blaesedalen Spatial Error Models
  ###
  
  # Define neighbourhood
  s2.bl.nb <- dnearneigh(s2.bl, d1 = 0, d2 = dist2)
  
  # Define spatial weight structure, equal weights
  s2.bl.lw <- nb2listw(s2.bl.nb, style = 'W', zero.policy = FALSE)
  # Define spatial wweight structure, inverse distance weights
  idw.s2.bl.lw <- nb2listwdist(s2.bl.nb, s2.bl, type="idw", style="U", 
                           alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)
  # SEM, equal weights
  sem.s2.bl.max <- errorsarlm(ndvi.max ~ snow.auc,
                              data = s2.bl, 
                              listw = s2.bl.lw,
                              zero.policy = FALSE)
  # SEM, IDW
  idw.sem.s2.bl.max <- errorsarlm(ndvi.max ~ snow.auc,
                              data = s2.bl, 
                              listw = idw.s2.bl.lw,
                              zero.policy = FALSE)
  
  
  ###
  # Kluane low Spatial Error Models
  ###
  
  # Define neighbourhood  
  s2.kl.nb <- dnearneigh(s2.kl, d1 = 0, dist2)

  # Define spatial weight structure, equal weights
  s2.kl.lw <- nb2listw(s2.kl.nb, style = 'W', zero.policy = TRUE)
  # Define spatial weight structure, IDW
  idw.s2.kl.lw <- nb2listwdist(s2.kl.nb, s2.kl, type="idw", style="U", 
                           alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)
  # SEM, equal weights
  sem.s2.kl.max <- errorsarlm(ndvi.max ~ snow.auc, 
                              data = s2.kl, 
                              listw = s2.kl.lw, 
                              zero.policy = TRUE)
  # SEM, IDW
  idw.sem.s2.kl.max <- errorsarlm(ndvi.max ~ snow.auc, 
                              data = s2.kl, 
                              listw = idw.s2.kl.lw, 
                              zero.policy = TRUE)
  
  
  ###
  # Kluane high
  ###
  
  # Define neighbourhood
  s2.kh.nb <- dnearneigh(s2.kh, d1 = 0, dist2)

  # Define spatial weight structure, equal weights
  s2.kh.lw <- nb2listw(s2.kh.nb, style = 'W', zero.policy = TRUE)
  # Define spatial weight structure, IDW
  idw.s2.kh.lw <- nb2listwdist(s2.kh.nb, s2.kh, type="idw", style="U", 
                           alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)
  # SEM, equal weights
  sem.s2.kh.max <- errorsarlm(ndvi.max ~ snow.auc, 
                              data = s2.kh, 
                              listw = s2.kh.lw,
                              zero.policy = TRUE)
  # SEM, IDW
  idw.sem.s2.kh.max <- errorsarlm(ndvi.max ~ snow.auc, 
                              data = s2.kh, 
                              listw = idw.s2.kh.lw,
                              zero.policy = TRUE)  
  
  # Examine model summaries ----
  
  # Blaesedalen
  stargazer(sem.s2.bl.max, idw.sem.s2.bl.max, type = 'html', 
            out = paste0('../../plots/nb-dist-sensitivity/idw-neighbour/model-output/bl/model-bl-maxndvi-', dist2, '.html'))
  
  # Kluane low
  stargazer(sem.s2.kl.max, idw.sem.s2.kl.max, type = 'html', 
            out = paste0('../../plots/nb-dist-sensitivity/idw-neighbour/model-output/kl/model-kl-maxndvi-', dist2, '.html'))
  
  # Kluane high
  stargazer(sem.s2.kh.max, idw.sem.s2.kh.max, type = 'html', 
            out = paste0('../../plots/nb-dist-sensitivity/idw-neighbour/model-output/bl/model-kh-maxndvi-', dist2, '.html'))
  
  
  # Extract the residuals from each model
  resid <- lm.s2.bl.max$residuals
  mor <- moran.mc(resid, s2.bl.lw, nsim = 999, zero.policy = FALSE)
  plot(mor)
  mx <- max(mor$res[1:999])
  mn <- min(mor$res)

# Run moran tests for models   
comparison_plot <- function(plot.name, data, lw, sem.model, idw.lw, idw.model) {  
  sem.moran <- moran.test(sem.model$residuals, lw)
  idw.sem.moran <- moran.test(idw.model$residuals, idw.lw)
  
  plot_grid(
    ggplot() +
      geom_point(data = data, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
      geom_point(data = data %>% mutate(fitted = fitted(sem.model)), aes(x = snow.auc, y = fitted), colour = 'green'),
    ggplot() + 
      geom_point(data = data, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
      geom_point(data = data %>% mutate(fitted = fitted(idw.model)), aes(x = snow.auc, y = fitted), colour = 'purple'),
    ggplot() +
      geom_point(data = data, aes(x = snow.auc, y = ndvi.max), colour = 'black', alpha = 0.1) +
      geom_point(data = data %>% mutate(fitted = fitted(sem.model)), aes(x = snow.auc, y = fitted), colour = 'green') +
      geom_point(data = data %>% mutate(fitted = fitted(idw.model)), aes(x = snow.auc, y = fitted), colour = 'purple'), 
    ggplot() +
      geom_line(data = lm.moran.data %>% filter(site == plot.name), aes(x = distance, y = moran, group = site), colour = 'red') +
      geom_vline(xintercept = dist2) +
      geom_hline(yintercept = mx) +
      geom_hline(yintercept = mn) +
      geom_hline(yintercept = sem.moran$estimate[[1]], colour = 'green', linewidth = 1) +
      geom_hline(yintercept = idw.sem.moran$estimate[[1]], colour = 'purple', linewidth = 1), 
    #geom_point(x = dist2, y = sem.moran$estimate[[1]]  , fill = 'black') +
    #geom_point(x = dist2, y = idw.sem.moran$estimate[[1]]  , fill = 'purple'),
    labels = c('Equal weights', 'IDW', 'Comparison', paste0("Moran's I, distance = ", dist2, 'm'))
  )
}    

bl.comp <- comparison_plot(plot.name = 'blaesedalen', 
                data = s2.bl, 
                lw = s2.bl.lw, 
                sem.model = sem.s2.bl.max,
                idw.lw = idw.s2.bl.lw,
                idw.model = idw.sem.s2.bl.max)

kl.comp <- comparison_plot(plot.name = 'kluane-low', 
                data = s2.kl, 
                lw = s2.kl.lw, 
                sem.model = sem.s2.kl.max, 
                idw.lw = idw.s2.kl.lw,
                idw.model = idw.sem.s2.kl.max)

kh.comp <- comparison_plot(plot.name = 'kluane-high', 
                data = s2.kh,
                lw = s2.kh.lw, 
                sem.model = sem.s2.kh.max,
                idw.lw = idw.s2.kh.lw, 
                idw.model = idw.sem.s2.kh.max)


# Save out comparison plots

# Blaesedalen
cowplot::save_plot(paste0('../../plots/nb-dist-sensitivity/idw-neighbour/plots/bl/bl-neighbourhood-invdist-', dist2, '-metres.png'), bl.comp, base_height = 140, base_width = 260, units = 'mm')
# Kluane low
cowplot::save_plot(paste0('../../plots/nb-dist-sensitivity/idw-neighbour/plots/kl/kl-neighbourhood-invdist-', dist2, '-metres.png'), kl.comp, base_height = 140, base_width = 260, units = 'mm')
# Kluane high
cowplot::save_plot(paste0('../../plots/nb-dist-sensitivity/idw-neighbour/plots/kh/kh-neighbourhood-invdist-', dist2, '-metres.png'), kh.comp, base_height = 140, base_width = 260, units = 'mm')
  
# Maps of observed and fitted data

# Blaesedalen
bl.maps <- plot_grid(
  ggplot() +
    geom_sf(data = s2.bl, aes(colour = ndvi.max, size = snow.auc)) +
    #scale_colour_viridis_c(),
    scale_colour_viridis_c(limits = c(0.2, 0.6)),
  ggplot() +
    geom_sf(data = s2.bl %>% mutate(
      fitted = fitted(sem.s2.bl.max)),
      mapping = aes(colour = fitted, size = snow.auc)) +
    #scale_colour_viridis_c(),
    scale_colour_viridis_c(limits = c(0.2, 0.6)),
  ggplot() +
    geom_sf(data = s2.bl %>% mutate(
      fitted = fitted(idw.sem.s2.bl.max)), 
      mapping = aes(colour = fitted, size = snow.auc)) +
    scale_colour_viridis_c(limits = c(0.2, 0.6)),
  labels = c('Observed', paste0(dist2, 'm ', 'Equal weights'), paste0(dist2, 'm ', 'IDW'))
)

cowplot::save_plot(paste0('../../plots/nb-dist-sensitivity/idw-neighbour/plots/bl/maps/bl-neighbourhood-invdist-', dist2, '-metres.png'), bl.maps, base_height = 140, base_width = 260, units = 'mm') 

  # Kluane Low
  kl.maps <- plot_grid(
    ggplot() +
      geom_sf(data = s2.kl, aes(colour = ndvi.max, size = snow.auc)) +
      #scale_colour_viridis_c(),
      scale_colour_viridis_c(limits = c(0.3, 0.6)),
    ggplot() +
      geom_sf(data = s2.kl %>% mutate(
        fitted = fitted(sem.s2.kl.max)),
        mapping = aes(colour = fitted, size = snow.auc)) +
      #scale_colour_viridis_c(),
      scale_colour_viridis_c(limits = c(0.3, 0.6)),
    ggplot() +
      geom_sf(data = s2.kl %>% mutate(
        fitted = fitted(idw.sem.s2.kl.max)), 
        mapping = aes(colour = fitted, size = snow.auc)) +
      scale_colour_viridis_c(limits = c(0.3, 0.6)),
    labels = c('Observed', paste0(dist2, 'm ', 'Equal weights'), paste0(dist2, 'm ', 'IDW'))
  )
  
  cowplot::save_plot(paste0('../../plots/nb-dist-sensitivity/idw-neighbour/plots/kl/maps/kl-neighbourhood-invdist-', dist2, '-metres.png'), kl.maps, base_height = 140, base_width = 260, units = 'mm') 
  
  # Kluane high
  kh.maps <- plot_grid(
    ggplot() +
      geom_sf(data = s2.kh, aes(colour = ndvi.max, size = snow.auc)) +
      #scale_colour_viridis_c(),
      scale_colour_viridis_c(limits = c(0.2, 0.6)),
    ggplot() +
      geom_sf(data = s2.kh %>% mutate(
        fitted = fitted(sem.s2.kh.max)),
        mapping = aes(colour = fitted, size = snow.auc)) +
      #scale_colour_viridis_c(),
      scale_colour_viridis_c(limits = c(0.2, 0.6)),
    ggplot() +
      geom_sf(data = s2.kh %>% mutate(
        fitted = fitted(idw.sem.s2.kh.max)), 
        mapping = aes(colour = fitted, size = snow.auc)) +
      scale_colour_viridis_c(limits = c(0.2, 0.6)),
    labels = c('Observed', paste0(dist2, 'm ', 'Equal weights'), paste0(dist2, 'm ', 'IDW'))
  )

  cowplot::save_plot(paste0('../../plots/nb-dist-sensitivity/idw-neighbour/plots/kh/maps/kh-neighbourhood-invdist-', dist2, '-metres.png'), kh.maps, base_height = 140, base_width = 260, units = 'mm') 
  
  # Investigate change in model summary values over increasing nb size, sem
  # Blaesedalen
  aic.bl <- AIC(sem.s2.bl.max) # aic
  ll.bl <- sem.s2.bl.max$LL[1,1] # log likelihood 
  sig2.bl <- sem.s2.bl.max$s2 # sigma^2  
  wald.bl <- Wald1.Sarlm(sem.s2.bl.max) # Run wald stat
  wald.s.bl <- wald.bl$statistic # Extract statistical value
  wald.p.bl <- wald.bl$p.value # Extract p value
  lr.bl <- LR1.Sarlm(sem.s2.bl.max) # Run likelihood ratio
  lr.s.bl <- lr.bl$statistic[1] # Extract statistic
  lr.p.bl <- lr.bl$p.value[1] # Extract p value
  
  # List the AIC values
  bl.new.row <- list(site = 'bl', weights = 'equal', dist = dist2, aic = aic.bl, loglike = ll.bl, 
                     sig2 = sig2.bl, walds = wald.s.bl, waldp = wald.p.bl, 
                     lograts = lr.s.bl, logratp = lr.p.bl)
  
  # Kluane high
  aic.kh <- AIC(sem.s2.kh.max) # aic
  ll.kh <- sem.s2.kh.max$LL[1,1] # log likelihood 
  sig2.kh <- sem.s2.kh.max$s2 # sigma^2  
  wald.kh <- Wald1.Sarlm(sem.s2.kh.max) # Run wald stat
  wald.s.kh <- wald.kh$statistic # Extract statistical value
  wald.p.kh <- wald.kh$p.value # Extract p value
  lr.kh <- LR1.Sarlm(sem.s2.kh.max) # Run likelihood ratio
  lr.s.kh <- lr.kh$statistic[1] # Extract statistic
  lr.p.kh <- lr.kh$p.value[1] # Extract p value
  
  # List the AIC values
  kh.new.row <- list(site = 'kh', weights = 'equal', dist = dist2, aic = aic.kh, loglike = ll.kh, 
                     sig2 = sig2.kh, walds = wald.s.kh, waldp = wald.p.kh, 
                     lograts = lr.s.kh, logratp = lr.p.kh)
  
  # Kluane low
  aic.kl <- AIC(sem.s2.kl.max) # aic
  ll.kl <- sem.s2.kl.max$LL[1,1] # log likelihood 
  sig2.kl <- sem.s2.kl.max$s2 # sigma^2  
  wald.kl <- Wald1.Sarlm(sem.s2.kl.max) # Run wald stat
  wald.s.kl <- wald.kl$statistic # Extract statistical value
  wald.p.kl <- wald.kl$p.value # Extract p value
  lr.kl <- LR1.Sarlm(sem.s2.kl.max) # Run likelihood ratio
  lr.s.kl <- lr.kl$statistic[1] # Extract statistic
  lr.p.kl <- lr.kl$p.value[1] # Extract p value
  
  # List the AIC values
  kl.new.row <- list(site = 'kl', weights = 'equal', dist = dist2, aic = aic.kl, loglike = ll.kl, 
                     sig2 = sig2.kl, walds = wald.s.kl, waldp = wald.p.kl, 
                     lograts = lr.s.kl, logratp = lr.p.kl)
  
  
  # Investigate change in model summary values over increasing nb size, idw weighting
  # Blaesedalen
  idw.aic.bl <- AIC(idw.sem.s2.bl.max) # aic
  idw.ll.bl <- idw.sem.s2.bl.max$LL[1,1] # log likelihood 
  idw.sig2.bl <- idw.sem.s2.bl.max$s2 # sigma^2  
  idw.wald.bl <- Wald1.Sarlm(idw.sem.s2.bl.max) # Run wald stat
  idw.wald.s.bl <- idw.wald.bl$statistic # Extract statistical value
  idw.wald.p.bl <- idw.wald.bl$p.value # Extract p value
  idw.lr.bl <- LR1.Sarlm(idw.sem.s2.bl.max) # Run likelihood ratio
  idw.lr.s.bl <- idw.lr.bl$statistic[1] # Extract statistic
  idw.lr.p.bl <- idw.lr.bl$p.value[1] # Extract p value
  
  # List the AIC values
  idw.bl.new.row <- list(site = 'bl', weights = 'idw', dist = dist2, aic = idw.aic.bl, loglike = idw.ll.bl, 
                     sig2 = idw.sig2.bl, walds = idw.wald.s.bl, waldp = idw.wald.p.bl, 
                     lograts = idw.lr.s.bl, logratp = idw.lr.p.bl)
  
  # Kluane high
  idw.aic.kh <- AIC(idw.sem.s2.kh.max) # aic
  idw.ll.kh <- idw.sem.s2.kh.max$LL[1,1] # log likelihood 
  idw.sig2.kh <- idw.sem.s2.kh.max$s2 # sigma^2  
  idw.wald.kh <- Wald1.Sarlm(idw.sem.s2.kh.max) # Run wald stat
  idw.wald.s.kh <- idw.wald.kh$statistic # Extract statistical value
  idw.wald.p.kh <- idw.wald.kh$p.value # Extract p value
  idw.lr.kh <- LR1.Sarlm(idw.sem.s2.kh.max) # Run likelihood ratio
  idw.lr.s.kh <- idw.lr.kh$statistic[1] # Extract statistic
  idw.lr.p.kh <- idw.lr.kh$p.value[1] # Extract p value
  
  # List the AIC values
  idw.kh.new.row <- list(site = 'kh', weights = 'idw', dist = dist2, aic = idw.aic.kh, loglike = idw.ll.kh, 
                     sig2 = idw.sig2.kh, walds = idw.wald.s.kh, waldp = idw.wald.p.kh, 
                     lograts = idw.lr.s.kh, logratp = idw.lr.p.kh)
  
  # Kluane low
  idw.aic.kl <- AIC(idw.sem.s2.kl.max) # aic
  idw.ll.kl <- idw.sem.s2.kl.max$LL[1,1] # log likelihood 
  idw.sig2.kl <- idw.sem.s2.kl.max$s2 # sigma^2  
  idw.wald.kl <- Wald1.Sarlm(idw.sem.s2.kl.max) # Run wald stat
  idw.wald.s.kl <- idw.wald.kl$statistic # Extract statistical value
  idw.wald.p.kl <- idw.wald.kl$p.value # Extract p value
  idw.lr.kl <- LR1.Sarlm(idw.sem.s2.kl.max) # Run likelihood ratio
  idw.lr.s.kl <- idw.lr.kl$statistic[1] # Extract statistic
  idw.lr.p.kl <- idw.lr.kl$p.value[1] # Extract p value
  
  # List the AIC values
  idw.kl.new.row <- list(site = 'kl', weights = 'idw', dist = dist2, aic = idw.aic.kl, loglike = idw.ll.kl, 
                     sig2 = idw.sig2.kl, walds = idw.wald.s.kl, waldp = idw.wald.p.kl, 
                     lograts = idw.lr.s.kl, logratp = idw.lr.p.kl)
  # Add to AIC df
  aic.df <- rbind(aic.df, bl.new.row, kh.new.row, kl.new.row, idw.bl.new.row, idw.kh.new.row, idw.kl.new.row)
}