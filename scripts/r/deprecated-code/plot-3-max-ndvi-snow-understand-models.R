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
library(tidyr)
library(tidyverse)

# Data ----
# Read in the data and reduce to unique id
# Blaesedalen
s2.bl <- read.csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

library(mgcv)
b <- gam(ndvi.max ~ s(snow.auc) + s(X,Y,k=400),method="REML",data=s2.bl)
plot(b)
par(mfrow=c(2,2))
gam.check(b)
b

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
  for (d in seq(10, 50, 10)) {
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
var.names <- c('site', 'dist', 'aic', 'loglike', 'sig2', 'walds', 'waldp', 'lograts', 'logratp')
aic.df <- data.frame(matrix(ncol = length(var.names), nrow = 0))
colnames(aic.df) <- var.names

# Set maximum distance of neighbourhood
dist2 <- 60

# For loop to iterate through increasing neighbourhood dist and produce plots
for (dist2 in seq(10, 30, 10)) {

  ###
  # Blaesedalen Spatial Error Models
  ###
  s2.bl.nb <- dnearneigh(s2.bl, d1 = 0, d2 = dist2)
  #s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
  # Redefine spatial weights for neighbourhoods
  #s2.bl.lw <- nb2listw(s2.bl.nb, style = 'W', zero.policy = FALSE)
  
  s2.bl.lw <- nb2listwdist(s2.bl.nb, s2.bl, type="idw", style="raw", 
                           alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)
  
  sem.s2.bl.max <- errorsarlm(ndvi.max ~ snow.auc,
                              data = s2.bl, 
                              listw = s2.bl.lw,
                              zero.policy = FALSE)
  
  
  ###
  # Kluane low Spatial Error Models
  ###
  
  s2.kl.nb <- dnearneigh(s2.kl, d1 = 0, dist2)
  #s2.kl.nb <- poly2nb(st_buffer(s2.kl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
  #s2.kl.lw <- nb2listw(s2.kl.nb, style = 'W', zero.policy = TRUE)
  s2.kl.lw <- nb2listwdist(s2.kl.nb, s2.kl, type="idw", style="raw", 
                           alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)
  
  sem.s2.kl.max <- errorsarlm(ndvi.max ~ snow.auc, 
                              data = s2.kl, 
                              listw = s2.kl.lw, 
                              zero.policy = TRUE)
  
  
  ###
  # Kluane high
  ###
  s2.kh.nb <- dnearneigh(s2.kh, d1 = 0, dist2)
  #s2.kh.nb <- poly2nb(st_buffer(s2.kh, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
  #s2.kh.lw <- nb2listw(s2.kh.nb, style = 'W', zero.policy = TRUE)
  s2.kh.lw <- nb2listwdist(s2.kh.nb, s2.kh, type="idw", style="raw", 
                           alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)
  
  sem.s2.kh.max <- errorsarlm(ndvi.max ~ snow.auc, 
                              data = s2.kh, 
                              listw = s2.kh.lw,
                              zero.policy = TRUE)
  
  
  # Examine model summaries ----
  #stargazer(sem.s2.bl.max, sem.s2.kl.max, sem.s2.kh.max, type = 'html',
  #          out = paste0('../../plots/moran-sensitivity/model-summaries/models-maxndvi-', dist2, '.html'))
  
  
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
  
  cowplot::save_plot(paste0('../../plots/nb-dist-sensitivity/neighbourhood-invdist-', dist2, '-metres.png'), out.plot, base_height = 140, base_width = 260, units = 'mm')
  
  
  # Maps of observed and fitted data
  # Kluane Low
  kl.maps <- plot_grid(
    ggplot() +
      geom_sf(data = s2.kl, aes(colour = ndvi.max, size = snow.auc)) +
      #scale_colour_viridis_c(),
      scale_colour_viridis_c(limits = c(0.1, 0.6)),
    ggplot() +
      geom_sf(data = s2.kl %>% mutate(
        fitted = fitted(sem.s2.kl.max)),
        mapping = aes(colour = fitted, size = snow.auc)) +
      #scale_colour_viridis_c(),
      scale_colour_viridis_c(limits = c(0.1, 0.6)),
    labels = c('KL', dist2)
  )
  
  #cowplot::save_plot(paste0('../../plots/moran-sensitivity/kl-maps/neighbourhood-', dist2, '-metres.png'), kl.maps, base_height = 140, base_width = 260, units = 'mm')
  
  # Kluane high
  kh.maps <- plot_grid(
    ggplot() +
      geom_sf(data = s2.kh, aes(colour = ndvi.max, size = snow.auc)) +
      #scale_colour_viridis_c(),
      scale_colour_viridis_c(limits = c(0.1, 0.6)),
    ggplot() +
      geom_sf(data = s2.kh %>% mutate(
        fitted = fitted(sem.s2.kh.max)),
        mapping = aes(colour = fitted, size = snow.auc)) +
      #scale_colour_viridis_c(),
      scale_colour_viridis_c(limits = c(0.1, 0.6)),
    labels = c('KH', dist2)
  )
  
  #cowplot::save_plot(paste0('../../plots/moran-sensitivity/kh-maps/neighbourhood-', dist2, '-metres.png'), kh.maps, base_height = 140, base_width = 260, units = 'mm')
  
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
      # scale_colour_viridis_c(),
      scale_colour_viridis_c(limits = c(0.1, 0.6)),
    labels = c('BL', dist2)
  )
  
  #cowplot::save_plot(paste0('../../plots/moran-sensitivity/bl-maps/neighbourhood-', dist2, '-metres.png'), bl.maps, base_height = 140, base_width = 260, units = 'mm')
  
  # Investigate change in model summary values over increasing nb size
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
  bl.new.row <- list(site = 'bl', dist = dist2, aic = aic.bl, loglike = ll.bl, 
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
  kh.new.row <- list(site = 'kh', dist = dist2, aic = aic.kh, loglike = ll.kh, 
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
  kl.new.row <- list(site = 'kl', dist = dist2, aic = aic.kl, loglike = ll.kl, 
                     sig2 = sig2.kl, walds = wald.s.kl, waldp = wald.p.kl, 
                     lograts = lr.s.kl, logratp = lr.p.kl)
  # Add to AIC df
  aic.df <- rbind(aic.df, bl.new.row, kh.new.row, kl.new.row)
}

# Save out the dataframe of model summary values
write_csv(aic.df, '../../plots/nb-dist-sensitivity/model-summary-output-ndvi-max-10-400m.csv')

stat.summary <- plot_grid(
  ggplot() +
    geom_line(data = aic.df, aes(x = dist, y = aic, group = site, colour = site)), 
  ggplot() +
    geom_line(data = aic.df, aes(x = dist, y = loglike, group = site, colour = site)),
  ggplot() +
    geom_line(data = aic.df, aes(x = dist, y = sig2, group = site, colour = site)),
  ggplot() +
    geom_line(data = aic.df, aes(x = dist, y = walds, group = site, colour = site)) +
    geom_point(data = aic.df, aes(x = dist, y = walds, colour = site, size = waldp)) +
    scale_size_continuous(breaks = c(0, 0.005, 0.01)),
  ggplot() +
    geom_line(data = aic.df, aes(x = dist, y = lograts, group = site, colour = site)) +
    geom_point(data = aic.df, aes(x = dist, y = lograts, colour = site, size = logratp)) +
    scale_size_continuous(breaks = c(0, 0.005, 0.01)), 
  sens.plot <- ggplot() +
    geom_line(data = sensit.all, aes(x = distance, y = moran, group = site, colour = site)), 
  labels = c('AIC', ' Log Likelihood', 'Sigma^2', 'Wald', 'Log ratio', "Moran's I"), 
  align = 'v'
)

cowplot::save_plot('../../plots/nb-dist-sensitivity/nb-sensitivity-model-summaries.png', stat.summary, base_height = 140, base_width = 260, units = 'mm')


stargazer(sem.s2.bl.max, type = 'text')

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
sens.plot

cowplot::save_plot(paste0('../../plots/moran-sensitivity/sensitivity-plot.png'), sens.plot, base_height = 140, base_width = 260, units = 'mm')

summary <- summary(sem.s2.bl.max)
summary$AIC_lm.model
AIC(sem.s2.bl.max)

br.tidy <- broom::tidy(sem.s2.bl.max)
br.tidy

stargazer(sem.s2.bl.max, type = 'text')
sum <- summary(sem.s2.bl.max)
sem.s2.bl.max$LL[1,1]
sem.s2.bl.max$s2
test <- sum$parameters[4]
test
AIC(sem.s2.bl.max)
LR

wald <- Wald1.Sarlm(sem.s2.bl.max)
wald$statistic
wald$p.value
lr <- LR1.Sarlm(sem.s2.bl.max)
lr$statistic[1]
lr$p.value[1]
