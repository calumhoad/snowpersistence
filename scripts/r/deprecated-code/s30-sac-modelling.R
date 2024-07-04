# Plot 3 for paper - NDVI max is influenced by snow persistence
# Calum Hoad, 1 Feb 2024

# This script specifically examines whether a combined spatial autoregressive model
#   (SAC, Durbin) better accounts for the spatial autocorrelation within the data.

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
library(sp)
library(gstat)
library(spdep)
library(spatialreg)
library(ggplot2)
library(cowplot)
library(colorspace)
library(pbapply)
library(broom)
library(kableExtra)
library(htmltools)

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

s30.bl <- read.csv('../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungrou

# Create linear models and test using Moran's I ----

# max.ndvi ~ snow.auc
lm.s30.bl.max <- lm(s30.bl$ndvi.max ~ s30.bl$snow.auc)
# max.ndvi.doy ~ snow.auc
lm.s30.bl.doy <- lm(s30.bl$ndvi.max.doy ~ s30.bl$snow.auc)

# Run Moran's I on the residuals of the linear models
# NDVI.max
lm.max.mt <- moran.test(residuals(lm.s30.bl.max), s30.bl.lw)
print(lm.max.mt)
# NDVI.max.doy
lm.doy.mt <- moran.test(residuals(lm.s30.bl.doy), s30.bl.lw)
print(lm.doy.mt)

# SAC Models, NDVI.max ----

###
# Blaesedalen SAC models
###

# Using the Matern semivariogram, we assume a neighbourhood of ~50 m
s30.bl.nb <- dnearneigh(s30.bl, d1 = 0, d2 = 50)
#s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s30.bl.lw <- nb2listw(s30.bl.nb, style = 'W', zero.policy = FALSE)

#s2.bl.lw <- nb2listwdist(s2.bl.nb, s2.bl, type="idw", style="U", 
#                         alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)

# SAC model
bl.sac <- sacsarlm(ndvi.max ~ snow.auc, # ndvi max is a function of snow.auc
          data = s30.bl,        # data is from s2.bl
          listw = s30.bl.lw,    # single set of list weights, because we don't have clear
          listw2 = NULL,        # enough documentation or justification to use two list weight ojects: this would only overcomplicate the model.
          Durbin = TRUE, 
          type = 'sac',
          method="eigen", 
          quiet=NULL, 
          zero.policy=NULL, 
          tol.solve=.Machine$double.eps,
          llprof=NULL, 
          interval1=NULL, interval2=NULL, 
          trs1=NULL, trs2=NULL,
          control = list())

bl.sem <- errorsarlm(ndvi.max ~ snow.auc, 
                            data = s30.bl, 
                            listw = s30.bl.lw, 
                            zero.policy = TRUE)


# For SAC
bl.sac.summary <- tidy(bl.sac) %>% kbl() %>% kable_classic(full_width = F, html_font = "Cambria")

# For SEM
bl.sem.summary <- tidy(bl.sem) %>% kbl() %>% kable_classic(full_width = F, html_font = "Cambria")

bl.sem.summary
bl.sac.summary

# Check for remaining spatial autocorrelation in the residuals
bl.sem.residuals <- residuals(bl.sem)
bl.sac.rediduals <- residuals(bl.sac)

bl.sem.moran <- moran.test(bl.sem.residuals, s30.bl.lw)
bl.sac.moran <- moran.test(bl.sac.rediduals, s30.bl.lw)

print(bl.sem.moran)
print(bl.sac.moran)


# Plot obs and preds
ggplot() +
  geom_point(data = s2.bl, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
  geom_point(data = s2.bl %>% mutate(fitted = fitted(bl.sac)), aes(x = snow.auc, y = fitted), colour = 'blue') +
  geom_point(data = s2.bl %>% mutate(fitted = fitted(bl.sem)), aes(x = snow.auc, y = fitted), colour = 'green')

ggplot(data = s2.bl %>% mutate(sem = fitted(bl.sem),
                               sac = fitted(bl.sac))) +
  geom_point(aes(x = sem, y = ndvi.max), colour = 'red') +
  geom_point(aes(x = sac, y = ndvi.max), colour = 'blue') +
  annotate("path", x = c(0, 0.5), y = c(0, 0.5), color = "red") # Add diagonal line


# SAC Models, NDVI.max.doy ----

###
# Blaesedalen SAC models
###

# Using the Matern semivariogram, we assume a neighbourhood of ~50 m
s30.bl.nb <- dnearneigh(s30.bl, d1 = 0, d2 = 40)
#s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s30.bl.lw <- nb2listw(s30.bl.nb, style = 'W', zero.policy = FALSE)

#s2.bl.lw <- nb2listwdist(s2.bl.nb, s2.bl, type="idw", style="U", 
#                         alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)


# SAC model
bl.sac <- sacsarlm(ndvi.max.doy ~ snow.auc, # ndvi max is a function of snow.auc
                   data = s30.bl,        # data is from s2.bl
                   listw = s30.bl.lw,    # single set of list weights, because we don't have clear
                   listw2 = NULL,        # enough documentation or justification to use two list weight ojects: this would only overcomplicate the model.
                   Durbin = TRUE, 
                   type = 'sac',
                   method="eigen", 
                   quiet=NULL, 
                   zero.policy=NULL, 
                   tol.solve=.Machine$double.eps,
                   llprof=NULL, 
                   interval1=NULL, interval2=NULL, 
                   trs1=NULL, trs2=NULL,
                   control = list())

bl.sem <- errorsarlm(ndvi.max.doy ~ snow.auc, 
                     data = s30.bl, 
                     listw = s30.bl.lw, 
                     zero.policy = TRUE)


# For SAC
bl.sac.summary <- tidy(bl.sac) %>% kbl() %>% kable_classic(full_width = F, html_font = "Cambria")

# For SEM
bl.sem.summary <- tidy(bl.sem) %>% kbl() %>% kable_classic(full_width = F, html_font = "Cambria")

bl.sem.summary
bl.sac.summary

# Check for remaining spatial autocorrelation in the residuals
bl.sem.residuals <- residuals(bl.sem)
bl.sac.rediduals <- residuals(bl.sac)

bl.sem.moran <- moran.test(bl.sem.residuals, s30.bl.lw)
bl.sac.moran <- moran.test(bl.sac.rediduals, s30.bl.lw)

print(bl.sem.moran)
print(bl.sac.moran)


# Plot obs and preds
ggplot() +
  geom_point(data = s2.bl, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
  geom_point(data = s2.bl %>% mutate(fitted = fitted(bl.sac)), aes(x = snow.auc, y = fitted), colour = 'blue') +
  geom_point(data = s2.bl %>% mutate(fitted = fitted(bl.sem)), aes(x = snow.auc, y = fitted), colour = 'green')

ggplot(data = s2.bl %>% mutate(sem = fitted(bl.sem),
                               sac = fitted(bl.sac))) +
  geom_point(aes(x = sem, y = ndvi.max), colour = 'red') +
  geom_point(aes(x = sac, y = ndvi.max), colour = 'blue') +
  annotate("path", x = c(0, 0.5), y = c(0, 0.5), color = "red") # Add diagonal line































# Plot semivariogram ----
# Fit variogram ndvi
vario_ndvi.max <- as_Spatial(s30.bl) %>% as("SpatialPointsDataFrame") %>%
  variogram(ndvi.max.doy ~ 1, data = ., cutoff = 130, width = 10)
vario_ndvi.max_fit <- fit.variogram(vario_ndvi.max, model = vgm(model = "Mat")) 
# Fit variogram residuals
vario_resid <- s30.bl %>% mutate(resid = lm.s30.bl.doy$residuals) %>%
  as_Spatial() %>% 
  as("SpatialPointsDataFrame") %>%
  variogram(resid ~ 1, data = ., cutoff = 130, width = 10)
vario_resid_fit <- fit.variogram(vario_resid, model = vgm(model = "Mat")) 

# Visualise results
plot <- plot_grid(
  ggplot(vario_ndvi.max) +
    geom_point(aes(x = dist, y = gamma)) +
    geom_line(aes(x = dist, y = gamma), 
              data = variogramLine(vario_ndvi.max_fit, 
                                   dist_vector = seq(10, 130, 10))) +
    geom_vline(xintercept = vario_ndvi.max_fit$range) +
    annotate("text", x = vario_ndvi.max_fit$range, y  = Inf, 
             label = paste0(" range = ", round(vario_ndvi.max_fit$range, 1), " m"),
             hjust = 0, vjust = 1.5) +
    #annotate("text", x = 30, y = Inf, 
    #         label = "Blaesedalen", 
    #         hjust = 0, vjust = 1.5) +
    #scale_x_continuous(limits = c(10, 200), breaks = seq(10,200,10)) +
    #labs(x = "lag distance (m)", y = "semivariance (gamma)") +
    theme_cowplot(),
  ggplot(vario_resid) +
    geom_point(aes(x = dist, y = gamma)) +
    geom_line(aes(x = dist, y = gamma), 
              data = variogramLine(vario_resid_fit, 
                                   dist_vector = seq(10, 130, 10))) +
    geom_vline(xintercept = vario_resid_fit$range) +
    annotate("text", x = vario_resid_fit$range, y  = Inf, 
             label = paste0(" range = ", round(vario_resid_fit$range, 1), " m"),
             hjust = 0, vjust = 1.5) +
    scale_x_continuous(limits = c(10, 130), breaks = seq(10, 130,10)) +
    labs(x = "lag distance (m)", y = "semivariance (gamma)") +
    theme_cowplot(),
  labels = c("a) S30 BL, ndvi.max.doy variogram", "b) residual variogram lm(ndvi.max.doy ~ snow.auc)"),
  hjust = 0,
  scale = 0.9,
  nrow = 2,
  ncol = 1
)

plot
cowplot::save_plot('../../plots/semivariograms/BL/s30-bl-ndvi-doy-spherical-semivariogram-130.png', plot, base_height = 140, base_width = 140, units = 'mm', bg = 'white')



# Output stats
save_html(bl.sac.summary, "../../data/statistical-output/SAC/bl-50-SAC.html")











