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


# Plot semivariogram ----
# Fit variogram ndvi
vario_ndvi.max <- as_Spatial(s2.kl) %>% as("SpatialPointsDataFrame") %>%
  variogram(ndvi.max.doy ~ 1, data = ., cutoff = 130, width = 10)
vario_ndvi.max_fit <- fit.variogram(vario_ndvi.max, model = vgm(model = "Sph")) 
# Fit variogram residuals
vario_resid <- s2.kl %>% mutate(resid = lm.s2.kl.doy$residuals) %>%
  as_Spatial() %>% 
  as("SpatialPointsDataFrame") %>%
  variogram(resid ~ 1, data = ., cutoff = 130, width = 10)
vario_resid_fit <- fit.variogram(vario_resid, model = vgm(model = "Sph")) 

# Visualise results
plot <- plot_grid(
  ggplot(vario_ndvi.max) +
    geom_point(aes(x = dist, y = gamma)) +
    geom_line(aes(x = dist, y = gamma), 
              data = variogramLine(vario_ndvi.max_fit, 
                                   dist_vector = seq(10,130,10))) +
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
                                   dist_vector = seq(10,130,10))) +
    geom_vline(xintercept = vario_resid_fit$range) +
    annotate("text", x = vario_resid_fit$range, y  = Inf, 
             label = paste0(" range = ", round(vario_resid_fit$range, 1), " m"),
             hjust = 0, vjust = 1.5) +
    scale_x_continuous(limits = c(10, 130), breaks = seq(10,130,10)) +
    labs(x = "lag distance (m)", y = "semivariance (gamma)") +
    theme_cowplot(),
  labels = c("a) KL, ndvi.doy variogram", "b) residual variogram lm(ndvi.max.doy ~ snow.auc)"),
  hjust = 0,
  scale = 0.9,
  nrow = 2,
  ncol = 1
)

plot
cowplot::save_plot('../../plots/semivariograms/KL/kl-ndvi-doy-spherical-semivariogram-130.png', plot, base_height = 140, base_width = 140, units = 'mm', bg = 'white')


# Models ----

###
# Blaesedalen SAC models
###

# Using the Matern semivariogram, we assume a neighbourhood of ~50 m
s2.bl.nb <- dnearneigh(s2.bl, d1 = 0, d2 = 50)
#s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s2.bl.lw <- nb2listw(s2.bl.nb, style = 'W', zero.policy = FALSE)

#s2.bl.lw <- nb2listwdist(s2.bl.nb, s2.bl, type="idw", style="U", 
#                         alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)

# SAC model
bl.sac <- sacsarlm(ndvi.max ~ snow.auc, # ndvi max is a function of snow.auc
          data = s2.bl,        # data is from s2.bl
          listw = s2.bl.lw,    # single set of list weights, because we don't have clear
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
                            data = s2.bl, 
                            listw = s2.bl.lw, 
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

bl.sem.moran <- moran.test(bl.sem.residuals, s2.bl.lw)
bl.sac.moran <- moran.test(bl.sac.rediduals, s2.bl.lw)

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


###
# Kluane Low SAC models
###

# Using the Matern semivariogram, we assume a neighbourhood of ~60 m
s2.kl.nb <- dnearneigh(s2.kl, d1 = 0, d2 = 60)
#s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s2.kl.lw <- nb2listw(s2.kl.nb, style = 'W', zero.policy = FALSE)

#s2.bl.lw <- nb2listwdist(s2.bl.nb, s2.bl, type="idw", style="U", 
#                         alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)

# SAC model
kl.sac <- sacsarlm(ndvi.max ~ snow.auc, # ndvi max is a function of snow.auc
                   data = s2.kl,        # data is from s2.bl
                   listw = s2.kl.lw,    # single set of list weights, because we don't have clear
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

kl.sem <- errorsarlm(ndvi.max ~ snow.auc, 
                     data = s2.kl, 
                     listw = s2.kl.lw, 
                     zero.policy = TRUE)


# For SAC
kl.sac.summary <- tidy(kl.sac) %>% kbl() %>% kable_classic(full_width = F, html_font = "Cambria")

# For SEM
kl.sem.summary <- tidy(kl.sem) %>% kbl() %>% kable_classic(full_width = F, html_font = "Cambria")

kl.sem.summary
kl.sac.summary

# Check for remaining spatial autocorrelation in the residuals
kl.sem.residuals <- residuals(kl.sem)
kl.sac.rediduals <- residuals(kl.sac)

kl.sem.moran <- moran.test(kl.sem.residuals, s2.kl.lw)
kl.sac.moran <- moran.test(kl.sac.rediduals, s2.kl.lw)

print(kl.sem.moran)
print(kl.sac.moran)

# Plot obs and preds
ggplot() +
  geom_point(data = s2.kl, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
  geom_point(data = s2.kl %>% mutate(fitted = fitted(kl.sac)), aes(x = snow.auc, y = fitted), colour = 'blue') +
  geom_point(data = s2.kl %>% mutate(fitted = fitted(kl.sem)), aes(x = snow.auc, y = fitted), colour = 'green')

ggplot(data = s2.bl %>% mutate(sem = fitted(bl.sem),
                               sac = fitted(bl.sac))) +
  geom_point(aes(x = sem, y = ndvi.max), colour = 'red') +
  geom_point(aes(x = sac, y = ndvi.max), colour = 'blue') +
  annotate("path", x = c(0, 0.5), y = c(0, 0.5), color = "red") # Add diagonal line


###
# Kluane High SAC models
###

# Using the Matern semivariogram, we assume a neighbourhood of ~40 m
s2.kh.nb <- dnearneigh(s2.kh, d1 = 0, d2 = 40)
#s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s2.kh.lw <- nb2listw(s2.kh.nb, style = 'W', zero.policy = FALSE)

#s2.bl.lw <- nb2listwdist(s2.bl.nb, s2.bl, type="idw", style="U", 
#                         alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)

# SAC model
kh.sac <- sacsarlm(ndvi.max ~ snow.auc, # ndvi max is a function of snow.auc
                   data = s2.kh,        # data is from s2.bl
                   listw = s2.kh.lw,    # single set of list weights, because we don't have clear
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

kh.sem <- errorsarlm(ndvi.max ~ snow.auc, 
                     data = s2.kh, 
                     listw = s2.kh.lw, 
                     zero.policy = TRUE)


# For SAC
kh.sac.summary <- tidy(kh.sac) %>% kbl() %>% kable_classic(full_width = F, html_font = "Cambria")

# For SEM
kh.sem.summary <- tidy(kh.sem) %>% kbl() %>% kable_classic(full_width = F, html_font = "Cambria")

kh.sem.summary
kh.sac.summary

# Check for remaining spatial autocorrelation in the residuals
kh.sem.residuals <- residuals(kh.sem)
kh.sac.rediduals <- residuals(kh.sac)

kh.sem.moran <- moran.test(kh.sem.residuals, s2.kh.lw)
kh.sac.moran <- moran.test(kh.sac.rediduals, s2.kh.lw)

print(kh.sem.moran)
print(kh.sac.moran)

# Plot obs and preds
ggplot() +
  geom_point(data = s2.kh, aes(x = snow.auc, y = ndvi.max), colour = 'red') +
  geom_point(data = s2.kh %>% mutate(fitted = fitted(kh.sac)), aes(x = snow.auc, y = fitted), colour = 'blue') +
  geom_point(data = s2.kh %>% mutate(fitted = fitted(kh.sem)), aes(x = snow.auc, y = fitted), colour = 'green')

ggplot(data = s2.bl %>% mutate(sem = fitted(bl.sem),
                               sac = fitted(bl.sac))) +
  geom_point(aes(x = sem, y = ndvi.max), colour = 'red') +
  geom_point(aes(x = sac, y = ndvi.max), colour = 'blue') +
  annotate("path", x = c(0, 0.5), y = c(0, 0.5), color = "red") # Add diagonal line



# Output stats
save_html(bl.sac.summary, "../../data/statistical-output/SAC/bl-50-SAC.html")











