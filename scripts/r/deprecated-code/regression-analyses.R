# Running spatial error models, using the spatialreg package
# Calum Hoad, 8th Jan 2024

# Conducting regresssion analysis of snow and vegetation metric data over two 
# focal sites (Blaesedalen, Greenland; Kluane, Yukon), while accounting for
# spatial autocorrelation.

# Libraries and options ----
library(ggplot2)
library(sf)
library(mapview)
library(spdep)
library(spatialreg)
library(stargazer)
library(tidyr)
library(dplyr)


# Turn off scientific notation
options(scipen = 999)


# Read in the data ----

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
  filter(row_number() ==7) %>%
  ungroup()
  
# ggplot(data = s2.bl, aes(x = snow.auc, y = ndvi.max)) +
#   geom_point(aes(color = snow.auc, size = snow.auc)) +
#   scale_color_viridis_c() +
#   geom_smooth(method = 'lm', aes(x = snow.auc, y = ndvi.max)) +
#   geom_line(data = prediction, aes(x = snow.auc, y = ndvi), color = 'red')


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

# Make a list of the linear model objects
lm.list <- list(lm.s2.bl.max, lm.s2.bl.doy,
             lm.s2.kl.max, lm.s2.kl.doy, 
             lm.s2.kh.max, lm.s2.kh.doy)

# Check models
stargazer(lm.s2.bl.max, lm.s2.bl.doy, type = 'text')


# Determining neighbours and assigning weights ----

# Note: 
#   Could determine neighbours here through:
#   - Adjacency
#   - Distance
#   - k-nearest

# Get neighbours for s2 pixels, distancen = 60 m
s2.bl.nb <- spdep::dnearneigh(s2.bl, d1 = 0, d2 = 60)
s2.kl.nb <- spdep::dnearneigh(s2.kl, d1 = 0, d2 = 60)
s2.kh.nb <- spdep::dnearneigh(s2.kh, d1 = 0, d2 = 60)

# Assign weights to the neighbours. Here the weights assigned are equal ('W')
s2.bl.lw <- nb2listw(s2.bl.nb, style = "W", zero.policy = TRUE)
s2.kl.lw <- nb2listw(s2.kl.nb, style = "W", zero.policy = TRUE)
s2.kh.lw <- nb2listw(s2.kh.nb, style = "W", zero.policy = TRUE)

# List of lw objects
lw.list <- list(s2.bl.lw, s2.bl.lw, 
                s2.kl.lw, s2.kl.lw, 
                s2.kh.lw, s2.kh.lw)


# Moran's Index test ----

# Moran's test, without zero policy as all samples have neighbours

# Blaesedalen
mt.lm.s2.bl.max <- lm.morantest(lm.s2.bl.max, s2.bl.lw)
mt.lm.s2.bl.doy <- lm.morantest(lm.s2.bl.doy, s2.bl.lw)
# kluane low
mt.lm.s2.kl.max <- lm.morantest(lm.s2.kl.max, s2.kl.lw)
mt.lm.s2.kl.doy <- lm.morantest(lm.s2.kl.doy, s2.kl.lw)
# Kluane high
mt.lm.s2.kh.max <- lm.morantest(lm.s2.kh.max, s2.kh.lw)
mt.lm.s2.kh.doy <- lm.morantest(lm.s2.kh.doy, s2.kh.lw)
  
# Plot the spatially lagged values using Moran plot
# Blaesedalen
moran.plot(lm.s2.bl.max$residuals, s2.bl.lw)
moran.plot(lm.s2.bl.doy$residuals, s2.bl.lw)
# Kluane low
moran.plot(lm.s2.kl.max$residuals, s2.kl.lw)
moran.plot(lm.s2.kl.doy$residuals, s2.kl.lw)
# Kluane high
moran.plot(lm.s2.kh.max$residuals, s2.kh.lw)
moran.plot(lm.s2.kh.doy$residuals, s2.kh.lw)


# Monte-carlo simulation to assess for significance
# Blaesedalen
mc.lm.s2.bl.max <- moran.mc(lm.s2.bl.max$residuals, 
                            s2.bl.lw,
                            nsim = 999)
mc.lm.s2.bl.doy <- moran.mc(lm.s2.bl.doy$residuals, 
                            s2.bl.lw,
                            nsim = 999)
# Kluane low
mc.lm.s2.kl.max <- moran.mc(lm.s2.kl.max$residuals, 
                            s2.kl.lw,
                            nsim = 999)
mc.lm.s2.kl.doy <- moran.mc(lm.s2.kl.doy$residuals, 
                            s2.kl.lw,
                            nsim = 999)
# Kluane high
mc.lm.s2.kh.max <- moran.mc(lm.s2.kh.max$residuals, 
                            s2.kh.lw,
                            nsim = 999)
mc.lm.s2.kh.doy <- moran.mc(lm.s2.kh.doy$residuals, 
                            s2.kh.lw,
                            nsim = 999)

# Plot the monte-carlo simulation
plot(mc.lm.s2.bl.max, main = '', las = 1)
plot(mc.lm.s2.bl.doy, main = '', las = 1)
plot(mc.lm.s2.kl.max, main = '', las = 1)
plot(mc.lm.s2.kl.doy, main = '', las = 1)
plot(mc.lm.s2.kh.max, main = '', las = 1)
plot(mc.lm.s2.kh.doy, main = '', las = 1)


# Sensitivity analysis ----

# Function to run sensitivity analysis with increasing distance of neighbourhood

# Loop d through a sequence ranging from 10 to 200m to check for sensitivity
sensit_check <- function(model, data) {
  # Empty vector for storage of Moran's I values
  moran.I <- c()
  for (d in seq(10, 200, 10)) {
    s2.sens.nb <- dnearneigh(data, d1 = 0, d2 = d)
    s2.sens.lw <- nb2listw(s2.sens.nb, style = 'W', zero.policy = TRUE)
    moran <- moran.mc(model$residuals, s2.sens.lw, 
                      nsim = 999, zero.policy = TRUE)
    moran.I <- c(moran.I, moran$statistic)
  }

  moran.I <- data.frame(moran = moran.I,
                        distance = seq(10, 200, 10))

  ggplot(data = moran.I, aes(x = distance, y = moran)) +
    geom_point(aes(y = moran)) +
    geom_line(aes(y = moran))
}

# Run function on each model
# Blaesedalen
sensit_check(lm.s2.bl.max, s2.bl)
sensit_check(lm.s2.bl.doy, s2.bl)
# Kluane low
sensit_check(lm.s2.kl.max, s2.kl)
sensit_check(lm.s2.kl.doy, s2.kl)
# Kluane high
sensit_check(lm.s2.kh.max, s2.kh)
sensit_check(lm.s2.kh.doy, s2.kh)


# Lagrange multiplier test, to check whether Spatial Lagg or Spatial Error appropriate ----

# Looking at the p-value of RLMerr and RLMlag for the model with p < 0.01, or lower p. 

# Blaesedalen
bl.max.lmt <- lm.LMtests(lm.s2.bl.max, s2.bl.lw, test = "all") # no indication
bl.doy.lmt <- lm.LMtests(lm.s2.bl.doy, s2.bl.lw, test = "all") # spatial error
# Kluane low
kl.max.lmt <- lm.LMtests(lm.s2.kl.max, s2.kl.lw, test = "all") # spatial error, though both sig.
kl.doy.lmt <- lm.LMtests(lm.s2.kl.doy, s2.kl.lw, test = "all") # spatial error
# Kluane high
kh.max.lmt <- lm.LMtests(lm.s2.kh.max, s2.kh.lw, test = "all") # no indication
kh.doy.lmt <- lm.LMtests(lm.s2.kh.doy, s2.kh.lw, test = "all") # spatial error, though both sig.


# The Lagrange multiplier test indicates for use of spatial error model over 
#   spatial lag model. Though some tests inconclusive, never recommends
#   spatial lag over spatial error.


######################
## GO TO SCRIPTS "PLOT-4-MAX-NDVI-DOY-SNOW.R" and "PLOT-3-MAX-NDVI-SNOW.R"
######################


# Code beyond here is duplicated in plotting scripts, and is more messy here.
# Will tidy in future.


# Spatial error model ----

# Query whether to use dnearneigh or poly2nb(data, queen = TRUE).

###
# Blaesedalen Spatial Error Models
###

# Run this at 60m
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

stargazer(sem.s2.bl.max, sem.s2.bl.doy, type = 'text', title = 'Blaesedalen')#, 
          out = '../../data/statistical-output/sem-blaesedalen.html')

# Generate predicted values based on model
to.predict <- seq(0, 25, 0.5)
to.predict

predicted.values <- data.frame(snow.auc = to.predict)

prediction <- predicted.values %>%
  mutate(ndvi = predict(sem.s2.bl.max, newdata = predicted.values))

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

stargazer(sem.s2.kl.max, sem.s2.kl.doy, type = 'html', title = 'Kluane low',
          out = '../../data/statistical-output/sem-kluane-low.html')

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

stargazer(sem.s2.kh.max, sem.s2.kh.doy, type = 'html', title = 'Kluane high',
          out = '../../data/statistical-output/sem-kluane-high.html')


check.lm <- moran.mc(lm.s2.bl.doy$residuals, s2.bl.lw, nsim = 999, zero.policy = TRUE)
check.sem <- moran.mc(sem.s2.bl.doy$residuals, s2.bl.lw, nsim = 999, zero.policy = TRUE)
plot(check.sem, main = '', las = 1)
