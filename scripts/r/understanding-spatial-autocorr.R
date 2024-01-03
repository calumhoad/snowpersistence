# Understanding the spdep and spatialreg packages
# Calum Hoad, 3rd Jan 2024

# Following the tutorial here:
#   https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab15_SpatialRegression.html

library(ggplot2)    # Graphics library
library(sf)         # Spatial data types and handling
library(mapview)    # Visualize spatial data
library(spdep)      # Diagnosing spatial dependence
library(spatialreg) # Spatial lag and spatial error model

# Get the example meuse dataset
library(sp)
data(meuse)

# Create sf from meuse
coords <- cbind(meuse$x, meuse$y)
meuse <- st_as_sf(meuse, coords = c("x", "y"))

# Linear model of the zinc as a function of elevation and distance
zinc.lm <- lm(log(zinc) ~ elev + sqrt(dist), data = meuse)
summary(zinc.lm)

# Get the residuals from the lm
meuse$residuals <- residuals(zinc.lm)
meuse$fitted <- fitted(zinc.lm)

# We can plot the location dataset using geom_sf() instead of geom_point()
ggplot(meuse, aes(col = residuals, size = residuals)) +
  geom_sf() +
  scale_color_gradient2()


# Find groups of neighbours ----
# Note:
#   When finding groups of neighbours in Earth Observation data, the spacing of 
#   the neighbours will be perfectly even (distance between raster cell centres,
#   or polygons. Should consider which method of finding neighbours, and which 
#   parameters are most appropriate.
meuse.nb <- spdep::dnearneigh(meuse, d1 = 0, d2 = 200)

# Assign weights to neighbours
meuse.lw <- nb2listw(meuse.nb, style = "W", zero.policy = TRUE)


# Moran's I test ----

# Run Moran's I on the meuse data
# Note:
#   The zero.policy parameter refers to observations which have no neighbours, 
#   which will never occur in EO data with regular cell spacing.
lm.morantest(zinc.lm, meuse.lw, zero.policy = T)

# Understanding Moran's I

# Using the information on which observation belongs to which neighborhood, 
# we can compute the average residuals of the neighborhood for each sample, 
# and plot them against the (global) residuals of the linear model.
Inc.lag <- lag.listw(meuse.lw, meuse$residuals, zero.policy = T)
plot(meuse$residuals, Inc.lag)

# Alternatively, can use the moran plot function
moran.plot(meuse$residuals, meuse.lw, zero.policy = TRUE)

# Moran's I is the slope of the line between the spatially lagged and observed variables
lm(Inc.lag ~ meuse$residuals)

# Monte-carlo simulation allows testing for significance of Moran's I
moran <- moran.mc(meuse$residuals, meuse.lw, nsim = 999, zero.policy = TRUE)
moran

# Plot Monte-carlo simulation
plot(moran, main="", las=1)


# Lagrange multiplier test ----

# The lagrange multiplier test indicates whether spatial lag or spatial error
# modelling is more appropriate to the given dataset

# Run the lagrange multiplier test
lm.LMtests(zinc.lm, meuse.lw, test="LMerr", zero.policy = T) # spatial error
lm.LMtests(zinc.lm, meuse.lw, test="LMlag", zero.policy = T) # spatial lag model


# Sensitivity analysis ----

# Spatial autocorrelation is scale dependent (size of neighbourhood affects result), 
# so should conduct a sensitivity analysis to determine the impact of neighbourhood 
# selection on the outcome. 

# empty vector for storage of Moran's I values
moran_I <- c()

# loop d through a sequence ranging from 50 to 2000
for (d in seq(50, 2000, 50)) {
  meuse.nb <- dnearneigh(meuse, d1 = 0, d2 = d)
  meuse.lw <- nb2listw(meuse.nb, style = "W", zero.policy = TRUE)
  moran <- moran.mc(zinc.lm$residuals, meuse.lw, nsim = 999, zero.policy = TRUE)
  moran_I <- c(moran_I, moran$statistic)
}

moran_I <- data.frame(moran = moran_I, 
                      distance = seq(50, 2000, 50))

ggplot(moran_I, aes(x = distance, y = moran)) + 
  geom_point() +
  geom_line()


# Spatial lag and spatial error models ----

# fit spatial lag moddel
zinc.slm <- lagsarlm(log(zinc) ~ elev + sqrt(dist), data = meuse, listw = meuse.lw, zero.policy = TRUE)

# fit spatial error model
zinc.sem <- errorsarlm(log(zinc) ~ elev + sqrt(dist), data = meuse, listw = meuse.lw, zero.policy = TRUE)

# Model comparison (simple)
AIC(zinc.lm, zinc.slm, zinc.sem)


# Model comparison, taking into account the scale sensitivity of the neighbourhood

# This step hasn't worked properly, why? ###

# create empty vectors for holding the Moran statistics
moran_I_lm <- c()
moran_I_slm <- c()
moran_I_sem <- c()

# loop through a distance vector d ranging from 50 to 2000 m
for (d in seq(50, 2000, 50)) {
  meuse.nb <- dnearneigh(meuse, d1 = 0, d2 = d)
  meuse.lw <- nb2listw(meuse.nb, style = "W", zero.policy = TRUE)
  
  moran_lm <- moran.mc(zinc.lm$residuals, meuse.lw, nsim = 999, zero.policy = TRUE)
  moran_I_lm <- c(moran_I_lm, moran_lm$statistic)
  
  moran_slm <- moran.mc(zinc.slm$residuals, meuse.lw, nsim = 999, zero.policy = TRUE)
  moran_I_slm <- c(moran_I_slm, moran_slm$statistic)
  
  moran_sem <- moran.mc(zinc.sem$residuals, meuse.lw, nsim = 999, zero.policy = TRUE)
  moran_I_sem <- c(moran_I_sem, moran_sem$statistic)
}


moran_I_lm  <- data.frame(moran_I = moran_I_lm, 
                          distance = seq(50, 2000, 50), 
                          model="Simple linear model")

moran_I_slm <- data.frame(moran_I = moran_I_slm, 
                          distance = seq(50, 2000, 50), 
                          model="Spatial lag model")

moran_I_sem <- data.frame(moran_I = moran_I_sem, 
                          distance = seq(50, 2000, 50), 
                          model="Spatial error model")

moran.df <- rbind(moran_I_lm, moran_I_slm, moran_I_sem)

ggplot(moran.df, aes(x = distance, y = moran_I, col=model)) + 
  geom_point() +
  geom_line()
