# Running spatial error models, using the spatialreg package
# Calum Hoad, 8th Jan 2024

# Conducting regresssion analysis of snow and vegetation metric data over two 
# focal sites (Blaesedalen, Greenland; Kluane, Yukon), while accounting for
# spatial autocorrelation.

library(ggplot2)
library(sf)
library(mapview)
library(spdep)
library(spatialreg)
library(stargazer)


# Turn off scientific notation
options(scipen = 999)

# Read in the data
s2 <- read.csv('../../data/sentinel-2/output/s2_all_models.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  drop_na() # drop rows where snow persistence is NA (edge of dataset)

# Plot the data
ggplot() + geom_sf(data = s2)

# Create linear models for later testing with Moran's I ----

# Maximum NDVI is a function of snow persistence
lm.ndvi.max_p <- lm(s2$ndvi.max_p ~ s2$snow.persist)
lm.ndvi.max_s <- lm(s2$ndvi.max_s ~ s2$snow.persist)
lm.ndvi.max_b <- lm(s2$ndvi.max_b ~ s2$snow.persist)

# Maximum NDVI day of year is a function of snow persistence 
lm.ndvi.max.doy_b <- lm(s2$ndvi.max.doy_b ~ s2$snow.persist)
lm.ndvi.max.doy_p <- lm(s2$ndvi.max.doy_p ~ s2$snow.persist)
lm.ndvi.max.doy_s <- lm(s2$ndvi.max.doy_s ~ s2$snow.persist)

# Check models
summary(lm.ndvi.max.doy_b)


# Determining neighbours and assigning weights ----

# Note: 
#   Could determine neighbours here through:
#   - Adjacency
#   - Distance
#   - k-nearest
#   It would seem to make sense in this use case to use distance, due to
#   arbitary pixel size?

# Get neighbours for s2 pixels
s2.nb <- spdep::dnearneigh(s2, d1 = 0, d2 = 10)

# Assign weights to the neighbours. Here the weights assigned are equal ('W')
s2.lw <- nb2listw(s2.nb, style = "W", zero.policy = TRUE)


# Moran's Index test ----

# Moran's test, without zero policy as all samples have neighbours
lm.morantest(lm.ndvi.max.doy_s, s2.lw)

# Plot the spatially lagged values using Moran plot
moran.plot(lm.ndvi.max.doy_s$residuals, s2.lw)

# Monte-carlo simulation to assess for significance
moran <- moran.mc(lm.ndvi.max.doy_s$residuals, 
                  s2.lw,
                  nsim = 999)

# Plot the monte-carlo simulation
plot(moran, main = '', las = 1)


# Lagrange multiplier test ----

# Test for whether Spatial error or spatial lag models is more appropriate
lm.LMtests(lm.ndvi.max.doy_s, s2.lw, test = 'LMerr') # Spatial error
lm.LMtests(lm.ndvi.max.doy_s, s2.lw, test = 'LMlag') # Spatial lag


# Sensitivity analysis ----

# Which model to conduct sensitivity analysis on?
model <- lm.ndvi.max

# Empty vector for storage of Moran's I values
moran.I <- c()

# Loop d through a sequence ranging from 10 to 100
for (d in seq(10, 200, 10)) {
  s2.sens.nb <- dnearneigh(s2, d1 = 0, d2 = d)
  s2.sens.lw <- nb2listw(s2.sens.nb, style = 'W', zero.policy = TRUE)
  moran <- moran.mc(model$residuals, s2.sens.lw, 
                    nsim = 999, zero.policy = TRUE)
  moran.I <- c(moran.I, moran$statistic)
}

moran.I <- data.frame(moran = moran.I,
                      distance = seq(10, 200, 10))

ggplot(moran.I, aes(x = distance, y = moran)) +
  geom_point() +
  geom_line()

help(seq)


# Spatial error model ----

# The Lagrange multiplier test does not give a clear preference of 
# spatial lag model or spatial error model. 

# In this case, it contextually makes most sense to use a spatial error model,
# as any autocorrelation is an artefact of no interest to the phenomenon in question.

# Redefine spatial neighbourhoods
s2.nb <- dnearneigh(s2, d1 = 0, d2 = 10)

# Redefine spatial weights for neighbourhoods
s2.lw <- nb2listw(s2.nb, style = 'W', zero.policy = TRUE)

sem.ndvi.max.doy_s <- errorsarlm(ndvi.max.doy_s ~ snow.persist,
                                 data = s2,
                                 listw = s2.lw, 
                                 zero.policy = TRUE)

stargazer(sem.ndvi.max.doy_s, type = 'text')
