# Get climatology for Kluane from CRU TS data
# Calum Hoad, 01/05/2024

# Libraries
library(ncdf4)
library(terra)
library(stars)
library(ncmeta)
library(tidyverse)
library(raster)
library(lubridate)

# Load the NetCDF files for precip
pre.2001.2010 <- brick("../../data/met-data/CRU/cru_ts4.07.2001.2010.pre.dat.nc",
                       varname = 'pre')

pre.2011.2020 <- brick("../../data/met-data/CRU/cru_ts4.07.2011.2020.pre.dat.nc",
                       varname = 'pre')

# Set the coords to extract data from the cdf files
coords <- data.frame(row.names = c('KL', 'BL'),
           lon = c(-138.40488817583486, -53.45637381046565),
           #lat = c(60.97933154339103, 69.29345972224118)) # Moved point northward, to avoid cell boundary.
           lat = c(61.1, 69.29345972224118))               # see 'check how close...' below

# Check location of coords is logical
plot(pre.2001.2010$X2001.01.16)
points(coords, pch = 16)

# Check how close to the edge of the grid cell Kluane is
kluane.area <- extent(-142, -134, 58, 62)
plot(crop(pre.2001.2010$X2001.01.16, kluane.area))
points(coords[1,])

# Extract the data
kl.2001.2010 <- data.frame(raster::extract(pre.2001.2010, coords[1,]))
kl.2011.2020 <- data.frame(raster::extract(pre.2011.2020, coords[1,]))
bl.2001.2010 <- data.frame(raster::extract(pre.2001.2010, coords[2,]))
bl.2011.2020 <- data.frame(raster::extract(pre.2011.2020, coords[2,]))

# Tidy the data ----
# KLUANE
kl.2001.2010 <- kl.2001.2010 %>%
  pivot_longer(cols = everything(),
               names_to = 'date', 
               values_to = 'pre')
kl.2011.2020 <- kl.2011.2020 %>%
  pivot_longer(cols = everything(), 
               names_to = 'date', 
               values_to = 'pre')

kl.data <- rbind(kl.2001.2010, kl.2011.2020) %>%
  mutate(date = ymd(gsub('X', '', date))) %>%
  mutate(year = year(date),
         month = month(date))

# BLAESEDALEN
bl.2001.2010 <- bl.2001.2010 %>%
  pivot_longer(cols = everything(),
               names_to = 'date', 
               values_to = 'pre')
bl.2011.2020 <- bl.2011.2020 %>%
  pivot_longer(cols = everything(), 
               names_to = 'date', 
               values_to = 'pre')

bl.data <- rbind(bl.2001.2010, bl.2011.2020) %>%
  mutate(date = ymd(gsub('X', '', date))) %>%
  mutate(year = year(date),
         month = month(date))

# Calculate annual precip ----
# KLUANE
kl.data.annual <- kl.data %>%
  group_by(year) %>%
  mutate(annual.pre = sum(pre)) %>%
  filter(row_number() == 1)

kl.annual.av.pre <- sum(kl.data.annual$annual.pre)/nrow(kl.data.annual)
kl.annual.av.pre

# BLAESEDALEN
bl.data.annual <- bl.data %>%
  group_by(year) %>%
  mutate(annual.pre = sum(pre)) %>%
  filter(row_number() ==1)

bl.annual.av.pre <- sum(bl.data.annual$annual.pre)/nrow(bl.data.annual)
bl.annual.av.pre


# Sanity check, plot out the data ----
ggplot() +
  geom_line(data = kl.data.annual, aes(x = year, y = annual.pre), colour = 'red') +
  geom_line(data = bl.data.annual, aes(x = year, y = annual.pre), colour = 'blue') +
  geom_hline(yintercept = kl.annual.av.pre, colour = 'red') +
  geom_hline(yintercept = bl.annual.av.pre, colour = 'blue')


# work out annual average precip in cold months, as proxy for snow ----

# Set 'snowy' months as November to February
month.start <- 11
month.end <- 2

# Blaesedalen, filter data to winter months
bl.winter <- bl.data %>%
  filter(month >= month.start | month <= month.end)

# Calculate annual winter precip
bl.winter.annual <- bl.winter %>%
  group_by(year) %>%
  mutate(annual.pre = sum(pre)) %>%
  filter(row_number() ==1)

# Get average
bl.winter.av.pre <- sum(bl.winter.annual$annual.pre)/nrow(bl.winter.annual)
bl.winter.av.pre

# Get fraction of 2000-2020 average precip which fell in winter
bl.winter.fraction <- (bl.winter.av.pre/bl.annual.av.pre)*100
bl.winter.fraction

# KLUANE, repeat as above
kl.winter <- kl.data %>%
  filter(month >= month.start | month <= month.end)

# Annual winter
kl.winter.annual <- kl.winter %>%
  group_by(year) %>%
  mutate(annual.pre = sum(pre)) %>%
  filter(row_number() == 1)

# Average winter 2000-2020
kl.winter.av.pre <- sum(kl.winter.annual$annual.pre)/nrow(kl.winter.annual)
kl.winter.av.pre

# Fraction of long term average which fell in winter
kl.winter.fraction <- (kl.winter.av.pre/kl.annual.av.pre)*100
kl.winter.fraction

# Check how close to the edge of the grid cell Kluane is
kluane.area <- extent(-142, -134, 58, 62)
plot(crop(pre.2001.2010$X2001.01.16, kluane.area))
points(coords[1,])
