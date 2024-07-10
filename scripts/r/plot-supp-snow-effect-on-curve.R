# Plot supp. snow effect on NDVI curve fit
# Calum Hoad, 9th July 2024

# Import necessary packages ----
library(terra)
library(dplyr)
library(rts)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(purrr)
library(broom)
library(viridis)
library(tidyverse)
library(bfast)
library(phenopix)
library(greenbrown)
library(greenbrown)
library(ggplot2)
library(cowplot)

# Load data
# Blasedalen
s2.bl.smooth <- read_csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv',
                         show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id)

s2.bl.snow <- read.csv('../../data/snow/snow-cover-10m-blaesedalen.csv') %>%
  dplyr::select(id, X2023.07.02, X2023.07.12, X2023.07.18, X2023.07.26)

s2.all <- left_join(s2.bl.smooth, s2.bl.snow, by = 'id')

no.snow <- s2.all %>% filter(X2023.07.02 == 0) # uav dates which never had snow
first.snow <- s2.all %>% filter(X2023.07.02 != 0 & X2023.07.12 == 0) # first uav date did have snow, second didn't
second.snow <- s2.all %>% filter(X2023.07.02 != 0 & X2023.07.12 != 0 & X2023.07.18 == 0) # third date has no snow, first two had snow
fourth.snow <- s2.all %>% filter(X2023.07.02 != 0 & X2023.07.12 != 0 
                                 & X2023.07.18 != 0 & X2023.07.26 == 0) # fourth date has no snow, first three did
fifth.snow <- s2.all %>% filter(X2023.07.02 != 0 & X2023.07.12 != 0 & X2023.07.18 != 0 &
                                X2023.07.26 != 0)
check <- fifth.snow %>% filter(id == 108)
plot_snow_curve <- function(data, uid, snowdate){
  filtered.data <- data %>% filter(id == uid)
  
  plot <- ggplot() +
            geom_point(data = filtered.data, aes(x = doy, y = ndvi), shape = 4, colour = 'black') +
            geom_line(data = filtered.data, aes(x = doy, y = ndvi.pred, group = id)) +
            geom_point(data = filtered.data, aes(x = ndvi.max.doy, y = ndvi.max)) +
            geom_vline(xintercept = yday(snowdate))
  return(plot)
}

plot_snow_curve(fifth.snow, uid = 108, '2023-07-26')
