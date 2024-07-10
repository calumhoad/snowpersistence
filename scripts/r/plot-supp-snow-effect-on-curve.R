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

# Function to create plots
plot_snow_curve <- function(data, uid, snowdate){
  filtered.data <- data %>% filter(id == uid)
  
  plot <- ggplot() +
            geom_point(data = filtered.data, aes(x = doy, y = ndvi), shape = 4, colour = 'black') +
            geom_line(data = filtered.data, aes(x = doy, y = ndvi.pred, group = id), colour = 'purple', linewidth = 1) +
            geom_point(data = filtered.data, aes(x = ndvi.max.doy, y = ndvi.max)) +
            geom_vline(xintercept = yday(snowdate)) +
            ylab('') +
            xlab('') +
            #coord_cartesian(ylim = c(0.1, 0.4)) +
            theme_cowplot()
  return(plot)
}

# Create plots
low.one <- plot_snow_curve(no.snow, uid = 167, '2023-07-02')
low.two <- plot_snow_curve(no.snow, uid = 837, '2023-07-02')
low.three <- plot_snow_curve(no.snow, uid = 942, '2023-07-02')
medium.one <- plot_snow_curve(second.snow, uid = 502, '2023-07-18')
medium.two <- plot_snow_curve(second.snow, uid = 1124, '2023-07-18')
medium.three <- plot_snow_curve(second.snow, uid = 252, '2023-07-18')
high.one <- plot_snow_curve(fifth.snow, uid = 108, '2023-07-26')
high.two <- plot_snow_curve(fifth.snow, uid = 332, '2023-07-26')
high.three <- plot_snow_curve(fifth.snow, uid = 1104, '2023-07-26')


combined <- plot_grid(low.one, low.two, low.three, 
          #medium.one, medium.two, medium.three, 
          high.one, high.two, high.three, 
          nrow = 2, ncol = 3)

cowplot::save_plot('../../plots/supplementary/snow-ndvi-curves.png', combined, 
                   base_width = 180, base_height = 180, units = 'mm',bg = 'white')
 