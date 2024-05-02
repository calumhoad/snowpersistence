# Plot 2 for paper, showing NDVI curves for each site
# Calum Hoad, 1 Feb 2024

library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(mapview)
library(spdep)
library(spatialreg)
library(stargazer)
library(cowplot)
library(terra)
library(magick)

# Data ----
# Read in the data and reduce to unique id
# Blaesedalen
s2.bl <- read.csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  drop_na(ndvi.max) %>%
  arrange(snow.auc)

# Kluane low
s2.kl <- read.csv('../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  group_by(id) %>%
  drop_na(ndvi.max)

# Kluane high
s2.kh <- read.csv('../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  group_by(id) %>%
  drop_na(ndvi.max)

# Plotting ----

# Function for plotting the data
plot_data <- function(data, rectangle.top, rug.alpha, line.width) {
  
  av.snow <- data %>% 
    filter(row_number() == 7) %>%
    group_by(ndvi.max.doy) %>%
    mutate(ndvi.max.doy.mean.snow = mean(snow.auc)) %>%
    filter(row_number() == 1) %>%
    ungroup()
  
  ggplot() +
    geom_rect(aes(xmin = (min(data$ndvi.max.doy) - 0.5), xmax = (max(data$ndvi.max.doy) + 0.5), 
                  ymin = -0.1, ymax = rectangle.top), fill = "grey", alpha = 0.2) +
    geom_rug(data = av.snow, sides = "t", inherit.aes = TRUE, length = unit(0.1, 'npc'), #position = 'jitter',
             aes(x = ndvi.max.doy, y = 0.05, color = ndvi.max.doy.mean.snow), linewidth = line.width, alpha = rug.alpha) +
    #geom_bar(data = av.snow, aes(x = ndvi.max.doy, y = 0.1, colour = ndvi.max.doy.mean.snow), alpha = 0.3) +
    geom_line(data = data, 
              aes(x = doy, y = ndvi.pred, group = id), colour = 'white',
              linewidth = 5) +
    geom_line(data = data %>% filter(snow.auc == 0), 
              aes(x = doy, y = ndvi.pred, group = id), 
              colour = '#e0f3db', alpha = 0.5, linewidth = 1) +
    geom_line(data = data %>% filter(snow.auc != 0), 
              aes(x = doy, y = ndvi.pred, group = id, colour = snow.auc, alpha = snow.auc),
              linewidth = 1) +
    #geom_segment(aes(x = earliest$ndvi.max.doy, y = earliest$ndvi.max, xend = earliest$ndvi.max.doy, yend = 0.75), color = "black", linetype = "solid") +
    #geom_segment(aes(x = highest$ndvi.max.doy, y = highest$ndvi.max, xend = highest$ndvi.max.doy, yend = 0.75), color = "black", linetype = "solid") + 
    #geom_segment(aes(x = latest$ndvi.max.doy, y = latest$ndvi.max, xend = latest$ndvi.max.doy, yend = 0.8), color = "red", linetype = "solid") +
    scale_color_gradientn(colors = c('#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081')) +
    scale_alpha_continuous(c(0.5, 1)) +# '#f7fcf0','#e0f3db', removed colors
    #scale_color_viridis_c(direction = -1) +
    #scale_color_distiller(palette = 'GnBu', type = 'seq', aesthetics =  'colour',
    #                     direction = 1) +
    labs( x = '', 
          y = '') +
    scale_x_continuous(breaks = c(200, 225, 250),
                       labels = c('200', '225', '250')) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6), 
                       labels = c('0', '0.2', '0.4', '0.6')) +
    # xlim(c(175, 260)) +
    #  ylim(c(0, 0.6)) +
    coord_cartesian(xlim = c(180, 250), ylim = c(0, 0.8)) + # 175, 260
    theme_cowplot() +
    theme(legend.position = 'none')
  
}


bl <- plot_data(s2.bl, line.width = 0.8, rug.alpha = 10, rectangle.top = 0.749)

kl <- plot_data(s2.kl, line.width = 0.8, rug.alpha = 10, rectangle.top = 0.749)
#kl
kh <- plot_data(s2.kh, line.width = 0.8, rug.alpha = 10, rectangle.top = 0.749)
#kh

combined.plots <- plot_grid(#logo.plot,
                            bl,
                            kl,
                            kh,
                            ncol = 3,
                            align = 'h',
                            labels = c('(a)', '(b)', '(c)'))


# output the plot

cowplot::save_plot('../../plots/figures/figure-2-r-output.png', combined.plots, base_height = 80, base_width = 180, units = 'mm', bg = 'white')
