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

# Functions for filtering the data
# Filter the data to obtain the upper and lower quartile of snow.auc
upper <- function(data) {
  data.no.zero <- data %>% filter(snow.auc != 0)
  upper.quartile <- quantile(data$snow.auc, probs = 0.75)
  data <- data %>% filter(snow.auc >= upper.quartile)
  return(data)
}

lower <- function(data) {
  lower.quartile <- quantile(data$snow.auc, probs = 0.25)
  data <- data %>% filter(snow.auc <= lower.quartile)
  #data <- data %>% filter(snow.auc == 0)
  return(data)
}

# Plotting ----

# Function for plotting the data
plot_data <- function(data, rectangle.top, rug.alpha, line.width) {
  
  av.snow <- data %>% 
    filter(row_number() == 7) %>%
    group_by(ndvi.max.doy) %>%
    mutate(ndvi.max.doy.mean.snow = mean(snow.auc)) %>%
    filter(row_number() == 1) %>%
    ungroup()
  
  # Get standard deviation and mean for peak NDVI and peak NDVI doy in each group
  lower.q.max = mean(lower(data)$ndvi.max)
  lower.q.doy = mean(lower(data)$ndvi.max.doy)
  sd.lower.q.max = sd(lower(data)$ndvi.max)
  sd.lower.q.doy = sd(lower(data)$ndvi.max.doy)
  
  upper.q.max = mean(upper(data)$ndvi.max)
  upper.q.doy = mean(upper(data)$ndvi.max.doy)
  sd.upper.q.max = sd(upper(data)$ndvi.max)
  sd.upper.q.doy = sd(upper(data)$ndvi.max.doy)
  
  
  
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
    annotate("point", x = upper.q.doy, y = upper.q.max, colour = "white", size = 3) +
    annotate("point", x = lower.q.doy, y = lower.q.max, colour = "white", size = 3) +
    geom_errorbar(aes(x = upper.q.doy, y = upper.q.max, ymin = upper.q.max + sd.upper.q.max, ymax = upper.q.max - sd.upper.q.max), linewidth = 2, colour = 'white' ) +
    geom_errorbar(aes(x = lower.q.doy, y = lower.q.max, ymin = lower.q.max + sd.lower.q.max, ymax = lower.q.max - sd.lower.q.max), linewidth = 2, colour = 'white' ) +
    geom_errorbarh(aes(y = upper.q.max, xmin = upper.q.doy + sd.upper.q.doy, xmax = upper.q.doy - sd.upper.q.doy), height = 0.005, linewidth = 2, colour = 'white' ) +
    geom_errorbarh(aes(y = lower.q.max, xmin = lower.q.doy + sd.lower.q.doy, xmax = lower.q.doy - sd.lower.q.doy), height = 0.005, linewidth = 2, colour = 'white' ) +
    annotate("point", x = lower.q.doy, y = lower.q.max, colour = "#a8ddb5", size = 2) +
    annotate("point", x = upper.q.doy, y = upper.q.max, colour = "#084081", size = 2) +
    geom_errorbar(aes(x = upper.q.doy, y = upper.q.max, ymin = upper.q.max + sd.upper.q.max, ymax = upper.q.max - sd.upper.q.max), linewidth = 0.5, colour = 'black' ) +
    geom_errorbar(aes(x = lower.q.doy, y = lower.q.max, ymin = lower.q.max + sd.lower.q.max, ymax = lower.q.max - sd.lower.q.max), linewidth = 0.5, colour = 'black' ) +
    geom_errorbarh(aes(y = upper.q.max, xmin = upper.q.doy + sd.upper.q.doy, xmax = upper.q.doy - sd.upper.q.doy), height = 0.005, linewidth = 0.5, colour = 'black' ) +
    geom_errorbarh(aes(y = lower.q.max, xmin = lower.q.doy + sd.lower.q.doy, xmax = lower.q.doy - sd.lower.q.doy), height = 0.005, linewidth = 0.5, colour = 'black' ) +
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
    scale_x_continuous(breaks = c(175, 200, 225, 250),
                       labels = c('175', '200', '225', '250')) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6), 
                       labels = c('0', '0.2', '0.4', '0.6')) +
    # xlim(c(175, 260)) +
    #  ylim(c(0, 0.6)) +
    coord_cartesian(xlim = c(174, 250), ylim = c(0, 0.8)) + # 175, 260
    theme_cowplot() +
    theme(legend.position = 'none')
  
}


bl <- plot_data(s2.bl, line.width = 0.8, rug.alpha = 10, rectangle.top = 0.749)

kl <- plot_data(s2.kl, line.width = 0.8, rug.alpha = 10, rectangle.top = 0.749)
#kl
kh <- plot_data(s2.kh, line.width = 0.8, rug.alpha = 10, rectangle.top = 0.749)
#kh

combined.plots <- plot_grid(#logo.plot,
                            kl,
                            kh,
                            bl,
                            ncol = 3,
                            align = 'h',
                            labels = c('(a)', '(b)', '(c)'))

combined.plots
# output the plot

cowplot::save_plot('../../plots/figures/figure-2-r-output-klkhbl-y07-x175-lowp25-highp75withzeros.png', combined.plots, base_height = 80, base_width = 180, units = 'mm', bg = 'white')
