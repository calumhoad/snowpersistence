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
  drop_na(ndvi.max)

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

data <- s2.bl

plot <- ggplot() +
  geom_line(data = data %>% arrange(desc(snow.auc)), 
            aes(x = doy, y = ndvi.pred, group = id, colour = snow.auc), 
            alpha = 0.1, linewidth = 1) +
  #scale_color_continuous()
  scale_color_distiller(palette = 'GnBu', type = 'seq', aesthetics =  'colour',
                        direction = 1) +
  labs( x = '', 
        y = '') +
  scale_x_continuous(breaks = c(200, 250),
                     labels = c('200', '250')) +
  # xlim(c(175, 260)) +
  #  ylim(c(0, 0.6)) +
  coord_cartesian(xlim = c(180, 250), ylim = c(0, 0.7)) + # 175, 260
  theme_cowplot()

plot

# Plotting ----

# Figure out the quantiles

# Function for plotting the data
plot_data <- function(data, snow.greater, snow.lesser) {
  
  # Sort the data
  data <- data %>%
    arrange(desc(snow.auc))
  
  tertile1 <- quantile(data$snow.auc, probs = 0.33)
  tertile2 <- quantile(data$snow.auc, probs = 0.66)
  
  quartile1 <- quantile(data$snow.auc, probs = 0.1)
  quartile3 <- quantile(data$snow.auc, probs = 0.9)
  
  # Use quantiles to standardise the data plotted for each site
  data.low <- data %>% filter(snow.auc <= quartile1)
  data.mid <- s2.bl %>% filter(snow.auc > quartile1) %>%
                        filter(snow.auc < quartile3)
  data.high <- data %>% filter(snow.auc >= quartile3)
  
  lower.mean <- data %>% filter(snow.auc <= quartile1) %>%
    group_by(doy) %>%
    mutate(mean.ndvi = mean(ndvi.pred)) %>%
    ungroup() %>%
    filter(id == first(unique(id)))
  
  upper.mean <- data %>% filter(snow.auc >= quartile3) %>%
    group_by(doy) %>%
    mutate(mean.ndvi = mean(ndvi.pred)) %>%
    ungroup() %>%
    filter(id == first(unique(id)))
  
  central.mean <- data %>% filter(snow.auc <= quartile3) %>%
    filter(snow.auc >= quartile1) %>%
    group_by(doy) %>%
    mutate(mean.ndvi = mean(ndvi.pred)) %>%
    ungroup() %>%
    filter(id == first(unique(id)))
  
  # Store min and max values
  #high.snow.max <- max(s2.data.high$ndvi.max.doy)
  #high.snow.min <- min(s2.data.high$ndvi.max.doy)
  #low.snow.max <- max(s2.data.low$ndvi.max.doy)
  #low.snow.min <- min(s2.data.low$ndvi.max.doy)
  
  # Store min and max values
  #high.snow.max.ndvi <- max(s2.data.high$ndvi.max)
  #high.snow.min.ndvi <- min(s2.data.high$ndvi.max)
  #low.snow.max.ndvi <- max(s2.data.low$ndvi.max)
  #low.snow.min.ndvi <- min(s2.data.low$ndvi.max)

  ggplot() +
    geom_line(data = data %>% arrange(desc(snow.auc)), 
              aes(x = doy, y = ndvi.pred, group = id, colour = snow.auc), 
              alpha = 0.5, linewidth = 1) +
    #scale_color_continuous()
    scale_color_distiller(palette = 'GnBu', type = 'seq', aesthetics =  'colour',
                          direction = 1) +
    labs( x = '', 
          y = '') +
    scale_x_continuous(breaks = c(200, 250),
                       labels = c('200', '250')) +
    # xlim(c(175, 260)) +
    #  ylim(c(0, 0.6)) +
    coord_cartesian(xlim = c(180, 250), ylim = c(0, 0.7)) + # 175, 260
    theme_cowplot()
  
  # Plot in the style of the conceptual plot from the beginning of this research proj.
  # ggplot() +
  #  # geom_line(data = data.low, aes(x = doy, y = ndvi.pred, group = id), colour = "#9ebc9fff", alpha = 0.1, size = 1) + #size = snow.auc)) +
  #   geom_line(data = data.mid, aes(x = doy, y = ndvi.pred, group = id), colour = "grey", alpha = 0.05, size = 1) +
  #   geom_line(data = data.high, aes(x = doy, y = ndvi.pred, group = id), colour = "#d08c6dff", alpha = 0.05, size = 1) +
  #   geom_line(data = data.low, aes(x = doy, y = ndvi.pred, group = id), colour = "#9ebc9fff", alpha = 0.05, size = 1) + 
  #   geom_line(data = lower.mean, aes(x = doy, y = mean.ndvi, group = id), colour = '#9ebc9fff', alpha = 1, size = 2) +
  #   geom_line(data = upper.mean, aes(x = doy, y = mean.ndvi, group = id), colour = "#d08c6dff", alpha = 1, size = 2) +
    #geom_line(data = central.mean, aes(x = doy, y = mean.ndvi, group = id), colour = "grey", alpha = 1, size = 2) +    
           #geom_line(data = data %>% filter(snow.auc > snow.greater), aes(x = doy, y = ndvi.pred, group = id), color = '#d08c6dff', alpha = 0.1, size = 1) +
    #geom_point(data = data, aes(x = ndvi.max.doy, y = ndvi.max, color = snow.auc), alpha = 0.5, size = 1) +
    #geom_point(data = data %>% filter(snow.auc > snow.greater), aes(x = ndvi.max.doy, y = ndvi.max), color = '#d08c6dff', alpha = 0.5, size = 1) +
    #geom_segment(aes(x = high.snow.min, xend = high.snow.max, y = 0.07, yend = 0.07), color = '#d08c6dff', size = 1) +
    #geom_segment(aes(x = low.snow.min, xend = low.snow.max, y = 0.08, yend = 0.08), color = '#9ebc9fff', size = 1) +
    #geom_segment(aes(y = high.snow.min.ndvi, yend = high.snow.max.ndvi, x = (low.snow.min -10), xend = (low.snow.min -10)), color = '#ffffff', size = 3) +
    #geom_segment(aes(y = low.snow.min.ndvi, yend = low.snow.max.ndvi, x = (low.snow.min - 10), xend = (low.snow.min - 10)), color = '#ffffff', size = 3) +
    #geom_segment(aes(y = high.snow.min.ndvi, yend = high.snow.max.ndvi, x = (low.snow.min -5), xend = (low.snow.min - 5)), color = '#d08c6dff', size = 1) +
    #geom_segment(aes(y = low.snow.min.ndvi, yend = low.snow.max.ndvi, x = (low.snow.min -5), xend = (low.snow.min - 5)), color = '#9ebc9fff', size = 1) +

    #scale_color_gradient(high = "#d08c6dff", low  = "#9ebc9fff", name = "Values") +
    # labs( x = '', 
    #       y = '') +
    # scale_x_continuous(breaks = c(200, 250),
    #                  labels = c('200', '250')) +
    # # xlim(c(175, 260)) +
    # #  ylim(c(0, 0.6)) +
    # coord_cartesian(xlim = c(180, 250), ylim = c(0, 0.7)) + # 175, 260
    # theme_cowplot()
}

bl <- plot_data(s2.bl, snow.greater = 10, snow.lesser = 10)
bl
kl <- plot_data(s2.kl, snow.greater = 2, snow.lesser = 2)
kl
kh <- plot_data(s2.kh, snow.greater = 2, snow.lesser = 2)
kh
logo <- ('../../illustration/both-lower-later-03.png')

logo.plot <- ggdraw() + draw_image(logo, scale = 0.4)

combined.plots <- plot_grid(#logo.plot,
                            bl,
                            kl,
                            kh,
                            ncol = 3,
                            align = 'h',
                            labels = c('(a)', '(b)', '(c)'))

combined.plots

# output the plot

cowplot::save_plot('../../plots/figures/figure-2v1.png', combined.plots, base_height = 80, base_width = 180, units = 'mm')
