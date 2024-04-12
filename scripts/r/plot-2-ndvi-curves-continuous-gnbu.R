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
# 
# av.snow <- data %>% 
#   filter(row_number() == 7) %>%
#   group_by(ndvi.max.doy) %>%
#   mutate(ndvi.max.doy.mean.snow = mean(snow.auc)) %>%
#   ungroup()
# 
# data <- s2.bl %>%
#   filter(row_number() == 7)
# 
# test <- data %>% group_by(id) %>%
#   filter(doy == ndvi.max.doy)
# 
# mind <- min(test$doy)
# maxd <- max(test$doy)
# 
# doys <- seq(mind, maxd, 1)
# 
# plot <- ggplot() +
#   geom_line(data = data %>% arrange(snow.auc), 
#             aes(x = doy, y = ndvi.pred, group = id, colour = snow.auc), 
#             alpha = 0.1, linewidth = 1) +
#   geom_rug(data = test, sides = "b", inherit.aes = FALSE, position = 'jitter',
#            aes(x = ndvi.max.doy, y = 0, color = ndvi.max.doy.mean.snow, linewidth = 0.1)) +
#   #geom_rug(data = data, sides = "l", inherit.aes = FALSE,
#   #         aes(y = ndvi.max, color = snow.auc, linewidth = 0.1)) +
#   #scale_color_continuous()
#   scale_color_distiller(palette = 'GnBu', type = 'seq', aesthetics =  'colour',
#                         direction = 1) +
#   labs( x = '', 
#         y = '') +
#   scale_x_continuous(breaks = c(200, 250),
#                      labels = c('200', '250')) +
#   # xlim(c(175, 260)) +
#   #  ylim(c(0, 0.6)) +
#   coord_cartesian(xlim = c(180, 250), ylim = c(0, 0.7)) + # 175, 260
#   theme_cowplot() +
#   theme(legend.position = 'none')
# 
# ggplot() +
#   geom_bin2d(data = data, aes(x = ndvi.max.doy, y = ndvi.max, colour = snow.auc), #alpha = ..count..), 
#              binwidth = c(1, 0.01)) +
#   theme_cowplot()
# 
# 
# # Preprocess data: Round ndvi to nearest 0.1 and doy to nearest 1
# data$ndvi_rounded <- as.factor(round(data$ndvi.max * 100) / 100)
# data$ndvi.max.doy <- as.factor(data$ndvi.max.doy)
# 
# # Calculate average snow.auc for each bin
# df_bins <- data %>%
#   group_by(ndvi_rounded, ndvi.max.doy) %>%
#   summarise(avg_snow_auc = mean(snow.auc, na.rm = TRUE))
# 
# # Create the heatmap
#   ggplot(df_bins, aes(x = ndvi.max.doy, y = ndvi_rounded, fill = avg_snow_auc)) +
#   geom_bin2d(binwidth = c(1, 0.01)) +
#   scale_fill_gradientn(colors = c('#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081')) +
#   scale_alpha_continuous(range = c(0.2, 1)) +
#   labs(x = "Day of Year (doy)", y = "NDVI", fill = "Average Snow AUC", alpha = "Count") +
#   theme_cowplot() +
#   theme(legend.position = 'none')


# Plotting ----

# Function for plotting the data
plot_data <- function(data) {

  #data <- data %>%
   # group_by(id) %>%
    #arrange(snow.auc)
  av.snow <- data %>% 
    filter(row_number() == 7) %>%
    group_by(ndvi.max.doy) %>%
    mutate(ndvi.max.doy.mean.snow = mean(snow.auc)) %>%
    ungroup()
  
  # earliest <- data %>% filter(row_number() == 7) %>% 
  #   ungroup() %>%
  #   filter(snow.auc == min(snow.auc)) %>%
  #   sample_n(1)
  # 
  # highest <- data %>% filter(row_number() == 7) %>%
  #   ungroup() %>%
  #   filter(snow.auc == max(snow.auc)) %>%
  #   sample_n(1)
  # 
  # latest <- data %>% filter(row_number() == 7) %>%
  #   ungroup() %>%
  #   filter(ndvi.max.doy == max(ndvi.max.doy))
  
  ggplot() +
    geom_rect(aes(xmin = (min(data$ndvi.max.doy) - 0.5), xmax = (max(data$ndvi.max.doy) + 0.5), 
                  ymin = -0.1, ymax = 0.751), fill = "grey", alpha = 0.2) +
    geom_rug(data = av.snow, sides = "t", inherit.aes = TRUE, length = unit(0.1, 'npc'), #position = 'jitter',
             aes(x = ndvi.max.doy, y = 0.05, color = ndvi.max.doy.mean.snow, linewidth = 0.01, alpha = 0.3)) +
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
  
  # Preprocess data: Round ndvi to nearest 0.1 and doy to nearest 1
  # data$ndvi_rounded <- as.factor(round(data$ndvi.max * 100) / 100)
  # data$ndvi.max.doy <- as.factor(data$ndvi.max.doy)
  
  # Calculate average snow.auc for each bin
  # df_bins <- data %>%
  #   group_by(ndvi_rounded, ndvi.max.doy) %>%
  #   summarise(avg_snow_auc = mean(snow.auc, na.rm = TRUE),
  #             n_peaks = n())
    
  
  # bottom <-   ggplot(df_bins, aes(x = ndvi.max.doy, y = ndvi_rounded, fill = avg_snow_auc)) +
  #   geom_bin2d(binwidth = c(1, 0.01)) +
  #   #geom_point(data = s2.bl, aes(x = ndvi.max.doy, y = ndvi.max, colour = snow.auc), alpha = 0.3, size = 2) + 
  #   scale_fill_gradientn(colors = c('#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081')) +
  #   scale_alpha_continuous(range = c(0.2, 1)) +
  #   #scale_x_discrete(limits = factor(210, 240)) +
  #   #scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
  #   labs(x = "", y = "", fill = "", alpha = "") +
  #   #coord_cartesian(xlim = c(180, 250), ylim = c(0, 0.6)) +
  #   theme_cowplot() +
  #   theme(legend.position = 'none') 
  # 
  # 
  # plot_grid(top, bottom, nrow = 2, align = 'v')
  # 
  # aligned_plots <- align_plots(top, bottom, align = "v", axis = "b")
  # 
  # # Arrange plots
  # plot_grid(aligned_plots[[1]], aligned_plots[[2]], ncol = 1, align = 'v')
}

bl <- plot_data(s2.bl)
#bl
kl <- plot_data(s2.kl)
#kl
kh <- plot_data(s2.kh)
#kh

bl

#logo.plot <- ggdraw() + draw_image(logo, scale = 0.4)

combined.plots <- plot_grid(#logo.plot,
                            bl,
                            kl,
                            kh,
                            ncol = 3,
                            align = 'h',
                            labels = c('(a)', '(b)', '(c)'))


max(s2.kh$ndvi.max.doy) - min(s2.k$ndvi.max.doy)

#combined.plots

# output the plot

cowplot::save_plot('../../plots/figures/colour-experimentation/gnbu-custom-a05-scale-alpha-l1-rug-and-greyarea-horizontal.png', combined.plots, base_height = 80, base_width = 180, units = 'mm', bg = 'white')
