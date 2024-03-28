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

av.snow <- data %>% 
  filter(row_number() == 7) %>%
  group_by(ndvi.max.doy) %>%
  mutate(ndvi.max.doy.mean.snow = mean(snow.auc)) %>%
  ungroup()

data <- s2.bl

test <- data %>% group_by(id) %>%
  filter(doy == ndvi.max.doy)

mind <- min(test$doy)
maxd <- max(test$doy)

doys <- seq(mind, maxd, 1)

for (i in doys) {
  data <- test %>% filter(doy == i)
  print(paste0('DOY:', i, 'average is:', ))
}

plot <- ggplot() +
  geom_line(data = data %>% arrange(snow.auc), 
            aes(x = doy, y = ndvi.pred, group = id, colour = snow.auc), 
            alpha = 0.1, linewidth = 1) +
  geom_rug(data = test, sides = "b", inherit.aes = FALSE, position = 'jitter',
           aes(x = ndvi.max.doy, y = 0, color = ndvi.max.doy.mean.snow, linewidth = 0.1)) +
  #geom_rug(data = data, sides = "l", inherit.aes = FALSE,
  #         aes(y = ndvi.max, color = snow.auc, linewidth = 0.1)) +
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
  theme_cowplot() +
  theme(legend.position = 'none')


plot

# Plotting ----

# Figure out the quantiles

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
  
  ggplot() +
    geom_rug(data = av.snow, sides = "t", inherit.aes = FALSE, length = unit(10, 'npc'), #position = 'jitter',
             aes(x = ndvi.max.doy, y = 0.05, color = ndvi.max.doy.mean.snow, linewidth = 0.1)) +
    geom_line(data = s2.bl, 
              aes(x = doy, y = ndvi.pred, group = id), colour = 'white',
              linewidth = 10) +
    geom_line(data = data %>% filter(snow.auc == 0), 
              aes(x = doy, y = ndvi.pred, group = id), 
              colour = '#e0f3db', alpha = 0.5, linewidth = 1) +
    geom_line(data = data %>% filter(snow.auc != 0), 
              aes(x = doy, y = ndvi.pred, group = id, colour = snow.auc, alpha = snow.auc),
              linewidth = 1) +
    scale_color_gradientn(colors = c('#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081')) +
    scale_alpha_continuous(c(0.5, 1)) +# '#f7fcf0','#e0f3db', removed colors
    #scale_color_viridis_c(direction = -1) +
    #scale_color_distiller(palette = 'GnBu', type = 'seq', aesthetics =  'colour',
    #                     direction = 1) +
    labs( x = '', 
          y = '') +
    scale_x_continuous(breaks = c(200, 250),
                       labels = c('200', '250')) +
    # xlim(c(175, 260)) +
    #  ylim(c(0, 0.6)) +
    coord_cartesian(xlim = c(180, 250), ylim = c(0, 0.7)) + # 175, 260
    theme_cowplot() +
    theme(legend.position = 'none')
  
}

bl <- plot_data(s2.bl)
#bl
kl <- plot_data(s2.kl)
#kl
kh <- plot_data(s2.kh)
#kh
#logo <- ('../../illustration/both-lower-later-03.png')

#logo.plot <- ggdraw() + draw_image(logo, scale = 0.4)

combined.plots <- plot_grid(#logo.plot,
                            bl,
                            kl,
                            kh,
                            ncol = 3,
                            align = 'h',
                            labels = c('(a)', '(b)', '(c)'))

#combined.plots

# output the plot

cowplot::save_plot('../../plots/figures/colour-experimentation/gnbu-custom-a05-scale-alpha-l1-withrug.png', combined.plots, base_height = 80, base_width = 180, units = 'mm', bg = 'white')
