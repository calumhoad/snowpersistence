# Plot 5b, S30 data across all plots and models
# Calum Hoad, 1 Feb 2024

library(ggplot2)
library(dplyr)
library(sf)
library(mapview)
library(spdep)
library(spatialreg)
library(stargazer)
library(cowplot)
library(tidyr)

# Data ----

# Blaesedalen, S30
s30.bl <- read.csv('../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  #filter(row_number() == 7) %>%
  ungroup() %>%
  drop_na(ndvi.max)

# Plot 2 comparison

# Function for plotting the data
plot_data <- function(data, snow.greater, snow.lesser) {
  
  # # Get the range of max doy from high and low snow
  # s2.data.low <- data %>% filter(snow.auc < snow.lesser)
  # s2.data.high <- data %>% filter(snow.auc > snow.greater)
  
  # Use quantiles to standardise the data plotted for each site
  s2.data.low <- data %>% filter(snow.auc <= quantile(data$snow.auc, probs = 0.25))
  s2.data.high <- data %>% filter(snow.auc >= quantile(data$snow.auc, probs = 0.75))
  
  # Store min and max values
  high.snow.max <- max(s2.data.high$ndvi.max.doy)
  high.snow.min <- min(s2.data.high$ndvi.max.doy)
  low.snow.max <- max(s2.data.low$ndvi.max.doy)
  low.snow.min <- min(s2.data.low$ndvi.max.doy)
  
  # Store min and max values
  high.snow.max.ndvi <- max(s2.data.high$ndvi.max)
  high.snow.min.ndvi <- min(s2.data.high$ndvi.max)
  low.snow.max.ndvi <- max(s2.data.low$ndvi.max)
  low.snow.min.ndvi <- min(s2.data.low$ndvi.max)
  
  # Plot in the style of the conceptual plot from the beginning of this research proj.
  ggplot() +
    geom_line(data = data %>% filter(snow.auc < snow.lesser), aes(x = doy, y = ndvi.pred, group = id), color = '#9ebc9fff', alpha = 0.4, size = 1) + #size = snow.auc)) +
    geom_line(data = data %>% filter(snow.auc > snow.greater), aes(x = doy, y = ndvi.pred, group = id), color = '#d08c6dff', alpha = 0.4, size = 1) +
    geom_point(data = data %>% filter(snow.auc < snow.lesser), aes(x = ndvi.max.doy, y = ndvi.max), color = '#9ebc9fff', alpha = 0.5, size = 1) +
    geom_point(data = data %>% filter(snow.auc > snow.greater), aes(x = ndvi.max.doy, y = ndvi.max), color = '#d08c6dff', alpha = 0.5, size = 1) +
    geom_segment(aes(x = high.snow.min, xend = high.snow.max, y = 0.07, yend = 0.07), color = '#d08c6dff', size = 1) +
    geom_segment(aes(x = low.snow.min, xend = low.snow.max, y = 0.08, yend = 0.08), color = '#9ebc9fff', size = 1) +
    #geom_segment(aes(y = high.snow.min.ndvi, yend = high.snow.max.ndvi, x = (low.snow.min -10), xend = (low.snow.min -10)), color = '#ffffff', size = 3) +
    #geom_segment(aes(y = low.snow.min.ndvi, yend = low.snow.max.ndvi, x = (low.snow.min - 10), xend = (low.snow.min - 10)), color = '#ffffff', size = 3) +
    #geom_segment(aes(y = high.snow.min.ndvi, yend = high.snow.max.ndvi, x = (low.snow.min -5), xend = (low.snow.min - 5)), color = '#d08c6dff', size = 1) +
    #geom_segment(aes(y = low.snow.min.ndvi, yend = low.snow.max.ndvi, x = (low.snow.min -5), xend = (low.snow.min - 5)), color = '#9ebc9fff', size = 1) +
    
    #scale_color_viridis_c() +
    labs( x = '', 
          y = '') +
    scale_x_continuous(breaks = c(200, 250),
                       labels = c('200', '250')) +
    # xlim(c(175, 260)) +
    #  ylim(c(0, 0.6)) +
    coord_cartesian(xlim = c(180, 250), ylim = c(0, 0.8)) + # 175, 260
    theme_cowplot()
}

plot.2 <- plot_data(s30.bl, snow.greater = 6, snow.lesser = 6)
plot.2


# Plot 3 
s30.bl.p3 <- read.csv('../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  drop_na(ndvi.max)

#s2.bl.nb <- dnearneigh(s2.bl, d1 = 0, d2 = 10)
s30.bl.nb <- poly2nb(st_buffer(s30.bl.p3, dist = 15, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s30.bl.lw <- nb2listw(s30.bl.nb, style = 'W', zero.policy = FALSE)

sem.s30.bl.max <- errorsarlm(ndvi.max ~ snow.auc,
                            data = s30.bl.p3, 
                            listw = s30.bl.lw,
                            zero.policy = FALSE)

sem.s30.bl.doy <- errorsarlm(ndvi.max.doy ~ snow.auc,
                            data = s30.bl.p3,
                            listw = s30.bl.lw, 
                            zero.policy = FALSE)

stargazer(sem.s30.bl.max, sem.s30.bl.doy, type = 'html', 
          out = '../../data/statistical-output/sem-s30-blaesedalen.html')

# Predictions
# Blaesedalen
bl.pred <- data.frame(snow.auc = seq(min(s30.bl.p3$snow.auc), 
                                     max(s30.bl.p3$snow.auc), 
                                     0.5))

bl.pred <- bl.pred %>% 
  mutate(ndvi = predict(sem.s30.bl.max, newdata = bl.pred))# %>%

# Plot 3
bl.coefficients <- sem.s30.bl.max$coefficients
bl.equation <- sprintf("y = %.3f + %.3f * x", bl.coefficients[1], bl.coefficients[2])

plot.3 <- ggplot() +
  geom_point(data = s30.bl.p3, aes(x = snow.auc, y = ndvi.max), alpha = 0.5, color = '#4984BF') +
  geom_line(data = bl.pred, aes(x = snow.auc, y = ndvi), color = 'black', linewidth = 1) +
  annotate("text", x = max(s30.bl$snow.auc), y = 0.9, label = bl.equation, hjust = 1, vjust = 1) +
  xlab('') +
  ylab('') +
  theme_cowplot()

plot.3


# Plot 4

# Predictions 
# Blaesedalen
bl.pred <- data.frame(snow.auc = seq(min(s30.bl.p3$snow.auc), 
                                     max(s30.bl.p3$snow.auc), 
                                     0.5))

bl.pred <- bl.pred %>% 
  mutate(doy = predict(sem.s30.bl.doy, newdata = bl.pred))# %>%


bl.coefficients <- sem.s30.bl.doy$coefficients
bl.equation <- sprintf("y = %.2f + %.2f * x", bl.coefficients[1], bl.coefficients[2])

plot.4 <- ggplot() +
  geom_point(data = s30.bl, aes(x = snow.auc, y = ndvi.max.doy), alpha = 0.5, color = '#4984BF') +
  geom_line(data = bl.pred, aes(x = snow.auc, y = doy), color = 'black', linewidth = 1) +
  annotate("text", x = max(s30.bl$snow.auc), y = max(s30.bl$ndvi.max.doy) - 13, label = bl.equation, hjust = 1, vjust = 1) +
  xlab('') +
  ylab('') +
  theme_cowplot()
plot.4

scatter <- plot_grid(plot.3, 
                     plot.4, 
                     ncol = 1, 
                     align = 'v')
combined <- plot_grid(plot.2, 
                      scatter,
                      #align = 'h',
                      ncol = 2)
combined

cowplot::save_plot('../../plots/figures/figure-5v2.png', combined, base_height = 80, base_width = 180, units = 'mm')
