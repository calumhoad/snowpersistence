# Create supplementary plots to compare smoothed spline, beck model fits
# Calum Hoad, 4th July 2024

# Import libs 
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

## Load data
# Blasedalen
s2.bl.smooth <- read_csv('../../data/ndvi/s2-bl-smooth.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id)
s2.bl.beck <- read_csv('../../data/ndvi/s2-bl-beck.csv',
                         show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id)

# Kluane low
s2.kl.smooth <- read_csv('../../data/ndvi/s2-kl-smooth.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id)
s2.kl.beck <- read_csv('../../data/ndvi/s2-kl-beck.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id)

# Kluane high
s2.kh.smooth <- read_csv('../../data/ndvi/s2-kh-smooth.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id)
s2.kh.beck <- read_csv('../../data/ndvi/s2-kh-beck.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id)

# Blaesedalen, S30
s30.bl.smooth <- read_csv("../../data/ndvi/s30-bl-smooth.csv",
                   show_col_types = FALSE) %>%
  st_as_sf(coords = c("X", "Y"), remove = F, crs = 32621) %>%
  group_by(id)

glimpse(s2.kl.beck)

# Plot out the data for comparison of model fit ----
# 100 random pixels overview
plot.data <- s2.kl.beck
beck <- s2.kh.beck 
smooth <- s2.kh.smooth 

ggplot(
  plot.data %>% filter(id %in% sample(unique(plot.data$id), 100)),
  aes(x = doy, colour = id, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred)) +
  theme_classic() +
  theme(legend.position = "none")

# 9 random pixels in detail
rand_id <- sample(plot.data %>% drop_na() %>% st_drop_geometry() %>% pull(id) %>% unique(), 9)

# For both models
model.comp <- ggplot() +
  # Beck
  geom_line(data = beck %>% filter(id %in% rand_id),
            aes(x = doy, y = ndvi.pred, color = 'green'), colour = 'green', linewidth = 1.2) +
  # Smoothed spline
  geom_line(data = smooth %>% filter(id %in% rand_id),
            aes(x = doy, y = ndvi.pred, color = 'purple'), colour = 'purple', linewidth = 0.8) +
  geom_point(data = smooth %>% filter(id %in% rand_id),
             aes(x = doy, y = ndvi, group = id), colour = 'black', alpha = 1, shape = 4, size = 2) +
  geom_point(data = beck %>% filter(id %in% rand_id),
             aes(x = ndvi.max.doy, y = ndvi.max), size = 2, color = "darkgreen") +
  geom_point(data = smooth %>% filter(id %in% rand_id),
             aes(x = ndvi.max.doy, y = ndvi.max), color = "navy", size = 2) +
  coord_cartesian(ylim = c(0, 0.6)) +
  facet_wrap(~id) +
  ylab('NDVI') +
  xlab('Day of Year') +
  theme_cowplot() +
  theme(legend.position = 'none')

model.comp

unique(s2.kh.smooth$id)
rand_id
# Save plots
cowplot::save_plot('../../plots/supplementary/smooth-beck-comp-kh-s2.png', model.comp, 
                   base_height = 140, base_width = 180, units = 'mm', 
                   bg = 'white')

# For both models
ggplot() +
  # Beck
  geom_point(data = s2.bl.beck %>% filter(id %in% rand_id),
             aes(x = doy, y = ndvi, group = id), colour = 'grey', alpha = 0.5) +
  geom_line(data = s2.bl.beck %>% filter(id %in% rand_id),
            aes(x = doy, y = ndvi.pred, color = 'red')) +
  geom_point(data = s2.bl.beck %>% filter(id %in% rand_id),
             aes(x = ndvi.max.doy, y = ndvi.max), color = "blue") +
  # Smoothed spline
  geom_point(data = s2.bl.smooth %>% filter(id %in% rand_id),
             aes(x = doy, y = ndvi, group = id)) +
  geom_line(data = s2.bl.smooth %>% filter(id %in% rand_id),
            aes(x = doy, y = ndvi.pred, color = 'blue')) +
  geom_point(data = s2.bl.smooth %>% filter(id %in% rand_id),
             aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id, labeller = label_parsed) +  # Parse label
  theme_cowplot() +
  theme(legend.position = "none")  # Remove legend

# Check effect of snow present in pixels on curve fitting ----

# Bring in the snow data
# Snow data
s2.bl.snow <- read.csv('../../data/snow/snow-cover-10m-blaesedalen.csv')

# Join data
s2.bl.joined <- left_join(s2.bl.smooth, s2.bl.snow, by = 'id') %>%
  drop_na(snow.auc) %>% # drop pixel ids for which we have no snow data
  # Create var for latest observed snow cover
  mutate(latest.snow = ifelse(X2023.07.26 != 0,'2023-07-26',
                              ifelse(X2023.07.18 != 0, '2023-07-18',
                                     ifelse(X2023.07.12 != 0, '2023-07-12',
                                            ifelse(X2023.07.02 != 0, '2023-07-02',
                                                   '2023-06-01'))))) # arbitrary date

snow.remains <- s2.bl.snow %>% filter(X2023.07.26 != 0) # 60
nrow(s2.bl.snow)
nrow(snow.remains)/nrow(s2.bl.snow)*100/1

# For one model
high.data <- s2.bl.joined %>% filter(latest.snow == '2023-07-26') # Which model?
low.data <- s2.bl.joined %>% filter(latest.snow == '2023-06-01')
# 9 random pixels in detail
high_rand_id <- sample(high.data %>% st_drop_geometry() %>% pull(id) %>% unique(), 3)
low_rand_id <- sample(low.data %>% st_drop_geometry() %>% pull(id) %>% unique(), 3)

lowplot <- ggplot(
  low.data %>% filter(id %in% low_rand_id),
  aes(x = doy, group = id)
) +
  geom_point(aes(y = ndvi), shape = 4) +
  geom_line(aes(y = ndvi.pred), colour = '#a8ddb5', linewidth = 1) +
  #geom_vline(xintercept = yday(plot.data$latest.snow), aes(color =  scale_color_viridis_c(snow.auc))) +
   #geom_segment(aes(x = yday(latest.snow) - 1, xend = yday(latest.snow) + 1, y = 0, yend = 0.5, color = snow.auc),
    #            size = 1) +
  scale_color_viridis_c() +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id, ncol = 3, nrow = 1) +
  scale_x_continuous(breaks = c(100, 150, 200, 250, 300, 350), 
                    labels = c('100', '150', '200', '250', '300', '350')) +
  ylab('NDVI') +
  xlab('Day of Year') +
  coord_cartesian(ylim = c(0, 0.6), xlim = c(100, 370)) +
  theme_cowplot()

highplot <- ggplot(
  high.data %>% filter(id %in% c(1105, 1196, 327)),
  aes(x = doy, group = id)
) +
  geom_vline(xintercept = yday('2023-07-26'), colour = 'grey', alpha = 0.5, linewidth = 2) +
  geom_point(aes(y = ndvi), shape = 4) +
  geom_line(aes(y = ndvi.pred), colour = '#084081', linewidth = 1) +
  #geom_segment(aes(x = yday(latest.snow) - 1, xend = yday(latest.snow) + 1, y = 0, yend = 0.5, color = snow.auc),
   #            size = 1) +
  scale_color_viridis_c() +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id, ncol = 3, nrow = 1) +
  scale_x_continuous(breaks = c(100, 150, 200, 250, 300, 350), 
                     labels = c('100', '150', '200', '250', '300', '350')) +
  ylab('NDVI') +
  xlab('Day of Year') +
  coord_cartesian(ylim = c(0, 0.6), xlim = c(100, 370)) +
  theme_cowplot()

combined <- plot_grid(lowplot, highplot, nrow = 2,
                      labels = c('(a)', '(b)'))
combined

cowplot::save_plot('../../plots/supplementary/snow-ndvi-curves.png', combined, 
                   base_width = 180, base_height = 180, units = 'mm',bg = 'white')

# For plotting, run curve fit again, but without the additional autumn synthetic 
# dates

source('a-ndvi-curve-fitting-functions.R')

# Blaesedalen
s2.bl <- read.csv('../../data/ndvi/s2-blaesedalen-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621)
# Kluane-low
s2.kl <- read.csv('../../data/ndvi/s2-kluane-low-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)
# Kluane-high
s2.kh <- read.csv('../../data/ndvi/s2-kluane-high-ndvi-ts-pt.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608)

# NASA HLS S30 ###
s30.bl <- read.csv('../../data/ndvi/s30-blaesedalen-ndvi-ts-pt.csv') %>%
  dplyr::select(-X2023.09.14, -X2023.09.22, -X2023.09.23, -X2023.10.03) %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621)

# Function for formatting the data, 
#   including setting all NDVI values < 0 to 0.
long_ndvi <- function(df) {
  df <- df %>%
    pivot_longer(!geometry & !id, names_to = 'doy', values_to = 'ndvi') %>%
    mutate(doy = sub('X', '', doy), 
           doy = sub('\\.', '-', doy), 
           doy = as_date(doy),
           doy = yday(doy)) %>%
    mutate(ndvi = ifelse(ndvi < 0, 0, ndvi)) %>%
    group_by(id)
}

# For S30, drop all data from Sept onwards due to quality error
s30.bl <- s30.bl %>% 
  select(-X2023.09.14, -X2023.09.22, -X2023.09.23, -X2023.10.03)

# Apply function
s2.bl <- long_ndvi(s2.bl)
s2.kl <- long_ndvi(s2.kl)
s2.kh <- long_ndvi(s2.kh)

s30.bl <- long_ndvi(s30.bl)
# Get rid of synthetic dates from dataset
s2.bl.asym <- s2.bl %>% filter(doy < 300)
s2.kl.asym <- s2.kl %>% filter(doy < 300)
s2.kh.asym <- s2.kh %>% filter(doy < 300)
s30.bl.asym <- s30.bl %>% filter(doy < 300)

# Generate curves for asymetric data
s2.bl.asym.curves <- s2.bl.asym %>%
  group_modify(~model_ndvi_smoothedspline(.x))
s2.kl.asym.curves <- s2.kl.asym %>%
  group_modify(~model_ndvi_smoothedspline(.x))
s2.kh.asym.curves <- s2.kh.asym %>%
  group_modify(~model_ndvi_smoothedspline(.x))
s30.bl.asym.curves <- s30.bl.asym %>%
  group_modify(~model_ndvi_smoothedspline(.x))

# Format asymetric data
s2.bl.asym.curves <- s2.bl.asym.curves %>%
  rename(ndvi.pred = 'ndvi.pred.doy.1')
s2.kl.asym.curves <- s2.kl.asym.curves %>%
  rename(ndvi.pred = 'ndvi.pred.doy.1')
s2.kh.asym.curves <- s2.kh.asym.curves %>%
  rename(ndvi.pred = 'ndvi.pred.doy.1')
s30.bl.asym.curves <- s30.bl.asym.curves %>%
  rename(ndvi.pred = 'ndvi.pred.doy.1')

# 9 random pixels in detail
plot.data <- s30.bl
rand_id <- sample(plot.data %>% st_drop_geometry() %>% pull(id) %>% unique(), 9)

# For one model
asym.data <- s30.bl.asym.curves %>% filter(id %in% rand_id)
smooth.data <- s30.bl.smooth %>% filter(id %in% rand_id)

# Which model?

extra.obs.plot <- ggplot() +
  geom_line(data = asym.data, aes(x = doy, y = ndvi.pred, group = id), linewidth = 1.2, colour = 'green') +
  geom_point(data = asym.data, aes(x = ndvi.max.doy, y = ndvi.max, group = id), size = 3, color = "green") +
  geom_point(data = smooth.data, aes(x = doy, y = ndvi, group = id), colour = 'black', shape = 4, alpha = 1) +
  geom_line(data = smooth.data, aes(x = doy, y = ndvi.pred, group = id), linewidth = 0.8, colour = 'purple') +
  geom_point(data = smooth.data, aes(x = ndvi.max.doy, y = ndvi.max, group = id), color = "purple") +
  geom_point(data = smooth.data %>% filter(doy > 300), aes(x = doy, y = ndvi, group = id), colour = 'red', shape = 4) +
  xlab('Day of Year') +
  ylab('NDVI') +
  scale_x_continuous(breaks = c(100, 200, 300),
                     labels = c('100', '200', '300')) +
  coord_cartesian(ylim = c(-0.1, 0.8)) +
  facet_wrap(~id) +
  theme_cowplot()

extra.obs.plot

cowplot::save_plot('../../plots/supplementary/s30-bl-extraobs-plot.png', extra.obs.plot, 
                   base_height = 180, base_width = 180, units = 'mm', 
                   bg = 'white')

# possible.dates <- c(329, 331) BL
# k.possible.dates <- c(315, 324) KL

