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

# Plot out the data for comparison of model fit ----
# 100 random pixels overview
plot.data <- s2.bl.smooth

ggplot(
  plot.data %>% filter(id %in% sample(unique(plot.data$id), 100)),
  aes(x = doy, colour = id, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred)) +
  theme_classic() +
  theme(legend.position = "none")

# 9 random pixels in detail
rand_id <- sample(plot.data %>% st_drop_geometry() %>% pull(id) %>% unique(), 9)

# For one model
plot.data <- s2.bl.smooth # Which model?

ggplot(
  plot.data %>% filter(id %in% rand_id),
  aes(x = doy, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred)) +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id) +
  theme_classic()

# For both models
model.comp <- ggplot() +
  # Beck
  geom_line(data = s2.bl.beck %>% filter(id %in% rand_id),
            aes(x = doy, y = ndvi.pred, color = 'green'), colour = 'green', linewidth = 1.2) +
  # Smoothed spline
  geom_line(data = s2.bl.smooth %>% filter(id %in% rand_id),
            aes(x = doy, y = ndvi.pred, color = 'purple'), colour = 'purple', linewidth = 0.8) +
  geom_point(data = s2.bl.beck %>% filter(id %in% rand_id),
             aes(x = doy, y = ndvi, group = id), colour = 'grey', alpha = 1) +
  geom_point(data = s2.bl.beck %>% filter(id %in% rand_id),
             aes(x = ndvi.max.doy, y = ndvi.max), size = 2, color = "darkgreen") +
  geom_point(data = s2.bl.smooth %>% filter(id %in% rand_id),
             aes(x = ndvi.max.doy, y = ndvi.max), color = "navy") +
  facet_wrap(~id) +
  theme_cowplot() +
  theme(legend.position = 'none')

model.comp

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
plot.data <- s2.bl.joined %>% filter(latest.snow == '2023-07-26') # Which model?

# 9 random pixels in detail
rand_id <- sample(plot.data %>% st_drop_geometry() %>% pull(id) %>% unique(), 9)

ggplot(
  plot.data %>% filter(id %in% rand_id),
  aes(x = doy, group = id)
) +
  geom_point(aes(y = ndvi)) +
  geom_line(aes(y = ndvi.pred)) +
  #geom_vline(xintercept = yday(plot.data$latest.snow), aes(color =  scale_color_viridis_c(snow.auc))) +
  # geom_segment(aes(x = yday(latest.snow) - 1, xend = yday(latest.snow) + 1, y = 0, yend = 0.5, color = snow.auc),
  #              size = 1) +
  scale_color_viridis_c() +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max), color = "red") +
  facet_wrap(~id) +
  theme_classic()

# For plotting, run curve fit again, but without the additional autumn synthetic 
# dates

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
  geom_point(data = smooth.data, aes(x = doy, y = ndvi, group = id), colour = 'grey', alpha = 0.5) +
  geom_line(data = smooth.data, aes(x = doy, y = ndvi.pred, group = id), linewidth = 0.8, colour = 'purple') +
  geom_point(data = smooth.data, aes(x = ndvi.max.doy, y = ndvi.max, group = id), color = "purple") +
  geom_point(data = smooth.data %>% filter(doy > 300), aes(x = doy, y = ndvi, group = id), colour = 'red') +
  facet_wrap(~id) +
  theme_classic()

extra.obs.plot

cowplot::save_plot('../../plots/figures/s2-bl-modelcomp-plot.png', model.comp, 
                   base_height = 180, base_width = 180, units = 'mm', 
                   bg = 'white')

# possible.dates <- c(329, 331) BL
# k.possible.dates <- c(315, 324) KL


