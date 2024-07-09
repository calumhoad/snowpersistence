# R script to generate binned scatterplots
# Jakob J. Assmann jakob.assmann@uzh.ch 1 May 2024,
# Modified by Calum Hoad, 08/05/2024

# Dependencies
library(tidyverse)
library(sf)
library(ggplot2)
library(cowplot)
library(pbapply)
library(INLA)
library(gt)

### Data Prep ##########

## Load data
# Blasedalen
s2.bl <- read_csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "Sentinel-2 BL", range = round(51.6/10), colour = '#4984BF',
         ymin = 0, ymax = 0.7, 
         dymin = 215, dymax = 250)

# Kluane low
s2.kl <- read_csv('../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "Sentinel-2 KL", range = round(61.9/10), colour = '#F5A40C', 
         ymin = 0.3, ymax = 0.8,
         dymin = 210, dymax = 240)

# Kluane high
s2.kh <- read_csv('../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "Sentinel-2 KH", range = round(38.5/10), colour = '#F23835', 
         ymin = 0.1, ymax = 0.75, 
         dymin = 210, dymax = 240)

# Blaesedalen, S30
s30.bl <- read_csv("../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv",
                   show_col_types = FALSE
) %>%
  st_as_sf(coords = c("X", "Y"), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "HLSS30 BL", range = round(110/30), colour = '#4984BF', 
         ymin = 0.4, ymax = 0.9, 
         dymin = 200, dymax = 250)

# Combine into one data object (keeping order of Calum's plots)
data_list <- list(s2.kl, s2.kh, s2.bl, s30.bl)



# Generate binned scatter plots with zeros in separate categories
hist_plots_max.ndvi <- map(data_list, function(site_data){
  # Set bin width
  bin_width <- 2
  if (unique(site_data$site) == "Sentinel-2 BL") bin_width <- 4
  # Remove geoms for easy handling
  site_data <- st_drop_geometry(site_data)
  # Separate and lable zeros
  zeros <- filter(site_data, snow.auc == 0)
  zeros$snow_auc_cat <- -1
  zeros$n <- nrow(zeros)
  zeros$scale_name <- "[0]"
  # Bin remainder and label
  site_data <- filter(site_data, snow.auc != 0) %>%
    mutate(snow_auc_cat = floor(snow.auc / bin_width))
  site_data <- site_data %>%
    group_by(snow_auc_cat) %>%
    tally() %>%
    mutate(scale_name = ifelse(snow_auc_cat == 1,
                               paste0("[", as.numeric(snow_auc_cat) * bin_width, " - ", (as.numeric(snow_auc_cat) + 1) * bin_width, "]"),
                               paste0("(", as.numeric(snow_auc_cat) * bin_width, " - ", (as.numeric(snow_auc_cat) + 1) * bin_width, "]"))) %>%
    full_join(site_data, .)
  # Merge zeros again and convert factors where needed
  site_data <- bind_rows(site_data, zeros) %>%
    arrange(snow_auc_cat) %>%
    mutate(snow_auc_cat = factor(snow_auc_cat))
  # Plot and return
  ggplot() +
    geom_jitter(aes(x = snow_auc_cat, y = ndvi.max, colour = snow_auc_cat), data = site_data) +
    geom_text(aes(x = snow_auc_cat, y = site_data$ymax[[1]] - 0.03, label = paste0("n = ", n)), vjust = 1.5, size = 4,
              angle = 45,
              data = filter(distinct(site_data, snow_auc_cat, n))) +
    scale_x_discrete(labels = unique(site_data$scale_name)) +
    labs(x = "Snow persistence\n(numerical range)", y = "peak NDVI", 
         title = paste0(unique(site_data$site), " ")) +
    coord_cartesian(ylim = c(site_data$ymin[[1]], site_data$ymax[[1]])) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_cowplot() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))
})

plot_grid(hist_plots_max.ndvi[[1]], hist_plots_max.ndvi[[2]],
          hist_plots_max.ndvi[[3]], hist_plots_max.ndvi[[4]],
          labels = c("(a)", "(b)", "(c)", "(d)"), 
          ncol = 2, nrow = 2) %>%
save_plot("../../plots/supplementary/binned_scatter_max.ndvi.png", .,
          units = 'mm', base_height = 180, base_width = 180, bg = 'white')
 
plot_grid(hist_plots_max.ndvi[1], hist_plots_max.ndvi[2],
          hist_plots_max.ndvi[3], hist_plots_max.ndvi[4],
          ncol = 2, nrow = 2,
          labels = c("(a)", "(b)", "(c)", "(d)"),
          bg = 'white') %>%
  save_plot("../../plots/supplementary/binned_scatter_max.ndvi.png", .,
            units = 'mm', base_height = 140, base_width = 180)

# Same for max.ndvi.doy
hist_plots_max.ndvi.doy <- map(data_list, function(site_data){
  # Set bin width
  bin_width <- 2
  if(unique(site_data$site) == "Sentinel-2 BL") bin_width <- 4
  # Remove geom
  site_data <- st_drop_geometry(site_data)
  # Separate and label zeros
  zeros <- filter(site_data, snow.auc == 0)
  zeros$snow_auc_cat <- -1
  zeros$n <- nrow(zeros)
  zeros$scale_name <- "[0]"
  # Bin remainder of the data and label
  site_data <- filter(site_data, snow.auc != 0) %>%
    mutate(snow_auc_cat = round(snow.auc / bin_width))
  site_data <- site_data %>%
    group_by(snow_auc_cat) %>%
    tally() %>%
    mutate(scale_name = ifelse(snow_auc_cat == 1,
                               paste0("[", as.numeric(snow_auc_cat) * bin_width, " - ", (as.numeric(snow_auc_cat) + 1) * bin_width, "]"),
                               paste0("(", as.numeric(snow_auc_cat) * bin_width, " - ", (as.numeric(snow_auc_cat) + 1) * bin_width, "]"))) %>%
    full_join(site_data, .)
  # Merge again with zeros
  site_data <- bind_rows(site_data, zeros) %>%
    arrange(snow_auc_cat) %>%
    mutate(snow_auc_cat = factor(snow_auc_cat))
  # Plot and return
  ggplot() +
    geom_jitter(aes(x = snow_auc_cat, y = ndvi.max.doy, colour = snow_auc_cat), data = site_data) +
    geom_text(aes(x = snow_auc_cat, y = site_data$dymax[[1]] - 4, label = paste0("n = ", n)), vjust = 1,
              size = 4, angle = 90,
              data = filter(distinct(site_data, snow_auc_cat, n))) +
    scale_x_discrete(labels = unique(site_data$scale_name)) +
    labs(x = "Snow persistence\n(numerical range)", y = "Peak NDVI DoY", 
         title = paste0(unique(site_data$site), "")) +
    coord_cartesian(ylim = c(site_data$dymin[[1]], site_data$dymax[[1]])) +
    theme_cowplot() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))
})

plot_grid(hist_plots_max.ndvi.doy[[1]], hist_plots_max.ndvi.doy[[2]], 
          hist_plots_max.ndvi.doy[[3]], hist_plots_max.ndvi.doy[[4]],
          labels = c("(a)", "(b)", "(c)", "(d)"), nrow = 2, ncol = 2)  %>%
  save_plot("../../plots/supplementary/binned_scatter_max.ndvi.doy.png", ., bg = "white",
            units = 'mm', base_height = 180, base_width = 180)
