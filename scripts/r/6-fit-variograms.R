# R script to fit final INLA models for Calum
# Updated 4th July 2024

# Dependencies
library(tidyverse)
library(sf)
library(ggplot2)
library(cowplot)
library(colorspace)
library(pbapply)
library(broom)
library(INLA)
library(rlang)
library(gstat)

source('b-modelling-helper-functions.R')

# Functions ----
## Helper function to generate INLA grid
get_inla_grid <- function(snow.auc_df, grid.size){
  # Calculate row and col_numbers
  snow.auc_df <- snow.auc_df %>% 
    # Calculate min X and max >
    mutate(min_x = min(.$X),
           max_y = max(.$Y)) %>%
    # Calczulate colum and row numbers for each cell (top-left to bottom-right corner)
    # We shift by one a R is 1- based
    mutate(col_number = 1 + ((X - min_x) / grid.size),
           row_number = 1 + ((max_y - Y) / grid.size)) %>%
    # Arrange dataset accordingly
    arrange(col_number, row_number)
  
  # Derrive spatial index in INLA fashion (running number)
  # Optain number of cols and rows
  n_col = max(snow.auc_df$col_number)
  n_row = max(snow.auc_df$row_number)
  # Expand grid and add index
  complete_grid <- expand.grid(1:n_col, 1:n_row) %>%
    # Rename columns
    select(col_number = Var1, row_number = Var2) %>%
    # Rearrange same as data frame above (top-left to bottom-right)
    arrange(col_number, row_number) %>%
    # Generate index ("node")
    mutate(node = 1:nrow(.))
  
  # Add node to kl_data, skipping those with no data (left_join)
  snow.auc_df <- snow.auc_df %>%
    left_join(complete_grid, by = c("col_number", "row_number"))
  
  # Return updated df
  return(snow.auc_df)
}

# Fit variogram for NDVI max doy
fit_variogram_doy <- function(input.data, var, grid.size, cutoff.dist, model.type,
                          ymin, ymax, yint) {
  
  vario <- as_Spatial(input.data) %>% as("SpatialPointsDataFrame") %>%
    variogram(ndvi.max.doy ~ 1, data = ., cutoff = cutoff.dist, width = grid.size)
  vario.max_fit <- fit.variogram(vario, model = vgm(model = model.type)) 
  
  
  # Visualise results
  vario.plot <- ggplot(data = vario) +
    geom_point(aes(x = dist, y = gamma)) +
    geom_line(aes(x = dist, y = gamma), 
              data = variogramLine(vario.max_fit, 
                                   dist_vector = seq(grid.size,
                                                     cutoff.dist,
                                                     grid.size))) +
    geom_vline(xintercept = vario.max_fit$range) +
    annotate("text", x = vario.max_fit$range, y  = Inf, 
             label = paste0(" range = ", round(vario.max_fit$range, 1), " m"),
             hjust = -0.1, vjust = 1.5, color = 'red') +
    scale_x_continuous(limits = c(grid.size, cutoff.dist +30), breaks = seq(grid.size,
                                                                        cutoff.dist + 30,
                                                                        grid.size*2)) +
    scale_y_continuous(breaks = c(seq(ymin, ymax, yint)), labels = (seq(ymin, ymax, yint))) +
    coord_cartesian(xlim = c(grid.size, cutoff.dist + 30), ylim = c(ymin, ymax+(0.1*ymax))) +
    labs(x = "lag distance (m)", y = "peak NDVI DoY") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(vario.plot)
}

# Fit variogram for NDVI max
fit_variogram_max <- function(input.data, var, grid.size, cutoff.dist, model.type,
                              ymin, ymax, yint) {
  
  vario <- as_Spatial(input.data) %>% as("SpatialPointsDataFrame") %>%
    variogram(ndvi.max ~ 1, data = ., cutoff = cutoff.dist, width = grid.size)
  vario.max_fit <- fit.variogram(vario, model = vgm(model = model.type)) 
  
  
  # Visualise results
  vario.plot <- ggplot(data = vario) +
    geom_point(aes(x = dist, y = gamma)) +
    geom_line(aes(x = dist, y = gamma), 
              data = variogramLine(vario.max_fit, 
                                   dist_vector = seq(grid.size,
                                                     cutoff.dist,
                                                     grid.size))) +
    geom_vline(xintercept = vario.max_fit$range) +
    annotate("text", x = vario.max_fit$range, y  = Inf, 
             label = paste0(" range = ", round(vario.max_fit$range, 1), " m"),
             hjust = -0.1, vjust = 1.5, color = 'red') +
    scale_x_continuous(limits = c(grid.size, cutoff.dist +30), breaks = seq(grid.size,
                                                                            cutoff.dist + 30,
                                                                            grid.size*2)) +
    scale_y_continuous(breaks = c(seq(ymin, ymax, yint)), labels = (seq(ymin, ymax, yint))) +
    coord_cartesian(xlim = c(grid.size, cutoff.dist + 30), ylim = c(ymin, ymax+(0.1*ymax))) +
    labs(x = "lag distance (m)", y = "peak NDVI") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(vario.plot)
}

# Fit variogram for snow persistence
fit_variogram_snow <- function(input.data, var, grid.size, cutoff.dist, model.type,
                              ymin, ymax, yint) {
  
  vario <- as_Spatial(input.data) %>% as("SpatialPointsDataFrame") %>%
    variogram(snow.auc ~ 1, data = ., cutoff = cutoff.dist, width = grid.size)
  vario.max_fit <- fit.variogram(vario, model = vgm(model = model.type)) 
  
  
  # Visualise results
  vario.plot <- ggplot(data = vario) +
    geom_point(aes(x = dist, y = gamma)) +
    geom_line(aes(x = dist, y = gamma), 
              data = variogramLine(vario.max_fit, 
                                   dist_vector = seq(grid.size,
                                                     cutoff.dist,
                                                     grid.size))) +
    geom_vline(xintercept = vario.max_fit$range) +
    annotate("text", x = vario.max_fit$range, y  = Inf, 
             label = paste0(" range = ", round(vario.max_fit$range, 1), " m"),
             hjust = -0.1, vjust = 1.5, color = 'red') +
    scale_x_continuous(limits = c(grid.size, cutoff.dist +30), breaks = seq(grid.size,
                                                                            cutoff.dist + 30,
                                                                            grid.size*2)) +
    scale_y_continuous(breaks = c(seq(ymin, ymax, yint)), labels = (seq(ymin, ymax, yint))) +
    coord_cartesian(xlim = c(grid.size, cutoff.dist + 30), ylim = c(ymin, ymax+(0.1*ymax))) +
    labs(x = "lag distance (m)", y = "Snow persistence") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(vario.plot)
}


# Get the diagonal of the datasets
get_diagonal <- function(input.data) {
  minx <- min(input.data$X)
  maxx <- max(input.data$X)
  miny <- min(input.data$Y)
  maxy <- max(input.data$Y)
  
  diag.dist <- sqrt((maxx - minx)^2 + (maxy - miny)^2)
  
  return(diag.dist)
}  

### Data Prep ----

## Load data
# Blasedalen
s2.bl <- read_csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Kluane low
s2.kl <- read_csv('../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Kluane high
s2.kh <- read_csv('../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Blaesedalen, S30
s30.bl <- read_csv('../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv',
                   show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()


# Add INLA grids to dfs
s2.bl <- get_inla_grid(s2.bl, grid.size = 10)
s2.kl <- get_inla_grid(s2.kl, grid.size = 10)
s2.kh <- get_inla_grid(s2.kh, grid.size = 10)
s30.bl <- get_inla_grid(s30.bl, grid.size = 30)


# Get range for NDVI max from matern variogram at each site

# Fit variograms, with cutoff set to 1/3 the diagonal of the dataset
bl.doy <- fit_variogram_doy(s2.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.bl)/3,
                    ymin = 0, ymax = 10, yint = 2) # 31.1
kl.doy <- fit_variogram_doy(s2.kl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kl)/3,
                    ymin = 0, ymax = 6, yint = 2) # 52.2
kh.doy <- fit_variogram_doy(s2.kh, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kh)/3,
                    ymin = 0, ymax = 6, yint = 2) # 30.6
bl.s30.doy <- fit_variogram_doy(s30.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s30.bl)/3,
                    ymin = 0, ymax = 30, yint = 5)# 21.8

# Fit variograms, with cutoff set to 1/3 the diagonal of the dataset
bl.max <- fit_variogram_max(s2.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.bl)/3,
                            ymin = 0, ymax = 0.006, yint = 0.001) # 64.2
kl.max <- fit_variogram_max(s2.kl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kl)/3,
                            ymin = 0, ymax = 0.003, yint = 0.001) # 61.7
kh.max <- fit_variogram_max(s2.kh, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kh)/3,
                            ymin = 0, ymax = 0.002, yint = 0.001) # 36.9
bl.s30.max <- fit_variogram_max(s30.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s30.bl)/3,
                                ymin = 0, ymax = 0.005, yint = 0.001)# 125

# Fit variograms, with cutoff set to 1/3 the diagonal of the dataset
bl.snow <- fit_variogram_snow(s2.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.bl)/3,
                            ymin = 0, ymax = 18, yint = 2) # 18.5
kl.snow <- fit_variogram_snow(s2.kl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kl)/3,
                            ymin = 0, ymax = 2, yint = 0.5) # 13.1
kh.snow <- fit_variogram_snow(s2.kh, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kh)/3,
                            ymin = 0, ymax = 1, yint = 0.2) # 13.5
bl.s30.snow <- fit_variogram_snow(s30.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s30.bl)/3,
                                ymin = 0, ymax = 8, yint = 2)# 19.2

# Put into cowplots
vgrams.max <- plot_grid(bl.max, kl.max, kh.max, bl.s30.max, nrow = 2, ncol = 2, labels = c('(a)', '(b)', '(c)', '(d)'))
vgrams.doy <- plot_grid(bl.doy, kl.doy, kh.doy, bl.s30.doy, nrow = 2, ncol = 2, labels = c('(a)', '(b)', '(c)', '(d)'))
vgrams.snow <- plot_grid(bl.snow, kl.snow, kh.snow, bl.s30.snow, nrow = 2, ncol = 2, labels = c('(a)', '(b)', '(c)', '(d)'))

vgrams.doy
vgrams.max
vgrams.snow

# Output variogram plots
cowplot::save_plot('../../plots/supplementary/variograms-mat-ndvi-max-doy.png', vgrams.doy, 
                   base_height = 140, base_width = 180, units = 'mm', 
                   bg = 'white')
cowplot::save_plot('../../plots/supplementary/variograms-mat-ndvi-max.png', vgrams.max, 
                   base_height = 140, base_width = 180, units = 'mm', 
                   bg = 'white')
cowplot::save_plot('../../plots/supplementary/variograms-snow-persistence.png', vgrams.snow, 
                   base_height = 140, base_width = 180, units = 'mm', 
                   bg = 'white')
