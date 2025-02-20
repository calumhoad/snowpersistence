# Script to generate the bones of plot 5
# Calum Hoad, 30/04/2024

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
library(stargazer)
library(ggnewscale)
library(tidyterra)
library(terra)

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

# Helper function to fit variograms
fit_variogram <- function(input.data, var, grid.size, cutoff.dist, model.type) {
  
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
             hjust = 0.5, vjust = 1.5, color = 'red') +
    scale_x_continuous(limits = c(grid.size, cutoff.dist), breaks = seq(grid.size,
                                                                        cutoff.dist,
                                                                        grid.size)) +
    labs(x = "lag distance (m)", y = "ndvi.max.doy") +
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

# Define helper function to fit matern model for various shape and scale 
# parameters as well as ranges
fit_matern_doy <- function(shape_param, scale_param, range_site, site_data, return_model = F){
  
  # Define max number of rows and cols
  n_col = max(site_data$col_number)
  n_row = max(site_data$row_number)
  
  # Specify priors for hyperparemeter
  log.range = list(initial = log(range_site), fixed = TRUE)
  hyperpar_matern = list(initial = 2, param = c(shape_param, scale_param))
  # Specify formula
  formula_matern <- ndvi.max.doy ~ snow.auc +
    f(node, 
      model = "matern2d", 
      nrow = n_row,
      ncol = n_col,
      hyper = list(range = log.range,
                   prec = hyperpar_matern))
  # Fit model
  model_matern <- inla(formula = formula_matern, 
                       data = site_data,
                       family = "gaussian",
                       control.predictor = list(compute = TRUE),
                       control.compute = list(dic = TRUE, 
                                              waic = TRUE, 
                                              cpo = TRUE, 
                                              return.marginals.predictor = TRUE),
                       keep = TRUE)
  
  # Return model?
  if(return_model) return(model_matern)
  
  # if not return summary
  return(
    data.frame(scale_param = scale_param,
               shape_param = shape_param, 
               intercept = round(model_matern$summary.fixed$mean[1],5),
               slope = round(model_matern$summary.fixed$mean[2],5),
               slope_q025 = round(model_matern$summary.fixed[2,3],5),
               slope_q975 = round(model_matern$summary.fixed[2,5],5),
               range = round(model_matern$summary.hyperpar$mean[3],2), 
               dic = model_matern$dic$dic,
               dic_eff_par = model_matern$dic$p.eff,
               dic_d_p = nrow(site_data)/model_matern$dic$p.eff,
               waic = model_matern$waic$waic,
               waic_eff_par = model_matern$waic$p.eff
    ))}

fit_matern_max <- function(shape_param, scale_param, range_site, site_data, return_model = F){
  
  # Define max number of rows and cols
  n_col = max(site_data$col_number)
  n_row = max(site_data$row_number)
  
  # Specify priors for hyperparemeter
  log.range = list(initial = log(range_site), fixed = TRUE)
  hyperpar_matern = list(initial = 2, param = c(shape_param, scale_param))
  # Specify formula
  formula_matern <- ndvi.max ~ snow.auc +
    f(node, 
      model = "matern2d", 
      nrow = n_row,
      ncol = n_col,
      hyper = list(range = log.range,
                   prec = hyperpar_matern))
  # Fit model
  model_matern <- inla(formula = formula_matern, 
                       data = site_data,
                       family = "gaussian",
                       control.predictor = list(compute = TRUE),
                       control.compute = list(dic = TRUE, 
                                              waic = TRUE, 
                                              cpo = TRUE, 
                                              return.marginals.predictor = TRUE),
                       keep = TRUE)
  
  # Return model?
  if(return_model) return(model_matern)
  
  # if not return summary
  return(
    data.frame(scale_param = scale_param,
               shape_param = shape_param, 
               intercept = round(model_matern$summary.fixed$mean[1],5),
               slope = round(model_matern$summary.fixed$mean[2],5),
               slope_q025 = round(model_matern$summary.fixed[2,3],5),
               slope_q975 = round(model_matern$summary.fixed[2,5],5),
               range = round(model_matern$summary.hyperpar$mean[3],2), 
               dic = model_matern$dic$dic,
               dic_eff_par = model_matern$dic$p.eff,
               dic_d_p = nrow(site_data)/model_matern$dic$p.eff,
               waic = model_matern$waic$waic,
               waic_eff_par = model_matern$waic$p.eff
    ))}

# Helper function to generate plots
plot_results <- function(model_matern, site_data, grid.size){
  # Calculate marginal predictions depending on snow.auc only
  site_data$preds <- model_matern$summary.fixed[1,1] + 
    site_data$snow.auc * model_matern$summary.fixed[2,1]
  
  # Calculate min and max values for 95 CI predictions to plot a ribbon
  preds_credible <- expand.grid(intercept = unlist(model_matern$summary.fixed[1, c(1,3,5)]), 
                                slope = unlist(model_matern$summary.fixed[2, c(1,3,5)])) %>%
    split(1:nrow(.)) %>%
    map(function(parameters) {
      data.frame(snow.auc =  seq((floor(min(site_data$snow.auc) * 10) / 10),
                                 (ceiling(max(site_data$snow.auc) * 10) / 10), 0.01),
                 ndvi_pred = parameters$intercept + seq((floor(min(site_data$snow.auc) * 10) / 10),
                                                        (ceiling(max(site_data$snow.auc) * 10) / 10), 0.01) * parameters$slope[1])
    }) %>% bind_rows() %>%
    group_by(snow.auc) %>% summarise(min_ndvi_pred = min(ndvi_pred), max_ndvi_pred = max(ndvi_pred))
  
  # Add model fitted values and residuals to data frame
  site_data$fitted <- model_matern$summary.fitted.values$mean
  site_data$residuals <- site_data$ndvi.max.doy - site_data$fitted
  
  # Plot model fit
  fit_plot <- ggplot(data = site_data) +
    geom_point(aes(x = snow.auc, y= ndvi.max.doy)) +
    geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                    ymax = max_ndvi_pred),
                data = preds_credible,
                fill = "blue", alpha = 0.3) +
    geom_line(aes(x = snow.auc, y = preds), colour = "blue") +
    theme_cowplot()
  
  # Plot predictions in space
  preds_in_space <- ggplot() +
    geom_sf(data = st_buffer(site_data, grid.size/2, endCapStyle = "SQUARE"), 
            aes(fill = preds)) +
    scale_fill_continuous_sequential("viridis", rev = F) +
    labs(fill = "predictions\nNDVI max") +
    theme_map()
  
  # Residual plotsukc
  resids_hist <- ggplot(site_data) +
    geom_histogram(aes(x = residuals), bins = 50) +
    labs(x = "residuals") +
    theme_cowplot()
  resids_space <- ggplot() +
    geom_sf(data = st_buffer(site_data, grid.size/2, endCapStyle = "SQUARE"), aes(fill = residuals)) +
    scale_fill_continuous_diverging() +
    theme_map()
  
  # plot the whole lot
  print(plot_grid(fit_plot, preds_in_space,
                  resids_hist, resids_space, nrow = 2))
  
  return(NULL)
}

# For generating scatterplots with INLA model fits
fit_plot_doy <- function(site_data,
                     colour.site, 
                     colour.darker, 
                     colour.lighter, 
                     colour.lightest, 
                     ymax, ymin, yint, 
                     xmax, xmin, xint){
  lm_plot <- ggplot(site_data, aes(x = snow.auc, y = ndvi.max.doy)) +
    geom_point(aes(x = snow.auc, y = ndvi.max.doy), colour = colour.lightest) +
    scale_y_continuous(breaks = c(seq(ymin, ymax, yint)), 
                       labels = c(as.character(seq(ymin, ymax, yint)))) +
    scale_x_continuous(breaks = c(seq(xmin, xmax, xint)), 
                       labels = c(as.character(seq(xmin, xmax, xint)))) +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
    theme_cowplot() +
      geom_smooth(method = "lm", formula = y ~ x, 
                  colour = colour.site,
                  fill = colour.site, 
                  linewidth = 1)
    #labs(x = "snow persistence", y = "ndvi.max.doy")#, 
    #title = paste0(unique(site_data$site), " (lm y ~ x)")) 
  return(lm_plot)
}

fit_plot_max <- function(site_data,
                     colour.site, 
                     colour.darker, 
                     colour.lighter, 
                     colour.lightest, 
                     ymax, ymin, yint,
                     xmax, xmin, xint){
  lm_plot <- ggplot(site_data, aes(x = snow.auc, y = ndvi.max)) +
    geom_point(aes(x = snow.auc, y = ndvi.max), colour = colour.lightest) +
    scale_y_continuous(breaks = c(seq(ymin, ymax, yint)), 
                       labels = c(as.character(seq(ymin, ymax, yint)))) +
    scale_x_continuous(breaks = c(seq(xmin, xmax, xint)), 
                       labels = c(as.character(seq(xmin, xmax, xint)))) +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
    theme_cowplot() +
      geom_smooth(method = "lm", formula = y ~ x, 
                  colour = colour.site,
                  fill = colour.site, 
                  linewidth = 1)
    #labs(x = "snow persistence", y = "ndvi.max.doy")#, 
    #title = paste0(unique(site_data$site), " (lm y ~ x)")) 
  return(lm_plot)
}
# Function to generate curves plot, as per fig 2
plot_data_fig2 <- function(data, rectangle.top, rug.alpha, line.width) {
  
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
    annotate("point", x = lower.q.doy, y = lower.q.max, colour = "black", size = 3) +
    #geom_errorbar(aes(x = upper.q.doy, y = upper.q.max, ymin = upper.q.max + sd.upper.q.max, ymax = upper.q.max - sd.upper.q.max), linewidth = 2, colour = 'white' ) +
    #geom_errorbar(aes(x = lower.q.doy, y = lower.q.max, ymin = lower.q.max + sd.lower.q.max, ymax = lower.q.max - sd.lower.q.max), linewidth = 2, colour = 'white' ) +
    #geom_errorbarh(aes(y = upper.q.max, xmin = upper.q.doy + sd.upper.q.doy, xmax = upper.q.doy - sd.upper.q.doy), height = 0.005, linewidth = 2, colour = 'white' ) +
    #geom_errorbarh(aes(y = lower.q.max, xmin = lower.q.doy + sd.lower.q.doy, xmax = lower.q.doy - sd.lower.q.doy), height = 0.005, linewidth = 2, colour = 'white' ) +
    annotate("point", x = lower.q.doy, y = lower.q.max, colour = "#a8ddb5", size = 2) +
    annotate("point", x = upper.q.doy, y = upper.q.max, colour = "#084081", size = 2) +
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
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), 
                       labels = c('0', '0.2', '0.4', '0.6', '0.8')) +
    # xlim(c(175, 260)) +
    #  ylim(c(0, 0.6)) +
    coord_cartesian(xlim = c(180, 250), ylim = c(0, 0.9)) + # 175, 260
    theme_cowplot() +
    theme(legend.position = 'none')
  
}

# Functions for filtering the data
# Filter the data to obtain the upper and lower quartile of snow.auc
upper <- function(data) {
  data.no.zero <- data %>% filter(snow.auc != 0)
  upper.quartile <- quantile(data.no.zero$snow.auc, probs = 0.9)
  data <- data %>% filter(snow.auc >= upper.quartile)
  return(data)
}

lower <- function(data) {
  lower.quartile <- quantile(data$snow.auc, probs = 0.1)
  data <- data %>% filter(snow.auc <= lower.quartile)
  #data <- data %>% filter(snow.auc == 0)
  return(data)
}


# Load in data----
# Blaesedalen, S30
s30.bl <- read_csv('../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv',
                   show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

s30.bl.full <- read_csv('../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv',
                                 show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  drop_na(ndvi.max)

# Sentinel-2 Blaesedalen data for map plots
s2.bl <- read_csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Generate plots ----
s30.bl_max_plot <- fit_plot_max(s30.bl, 
                            colour.site = '#4984BF', 
                            colour.darker = '#2E5277',
                            colour.lighter = '#9BB2DA',
                            colour.lightest = '#BECBE7', 
                            ymax = 0.8, ymin = 0.4, yint = 0.1, 
                            xmax = 12, xmin = 0, xint = 2)

s30.bl_doy_plot <- fit_plot_doy(s30.bl,
                           colour.site = '#4984BF', 
                           colour.darker = '#2E5277',
                           colour.lighter = '#9BB2DA',
                           colour.lightest = '#BECBE7', 
                           ymax = 240, ymin = 210, yint = 10,
                           xmax = 12, xmin = 0, xint = 2) 
s30.bl_max_plot
s30.bl_doy_plot

s30.fig2 <- plot_data_fig2(s30.bl.full, rectangle.top = 0.845, rug.alpha = 10, line.width = 1.5)
s30.fig2

# Maps
s30.doy.map <- ggplot() +
  geom_sf(data = st_buffer(s30.bl, 15, endCapStyle = "SQUARE"), aes(fill = ndvi.max.doy)) +
  scale_fill_viridis_c() +
  theme_void() +
  theme(legend.position = 'None')

s2.doy.map <- ggplot() +
  geom_sf(data = st_buffer(s2.bl, 5, endCapStyle = "SQUARE"), aes(fill = ndvi.max.doy)) +
  scale_fill_viridis_c() +
  theme_void() +
  theme(legend.position = 'None')

s2.snow.map <- ggplot() +
  geom_sf(data = st_buffer(s2.bl, 5, endCapStyle = "SQUARE"), aes(fill = snow.auc)) +
  scale_fill_gradientn(colors = c('#ccebc5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081')) +
  theme_void() +
  theme(legend.position = 'None')
  
            

# Combine plots
maps <- plot_grid(s2.doy.map, s2.snow.map, s30.doy.map, ncol = 3, align = 'h')
maps

fits <- plot_grid(s30.bl_max_plot, s30.bl_doy_plot, nrow = 2, align = 'v')

bottom <- plot_grid(s30.fig2, fits, ncol = 2)

full.fig <- plot_grid(maps, bottom, nrow = 2)

full.fig

# Save plots
cowplot::save_plot('../../plots/figures/figure-5-lm-correctaxis-points.png', full.fig, 
                   base_height = 200, base_width = 180, units = 'mm',#)#, 
                   bg = 'white')

