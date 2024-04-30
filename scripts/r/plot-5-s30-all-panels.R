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
fit_plot_doy <- function(model_matern, site_data, grid.size, 
                     colour.site, 
                     colour.darker, 
                     colour.lighter,
                     colour.lightest, 
                     ymax, ymin){
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
  site_data$residuals <- site_data$ndvi.max - site_data$fitted
  
  # Plot model fit
  fit_plot <- ggplot(data = site_data) +
    geom_point(aes(x = snow.auc, y = ndvi.max.doy), colour = colour.lightest) +
    #geom_point(data = site_data %>% distinct(snow.auc, ndvi.max.doy), aes(x = snow.auc, y= ndvi.max.doy), colour = 'black', alpha = 0.1) +
    geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                    ymax = max_ndvi_pred),
                data = preds_credible,
                fill = colour.site, alpha = 0.3) +
    geom_line(aes(x = snow.auc, y = preds), colour = colour.site, linewidth = 1) +
    xlab('') +
    coord_cartesian(ylim = c(ymin, ymax)) +
    theme_cowplot()
  
  return(fit_plot)
}

fit_plot_max <- function(model_matern, site_data, grid.size, 
                         colour.site, 
                         colour.darker, 
                         colour.lighter,
                         colour.lightest, 
                         ymax, ymin){
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
  site_data$residuals <- site_data$ndvi.max - site_data$fitted
  
  # Plot model fit
  fit_plot <- ggplot(data = site_data) +
    geom_point(aes(x = snow.auc, y = ndvi.max), colour = colour.lightest) +
    #geom_point(data = site_data %>% distinct(snow.auc, ndvi.max.doy), aes(x = snow.auc, y= ndvi.max.doy), colour = 'black', alpha = 0.1) +
    geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                    ymax = max_ndvi_pred),
                data = preds_credible,
                fill = colour.site, alpha = 0.3) +
    geom_line(aes(x = snow.auc, y = preds), colour = colour.site, linewidth = 1) +
    xlab('') +
    coord_cartesian(ylim = c(ymin, ymax)) +
    theme_cowplot()
  
  return(fit_plot)
}

# Function to generate curves plot, as per fig 2
plot_data_fig2 <- function(data) {
  
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
    coord_cartesian(xlim = c(180, 250), ylim = c(0, 0.9)) + # 175, 260
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
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  group_by(id) %>%
  drop_na(ndvi.max)

# Add INLA grid to df
s30.bl <- get_inla_grid(s30.bl, grid.size = 30)

# Set variogram ranges for ndvi.max and ndvi.max.doy
range.max.doy <- round(21.8/30)
range.max <- round(125/30)

# Fit models ----
s30.bl_max_fit <- fit_matern_max(scale_param = 0.01,
                             shape_param = 1,
                             range_site = range.max, 
                             site_data = s30.bl,
                             return_model = T)

s30.bl_doy_fit <- fit_matern_doy(scale_param = 0.01, 
                                 shape_param = 1, 
                                 range_site = range.max.doy,
                                 site_data = s30.bl, 
                                 return_model = T)

# Generate plots ----
s30.bl_max_plot <- fit_plot_max(s30.bl_max_fit, s30.bl, 10, 
                            colour.site = '#4984BF', 
                            colour.darker = '#2E5277',
                            colour.lighter = '#9BB2DA',
                            colour.lightest = '#BECBE7', 
                            ymax = 0.8, ymin = 0.45)

s30.bl_doy_plot <- fit_plot_doy(s30.bl_doy_fit, s30.bl, 10, 
                           colour.site = '#4984BF', 
                           colour.darker = '#2E5277',
                           colour.lighter = '#9BB2DA',
                           colour.lightest = '#BECBE7', 
                           ymax = 240, ymin = 210) 
s30.bl_max_plot
s30.bl_doy_plot

s30.fig2 <- plot_data_fig2(s30.bl.full)
s30.fig2
