# Functions for INLA and OLS model fitting in script 7
# Calum Hoad, adapted from Jakob Assmann, 4th July 2024

# Dependencies ----
library(tidyverse)
library(sf)
library(ggplot2)
library(cowplot)
library(pbapply)
library(INLA)
library(gt)
library(stargazer)
library(gtsummary)
library(broom)
library(gstat)

# Helper functions ----

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
  
  # Derive spatial index in INLA fashion (running number)
  # Obtain number of cols and rows
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
fit_variogram_max <- function(input.data, var, grid.size, cutoff.dist, model.type) {
  
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
             hjust = 0.5, vjust = 1.5, color = 'red') +
    scale_x_continuous(limits = c(grid.size, cutoff.dist), breaks = seq(grid.size,
                                                                        cutoff.dist,
                                                                        grid.size)) +
    labs(x = "lag distance (m)", y = "ndvi.max") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(vario.plot)
}

fit_variogram_doy <- function(input.data, var, grid.size, cutoff.dist, model.type) {
  
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
    labs(x = "lag distance (m)", y = "ndvi.max") +
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

summarise_matern <- function(model_matern, site.name, sensor.name) {
  
  fixed_effects <- data.frame(
    Parameter = rownames(model_matern$summary.fixed),
    Mean = model_matern$summary.fixed[, "mean"],
    SD = model_matern$summary.fixed[, "sd"],
    Q025 = model_matern$summary.fixed[, "0.025quant"],
    #Q50 = model_matern$summary.fixed[, "0.5quant"],
    Q975 = model_matern$summary.fixed[, "0.975quant"],
    Mode = model_matern$summary.fixed[, "mode"], 
    site = site.name, 
    sensor = sensor.name)
  
  # Extract and summarize hyperparameters
  hyperparameters <- data.frame(
    Parameter = rownames(model_matern$summary.hyperpar),
    Mean = model_matern$summary.hyperpar[, "mean"],
    SD = model_matern$summary.hyperpar[, "sd"],
    Q025 = model_matern$summary.hyperpar[, "0.025quant"],
    Q50 = model_matern$summary.hyperpar[, "0.5quant"],
    Q975 = model_matern$summary.hyperpar[, "0.975quant"],
    Mode = model_matern$summary.hyperpar[, "mode"])
  
  return(fixed_effects)
}