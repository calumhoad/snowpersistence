# R script to test INLA model sensitivity
# Jakob Assmann, modified by Calum Hoad July 5th 2024

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
fit_matern <- function(shape_param, scale_param, range_site, site_data, return_model = F){
  
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
plot_results <- function(model_matern, site_data, grid.size, 
                         symin, symax, syint, scountmax, scountint){
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
    geom_point(aes(x = snow.auc, y= ndvi.max)) +
    geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                    ymax = max_ndvi_pred),
                data = preds_credible,
                fill = "blue", alpha = 0.3) +
    xlab("Snow persistence") +
    ylab("peak NDVI") +
    scale_y_continuous(breaks = c(seq(symin, symax-syint, syint))) + 
    coord_cartesian(ylim = c(symin, symax)) +
    geom_line(aes(x = snow.auc, y = preds), colour = "blue") +
    theme_cowplot()
  
  # Plot predictions in space
  preds_in_space <- ggplot() +
    geom_sf(data = st_buffer(site_data, grid.size/2, endCapStyle = "SQUARE"), 
            aes(fill = preds)) +
    scale_fill_continuous_sequential("viridis", rev = F) +
    labs(fill = "Predictions\npeak NDVI") +
    theme_map()
  
  # Residual plotsukc
  resids_hist <- ggplot(site_data) +
    geom_histogram(aes(x = residuals), bins = 50) +
    labs(x = "Residuals") +
    labs(y = 'Count') +
    scale_y_continuous(breaks = c(seq(0, scountmax-scountint, scountint)),
                       labels = c(as.character(seq(0, scountmax-scountint, scountint)))) +
    coord_cartesian(ylim = c(0, scountmax)) +
    theme_cowplot()
  resids_space <- ggplot() +
    geom_sf(data = st_buffer(site_data, grid.size/2, endCapStyle = "SQUARE"), aes(fill = residuals)) +
    scale_fill_continuous_diverging() +
    labs(fill = "Residuals") +
    theme_map()
  
  # plot the whole lot
  whole.plot <- plot_grid(fit_plot, preds_in_space,
                  resids_hist, resids_space, nrow = 2, ncol = 2, 
                  labels = c('(a)', '(b)', '(c)', '(d)'))
  
  return(whole.plot)
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
fit_variogram(s2.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.bl)/3) # 64.2
fit_variogram(s2.kl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kl)/3) # 61.7
fit_variogram(s2.kh, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kh)/3) # 36.9
fit_variogram(s30.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s30.bl)/3)# 125


# Set empirical ranges from variograms above
ranges_df <- data.frame(site = c("s2.bl",
                                 "s2.kl",
                                 "s2.kh",
                                 "s30.bl"),
                        range = c(round(64.2/10),
                                  round(61.7/10),
                                  round(36.9/10), 
                                  round(125/30))
)


### Test hyperparameters ----

# Define a set of shape parameters to test
shape_param <- c(1,10,20,50,75,100,150,400,500)

# Define a set of scale parameters to test
scale_param <- c(0.1,0.01,0.001,0.0001)


# Iteratively fit models to test for sensitivity to shape and scale

# Fit models for s2.bl
s2.bl.sensitivity <- map(scale_param, function(x){
  cat("Scale param = ", x, "\n")
  pblapply(shape_param, fit_matern, scale_param = x, 
           range_site = ranges_df$range[1], 
           site_data = s2.bl) %>% bind_rows()
}) %>% bind_rows()
write_csv(s2.bl.sensitivity, "../../data/statistical-output/s2-BL-ndvi-max-prior-sensitivity.csv")

# => This model is somehwat sensitive to the hyperparameter choice!
# Optimum paprameter values will estimate the range well
# We're aiming for a range around 5 cells -> choose scale of 0.1 and shape of
# 100, there seems to be the least sensitivity for the range in that area. 

# Fit models for s2.kl
s2.kl.sensitivity <- map(scale_param, function(x){
  cat("Scale param = ", x, "\n")
  pblapply(shape_param, fit_matern, scale_param = x, 
           range_site = ranges_df$range[2], 
           site_data = s2.kl) %>% bind_rows()
}) %>% bind_rows()
write_csv(s2.kl.sensitivity, "../../data/statistical-output/s2-KL-ndvi-max-prior-sensitivity.csv")
# This model is pretty consistent in the effect estimation across hyperparameters
# We aim for a range of around 6 cells, this seems to be most stable around  
# scale = 0.01 and shape = 75

# Fit models for s2.kh
s2.kh.sensitivity <- map(scale_param, function(x){
  cat("Scale param = ", x, "\n")
  pblapply(shape_param, fit_matern, scale_param = x, 
           range_site = ranges_df$range[3], 
           site_data = s2.kh) %>% bind_rows()
}) %>% bind_rows() 
write_csv(s2.kh.sensitivity, "../../data/statistical-output/s2-KH-ndvi-max-prior-sensitivity.csv")
# Again fairly stable, we're aimaing for a range of 4, so scale = 0.1 and shape =
# 400 seems to be best fitting

# Fit models for s30.bl
s30.bl.sensitivity <- map(scale_param, function(x){
  cat("Scale param = ", x, "\n")
  pblapply(shape_param, fit_matern, scale_param = x, 
           range_site = ranges_df$range[4], 
           site_data = s30.bl) %>% bind_rows()
}) %>% bind_rows() 
write_csv(s30.bl.sensitivity, "../../data/statistical-output/s30-BL-ndvi-max-prior-sensitivity.csv")
# => Very consitent reuslts, aim for 11 cells -> scale 0.1, shape 75

#### Now remove all model data on the hard drive!!! (do that manually for safety)
#list.dirs() 

# Fit final models ----
s2.bl_fit <- fit_matern(scale_param = 0.1,
                        shape_param = 1,
                        range_site = ranges_df$range[1], 
                        site_data = s2.bl,
                        return_model = T)
s2.kl_fit <- fit_matern(scale_param = 0.01,
                        shape_param = 1,
                        range_site = ranges_df$range[2], 
                        site_data = s2.kl,
                        return_model = T)
s2.kh_fit <- fit_matern(scale_param = 0.01,
                        shape_param = 1,
                        range_site = ranges_df$range[3], 
                        site_data = s2.kh,
                        return_model = T)
s30.bl_fit <- fit_matern(scale_param = 0.01,
                         shape_param = 1,
                         range_site = ranges_df$range[4], 
                         site_data = s30.bl,
                         return_model = T)

### Visualise results

# Plot results for all models
bl.plot <- plot_results(s2.bl_fit, s2.bl, 10, symin = 0.1, symax = 0.7, syint = 0.1, 
                        scountmax = 80, scountint = 20)
kl.plot <- plot_results(s2.kl_fit, s2.kl, 10, symin = 0.3, symax = 0.8, syint = 0.1, 
                        scountmax = 70, scountint = 10)
kh.plot <- plot_results(s2.kh_fit, s2.kh, 10, symin = 0.1, symax = 0.7, syint = 0.1, 
                        scountmax = 160, scountint = 20)
bl.s30.plot <- plot_results(s30.bl_fit, s30.bl, 30, symin = 0.4, symax = 1,
                            syint = 0.1, scountmax = 8, scountint = 2)

# Output plots
cowplot::save_plot('../../plots/supplementary/inla-bl-residuals.png', bl.plot,
                   base_height = 180, base_width = 180, units = 'mm', bg = 'white')
cowplot::save_plot('../../plots/supplementary/inla-kl-residuals.png', kl.plot,
                   base_height = 180, base_width = 180, units = 'mm', bg = 'white')
cowplot::save_plot('../../plots/supplementary/inla-kh-residuals.png', kh.plot,
                   base_height = 180, base_width = 180, units = 'mm', bg = 'white')
cowplot::save_plot('../../plots/supplementary/inla-bl30-residuals.png', bl.s30.plot,
                   base_height = 180, base_width = 180, units = 'mm', bg = 'white')

### Alternative for S2 BL: fit a break point model (manually identifited as snow.auc = 5)

# Split data into two
s2.bl_snow_auc_lte5 <- s2.bl %>% filter(snow.auc <= 5) 
s2.bl_snow_auc_gt5 <- s2.bl %>% filter(snow.auc > 5) 

# Fit two models lte 5 first with out any modifications
s2.bl_fit_lte5 <- fit_matern(scale_param = 0.01,
                             shape_param = 1,
                             range_site = ranges_df$range[1], 
                             site_data = s2.bl_snow_auc_lte5,
                             return_model = T) 
summary(s2.bl_fit_lte5)

# gt5 next with fixed intercept prior from s2.bl_fit_lte5 model

# Calculate prediction for 5 for first model
intercept_gt5 <- s2.bl_fit_lte5$summary.fixed[1,1] + 5 * s2.bl_fit_lte5$summary.fixed[2,1]

# Remove intercept from data and remove offest from 0 on x axis
s2.bl_snow_auc_gt5$snow.auc <- s2.bl_snow_auc_gt5$snow.auc - 5

# Define max number of rows and cols
n_col = max(s2.bl$col_number)
n_row = max(s2.bl$row_number)

# Specify priors for hyperparemeter
log.range = list(initial = log(ranges_df$range[1]), fixed = TRUE)
hyperpar_matern = list(initial = 2, param = c(1,  0.01))

# Specify formula (no intercept)
formula_matern <- ndvi.max ~ snow.auc +
  f(node, 
    model = "matern2d", 
    nrow = n_row,
    ncol = n_col,
    hyper = list(range = log.range,
                 prec = hyperpar_matern))
# Fit model
s2.bl_fit_gt5 <- inla(formula = formula_matern, 
                      data = s2.bl_snow_auc_gt5,
                      family = "gaussian",
                      # Define prior for intercept with very high precision value
                      control.fixed = list(mean.intercept = intercept_gt5,
                                           prec.intercept = 0.1),
                      control.predictor = list(compute = TRUE),
                      control.compute = list(dic = TRUE, 
                                             waic = TRUE, 
                                             cpo = TRUE, 
                                             return.marginals.predictor = TRUE),
                      keep = TRUE)

# Get the summary
summary(s2.bl_fit_gt5)  

# Calculate marginal predictions depending on snow.auc only
s2.bl_snow_auc_lte5$preds <- s2.bl_fit_lte5$summary.fixed[1,1] + 
  s2.bl_snow_auc_lte5$snow.auc * s2.bl_fit_lte5$summary.fixed[2,1]
s2.bl_snow_auc_gt5$preds <- s2.bl_fit_gt5$summary.fixed[1,1] + 
  s2.bl_snow_auc_gt5$snow.auc * s2.bl_fit_gt5$summary.fixed[2,1]

# Calculate min and max values for 95 CI predictions to plot a ribbon
preds_credible_lte5 <- expand.grid(intercept = unlist(s2.bl_fit_lte5$summary.fixed[1, c(1,3,5)]), 
                                   slope = unlist(s2.bl_fit_lte5$summary.fixed[2, c(1,3,5)])) %>%
  split(1:nrow(.)) %>%
  map(function(parameters) {
    data.frame(snow.auc =  seq((floor(min(s2.bl_snow_auc_lte5$snow.auc) * 10) / 10),
                               (ceiling(max(s2.bl_snow_auc_lte5$snow.auc) * 10) / 10), 0.01),
               ndvi_pred = parameters$intercept[1] + seq((floor(min(s2.bl_snow_auc_lte5$snow.auc) * 10) / 10),
                                                         (ceiling(max(s2.bl_snow_auc_lte5$snow.auc) * 10) / 10), 0.01) * parameters$slope[1])
  }) %>% bind_rows() %>%
  group_by(snow.auc) %>% summarise(min_ndvi_pred = min(ndvi_pred), max_ndvi_pred = max(ndvi_pred))
preds_credible_gt5 <- expand.grid(intercept = unlist(s2.bl_fit_gt5$summary.fixed[1, c(1,3,5)]), 
                                  slope = unlist(s2.bl_fit_gt5$summary.fixed[2, c(1,3,5)])) %>%
  split(1:nrow(.)) %>%
  map(function(parameters) {
    data.frame(snow.auc =  seq((floor(min(s2.bl_snow_auc_gt5$snow.auc) * 10) / 10),
                               (ceiling(max(s2.bl_snow_auc_gt5$snow.auc) * 10) / 10), 0.01),
               ndvi_pred = parameters$intercept[1] + seq((floor(min(s2.bl_snow_auc_gt5$snow.auc) * 10) / 10),
                                                         (ceiling(max(s2.bl_snow_auc_gt5$snow.auc) * 10) / 10), 0.01) * parameters$slope[1])
  }) %>% bind_rows() %>%
  group_by(snow.auc) %>% summarise(min_ndvi_pred = min(ndvi_pred), max_ndvi_pred = max(ndvi_pred)) %>% 
  mutate(snow.auc = snow.auc + 5)

# Add model fitted values and residuals to data frame
s2.bl_snow_auc_lte5$fitted <- s2.bl_fit_lte5$summary.fitted.values$mean
s2.bl_snow_auc_lte5$residuals <- s2.bl_snow_auc_lte5$ndvi.max - s2.bl_snow_auc_lte5$fitted
s2.bl_snow_auc_gt5$fitted <- s2.bl_fit_gt5$summary.fitted.values$mean 
s2.bl_snow_auc_gt5$residuals <- s2.bl_snow_auc_gt5$ndvi.max - s2.bl_snow_auc_gt5$fitted

# Re-add intercept to ndvi.max for gt5 and shift data on x axis
s2.bl_snow_auc_gt5$snow.auc <- s2.bl_snow_auc_gt5$snow.auc + 5


  # Plot model fits
bl_plot <- ggplot() +
  geom_point(data = s2.bl_snow_auc_lte5, aes(x = snow.auc, y= ndvi.max), colour = '#BECBE7') +
  geom_point(data = s2.bl_snow_auc_gt5, aes(x = snow.auc, y= ndvi.max), colour ='#BECBE7') +
  geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                  ymax = max_ndvi_pred),
              data = preds_credible_lte5,
              fill = "#4984BF", alpha = 0.5) +
  geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                  ymax = max_ndvi_pred),
              data = preds_credible_gt5,
              fill = '#4984BF', alpha = 0.5) +
  geom_line(aes(x = snow.auc, y = preds), colour = '#4984BF', data = s2.bl_snow_auc_lte5) +
  geom_line(aes(x = snow.auc, y = preds), colour = '#4984BF', data = s2.bl_snow_auc_gt5) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6), 
                     labels = c('0.0', '0.2', '0.4', '0.6')) +
  coord_cartesian(ylim = c(0, 0.7), xlim = c(0, 25)) +
  ylab('peak NDVI') +
  xlab('Snow persistence') +
  theme_cowplot()
bl_plot

# Residual plots
hist_lte5 <- ggplot(s2.bl_snow_auc_lte5) +
  geom_histogram(aes(x = residuals), bins = 50) +
  labs(x = "Residuals\n(Pixels less than or\nequal to the snow\npersistence value of 5)") +
  ylab('Count') +
  coord_cartesian( ylim = c(0, 100), xlim = c(-0.002, 0.002)) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80), 
                     labels = c('0', '20', '40', '60', '80')) +
  theme_cowplot()
hist_gt5 <- ggplot(s2.bl_snow_auc_gt5) +
  geom_histogram(aes(x = residuals), bins = 50) +
  labs(x = "Residuals\n(Pixels greater than the snow\npersistence value of 5)") +
  ylab('Count') +
  scale_y_continuous(breaks = c(0, 4, 8, 12), 
                     labels = c('0', '4', '8', '12')) +
  coord_cartesian(ylim = c(0, 14), xlim = c(-0.002, 0.002)) +
  theme_cowplot()
plot_grid(hist_lte5, hist_gt5, ncol = 2, labels = c('lte5', 'gt5'))

resids_space <- ggplot() +
  geom_sf(data = st_buffer(s2.bl_snow_auc_lte5, 5, endCapStyle = "SQUARE"), aes(fill = residuals)) +
  geom_sf(data = st_buffer(s2.bl_snow_auc_gt5, 5, endCapStyle = "SQUARE"), aes(fill = residuals)) +
  scale_fill_continuous_diverging() +
  labs(fill = 'Residuals') +
  theme_map()
resids_space

# Plot predictions in space
preds_in_space <- ggplot() +
  geom_sf(data = st_buffer(bind_rows(s2.bl_snow_auc_lte5, s2.bl_snow_auc_gt5), 5, endCapStyle = "SQUARE"), 
          aes(fill = preds)) +
  scale_fill_continuous_sequential("viridis", rev = F) +
  labs(fill = "Predictions\nNDVI max") +
  theme_map()
data_in_space <- ggplot() +
  geom_sf(data = st_buffer(s2.bl, 5, endCapStyle = "SQUARE"), 
          aes(fill = ndvi.max)) +
  scale_fill_continuous_sequential("viridis", rev = F) +
  labs(fill = "true\nNDVI max") +
  theme_map()
plot_grid(preds_in_space, data_in_space)

middle <- cowplot::plot_grid(hist_lte5, hist_gt5, nrow = 1,
                   ncol = 2, align = 'h', labels = c('(b)', '(c)'))
bottom <- cowplot::plot_grid(resids_space, preds_in_space, nrow = 1, ncol = 2, 
                             labels = c('(d)', '(e)'))
combined <- cowplot::plot_grid(bl_plot, middle, bottom, nrow = 3, labels = c('(a)', '', ''))
combined

cowplot::save_plot('../../plots/supplementary/inla-bl-break-residuals.png', combined,
                   base_height = 180, base_width = 180, units = 'mm', bg = 'white')
# Final plots ----
fit_plot <- function(model_matern, site_data, grid.size, 
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
    geom_point(aes(x = snow.auc, y= ndvi.max), colour = colour.lightest) +
    geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                    ymax = max_ndvi_pred),
                data = preds_credible,
                fill = colour.site, alpha = 0.3) +
    geom_line(aes(x = snow.auc, y = preds), colour = colour.site, linewidth = 1) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    theme_cowplot()
  
  return(fit_plot)
}


kl <- fit_plot(s2.kl_fit, s2.kl, 10, 
               colour.site = '#F5A40C', 
               colour.darker = '#946606',
               colour.lighter = '#FBCA7F',
               colour.lightest = '#FDDCAC', 
               ymax = 0.7, ymin = 0.4)
kh <- fit_plot(s2.kh_fit, s2.kh, 10, 
               colour.site = '#F23835', 
               colour.darker = '#8D271E',
               colour.lighter = '#F29580',
               colour.lightest = '#F8BBAA', 
               ymax = 0.6, ymin = 0.15)

top <- plot_grid(kl, kh, ncol = 2, align = 'h')
bottom <- plot_grid(bl_plot, bl_plot, ncol = 2, align = 'h')
combined <- plot_grid(top, bottom, nrow = 2, align = 'v')
combined

# Save plots
cowplot::save_plot('../../plots/figures/figure-3-inla-v2.png', combined, 
                   base_height = 140, base_width = 180, units = 'mm', 
                   bg = 'white')

# S30 plots
s30.bl.plot <- fit_plot(s30.bl_fit, s30.bl, 10, 
                        colour.site = '#4984BF', 
                        colour.darker = '#2E5277',
                        colour.lighter = '#9BB2DA',
                        colour.lightest = '#BECBE7', 
                        ymax = 0.8, ymin = 0.45)
s30.bl.plot

cowplot::save_plot('../../plots/figures/figure-5-ndvimax.png', s30.bl.plot, 
                   base_height = 70, base_width = 90, units = 'mm', 
                   bg = 'white')


# Get summary from S2 BL break-point model
model_matern <- s2.bl_fit_gt5
site_data <- s2.bl_snow_auc_gt5

data.frame(scale_param = 0.01,
           shape_param = 1, 
           intercept = round(model_matern$summary.fixed$mean[1],5),
           slope = round(model_matern$summary.fixed$mean[2],5),
           slope_q025 = round(model_matern$summary.fixed[2,3],5),
           slope_q975 = round(model_matern$summary.fixed[2,5],5),
           range = round(model_matern$summary.hyperpar$mean[3],2), 
           dic = model_matern$dic$dic,
           dic_eff_par = model_matern$dic$p.eff,
           dic_d_p = nrow(site_data)/model_matern$dic$p.eff,
           waic = model_matern$waic$waic,
           waic_eff_par = model_matern$waic$p.eff)
