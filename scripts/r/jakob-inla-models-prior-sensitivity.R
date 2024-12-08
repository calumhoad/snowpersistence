# R script to fit final INLA models for Calum
# Dependencies
library(tidyverse)
library(sf)
library(ggplot2)
library(cowplot)
library(colorspace)
library(pbapply)
library(broom)
library(INLA)

### Data Prep ##########

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

# Helper function to generate INLA grid
get_inla_grid <- function(snow.auc_df){
  # Calculate row and col_numbers
  snow.auc_df <- snow.auc_df %>% 
    # Calculate min X and max >
    mutate(min_x = min(.$X),
           max_y = max(.$Y)) %>%
    # Calczulate colum and row numbers for each cell (top-left to bottom-right corner)
    # We shift by one a R is 1- based
    mutate(col_number = 1 + ((X - min_x) / 10),
           row_number = 1 + ((max_y - Y) / 10)) %>%
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

# Add INLA grids to dfs
s2.bl <- get_inla_grid(s2.bl)
s2.kl <- get_inla_grid(s2.kl)
s2.kh <- get_inla_grid(s2.kh)
s30.bl <- get_inla_grid(s30.bl)

# # Quick and dirty subsampling for blasedalen site
# hist(s2.bl$snow.auc, breaks = seq(-0.5, 23.5,1))
# s2.bl$bin <- round(s2.bl$snow.auc)
# group_by(s2.bl, bin) %>% st_drop_geometry() %>% tally() 
# group_by(s2.bl, bin) %>% st_drop_geometry() %>% tally() %>% summarize(mean(n))
# s2.bl <- group_by(s2.bl, bin) %>% group_map(function(x, ...) {
#   if(nrow(x) >= 36) x <- sample_n(x, 36)
#   return(x)
# }) %>% bind_rows()
# s2.bl %>% mutate(bin = round(snow.auc)) %>% 
#   group_by(bin) %>% st_drop_geometry() %>% tally() 

# Set emperical ranges (from Calum's variograms)
ranges_df <- data.frame(site = c("s2.bl",
                    "s2.kl",
                    "s2.kh",
                    "s30.bl"),
           range = c(round(51.6/10),
                     round(61.9/10),
                     round(38.5/10), 
                     round(110/10))
           )

### Test hyperparameters
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

# Define a set of shape parameters to test
shape_param <- c(1,10,20,50,75,100,150,400,500)

# Define a set of scale parameters to test
scale_param <- c(0.1,0.01,0.001,0.0001)

# Fit models for s2.bl
map(scale_param, function(x){
  cat("Scale param = ", x, "\n")
  pblapply(shape_param, fit_matern, scale_param = x, 
           range_site = ranges_df$range[1], 
           site_data = s2.bl) %>% bind_rows()
}) %>% bind_rows()
# => This model is somehwat sensitive to the hyperparameter choice!
# Optimum paprameter values will estimate the range well
# We're aiming for a range around 5 cells -> choose scale of 0.1 and shape of
# 100, there seems to be the least sensitivity for the range in that area. 

# Fit models for s2.kl
map(scale_param, function(x){
  cat("Scale param = ", x, "\n")
  pblapply(shape_param, fit_matern, scale_param = x, 
           range_site = ranges_df$range[2], 
           site_data = s2.kl) %>% bind_rows()
}) %>% bind_rows() 
# This model is pretty consistent in the effect estimation across hyperparameters
# We aim for a range of around 6 cells, this seems to be most stable around  
# scale = 0.01 and shape = 75

# Fit models for s2.kh
map(scale_param, function(x){
  cat("Scale param = ", x, "\n")
  pblapply(shape_param, fit_matern, scale_param = x, 
           range_site = ranges_df$range[3], 
           site_data = s2.kh) %>% bind_rows()
}) %>% bind_rows() 
# Again fairly stable, we're aimaing for a range of 4, so scale = 0.1 and shape =
# 400 seems to be best fitting

# Fit models for s30.bl
map(scale_param, function(x){
  cat("Scale param = ", x, "\n")
  pblapply(shape_param, fit_matern, scale_param = x, 
           range_site = ranges_df$range[4], 
           site_data = s30.bl) %>% bind_rows()
}) %>% bind_rows() 
# => Very consitent reuslts, aim for 11 cells -> scale 0.1, shape 75

#### Now remove all model data on the hard drive!!! (do that manually for safety)
#list.dirs() 

# Fit final models
s2.bl_fit <- fit_matern(scale_param = 0.01,
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

# Helper function to generate plots
plot_results <- function(model_matern, site_data){
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
    geom_line(aes(x = snow.auc, y = preds), colour = "blue") +
    theme_cowplot()
  
  # Plot predictions in space
  preds_in_space <- ggplot() +
    geom_sf(data = st_buffer(site_data, 5, endCapStyle = "SQUARE"), 
            aes(fill = preds)) +
    scale_fill_continuous_sequential("viridis", rev = F) +
    labs(fill = "predictions\nNDVI max") +
    theme_map()
  
  # Residual plots
  resids_hist <- ggplot(site_data) +
    geom_histogram(aes(x = residuals), bins = 50) +
    labs(x = "residuals") +
    theme_cowplot()
  resids_space <- ggplot() +
    geom_sf(data = st_buffer(site_data, 5, endCapStyle = "SQUARE"), aes(fill = residuals)) +
    scale_fill_continuous_diverging() +
    theme_map()
  
  # plot the whole lot
  print(plot_grid(fit_plot, preds_in_space,
                  resids_hist, resids_space, nrow = 2))
  
  return(NULL)
}

# Plot results for all models
plot_results(s2.bl_fit, s2.bl)
plot_results(s2.kl_fit, s2.kl)
plot_results(s2.kh_fit, s2.kh)
plot_results(s30.bl_fit, s30.bl)


#### 3rd degree polynom for BL data
# see here: https://becarioprecario.bitbucket.io/inla-gitbook/ch-smoothing.html
site_data <- s2.bl

# Define max number of rows and cols
n_col = max(site_data$col_number)
n_row = max(site_data$row_number)

# Specify priors for hyperparemeter
log.range = list(initial = log(ranges_df$range[1]), fixed = TRUE)
hyperpar_matern = list(initial = 2, param = c(1,  0.01))

# Specify formula
formula_matern <- ndvi.max ~snow.auc +  I(snow.auc^2) + I(snow.auc^3) +
  f(node,
    model = "matern2d",
    nrow = n_row,
    ncol = n_col,
    hyper = list(range = log.range,
                 prec = hyperpar_matern))

# Model works a lot better without spatial effect... but then what's the point?
# formula_matern <- ndvi.max ~snow.auc +  I(snow.auc^2) + I(snow.auc^3)

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
summary(model_matern)

# Calculate marginal predictions depending on snow.auc only
site_data$preds <- model_matern$summary.fixed[1, 1] +
  site_data$snow.auc * model_matern$summary.fixed[2, 1] +
  site_data$snow.auc^2 * model_matern$summary.fixed[3, 1] +
  site_data$snow.auc^3 * model_matern$summary.fixed[4,1]

# Calculate min and max values for 95 CI predictions to plot a ribbon
preds_credible <- expand.grid(intercept = unlist(model_matern$summary.fixed[1, c(1,3,5)]), 
                              linear = unlist(model_matern$summary.fixed[2, c(1,3,5)]),
                              square = unlist(model_matern$summary.fixed[3, c(1,3,5)]),
                              cube = unlist(model_matern$summary.fixed[4, c(1,3,5)])) %>%
  split(1:nrow(.)) %>%
  map(function(parameters) {
    snow.auc_new_data = seq((floor(min(site_data$snow.auc) * 10) / 10),
                            (ceiling(max(site_data$snow.auc) * 10) / 10), 0.01)
    data.frame(
      snow.auc = snow.auc_new_data,
      ndvi_pred = parameters$intercept + snow.auc_new_data * parameters$linear +
        snow.auc_new_data^2 * parameters$square + snow.auc_new_data^3 * parameters$cube
    )
  }) %>% bind_rows() %>%
  group_by(snow.auc) %>% summarise(min_ndvi_pred = min(ndvi_pred), max_ndvi_pred = max(ndvi_pred))

# Add model fitted values and residuals to data frame
site_data$fitted <- model_matern$summary.fitted.values$mean
site_data$residuals <- site_data$ndvi.max - site_data$fitted

# Plot model fit
fit_plot <- ggplot(data = site_data) +
  geom_point(aes(x = snow.auc, y= ndvi.max)) +
  geom_point(aes(x = snow.auc, y= fitted), colour = "red", alpha = 0.1) +
  geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                  ymax = max_ndvi_pred),
              data = preds_credible,
              fill = "blue", alpha = 0.3) +
  geom_line(aes(x = snow.auc, y = preds), colour = "blue") +
  theme_cowplot()

# Plot predictions in space
preds_in_space <- ggplot() +
  geom_sf(data = st_buffer(site_data, 5, endCapStyle = "SQUARE"), 
          aes(fill = preds)) +
  scale_fill_continuous_sequential("viridis", rev = F) +
  labs(fill = "predictions\nNDVI max") +
  theme_map()

# Residual plots
resids_hist <- ggplot(site_data) +
  geom_histogram(aes(x = residuals), bins = 50) +
  labs(x = "residuals") +
  theme_cowplot()
resids_space <- ggplot() +
  geom_sf(data = st_buffer(site_data, 5, endCapStyle = "SQUARE"), aes(fill = residuals)) +
  scale_fill_continuous_diverging() +
  theme_map()

# plot the whole lot
print(plot_grid(fit_plot, preds_in_space,
                resids_hist, resids_space, nrow = 2))


### Other alternative - fit a break point model (manually identifited as snow.auc = 5)

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

# Specify formula (no itnercept)
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
fit_plot <- ggplot() +
  geom_point(aes(x = snow.auc, y= ndvi.max), data = s2.bl_snow_auc_lte5) +
  geom_point(aes(x = snow.auc, y= ndvi.max), data = s2.bl_snow_auc_gt5) +
  geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                  ymax = max_ndvi_pred),
              data = preds_credible_lte5,
              fill = "blue", alpha = 0.3) +
  geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                  ymax = max_ndvi_pred),
              data = preds_credible_gt5,
              fill = "blue", alpha = 0.3) +
  geom_line(aes(x = snow.auc, y = preds), colour = "blue", data = s2.bl_snow_auc_lte5) +
  geom_line(aes(x = snow.auc, y = preds), colour = "blue", data = s2.bl_snow_auc_gt5) +
  theme_cowplot()
fit_plot

# Residual plots
ggplot(s2.bl_snow_auc_lte5) +
  geom_histogram(aes(x = residuals), bins = 50) +
  labs(x = "residuals") +
  theme_cowplot()
ggplot(s2.bl_snow_auc_gt5) +
  geom_histogram(aes(x = residuals), bins = 50) +
  labs(x = "residuals") +
  theme_cowplot()

# Plot predictions in space
preds_in_space <- ggplot() +
  geom_sf(data = st_buffer(bind_rows(s2.bl_snow_auc_lte5, s2.bl_snow_auc_gt5), 5, endCapStyle = "SQUARE"), 
          aes(fill = preds)) +
  scale_fill_continuous_sequential("viridis", rev = F) +
  labs(fill = "predictions\nNDVI max") +
  theme_map()
data_in_space <- ggplot() +
  geom_sf(data = st_buffer(s2.bl, 5, endCapStyle = "SQUARE"), 
          aes(fill = ndvi.max)) +
  scale_fill_continuous_sequential("viridis", rev = F) +
  labs(fill = "true\nNDVI max") +
  theme_map()
plot_grid(preds_in_space, data_in_space)
