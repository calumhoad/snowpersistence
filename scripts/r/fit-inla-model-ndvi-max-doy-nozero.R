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
fit_matern <- function(shape_param, scale_param, range_site, site_data, return_model = F){
  
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
fit_variogram(s2.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.bl)/3) # 31.3
fit_variogram(s2.kl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kl)/3) # 52.2
fit_variogram(s2.kh, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kh)/3) # 30.6
fit_variogram(s30.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s30.bl)/3)# 21.8


# Set empirical ranges from variograms above
ranges_df <- data.frame(site = c("s2.bl",
                                 "s2.kl",
                                 "s2.kh",
                                 "s30.bl"),
                        range = c(round(31.3/10),
                                  round(52.2/10),
                                  round(30.6/10), 
                                  round(21.8/30))
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
write_csv(s2.bl.sensitivity, "../../data/statistical-output/s2-BL-ndvi-doy-prior-sensitivity.csv")

# This model is insensitive to hyperparameter choice

# Fit models for s2.kl
s2.kl.sensitivity <- map(scale_param, function(x){
  cat("Scale param = ", x, "\n")
  pblapply(shape_param, fit_matern, scale_param = x, 
           range_site = ranges_df$range[2], 
           site_data = s2.kl) %>% bind_rows()
}) %>% bind_rows()
write_csv(s2.kl.sensitivity, "../../data/statistical-output/s2-KL-ndvi-doy-prior-sensitivity.csv")
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
write_csv(s2.kh.sensitivity, "../../data/statistical-output/s2-KH-ndvi-doy-prior-sensitivity.csv")
# Again fairly stable, we're aimaing for a range of 4, so scale = 0.1 and shape =
# 400 seems to be best fitting

# Fit models for s30.bl
s30.bl.sensitivity <- map(scale_param, function(x){
  cat("Scale param = ", x, "\n")
  pblapply(shape_param, fit_matern, scale_param = x, 
           range_site = ranges_df$range[4], 
           site_data = s30.bl) %>% bind_rows()
}) %>% bind_rows() 
write_csv(s30.bl.sensitivity, "../../data/statistical-output/s30-BL-ndvi-doy-prior-sensitivity.csv")
# => Very consitent reuslts, aim for 11 cells -> scale 0.1, shape 75

#### Now remove all model data on the hard drive!!! (do that manually for safety)
#list.dirs() 

# Try removing 0 values from the data, see effect on models
s2.bl.gt0 <- s2.bl %>% filter(snow.auc > 0)

# Fit final models ----
s2.bl_fit <- fit_matern(scale_param = 0.01,
                        shape_param = 1,
                        range_site = ranges_df$range[1], 
                        site_data = s2.bl.gt0,
                        return_model = T)
s2.kl_fit <- fit_matern(scale_param = 0.01,
                        shape_param = 1,
                        range_site = ranges_df$range[2], 
                        site_data = s2.kl %>% filter(snow.auc > 0),
                        return_model = T)
s2.kh_fit <- fit_matern(scale_param = 0.01,
                        shape_param = 1,
                        range_site = ranges_df$range[3], 
                        site_data = s2.kh %>% filter(snow.auc > 0),
                        return_model = T)
s30.bl_fit <- fit_matern(scale_param = 0.01,
                         shape_param = 1,
                         range_site = ranges_df$range[4], 
                         site_data = s30.bl %>% filter(snow.auc > 0),
                         return_model = T)

### Visualise results

# Plot results for all models
plot_results(s2.bl_fit, s2.bl.gt0, 10)
plot_results(s2.kl_fit, s2.kl %>% filter(snow.auc > 0), 10)
plot_results(s2.kh_fit, s2.kh %>% filter(snow.auc > 0), 10)
plot_results(s30.bl_fit, s30.bl %>% filter(snow.auc > 0), 30)


# S2 KH has a very low intercept compared with the data,
# try thinning out the 0 snow.auc values to reduce their weight on the model

# # Quick and dirty subsampling for KH site
hist(s2.kh.thin$snow.auc, breaks = seq(-0.5, 23.5,1))
s2.kh$bin <- round(s2.kh$snow.auc)
group_by(s2.kh, bin) %>% st_drop_geometry() %>% tally()
group_by(s2.kh, bin) %>% st_drop_geometry() %>% tally() %>% summarize(mean(n))
s2.kh.thin <- group_by(s2.kh, bin) %>% group_map(function(x, ...) {
  if(nrow(x) >= 20) x <- sample_n(x, 350)
  return(x)
}) %>% bind_rows()
s2.kh.thin %>% mutate(bin = round(snow.auc)) %>%
  group_by(bin) %>% st_drop_geometry() %>% tally()


s2.kh.thin <- s2.kh


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
    geom_point(aes(x = snow.auc, y= ndvi.max.doy), colour = colour.lightest) +
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

bl <- fit_plot(s2.bl_fit, s2.bl %>% filter(snow.auc > 0), 10, 
               colour.site = '#4984BF', 
               colour.darker = '#2E5277',
               colour.lighter = '#9BB2DA',
               colour.lightest = '#BECBE7', 
               ymax = 240, ymin = 215)
kl <- fit_plot(s2.kl_fit, s2.kl %>% filter(snow.auc > 0), 10, 
               colour.site = '#F5A40C', 
               colour.darker = '#946606',
               colour.lighter = '#FBCA7F',
               colour.lightest = '#FDDCAC', 
               ymax = 235, ymin = 210)
kh <- fit_plot(s2.kh_fit, s2.kh %>% filter(snow.auc > 0), 10, 
               colour.site = '#F23835', 
               colour.darker = '#8D271E',
               colour.lighter = '#F29580',
               colour.lightest = '#F8BBAA', 
               ymax = 235, ymin = 210)

# Create panel plot
top <- plot_grid(kl, kh, ncol = 2, align = 'h')
bottom <- plot_grid(bl, time.plot, ncol = 2, align = 'h')
combined <- plot_grid(top, bottom, nrow = 2, align = 'v')

combined.2 <- plot_grid(kl, kh, bl, time.plot, ncol = 2, nrow = 2, align = 'hv')
combined.2
#Check the plot
combined 

# Save plots
cowplot::save_plot('../../plots/figures/figure-4-inla.png', combined.2, 
                   base_height = 140, base_width = 180, units = 'mm', 
                   bg = 'white')

# S30 plots
s30.bl.plot <- fit_plot(s30.bl_fit, s30.bl, 10, 
                        colour.site = '#4984BF', 
                        colour.darker = '#2E5277',
                        colour.lighter = '#9BB2DA',
                        colour.lightest = '#BECBE7', 
                        ymax = 240, ymin = 210)
s30.bl.plot

cowplot::save_plot('../../plots/figures/figure-5-ndvimaxdoy.png', s30.bl.plot, 
                   base_height = 70, base_width = 90, units = 'mm', 
                   bg = 'white')

# Create panel (d) visualising imagery dates
bl.dates <- c('2023-07-02', '2023-07-12', '2023-07-18', '2023-07-26')
kl.dates <- c('2022-06-29', '2022-07-05', '2022-07-18', '2022-08-01', '2022-08-14')
kh.dates <- c('2022-07-09', '2022-07-19', '2022-07-29', '2022-08-04', '2022-08-13')

bl.dates.df <- tibble(bl.dates) %>% mutate(site = 'BL') %>% rename(date ='bl.dates')
kl.dates.df <- tibble(kl.dates) %>% mutate(site = 'KL') %>% rename(date = 'kl.dates')
kh.dates.df <- tibble(kh.dates) %>% mutate(site = 'KH') %>% rename(date = 'kh.dates')

imagery.dates <- rbind(bl.dates.df, kl.dates.df, kh.dates.df) %>%
  mutate(date = date(date)) %>%
  mutate(date.no_year = format(date, '%m-%d'))

time.plot <- ggplot() +
  geom_vline(xintercept = 182, linetype = 'solid', alpha = 0.2) +
  geom_vline(xintercept = 213, linetype = 'solid', alpha = 0.2) +
  geom_line(data = imagery.dates, aes(x = yday(date), y = site, group = site, colour = site),
            linewidth = 2) +
  geom_point(data = imagery.dates, aes(x = yday(date), y = site, group = site, colour = site), 
             size = 4) +
  scale_color_manual(values = c('#4984BF', '#F5A40C', '#F23835'), breaks = c('BL', 'KL', 'KH')) +
  new_scale_colour() +
  geom_point(data = imagery.dates, aes(x = yday(date), y = site, group = site, colour = site), 
             size = 2) +
  scale_color_manual(values = c('#BECBE7', '#FBCA7F', '#F29580'), breaks = c('BL', 'KL', 'KH')) +
  scale_x_continuous(breaks = c(182, 189, 196, 203, 210, 213, 217, 224),
                    labels = c('1\nJuly', '8', '15', '22', '29',
                               '\nAugust', '5', '12')) +
  ylab('') +
  xlab('') +
  theme_cowplot() +
  theme(legend.position = 'none', 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) 

bl <- fit_plot(s2.bl_fit, s2.bl, 10, 
               colour.site = '#4984BF', 
               colour.darker = '#2E5277',
               colour.lighter = '#9BB2DA',
               colour.lightest = '#BECBE7', 
               ymax = 240, ymin = 215)
kl <- fit_plot(s2.kl_fit, s2.kl, 10, 
               colour.site = '#F5A40C', 
               colour.darker = '#946606',
               colour.lighter = '#FBCA7F',
               colour.lightest = '#FDDCAC', 
               ymax = 235, ymin = 210)
kh <- fit_plot(s2.kh_fit, s2.kh, 10, 
               colour.site = '#F23835', 
               colour.darker = '#8D271E',
               colour.lighter = '#F29580',
               colour.lightest = '#F8BBAA', 
               ymax = 235, ymin = 210)