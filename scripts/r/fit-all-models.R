# Fit all models, output tables
# Calum Hoad, calum.hoad@ed.ac.uk, 09/05/2024

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
# Data prep ----

## Load data
# Blasedalen
s2.bl <- read_csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "S2 BL", range = round(51.6/10), colour = '#4984BF')

# Kluane low
s2.kl <- read_csv('../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "S2 KL", range = round(61.9/10), colour = '#F5A40C')

# Kluane high
s2.kh <- read_csv('../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv',
                  show_col_types = FALSE) %>%
  st_as_sf(coords = c('X', 'Y'), remove = F, crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "S2 KH", range = round(38.5/10), colour = '#F23835')

# Blaesedalen, S30
s30.bl <- read_csv("../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv",
                   show_col_types = FALSE
) %>%
  st_as_sf(coords = c("X", "Y"), remove = F, crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() %>%
  mutate(site = "S30 BL", range = round(110/30), colour = '#4984BF')

# Combine into one data object (keeping order of Calum's plots)
data_list <- list(s2.kl, s2.kh, s2.bl, s30.bl)


# Spearmans rho ----

# Simple correlation tables (from Jakob's script)
bind_rows(
  map(data_list, function(site_data) {
    tibble(
      site_name = unique(site_data$site),
      var_x = "snow.auc",
      var_y = "ndvi.max",
      spearmans_rho = round(cor(site_data$snow.auc, site_data$ndvi.max, method = "spearman"), 3),
      p_value = round(cor.test(site_data$snow.auc, site_data$ndvi.max, method = "spearman")$p.value, 3)
    )
  }),
  map(data_list, function(site_data) {
    tibble(
      site_name = unique(site_data$site),
      var_x = "snow.auc",
      var_y = "ndvi.max.doy",
      spearmans_rho = round(cor(site_data$snow.auc, site_data$ndvi.max.doy, method = "spearman"), 3),
      p_value = round(cor.test(site_data$snow.auc, site_data$ndvi.max.doy, method = "spearman")$p.value, 3)
    )
  })
) %>% gt() %>%
  gtsave("../../plots/stats-tables/spearmans_table.png")


# Linear models, OLS regression ----

# peak-NDVI magnitude
lm.bl.max <- lm(data = s2.bl, formula = ndvi.max ~ snow.auc)
lm.kl.max <- lm(data = s2.kl, formula = ndvi.max ~ snow.auc)
lm.kh.max <- lm(data = s2.kh, formula = ndvi.max ~ snow.auc)
lm.s30bl.max <- lm(data = s30.bl, formula = ndvi.max ~ snow.auc)

stargazer(lm.kl.max, lm.kh.max, lm.bl.max, lm.s30bl.max, type = 'html',
          column.labels = c("Kluane Low", "Kluane High", "Blaesedalen, S2", "Blaesedalen, S30"), 
          out = '../../plots/stats-tables/linear-models-max.html')

# peak-NDVI magnitude, KH (y = ln(x + 1))
ln.kh.max <- lm(data = s2.kh, formula = ndvi.max ~ log(snow.auc + 1))

stargazer(ln.kh.max, type = 'html',
          column.labels = "Kluane High", 
          out = '../../plots/stats-tables/linear-models-kh-ln-max.html')

# peak-NDVI timing
lm.bl.doy <- lm(data = s2.bl, formula = ndvi.max.doy ~ snow.auc)
lm.kl.doy <- lm(data = s2.kl, formula = ndvi.max.doy ~ snow.auc)
lm.kh.doy <- lm(data = s2.kh, formula = ndvi.max.doy ~ snow.auc)
lm.s30bl.doy <- lm(data = s30.bl, formula = ndvi.max.doy ~ snow.auc)

stargazer(lm.kl.doy, lm.kh.doy, lm.bl.doy, lm.s30bl.doy, type = 'html',
          column.labels = c("Kluane Low", "Kluane High", "Blaesedalen, S2", "Blaesedalen, S30"), 
          out = '../../plots/stats-tables/linear-models-doy.html')

# peak-NDVI timing, KH (y = ln(x + 1))
ln.kh.doy <- lm(data = s2.kh, formula = ndvi.max.doy ~ log(snow.auc + 1))

stargazer(ln.kh.doy, type = 'html',
          column.labels = "Kluane High", 
          out = '../../plots/stats-tables/linear-models-kh-ln-doy.html')


# INLA Matern 2D models ----

# Add INLA grids to dfs
s2.bl <- get_inla_grid(s2.bl, grid.size = 10)
s2.kl <- get_inla_grid(s2.kl, grid.size = 10)
s2.kh <- get_inla_grid(s2.kh, grid.size = 10)
s30.bl <- get_inla_grid(s30.bl, grid.size = 30)


# peak-NDVI magnitude ###
# Get range for NDVI max from matern variogram at each site

# Fit variograms, with cutoff set to 1/3 the diagonal of the dataset
fit_variogram_max(s2.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.bl)/3) # 64.2
fit_variogram_max(s2.kl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kl)/3) # 61.7
fit_variogram_max(s2.kh, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kh)/3) # 36.9
fit_variogram_max(s30.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s30.bl)/3)# 125

# Set empirical ranges from variograms above
max_ranges_df <- data.frame(site = c("s2.bl",
                                 "s2.kl",
                                 "s2.kh",
                                 "s30.bl"),
                        range = c(round(64.2/10),
                                  round(61.7/10),
                                  round(36.9/10), 
                                  round(125/30))
)

# Fit models
max.s2.bl_fit <- fit_matern_max(scale_param = 0.01,
                            shape_param = 1,
                            range_site = max_ranges_df$range[1], 
                            site_data = s2.bl,
                            return_model = T)
max.s2.kl_fit <- fit_matern_max(scale_param = 0.01,
                            shape_param = 1,
                            range_site = max_ranges_df$range[2], 
                            site_data = s2.kl,
                            return_model = T)
max.s2.kh_fit <- fit_matern_max(scale_param = 0.01,
                            shape_param = 1,
                            range_site = max_ranges_df$range[3], 
                            site_data = s2.kh,
                            return_model = T)
max.s30.bl_fit <- fit_matern_max(scale_param = 0.01,
                             shape_param = 1,
                             range_site = max_ranges_df$range[4], 
                             site_data = s30.bl,
                             return_model = T)


# Extract and summarise fixed effects
max.s2.bl.mat <- summarise_matern(max.s2.bl_fit, site.name = "Blaesedalen, Sentinel-2", sensor.name = "Sentinel-2")
max.s2.kl.mat <- summarise_matern(max.s2.kl_fit, site.name = "Kluane Low", sensor.name = "Sentinel-2")
max.s2.kh.mat <- summarise_matern(max.s2.kh_fit, site.name = "Kluane High", sensor.name = "Sentinel-2")
max.s30.bl.mat <- summarise_matern(max.s30.bl_fit, site.name = "Blaesedalen, HLS S30", sensor.name = "HLS S30")

max.all.mat <- rbind(max.s2.kl.mat,
                     max.s2.kh.mat,
                     max.s2.bl.mat,
                     max.s30.bl.mat)


max.summary <- gt(max.all.mat %>% group_by(site) %>% select(-sensor)) %>%
  fmt_number(
    columns = vars(Mean, SD, Q025, Q975, Mode),
    decimals = 5, # Set the number of decimal places here
    use_seps = FALSE
  ) %>%
  cols_label(
    Parameter = "Parameter",
    Mean = "Effect size",
    SD = "SD",
    Q025 = "Lower 95% CI",
    Q975 =  "Upper 95% CI",
    Mode = "Mode"
  ) %>%
  gtsave(filename = "../../plots/stats-tables/inla-max-summaries.html")



# peak-NDVI timing ###
# Fit variograms, with cutoff set to 1/3 the diagonal of the dataset
fit_variogram_doy(s2.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.bl)/3) # 31.1
fit_variogram_doy(s2.kl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kl)/3) # 52.2
fit_variogram_doy(s2.kh, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s2.kh)/3) # 30.6
fit_variogram_doy(s30.bl, model.type = 'Mat', grid.size = 10, cutoff.dist = get_diagonal(s30.bl)/3)# 21.8

# Set empirical ranges from variograms above
doy_ranges_df <- data.frame(site = c("s2.bl",
                                 "s2.kl",
                                 "s2.kh",
                                 "s30.bl"),
                        range = c(round(31.1/10),
                                  round(52.2/10),
                                  round(30.6/10), 
                                  round(21.8/30))
)

# Fit models
doy.s2.bl_fit <- fit_matern_doy(scale_param = 0.01,
                                shape_param = 1,
                                range_site = doy_ranges_df$range[1], 
                                site_data = s2.bl,
                                return_model = T)
doy.s2.kl_fit <- fit_matern_doy(scale_param = 0.01,
                                shape_param = 1,
                                range_site = doy_ranges_df$range[2], 
                                site_data = s2.kl,
                                return_model = T)
doy.s2.kh_fit <- fit_matern_doy(scale_param = 0.01,
                                shape_param = 1,
                                range_site = doy_ranges_df$range[3], 
                                site_data = s2.kh,
                                return_model = T)
doy.s30.bl_fit <- fit_matern_doy(scale_param = 0.01,
                                 shape_param = 1,
                                 range_site = doy_ranges_df$range[4], 
                                 site_data = s30.bl,
                                 return_model = T)

summary(doy.s2.bl_fit)




# Extract and summarise fixed effects
doy.s2.bl.mat <- summarise_matern(doy.s2.bl_fit, site.name = "Blaesedalen, Sentinel-2", sensor.name = "Sentinel-2")
doy.s2.kl.mat <- summarise_matern(doy.s2.kl_fit, site.name = "Kluane Low", sensor.name = "Sentinel-2")
doy.s2.kh.mat <- summarise_matern(doy.s2.kh_fit, site.name = "Kluane High", sensor.name = "Sentinel-2")
doy.s30.bl.mat <- summarise_matern(doy.s30.bl_fit, site.name = "Blaesedalen, HLS S30", sensor.name = "HLS S30")

doy.all.mat <- rbind(doy.s2.kl.mat,
                     doy.s2.kh.mat,
                     doy.s2.bl.mat,
                     doy.s30.bl.mat)


doy.summary <- gt(doy.all.mat %>% group_by(site) %>% select(-sensor)) %>%
  fmt_number(
    columns = vars(Mean, SD, Q025, Q975, Mode),
    decimals = 5, # Set the number of decimal places here
    use_seps = FALSE
  ) %>%
  cols_label(
    Parameter = "Parameter",
    Mean = "Effect size",
    SD = "SD",
    Q025 = "Lower 95% CI",
    Q975 =  "Upper 95% CI",
    Mode = "Mode"
  ) %>%
  gtsave(filename = "../../plots/stats-tables/inla-doy-summaries.html")






# Fit one final model, adding a breakpoint at snow.auc == 5 for S2 BL peak-NDVI magnitude
# Break point was manually identified

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


# Extract and summarise fixed effects
max.s2.bl.lte5.matbreak <- summarise_matern(s2.bl_fit_lte5, site.name = "Blaesedalen, Sentinel-2, snow persistence <= 5", sensor.name = "Sentinel-2")
max.s2.bl.gt5.matbreak <- summarise_matern(s2.bl_fit_gt5, site.name = "Blaesedalen, Sentinel-2, snow persistence > 5", sensor.name = "Sentinel-2")

max.s2.bl.breakpoint <- rbind(max.s2.bl.lte5.matbreak,
                              max.s2.bl.gt5.matbreak)


max.bl.break.summary <- gt(max.s2.bl.breakpoint %>% group_by(site) %>% select(-sensor)) %>%
  fmt_number(
    columns = vars(Mean, SD, Q025, Q975, Mode),
    decimals = 5, # Set the number of decimal places here
    use_seps = FALSE
  ) %>%
  cols_label(
    Parameter = "Parameter",
    Mean = "Effect size",
    SD = "SD",
    Q025 = "Lower 95% CI",
    Q975 =  "Upper 95% CI",
    Mode = "Mode"
  ) %>%
  gtsave(filename = "../../plots/stats-tables/inla-bl-max-break-summaries.html")





















