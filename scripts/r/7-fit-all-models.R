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

source('b-modelling-helper-functions.R')

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





















