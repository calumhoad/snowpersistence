# R script to fit final INLA models for Calum
# Jakob J. Assmann jakob.assmann@uzh.ch 1 May 2024

# Dependencies
library(tidyverse)
library(sf)
library(ggplot2)
library(cowplot)
library(pbapply)
library(inla)
library(gt)


# Set wd
setwd("C:/Users/jakob/repositories/c1-analyses/scripts/jakob/")

### Data Prep ##########

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
data_list <- list(s2.bl, s2.kh, s2.kl, s30.bl)

# Simple correlation tables
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
gtsave("spearmans_table.png")

# Generate binned scatter plots with zeros in separate categories
hist_plots_max.ndvi <- map(data_list, function(site_data){
  # Set bin width
  bin_width <- 2
  if (unique(site_data$site) == "S2 BL") bin_width <- 4
  # Remove geoms for easy handling
  site_data <- st_drop_geometry(site_data)
  # Separate and lable zeros
  zeros <- filter(site_data, snow.auc == 0)
  zeros$snow_auc_cat <- -1
  zeros$n <- nrow(zeros)
  zeros$scale_name <- "0"
  # Bin remainder and label
  site_data <- filter(site_data, snow.auc != 0) %>%
    mutate(snow_auc_cat = round(snow.auc / bin_width))
  site_data <- site_data %>%
    group_by(snow_auc_cat) %>%
    tally() %>%
    mutate(scale_name = paste0(
      snow_auc_cat, "-",
      as.numeric(snow_auc_cat) + bin_width
    )) %>%
    full_join(site_data, .)
  # Merge zeros again and convert factors where needed
  site_data <- bind_rows(site_data, zeros) %>%
      arrange(snow_auc_cat) %>%
      mutate(snow_auc_cat = factor(snow_auc_cat))
  # Plot and return
  ggplot() +
    geom_jitter(aes(x = snow_auc_cat, y = ndvi.max, colour = snow_auc_cat), data = site_data) +
    geom_text(aes(x = snow_auc_cat, y = Inf, label = paste0("n = ", n)), vjust = 1.5, 
                  data = filter(distinct(site_data, snow_auc_cat, n))) +
    scale_x_discrete(labels = unique(site_data$scale_name)) +
        labs(x = "snow persistance (numerical range)", y = "ndvi.max", 
              title = paste0(unique(site_data$site), " (lm, no 0s)")) +
    theme_cowplot() +
    theme(legend.position = "none")
})
plot_grid(plotlist = hist_plots_max.ndvi, 
          labels = paste0("(", letters[1:4], ")"),
          label_size = 18)  %>%
    save_plot("binned_scatter_max.ndvi.png", ., nrow = 2, ncol = 2, bg = "white")

# Same for max.ndvi.doy
hist_plots_max.ndvi.doy <- map(data_list, function(site_data){
  # Set bin width
  bin_width <- 2
  if(unique(site_data$site) == "S2 BL") bin_width <- 4
  # Remove geom
  site_data <- st_drop_geometry(site_data)
  # Separate and label zeros
  zeros <- filter(site_data, snow.auc == 0)
  zeros$snow_auc_cat <- -1
  zeros$n <- nrow(zeros)
  zeros$scale_name <- "0"
  # Bin remainder of the data and label
  site_data <- filter(site_data, snow.auc != 0) %>%
    mutate(snow_auc_cat = round(snow.auc / bin_width))
  site_data <- site_data %>%
    group_by(snow_auc_cat) %>%
    tally() %>%
    mutate(scale_name = paste0(
      snow_auc_cat, "-",
      as.numeric(snow_auc_cat) + bin_width
    )) %>%
    full_join(site_data, .)
  # Merge again with zeros
  site_data <- bind_rows(site_data, zeros) %>%
      arrange(snow_auc_cat) %>%
      mutate(snow_auc_cat = factor(snow_auc_cat))
  # Plot and return
  ggplot() +
    geom_jitter(aes(x = snow_auc_cat, y = ndvi.max.doy, colour = snow_auc_cat), data = site_data) +
    geom_text(aes(x = snow_auc_cat, y = Inf, label = paste0("n = ", n)), vjust = 1.5, 
                  data = filter(distinct(site_data, snow_auc_cat, n))) +
    scale_x_discrete(labels = unique(site_data$scale_name)) +
        labs(x = "snow persistance (numerical range)", y = "ndvi.max.doy", 
              title = paste0(unique(site_data$site), " (lm, no 0s)")) +
    theme_cowplot() +
    theme(legend.position = "none")
})
plot_grid(plotlist = hist_plots_max.ndvi.doy, 
          labels = paste0("(", letters[1:4], ")"),
          label_size = 18)  %>%
    save_plot("binned_scatter_max.ndvi.doy.png", ., nrow = 2, ncol = 2, bg = "white")

# Generate lm plots with zeros for nvdi.max
ndvi_max_plots <- map(data_list, function(site_data){
  ggplot(site_data, aes(x = snow.auc, y = ndvi.max)) +
        geom_point(colour = unique(site_data$colour)) +
        geom_smooth(method = "lm", formula = y ~ x, 
                    colour = unique(site_data$colour),
                    fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max", 
              title = paste0(unique(site_data$site), " (lm)")) +
        theme_cowplot()
})
plot_grid(plotlist = ndvi_max_plots, 
          labels = paste0("(", letters[1:4], ")"),
          label_size = 18)  %>%
    save_plot("ndvi_max_plots_lm.png", ., nrow = 2, ncol = 2, bg = "white")

# Generate lm plots with zeros for nvdi.max
max_doy_plots <- map(data_list, function(site_data){
    ggplot(site_data, aes(x = snow.auc, y = ndvi.max.doy)) +
        geom_point(colour = unique(site_data$colour)) +
        geom_smooth(method = "lm", formula = y ~ x, 
                    colour = unique(site_data$colour),
                    fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max", 
            title = paste0(unique(site_data$site), " (lm)")) +
        theme_cowplot()
})
plot_grid(plotlist = max_doy_plots, 
          labels = paste0("(", letters[1:4], ")"),
          label_size = 18) %>%
    save_plot("max_doy_plots_lm.png", ., nrow = 2, ncol = 2, bg = "white")


# Generate lm plots without zeros for nvdi.max and log link of kh
ndvi_max_plots_ln <- map(data_list, function(site_data){
    lm_plot <- ggplot(site_data, aes(x = snow.auc, y = ndvi.max)) +
        geom_point(colour = unique(site_data$colour)) +
        theme_cowplot()
    if(unique(site_data$site) == "S2 KH") {
      lm_plot <- lm_plot +
        geom_smooth(method = "lm", formula = y ~ log(x + 1), 
                    colour = unique(site_data$colour),
                   fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max", 
             title = paste0(unique(site_data$site), " (lm y ~ ln(x + 1), no 0s)")) 
    } else {
      lm_plot <- lm_plot +
        geom_smooth(method = "lm", formula = y ~ x, 
                    colour = unique(site_data$colour),
                    fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max", 
            title = paste0(unique(site_data$site), " (lm y ~ x, no 0s)")) 
    }
    return(lm_plot)
})
plot_grid(plotlist = ndvi_max_plots_ln, 
          labels = paste0("(", letters[1:4], ")"),
          label_size = 18)  %>%
    save_plot("ndvi_max_plots_lm_kh-ln.png", ., nrow = 2, ncol = 2, bg = "white")

# Generate lm plots without zeros for nvdi.max and log link of kh
max_doy_plots_ln <- map(data_list, function(site_data){
    lm_plot <- ggplot(site_data, aes(x = snow.auc, y = ndvi.max.doy)) +
        geom_point(colour = unique(site_data$colour)) +
        theme_cowplot()
    if(unique(site_data$site) == "S2 KH") {
      lm_plot <- lm_plot +
        geom_smooth(method = "lm", formula = y ~ log(x + 1), 
                    colour = unique(site_data$colour),
                   fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max.doy", 
             title = paste0(unique(site_data$site), " (lm y ~ ln(x + 1), no 0s)")) 
    } else {
      lm_plot <- lm_plot +
        geom_smooth(method = "lm", formula = y ~ x, 
                    colour = unique(site_data$colour),
                    fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max.doy", 
            title = paste0(unique(site_data$site), " (lm y ~ x, no 0s)")) 
    }
    return(lm_plot)
})
plot_grid(plotlist = max_doy_plots_ln, 
          labels = paste0("(", letters[1:4], ")"),
          label_size = 18)  %>%
    save_plot("max_doy_plots_lm_kh-ln.png", ., nrow = 2, ncol = 2, bg = "white")



# Generate lm plots without zeros for nvdi.max
ndvi_max_plots_no <- map(data_list, function(site_data){
    site_data <- filter(site_data, snow.auc != 0)
    if(unique(site_data$site) == "S2 BL") site_data <- filter(site_data, snow.auc > 5)
    ggplot(site_data, aes(x = snow.auc, y = ndvi.max)) +
        geom_point(colour = unique(site_data$colour)) +
        geom_smooth(method = "lm", formula = y ~ x, 
                    colour = unique(site_data$colour),
                    fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max", 
              title = paste0(unique(site_data$site), " (lm, no 0s)")) +
        theme_cowplot()
})
plot_grid(plotlist = ndvi_max_plots_no, 
          labels = paste0("(", letters[1:4], ")"),
          label_size = 18)  %>%
    save_plot("no_zero_ndvi_max_plots_lm.png", ., nrow = 2, ncol = 2, bg = "white")

# Generate lm plots without zeros for nvdi.max
max_doy_plots_no <- map(data_list, function(site_data){
    site_data <- filter(site_data, snow.auc != 0)
    if(unique(site_data$site) == "S2 BL") site_data <- filter(site_data, snow.auc > 5)
    ggplot(site_data, aes(x = snow.auc, y = ndvi.max.doy)) +
        geom_point(colour = unique(site_data$colour)) +
        geom_smooth(method = "lm", formula = y ~ x, 
                    colour = unique(site_data$colour),
                    fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max", 
            title = paste0(unique(site_data$site), " (lm, no 0s)")) +
        theme_cowplot()
})
plot_grid(plotlist = max_doy_plots_no, 
          labels = paste0("(", letters[1:4], ")"),
          label_size = 18) %>%
    save_plot("no_zero_max_doy_plots_lm.png", ., nrow = 2, ncol = 2, bg = "white")

# Generate lm plots without zeros for nvdi.max and log link of kh
ndvi_max_plots_ln_no <- map(data_list, function(site_data){
    site_data <- filter(site_data, snow.auc != 0)
    if(unique(site_data$site) == "S2 BL") site_data <- filter(site_data, snow.auc > 5)
    lm_plot <- ggplot(site_data, aes(x = snow.auc, y = ndvi.max)) +
        geom_point(colour = unique(site_data$colour)) +
        theme_cowplot()
    if(unique(site_data$site) == "S2 KH") {
      lm_plot <- lm_plot +
        geom_smooth(method = "lm", formula = y ~ log(x), 
                    colour = unique(site_data$colour),
                   fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max", 
             title = paste0(unique(site_data$site), " (lm y ~ ln(x), no 0s)")) 
    } else {
      lm_plot <- lm_plot +
        geom_smooth(method = "lm", formula = y ~ x, 
                    colour = unique(site_data$colour),
                    fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max", 
            title = paste0(unique(site_data$site), " (lm y ~ x, no 0s)")) 
    }
    return(lm_plot)
})
plot_grid(plotlist = ndvi_max_plots_ln_no, 
          labels = paste0("(", letters[1:4], ")"),
          label_size = 18)  %>%
    save_plot("no_zero_ndvi_max_plots_lm_kh-ln.png", ., nrow = 2, ncol = 2, bg = "white")

# Generate lm plots without zeros for nvdi.max and log link of kh
max_doy_plots_ln_no <- map(data_list, function(site_data){
    site_data <- filter(site_data, snow.auc != 0)
    if(unique(site_data$site) == "S2 BL") site_data <- filter(site_data, snow.auc > 5)
    lm_plot <- ggplot(site_data, aes(x = snow.auc, y = ndvi.max.doy)) +
        geom_point(colour = unique(site_data$colour)) +
        theme_cowplot()
    if(unique(site_data$site) == "S2 KH") {
      lm_plot <- lm_plot +
        geom_smooth(method = "lm", formula = y ~ log(x), 
                    colour = unique(site_data$colour),
                   fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max.doy", 
             title = paste0(unique(site_data$site), " (lm y ~ ln(x), no 0s)")) 
    } else {
      lm_plot <- lm_plot +
        geom_smooth(method = "lm", formula = y ~ x, 
                    colour = unique(site_data$colour),
                    fill = unique(site_data$colour)) +
        labs(x = "snow persistence", y = "ndvi.max.doy", 
            title = paste0(unique(site_data$site), " (lm y ~ x, no 0s)")) 
    }
    return(lm_plot)
})
plot_grid(plotlist = max_doy_plots_ln_np, 
          labels = paste0("(", letters[1:4], ")"),
          label_size = 18)  %>%
    save_plot("no_zero_max_doy_plots_lm_kh-ln.png", ., nrow = 2, ncol = 2, bg = "white")

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

# Helper function to fit matern model
fit_matern_ndvi.max <- function(site_data) {
    # Get range
    range_site <- unique(site_data$range)
    # Define max number of rows and cols
    n_col = max(site_data$col_number)
    n_row = max(site_data$row_number)

    # Specify priors for hyperparemeter
    log.range = list(initial = log(range_site), fixed = TRUE)
    hyperpar_matern = list(initial = 2, param = c(1, 0.01))
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
    # Return model
    return(model_matern)
}

# Helper function to fit matern model
fit_matern_ndvi.max.doy <- function(site_data) {
    # Get range
    range_site <- unique(site_data$range)
    # Define max number of rows and cols
    n_col = max(site_data$col_number)
    n_row = max(site_data$row_number)

    # Specify priors for hyperparemeter
    log.range = list(initial = log(range_site), fixed = TRUE)
    hyperpar_matern = list(initial = 2, param = c(1, 0.01))
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
    # Return model
    return(model_matern)
}

# Helper function to generate plots
plot_results_ndvi.max <- function(model_matern, site_data){
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
    geom_point(aes(x = snow.auc, y= ndvi.max),
                colour = unique(site_data$colour)) +
    geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                    ymax = max_ndvi_pred),
                data = preds_credible,
                fill = unique(site_data$colour), alpha = 0.3) +
    geom_line(aes(x = snow.auc, y = preds), colour = unique(site_data$colour)) +
    labs(x = "snow persistence", y = "ndvi.max.doy", 
            title = paste0(unique(site_data$site), " (INLA_matern2d y ~ x, no 0s)")) +
    theme_cowplot()

  return(fit_plot)
}

# Helper function to generate plots
plot_results_ndvi.max.doy <- function(model_matern, site_data){
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
    geom_point(aes(x = snow.auc, y= ndvi.max.doy), 
              colour = unique(site_data$colour)) +
    geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                    ymax = max_ndvi_pred),
                data = preds_credible,
                fill = unique(site_data$colour), alpha = 0.3) +
    geom_line(aes(x = snow.auc, y = preds), colour = unique(site_data$colour)) +
    labs(x = "snow persistence", y = "ndvi.max.doy", 
            title = paste0(unique(site_data$site), " (INLA_matern2d y ~ x, no 0s)")) +
    theme_cowplot()

  return(fit_plot)
}

# Generate R INLA plots without zeros for nvdi.max
ndvi_max_plots_inla <- map(data_list, function(site_data){
    # Add INLA grid
    site_data <- get_inla_grid(site_data)
    # Filter 0 values
    site_data <- filter(site_data, snow.auc != 0)
    if(unique(site_data$site) == "S2 BL") site_data <- filter(site_data, snow.auc > 5)
    # Fit model
    site_model <- fit_matern_ndvi.max(site_data)
    # Print model summary
    print(summary(site_model))
    # Generate plot and return
    plot_results_ndvi.max(site_model, site_data)
})

plot_grid(plotlist = ndvi_max_plots_inla) %>%
    save_plot("no_zero_ndvi_max_plots_INLA.png", ., nrow = 2, ncol = 2, bg = "white")

# Generate R INLA plots without zeros for ndvi.max.doy
max_doy_plots_inla <- map(data_list, function(site_data){
    # Add INLA grid
    site_data <- get_inla_grid(site_data)
    # Filter 0 values
    site_data <- filter(site_data, snow.auc != 0)
    if(unique(site_data$site) == "S2 BL") site_data <- filter(site_data, snow.auc > 5)
    # Fit model
    site_model <- fit_matern_ndvi.max.doy(site_data)
        # Print model summary
    print(summary(site_model))
    # Generate plot and return
    plot_results_ndvi.max.doy(site_model, site_data)
})

plot_grid(plotlist = max_doy_plots_inla) %>%
    save_plot("no_zero_max_doy_plots_INLA.png", ., nrow = 2, ncol = 2, bg = "white")

# Finally INLA matern2d + ln for KH

# Helper function to fit matern model
fit_matern_ndvi.max_ln <- function(site_data) {
    # Get range
    range_site <- unique(site_data$range)
    # Define max number of rows and cols
    n_col = max(site_data$col_number)
    n_row = max(site_data$row_number)

    # Specify priors for hyperparemeter
    log.range = list(initial = log(range_site), fixed = TRUE)
    hyperpar_matern = list(initial = 2, param = c(1, 0.01))
    # Specify formula
    formula_matern <- ndvi.max ~ log(snow.auc) +
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
    # Return model
    return(model_matern)
}

# Helper function to fit matern model
fit_matern_ndvi.max.doy_ln <- function(site_data) {
    # Get range
    range_site <- unique(site_data$range)
    # Define max number of rows and cols
    n_col = max(site_data$col_number)
    n_row = max(site_data$row_number)

    # Specify priors for hyperparemeter
    log.range = list(initial = log(range_site), fixed = TRUE)
    hyperpar_matern = list(initial = 2, param = c(1, 0.01))
    # Specify formula
    formula_matern <- ndvi.max.doy ~ log(snow.auc) +
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
    # Return model
    return(model_matern)
}

# Helper function to generate plots
plot_results_ndvi.max_ln <- function(model_matern, site_data){
  # Calculate marginal predictions depending on snow.auc only
  site_data$preds <- model_matern$summary.fixed[1, 1] +
    log(site_data$snow.auc) * model_matern$summary.fixed[2, 1]
  
  # Calculate min and max values for 95 CI predictions to plot a ribbon
  preds_credible <- expand.grid(intercept = unlist(model_matern$summary.fixed[1, c(1,3,5)]), 
              slope = unlist(model_matern$summary.fixed[2, c(1,3,5)])) %>%
    split(1:nrow(.)) %>%
    map(function(parameters) {
      data.frame(snow.auc =  seq((floor(min(site_data$snow.auc) * 10) / 10),
                                 (ceiling(max(site_data$snow.auc) * 10) / 10), 0.01),
                 ndvi_pred = parameters$intercept + log(seq((floor(min(site_data$snow.auc) * 10) / 10),
                                                        (ceiling(max(site_data$snow.auc) * 10) / 10), 0.01)) * parameters$slope[1])
      }) %>% bind_rows() %>%
    group_by(snow.auc) %>% summarise(min_ndvi_pred = min(ndvi_pred), max_ndvi_pred = max(ndvi_pred))
  
  # Add model fitted values and residuals to data frame
  site_data$fitted <- model_matern$summary.fitted.values$mean
  site_data$residuals <- site_data$ndvi.max - site_data$fitted
  
  # Plot model fit
  fit_plot <- ggplot(data = site_data) +
    geom_point(aes(x = snow.auc, y= ndvi.max),
                colour = unique(site_data$colour)) +
    geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                    ymax = max_ndvi_pred),
                data = preds_credible,
                fill = unique(site_data$colour), alpha = 0.3) +
    geom_line(aes(x = snow.auc, y = preds), colour = unique(site_data$colour)) +
    labs(x = "snow persistence", y = "ndvi.max.doy", 
            title = paste0(unique(site_data$site), " (INLA_matern2d y ~ ln(x), no 0s)")) +
    theme_cowplot()

  return(fit_plot)
}

# Helper function to generate plots
plot_results_ndvi.max.doy_ln <- function(model_matern, site_data){
  # Calculate marginal predictions depending on snow.auc only
  site_data$preds <- model_matern$summary.fixed[1, 1] +
    log(site_data$snow.auc) * model_matern$summary.fixed[2, 1]
  
  # Calculate min and max values for 95 CI predictions to plot a ribbon
  preds_credible <- expand.grid(intercept = unlist(model_matern$summary.fixed[1, c(1,3,5)]), 
              slope = unlist(model_matern$summary.fixed[2, c(1,3,5)])) %>%
    split(1:nrow(.)) %>%
    map(function(parameters) {
      data.frame(snow.auc =  seq((floor(min(site_data$snow.auc) * 10) / 10),
                                 (ceiling(max(site_data$snow.auc) * 10) / 10), 0.01),
                 ndvi_pred = parameters$intercept + log(seq((floor(min(site_data$snow.auc) * 10) / 10),
                                                        (ceiling(max(site_data$snow.auc) * 10) / 10), 0.01)) * parameters$slope[1])
      }) %>% bind_rows() %>%
    group_by(snow.auc) %>% summarise(min_ndvi_pred = min(ndvi_pred), max_ndvi_pred = max(ndvi_pred))
  
  # Add model fitted values and residuals to data frame
  site_data$fitted <- model_matern$summary.fitted.values$mean
  site_data$residuals <- site_data$ndvi.max - site_data$fitted
  
  # Plot model fit
  fit_plot <- ggplot(data = site_data) +
    geom_point(aes(x = snow.auc, y= ndvi.max.doy), 
              colour = unique(site_data$colour)) +
    geom_ribbon(aes(x = snow.auc, ymin = min_ndvi_pred,
                    ymax = max_ndvi_pred),
                data = preds_credible,
                fill = unique(site_data$colour), alpha = 0.3) +
    geom_line(aes(x = snow.auc, y = preds), colour = unique(site_data$colour)) +
    labs(x = "snow persistence", y = "ndvi.max.doy", 
            title = paste0(unique(site_data$site), " (INLA_matern2d y ~ ln(x), no 0s)")) +
    theme_cowplot()

  return(fit_plot)
}

# Generate R INLA plots without zeros for nvdi.max
ndvi_max_plots_inla_ln <- map(data_list, function(site_data){
    # Add INLA grid
    site_data <- get_inla_grid(site_data)
    # Filter 0 values
    site_data <- filter(site_data, snow.auc != 0)
    if(unique(site_data$site) == "S2 BL") site_data <- filter(site_data, snow.auc > 5)
    # For KL do ln fit
    if(unique(site_data$site) == "S2 KH"){
    # Fit model
    site_model <- fit_matern_ndvi.max_ln(site_data)
    # Print model summary
    print(summary(site_model))
    # Generate plot and return
    return(plot_results_ndvi.max_ln(site_model, site_data))
    } else{
          # Fit model
    site_model <- fit_matern_ndvi.max(site_data)
    # Print model summary
    print(summary(site_model))
    # Generate plot and return
    return(plot_results_ndvi.max(site_model, site_data))
    }
})

plot_grid(plotlist = ndvi_max_plots_inla_ln) %>%
    save_plot("no_zero_ndvi_max_plots_INLA_ln.png", ., nrow = 2, ncol = 2, bg = "white")

# Generate R INLA plots without zeros for ndvi.max.doy
max_doy_plots_inla_ln <- map(data_list, function(site_data){
    # Add INLA grid
    site_data <- get_inla_grid(site_data)
    # Filter 0 values
    site_data <- filter(site_data, snow.auc != 0)
    if(unique(site_data$site) == "S2 BL") site_data <- filter(site_data, snow.auc > 5)
        # For KL do ln fit
    if(unique(site_data$site) == "S2 KH"){
      # Fit model
    site_model <- fit_matern_ndvi.max.doy_ln(site_data)
        # Print model summary
    print(summary(site_model))
    # Generate plot and return
    return(plot_results_ndvi.max.doy_ln(site_model, site_data))
    } else {
         site_model <- fit_matern_ndvi.max.doy(site_data)
        # Print model summary
    print(summary(site_model))
    # Generate plot and return
    return(plot_results_ndvi.max.doy(site_model, site_data))
    }
})

plot_grid(plotlist = max_doy_plots_inla_ln) %>%
    save_plot("no_zero_max_doy_plots_INLA_ln.png", ., nrow = 2, ncol = 2, bg = "white")
