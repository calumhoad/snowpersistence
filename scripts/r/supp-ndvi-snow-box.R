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

# Data load
bl <- read_csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv')
kl <- read_csv('../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv')
kh <- read_csv('../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv')
bls30 <- read_csv('../../data/combined-ndvi-snow/s30-bl-smooth-joined.csv')


# Get the upper, lower and mid values for peak NDVI
categorise_ndvi <- function(data) {
  data.filtered <- data %>% 
    group_by(id) %>% 
    filter(row_number() == 10 & !is.na(ndvi.max) & snow.auc != 0)
  
  data.cat <- data.filtered %>% 
    mutate(category = factor(ifelse(ndvi.max < 0.2, '[0 - 0.2)',
                                    ifelse(ndvi.max >= 0.2 & ndvi.max < 0.4, '[0.2 - 0.4)',
                                           ifelse(ndvi.max >= 0.4 & ndvi.max < 0.6, '[0.4 - 0.6)',
                                                  '[0.6 - 1)'))),
                             levels = c('[0 - 0.2)', '[0.2 - 0.4)', '[0.4 - 0.6)', '[0.6 - 1)')))
  return(data.cat)
}


bl <- categorise_ndvi(bl)
kl <- categorise_ndvi(kl)
kh <- categorise_ndvi(kh)
bls30 <- categorise_ndvi(bls30)

# Plotting function
ndvi_box <- function(data, label.step, label.max) {
  out.plot <- ggplot() +
    geom_boxplot(data = data, aes(x = category, y = snow.auc, fill = category)) +
    scale_fill_manual(values = c('[0 - 0.2)' = "#9AE091", '[0.2 - 0.4)' = "#54E841", '[0.4 - 0.6)' = "#399C2C", '[0.6 - 1)' = "#205819"),
                      guide = 'none') +
    scale_y_continuous(breaks = c(seq(0, label.max, label.step)),
                       labels = c(as.character(seq(0, label.max, label.step)))) +
    xlab("peak NDVI bracket") +
    ylab("Snow persistence") +
    coord_cartesian(ylim = c(0, label.max + 1)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


bl.p <- ndvi_box(bl, label.max = 25, label.step = 5)
kl.p <- ndvi_box(kl, label.max = 15, label.step = 5)
kh.p <- ndvi_box(kh, label.max = 15, label.step = 5)
bl30.p <- ndvi_box(bls30, label.max = 15, label.step = 5)

combined <- plot_grid(kl.p, kh.p, bl.p, bl30.p, nrow = 2, ncol = 2)#, labels = c('KL', 'KH', 'BL', 'HLSS30 BL'))

cowplot::save_plot('../../plots/supplementary/veg-boxplots.png', combined, 
                   base_height = 180, base_width = 140, bg = 'white', units = 'mm')

