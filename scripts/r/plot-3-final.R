# Plot 3 (peak-NDVI magnitude ~ snow persistence)
# Calum Hoad, 08/05/2024
# Code adapted from Jakob J. Assmann

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


# plotting ----

# Generate lm plots without zeros for nvdi.max and log link of kh
fit_plot <- function(site_data,
                     colour.site, 
                     colour.darker, 
                     colour.lighter, 
                     colour.lightest, 
                     ymax, ymin,
                     xmax, xmin, xint){
  lm_plot <- ggplot(site_data, aes(x = snow.auc, y = ndvi.max)) +
    geom_point(aes(x = snow.auc, y = ndvi.max), colour = colour.lightest) +
    scale_y_continuous(breaks = c(seq(0, ymax, 0.1)), 
                       labels = c(as.character(seq(0, ymax, 0.1)))) +
    scale_x_continuous(breaks = c(seq(xmin, xmax, xint)), 
                       labels = c(as.character(seq(xmin, xmax, xint)))) +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
    theme_cowplot()
  if(unique(site_data$site) == "S2 KH") {
    lm_plot <- lm_plot +
      geom_smooth(method = "lm", formula = y ~ log(x + 1), 
                  colour = colour.site,
                  fill = colour.site, 
                  linewidth = 1)# +
    #labs(x = "snow persistence")#, y = "ndvi.max.doy", 
    #title = paste0(unique(site_data$site), " (lm y ~ ln(x + 1))")) 
  } else {
    lm_plot <- lm_plot +
      geom_smooth(method = "lm", formula = y ~ x, 
                  colour = colour.site,
                  fill = colour.site, 
                  linewidth = 1) #+
    #labs(x = "snow persistence", y = "ndvi.max.doy")#, 
    #title = paste0(unique(site_data$site), " (lm y ~ x)")) 
  }
  return(lm_plot)
}

bl <- fit_plot(s2.bl, 
               colour.site = '#4984BF', 
               colour.darker = '#2E5277',
               colour.lighter = '#9BB2DA',
               colour.lightest = '#BECBE7', 
               ymax = 0.6, ymin = 0.1,
               xmin = 0, xmax = 25, xint = 5)
kl <- fit_plot(s2.kl, 
               colour.site = '#F5A40C', 
               colour.darker = '#946606',
               colour.lighter = '#FBCA7F',
               colour.lightest = '#FDDCAC', 
               ymax = 0.7, ymin = 0.4,
               xmin = 0, xmax = 15, xint = 5)
kh <- fit_plot(s2.kh, 
               colour.site = '#F23835', 
               colour.darker = '#8D271E',
               colour.lighter = '#F29580',
               colour.lightest = '#F8BBAA', 
               ymax = 0.6, ymin = 0.1,
               xmin = 0, xmax = 15, xint = 5)


# Plot out the full panel
combined.2 <- plot_grid(kl, kh, bl, ncol = 2, nrow = 2, align = 'hv')
combined.2
#Check the plot
combined 

# Save plots
cowplot::save_plot('../../plots/figures/figure-3-final-axiscorrect.png', combined.2, 
                   base_height = 140, base_width = 180, units = 'mm', 
                   bg = 'white')
