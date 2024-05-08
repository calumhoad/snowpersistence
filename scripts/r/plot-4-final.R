# Plot 3 (peak-NDVI magnitude ~ snow persistence)
# Calum Hoad, 08/05/2024
# Code adapted from Jakob J. Assmann

# Dependencies
library(tidyverse)
library(sf)
library(ggplot2)
library(cowplot)
library(pbapply)
library(INLA)
library(gt)
library(ggnewscale)

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
data_list <- list(s2.kl, s2.kh, s2.bl)


# plotting ----

# Generate lm plots without zeros for nvdi.max and log link of kh
fit_plot <- function(site_data,
                             colour.site, 
                             colour.darker, 
                             colour.lighter, 
                             colour.lightest, 
                             ymax, ymin){
  lm_plot <- ggplot(site_data, aes(x = snow.auc, y = ndvi.max.doy)) +
    geom_point(aes(x = snow.auc, y = ndvi.max.doy), colour = colour.lightest) +
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
                 ymax = 240, ymin = 215)
kl <- fit_plot(s2.kl, 
               colour.site = '#F5A40C', 
               colour.darker = '#946606',
               colour.lighter = '#FBCA7F',
               colour.lightest = '#FDDCAC', 
               ymax = 235, ymin = 210)
kh <- fit_plot(s2.kh, 
               colour.site = '#F23835', 
               colour.darker = '#8D271E',
               colour.lighter = '#F29580',
               colour.lightest = '#F8BBAA', 
               ymax = 235, ymin = 210)


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


# Plot out the full panel
combined.2 <- plot_grid(kl, kh, bl, time.plot, ncol = 2, nrow = 2, align = 'hv')
combined.2
#Check the plot
combined 

# Save plots
cowplot::save_plot('../../plots/figures/figure-4-final.png', combined.2, 
                   base_height = 140, base_width = 180, units = 'mm', 
                   bg = 'white')
