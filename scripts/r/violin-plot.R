library(ggplot2)
library(tidyverse)
library(lubridate)
library(cowplot)

# Load in the data
bl.snow <- read_csv('../../data/snow/snow-cover-10m-blaesedalen.csv') %>%
  mutate(site = 'BL') %>%
  select(-X & -Y) %>%
  pivot_longer(!id & !snow.auc & !site & !snow.av, names_to = 'date', values_to = 'snow.pcnt')


kl.snow <- read_csv('../../data/snow/snow-cover-10m-kluane-low.csv') %>%
  mutate(site = 'KL') %>%
  select(-X & -Y) %>%
  pivot_longer(!id & !snow.auc & !site & !snow.av, names_to = 'date', values_to = 'snow.pcnt')

kh.snow <- read_csv('../../data/snow/snow-cover-10m-kluane-high.csv') %>%
  mutate(site = 'KH') %>%
  select(-X & -Y) %>%
  pivot_longer(!id & !snow.auc & !site & !snow.av, names_to = 'date', values_to = 'snow.pcnt')

snow.all <- rbind(bl.snow, kl.snow, kh.snow)


# Plot the data
plot_grid(
  ggplot(snow.all %>% filter(site == 'BL') %>% filter(snow.pcnt != 0), aes(x = yday(date), y = snow.pcnt, group = date), fill = 'purple') +
    geom_violin(trim = TRUE) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    labs(x = "Day of Year", y = "Snow cover percentage") +
    theme(legend.position = "bottom"),
  ggplot(snow.all %>% filter(site == 'KL'), aes(x = yday(date), y = snow.pcnt, group = date), fill = 'blue') +
    geom_violin(trim = TRUE) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    labs(x = "Day of Year", y = "Snow cover percentage") +
    theme(legend.position = "bottom"),
  ggplot(snow.all %>% filter(site == 'KH'), aes(x = yday(date), y = snow.pcnt, group = date), fill = 'red') +
    geom_violin(trim = TRUE) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    labs(x = "Day of Year", y = "Snow cover percentage") +
    theme(legend.position = "bottom"),
  nrow = 3
)
  

