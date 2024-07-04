# Plot 2 for paper - NDVI max is influenced by snow persistence
# Calum Hoad, 1 Feb 2024

library(ggplot2)
library(dplyr)
library(sf)
library(mapview)
library(spdep)
library(spatialreg)
library(stargazer)
library(cowplot)

# Data ----
# Read in the data and reduce to unique id
# Blaesedalen
s2.bl <- read.csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Kluane low
s2.kl <- read.csv('../../data/combined-ndvi-snow/s2-kl-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Kluane high
s2.kh <- read.csv('../../data/combined-ndvi-snow/s2-kh-smooth-joined.csv') %>%
  st_as_sf(coords = c('X', 'Y'), crs = 32608) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup()

# Models ----

###
# Blaesedalen Spatial Error Models
###

s2.bl.nb.10 <- dnearneigh(s2.bl, d1 = 0, d2 = 10)
#s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s2.bl.lw.10 <- nb2listw(s2.bl.nb.10, style = 'W', zero.policy = FALSE)

sem.s2.bl.max.10 <- errorsarlm(ndvi.max ~ snow.auc,
                            data = s2.bl, 
                            listw = s2.bl.lw.10,
                            zero.policy = FALSE)

sem.s2.bl.doy.10 <- errorsarlm(ndvi.max.doy ~ snow.auc,
                            data = s2.bl,
                            listw = s2.bl.lw.10, 
                            zero.policy = FALSE)

s2.bl.nb.60 <- dnearneigh(s2.bl, d1 = 0, d2 = 60)
#s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s2.bl.lw.60 <- nb2listw(s2.bl.nb.60, style = 'W', zero.policy = FALSE)

sem.s2.bl.max.60 <- errorsarlm(ndvi.max ~ snow.auc,
                            data = s2.bl, 
                            listw = s2.bl.lw.60,
                            zero.policy = FALSE)

sem.s2.bl.doy.60 <- errorsarlm(ndvi.max.doy ~ snow.auc,
                            data = s2.bl,
                            listw = s2.bl.lw.60, 
                            zero.policy = FALSE)

s2.bl.nb.100 <- dnearneigh(s2.bl, d1 = 0, d2 = 100)
#s2.bl.nb <- poly2nb(st_buffer(s2.bl, dist = 5, endCapStyle = "SQUARE"), queen = TRUE)
# Redefine spatial weights for neighbourhoods
s2.bl.lw.100 <- nb2listw(s2.bl.nb.100, style = 'W', zero.policy = FALSE)

sem.s2.bl.max.100 <- errorsarlm(ndvi.max ~ snow.auc,
                            data = s2.bl, 
                            listw = s2.bl.lw.100,
                            zero.policy = FALSE)

sem.s2.bl.doy.100 <- errorsarlm(ndvi.max.doy ~ snow.auc,
                            data = s2.bl,
                            listw = s2.bl.lw.100, 
                            zero.policy = FALSE)

stargazer(sem.s2.bl.max.60, type = 'text')

stargazer(sem.s2.bl.max.10, sem.s2.bl.doy.10, 
          sem.s2.bl.max.60, sem.s2.bl.doy.60, type = 'html', 
          out = '../../data/statistical-output/sem-blaesedalen-nb-dist-10-60.html')#,
          
stargazer(sem.s2.bl.max.100, sem.s2.bl.doy.100, type = 'html',
          out = '../../data/statistical-output/sem-blaesedalen-nb-dist-100.html') #sem.s2.bl.doy.100, type = 'text')#, 
        


# Predictions ----
# Generate predicted values based on model

# Blaesedalen
bl.pred.10 <- data.frame(snow.auc = seq(min(s2.bl$snow.auc), 
                                     max(s2.bl$snow.auc), 
                                     0.5))

bl.pred.10 <- bl.pred %>% 
  mutate(ndvi = predict(sem.s2.bl.max.10, newdata = bl.pred))# %>%


bl.pred.60 <- data.frame(snow.auc = seq(min(s2.bl$snow.auc), 
                                        max(s2.bl$snow.auc), 
                                        0.5))

bl.pred.60 <- bl.pred %>% 
  mutate(ndvi = predict(sem.s2.bl.max.60, newdata = bl.pred))# %>% 

bl.pred.100 <- data.frame(snow.auc = seq(min(s2.bl$snow.auc), 
                                        max(s2.bl$snow.auc), 
                                        0.5))

bl.pred.100 <- bl.pred %>% 
  mutate(ndvi = predict(sem.s2.bl.max.100, newdata = bl.pred))# %>% 


# Plot ----
bl.coefficients.10 <- sem.s2.bl.max.10$coefficients
bl.equation.10 <- sprintf("10m, y = %.3f + %.3f * x", bl.coefficients.10[1], bl.coefficients.10[2])

bl.coefficients.60 <- sem.s2.bl.max.60$coefficients
bl.equation.60 <- sprintf("60m, y = %.3f + %.3f * x", bl.coefficients.60[1], bl.coefficients.60[2])

bl.coefficients.100 <- sem.s2.bl.max.100$coefficients
bl.equation.100 <- sprintf("100m, y = %.3f + %.3f * x", bl.coefficients.100[1], bl.coefficients.100[2])

bl <- ggplot() +
  geom_point(data = s2.bl, aes(x = snow.auc, y = ndvi.max), alpha = 0.5, color = '#4984BF') +
  geom_line(data = bl.pred.10, aes(x = snow.auc, y = ndvi), color = 'black', linewidth = 1) +
  geom_line(data = bl.pred.60, aes(x = snow.auc, y = ndvi), color = 'red', linewidth = 1) +
  geom_line(data = bl.pred.100, aes(x = snow.auc, y = ndvi), color = 'purple', linewidth = 1) +
  annotate("text", x = max(s2.bl$snow.auc), y = 0.4, label = bl.equation.10, hjust = 1, vjust = 1, color = 'black') +
  annotate("text", x = max(s2.bl$snow.auc), y = 0.3, label = bl.equation.60, hjust = 1, vjust = 1, color = 'red') +
  annotate("text", x = max(s2.bl$snow.auc), y = -0.1, label = bl.equation.100, hjust = 1, vjust = 1, color = 'purple') +
  xlab('') +
  ylab('') +
  theme_cowplot()
bl
kl.coefficients <- sem.s2.kl.max$coefficients
kl.equation <- sprintf("y = %.3f + %.3f * x", kl.coefficients[1], kl.coefficients[2])

kl <- ggplot() +
  geom_point(data = s2.kl, aes(x = snow.auc, y = ndvi.max), alpha = 0.5, color = '#F5A40C') +
  geom_line(data = kl.pred, aes(x = snow.auc, y = ndvi), color = 'black', linewidth = 1) +
  annotate("text", x = max(s2.kl$snow.auc), y = 0.6, label = kl.equation, hjust = 1, vjust = 1) + 
   xlab('') +
  ylab('') +
  theme_cowplot()

kh.coefficients <- sem.s2.kh.max$coefficients
kh.equation <- sprintf("y = %.3f + %.3f * x", kh.coefficients[1], kh.coefficients[2])

kh <- ggplot() +
  geom_point(data = s2.kh, aes(x = snow.auc, y = ndvi.max), alpha = 0.5, color = '#F23835') +
  geom_line(data = kh.pred, aes(x = snow.auc, y = ndvi), color = 'black', linewidth = 1) +
  annotate("text", x = max(s2.kh$snow.auc), y = 0.6, label = kh.equation, hjust = 1, vjust = 1) +  
  xlab('') +
  ylab('') +
  theme_cowplot()

kluane.plot <- plot_grid(kl, kh, 
                         ncol = 2, 
                         align = 'h')

combined.plots <- plot_grid(bl,
                            kluane.plot,
                            nrow = 2, 
                            align = 'h', 
                            labels = 'AUTO')
combined.plots

