# NASA HLS 30m NDVI metrics
# Calum Hoad, 29 Jan 2024

# Extract NDVI time series from NASA HLS S30 data, over Kluane
# High, Kluane Low, and Blaesedalen plot locations.

# Import necessary libraries
library(terra)
library(dplyr)
library(rts)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(purrr)
library(broom)
library(viridis)
library(pbapply)
library(tidyterra)
library(cowplot)

###
# NASA HLS S30
###

# Blaesedalen ###

# Create a list of the S30 files for each S30 scene
d20230406 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-04-06/', full.names = TRUE)
d20230501 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-05-01/', full.names = TRUE)
d20230516 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-05-16/', full.names = TRUE)
d20230522 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-05-22/', full.names = TRUE)
d20230608 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-06-08/', full.names = TRUE)
d20230626 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-06-26/', full.names = TRUE)
d20230708 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-07-08/', full.names = TRUE)
d20230726 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-07-26/', full.names = TRUE)
d20230729 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-07-29/', full.names = TRUE)
d20230807 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-08-07/', full.names = TRUE)
d20230808 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-08-08/', full.names = TRUE)
d20230914 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-09-14/', full.names = TRUE)
d20230922 <- list.files('../../data/nasa-hls/blaesedalen/s30/X2023-09-22/', full.names = TRUE)
d20230923 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-09-23/', full.names = TRUE)
d20231003 <- list.files('../../data/nasa-hls/blaesedalen/s30/2023-10-03/', full.names = TRUE)

# List of imagery dates, for later use
bl.dates <- c('2023-04-06',
           '2023-05-01', 
           '2023-05-16',
           '2023-05-22', 
           '2023-06-08', 
           '2023-06-26', 
           '2023-07-08', 
           '2023-07-26', 
           '2023-07-29', 
           '2023-08-07', 
           '2023-08-08', 
           '2023-09-14',
           '2023-09-22',
           '2023-09-23', 
           '2023-10-03')

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s30.bl.data <- list(d20230406,
                 d20230501,
                 d20230516,
                 d20230522,
                 d20230608,
                 d20230626,
                 d20230708,
                 d20230726,
                 d20230729,
                 d20230807,
                 d20230808,
                 d20230914,
                 d20230922,
                 d20230923,
                 d20231003)

# Kluane ###

# Create a list of the S30 files for each S30 scene
d20220404 <- list.files('../../data/nasa-hls/kluane/s30/2022-04-04/', full.names = TRUE)
d20220512 <- list.files('../../data/nasa-hls/kluane/s30/2022-05-12/', full.names = TRUE)
d20220529 <- list.files('../../data/nasa-hls/kluane/s30/2022-05-29/', full.names = TRUE)
d20220608 <- list.files('../../data/nasa-hls/kluane/s30/2022-06-08/', full.names = TRUE)
d20220708 <- list.files('../../data/nasa-hls/kluane/s30/2022-07-08/', full.names = TRUE)
d20220721 <- list.files('../../data/nasa-hls/kluane/s30/2022-07-21/', full.names = TRUE)
d20220726 <- list.files('../../data/nasa-hls/kluane/s30/2022-07-26/', full.names = TRUE)
d20220812 <- list.files('../../data/nasa-hls/kluane/s30/2022-08-12/', full.names = TRUE)
d20220916 <- list.files('../../data/nasa-hls/kluane/s30/2022-09-16/', full.names = TRUE)
d20220924 <- list.files('../../data/nasa-hls/kluane/s30/2022-09-24/', full.names = TRUE)
d20220929 <- list.files('../../data/nasa-hls/kluane/s30/2022-09-29/', full.names = TRUE)
d20221006 <- list.files('../../data/nasa-hls/kluane/s30/2022-10-06/', full.names = TRUE)


# List of imagery dates, for later use
k.dates <- c('2022-04-04',
           '2022-05-12',
           '2022-05-29',
           '2022-06-08',
           '2022-07-08',
           '2022-07-21',
           '2022-07-26',
           '2022-08-12',
           '2022-09-16',
           '2022-09-24',
           '2022-09-29',
           '2022-10-06')

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s30.k.data <- list(d20220404,
                d20220512,
                d20220529,
                d20220608,
                d20220708,
                d20220721,
                d20220726,
                d20220812,
                d20220916,
                d20220924,
                d20220929,
                d20221006)


# Import UAV data for each site, to give area for NDVI time series extraction

# Blaesedalen
bl.uav <- rast('../../data/uav/orthomosaics/m3m/5cm/2023-07-02-5cm-clipped.tif')

# # Kluane, low, get uav imagery
kl.uav <- rast('../../data/uav/orthomosaics/maia/kluane-low/5cm/2022-07-05-5cm-clipped.tif')
# 
# # Kluane, high, get UAV imagery
kh.uav <- rast('../../data/uav/orthomosaics/maia/kluane-high/5cm/2022-07-09-5cm-clipped.tif')


# Create function for stacking rasters from lists of filepaths, then cropping to extent of UAV imagery
import_s30 <- function(x, plot) {
  x <- rast(x[1:13]) %>%
    crop(plot)
}

# Import the data
s30.bl.import <- pblapply(s30.bl.data, import_s30, plot = bl.uav) # Blaesedalen
s30.kl.import <- pblapply(s30.k.data, import_s30, plot = kl.uav) # Kluane-low
s30.kh.import <- pblapply(s30.k.data, import_s30, plot = kh.uav) # Kluane-high

# Check import function works by plotting rasters from list
ggplot() +
  geom_spatraster_rgb(data = s30.bl.import[[13]], r = 4, g = 3, b = 2, 
                      max_col_value = 0.2, 
                      interpolate = FALSE)

# Function to calculate NDVI (B8A = narrow NIR, B4 = red)
s30_ndvi <- function(x) {
  x <- (x$NIR_Narrow-x$Red)/(x$NIR_Narrow+x$Red)
}

# Apply function to calculate S2 NDVI
s30.bl.ndvi <- pblapply(s30.bl.import, s30_ndvi)
s30.kl.ndvi <- pblapply(s30.kl.import, s30_ndvi)
s30.kh.ndvi <- pblapply(s30.kh.import, s30_ndvi)

# Stack NDVI rasters into single spatRast
s30.bl.ndvi <- rast(s30.bl.ndvi)
s30.kl.ndvi <- rast(s30.kl.ndvi)
s30.kh.ndvi <- rast(s30.kh.ndvi)

# Name spatRaster layers with dates
names(s30.bl.ndvi) <- bl.dates
names(s30.kl.ndvi) <- k.dates
names(s30.kh.ndvi) <- k.dates

# Plot to check
ggplot() +
  geom_spatraster(data = s30.bl.ndvi) +
  scale_fill_viridis_c() +
  facet_wrap(~lyr)

# Extract raster time series to points
s30.bl.ndvi.points <- st_as_sf(as.points(s30.bl.ndvi, values = TRUE)) %>%
  mutate(id = row_number())

s30.kl.ndvi.points <- st_as_sf(as.points(s30.kl.ndvi, values = TRUE)) %>%
  mutate(id = row_number())

s30.kh.ndvi.points <- st_as_sf(as.points(s30.kh.ndvi, values = TRUE)) %>%
  mutate(id = row_number())


# Synthesise 0 values for days late in year, due to assymetry of the datasets ----

# Using the same dates as for the S2 data

# Blaesedalen ###
origin <- '2023-01-01'
d1 <- as.character(as_date(290, origin = origin)) # S30 data not available, but S2 data available, and 0 value predominant
d2 <- as.character(as_date(329, origin = origin)) # as per S2 synthetic imagery dates, and d3 and d4 the same
d3 <- as.character(as_date(331, origin = origin))
d4 <- as.character(as_date(344, origin = origin))

# Assign 0 NDVI value to the synthetic imagery dates generated above
s30.bl.ndvi.points <- s30.bl.ndvi.points %>%
  mutate(!!d1 := 0, 
         !!d2 := 0, 
         !!d3 := 0,
         !!d4 := 0)

# # Kluane ###
# origin <- '2022-01-01'
# d1 <- as.character(as_date(315, origin = origin)) 
# d2 <- as.character(as_date(324, origin = origin)) 
# d3 <- as.character(as_date(346, origin = origin))
# 
# # Assign 0 NDVI value to the synthetic imagery dates generated above
# s30.kl.ndvi.points <- s30.kl.ndvi.points %>%
#   mutate(!!d1 := 0,
#          !!d2 := 0, 
#          !!d3 := 0)
# s30.kh.ndvi.points <- s30.kl.ndvi.points %>%
#   mutate(!!d1 := 0,
#          !!d2 := 0, 
#          !!d3 := 0)

# Checking the data is logical ----
test <- s30.bl.ndvi.points %>% 
  pivot_longer(!id & !geometry, names_to = 'date', values_to = 'ndvi') %>%
  group_by(id)

# Examine the NDVI error in September and October 
whole.ts.plot <- ggplot() +
                 geom_vline(xintercept = 274, colour = 'red') +
                 geom_vline(xintercept = 244, colour = 'red') +
                 geom_line(data = test, aes(x = yday(date), y = ndvi, group = id), alpha = 0.1, linewidth = 0.5) +
                 xlab('Day of Year') +
                 ylab('NDVI (HLSS30)') +
                 theme_cowplot()
# Create another plot which focuses on the part of the year with apparent error
focus.ts.plot <- ggplot() +
                 geom_hline(yintercept = 1, colour = 'blue', linewidth = 1, alpha = 0.5, linetype = 'dashed') +
                 geom_vline(xintercept = 274, colour = 'red') +
                 geom_vline(xintercept = 244, colour = 'red') +
                 geom_line(data = test %>% filter(yday(date) > 210), 
                           aes(x = yday(date), y = ndvi, group = id), 
                           alpha = 0.1, linewidth = 0.5) + 
                 annotate("text", x = 275, y = 1.5, label = 'October 1st', hjust = 0, angle = 90) +
                 annotate("text", x = 245, y = 1.5, label = 'September 1st', hjust = 0, angle = 90) +
                 annotate("text", x = 245, y = 1.1,
                          label = "Standard maximum of\nthe NDVI = 1", hjust = 0) +
                 xlab('Day of Year') +
                 ylab('NDVI (HLSS30)') +
                 coord_cartesian(xlim = c(240, 280), ylim = c(0, 2.3)) +
                 theme_cowplot()

# Compare with NDVI derived from Sentinel-2
s2.ndvi <- read_csv('../../data/ndvi/s2-blaesedalen-ndvi-ts-pt.csv') %>%
  pivot_longer(!id & !X & !Y, names_to = 'date', values_to = 'ndvi') %>%
  group_by(id)

s2.ts.plot <- ggplot() +
  geom_vline(xintercept = 274, colour = 'red') +
  geom_vline(xintercept = 244, colour = 'red') +
  geom_line(data = s2.ndvi, aes(x = yday(date), y = ndvi, group = id), alpha = 0.1, linewidth = 0.5) +
  xlab('Day of Year') +
  ylab('NDVI (Sentinel-2)') +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75),
                     labels = c('0', '0.25', '0.50', '0.75')) +
  coord_cartesian(xlim = c(100, 280), ylim = c(0, 1)) +
  theme_cowplot()

# Combine the whole ts and the focus 
one <- cowplot::plot_grid(whole.ts.plot, s2.ts.plot, nrow = 2, align = 'v', labels = c('(a)', '(b)'))
two <- cowplot::plot_grid(one, focus.ts.plot, ncol = 2, labels = c('', '(c)'))
two

# Save out the plot
cowplot::save_plot('../../plots/figures/supp-fig-hlss30-ndvi-erratic.png', two, 
                   base_height = 180, base_width = 180, units = 'mm', 
                   bg = 'white')

# Outputs ---- 

# Wide format
st_write(st_as_sf(s30.bl.ndvi.points),  '../../data/ndvi/s30-blaesedalen-ndvi-ts-pt.csv', 
         layer_options = "GEOMETRY=AS_XY")
st_write(st_as_sf(s30.kl.ndvi.points), '../../data/ndvi/s30-kluanelow-ndvi-ts-pt.csv',
         layer_options = "GEOMETRY=AS_XY")
st_write(st_as_sf(s30.kh.ndvi.points), '../../data/ndvi/s30-kluanehigh-ndvi-ts-pt.csv',
         layer_options = "GEOMETRY=AS_XY")

# Long format
s30.modelled.export.long <- s30.bl.ndvi.points %>%
  pivot_longer(!id & !geometry, names_to = 'doy', values_to = 'ndvi')
rename(doy.obs = 'doy', 
       ndvi.obs = 'ndvi', 
       ndvi.pred = 'ndvi.pred.doy.1') %>%
  select(-ndvi.max.date)

st_write(st_as_sf(s30.modelled.export.long),  '../../data/nasa-hls/output/s30_kluane-high-modelled_smoothed_spline_point_long.csv', 
         layer_options = "GEOMETRY=AS_XY")



# Investigate what is going on with weird NDVI values? ----

# Get Fmask for all of the rasters
get_fmask <- function(x) {
  x <- crop(rast(x[14]), bl.uav)
}

fmask.ts <- pblapply(s30.bl.data, get_fmask)

names(fmask.ts) <- bl.dates

toBinary <- function(i){
  valuesbin = paste0(as.integer(intToBits(i)[1:8]), collapse = "")
  valuesbin = unlist(strsplit(valuesbin, split = ""))
  as.numeric(valuesbin) 
}

toB <- function(x) {
  t(sapply(x, toBinary))
}

fmask.ts <- pblapply(fmask.ts, terra::app, toB)

rename_bands <- function(x, bands) {
  names(x) <- bands
  update(x, names = TRUE)
}

bits <- c('Cirrus', 'Cloud', 'Adjacent C', 'Shadow', 'Snow/ice', 'Water', 'Aerosol1', 'Aerosol2')
fmask.ts <- pblapply(fmask.ts, rename_bands, bits)

# Plot as a large grid
ggplot() +
  geom_spatraster(data = fmask.ts[[1]]) +
  scale_fill_viridis_c() +
  facet_grid(~lyr)

# Iterate over the list of spatraster objects
plots <- list()
for (i in seq_along(fmask.ts)) {
  # Plot each spatraster and store the ggplot object in the list
  plot <- ggplot() +
    geom_spatraster(data = fmask.ts[[i]]) +
    scale_fill_viridis_c() +
    facet_grid(~lyr) +
    theme_void() +
    guides(fill = FALSE) +
    labs(title = '')
  
  plots[[i]] <- plot
}

plots[[2]]
bl.dates
# Combine all ggplot objects in the list into a single plot
combined_plot <- cowplot::plot_grid(plotlist = plots, nrow = 14, ncol = 1,
                                    labels = bl.dates,
                                    label_x = 0.34, label_y = 0.5)

combined_plot
# Save plots
cowplot::save_plot('../../plots/figures/figure-sup-hlss30-fmask.png', combined_plot, 
                   base_height = 180, base_width = 180, units = 'mm', 
                   bg = 'white')


# Print or save the combined plot
print(combined_plot)
  
rast(fmask.ts)

crs(s30.bl.import[[14]])

d2023
test <- crop(rast(c(d20230922[[4]], d20230922[[13]])), bl.uav)
test$NDVI <- (test$NIR_Narrow-test$Red)/(test$NIR_Narrow+test$Red)
ggplot() +
  geom_spatraster(data = test) +
  scale_color_viridis_c() +
  facet_wrap(~lyr)

# Try plotting histogram of values 
s30.bl.import[[1]]$NIR_Narrow

# Get only NIR_Narrow or Red
select_band <- function(x) {
  x <- x$NIR_Narrow
}

test <- pblapply(s30.bl.import, select_band)

test <- rast(test)

names(test) <- bl.dates

plot(test)
ggplot() +
  geom_bar(data = s30.bl.import[[1]]$NIR_Narrow)

ndvi <- pblapply(s30.bl.import, s30_ndvi)
ndvi <- rast(ndvi)
names(ndvi) <- bl.dates
plot(ndvi)

terra::hist(s30.bl.import[[1]])


# Load your Fmask raster
fmask_raster <- crop(rast('../../data/nasa-hls/blaesedalen/s30/2023-09-23/HLS.S30.T21WXT.2023266T151941.v2.0.Fmask.tif'), bl.uav)
fmask_raster
# Function to get FMask values
toBinary <- function(i){
  valuesbin = paste0(as.integer(rev(intToBits(i)[1:8])), collapse = "")
  valuesbin = unlist(strsplit(valuesbin, split = ""))
  as.numeric(valuesbin) 
}

toB <- function(x) {
  t(sapply(x, toBinary))
}

r <- fmask_raster
a <- terra::app(r, toB)
a

ggplot() +
  geom_spatraster(data = a) +
  scale_fill_viridis_c() +
  facet_grid(~lyr)
# Apply the mask to your Fmask raster
masked_fmask <- mask_quality_bits(fmask_raster)

# Apply the mask to your Fmask raster
masked_fmask <- mask_quality_bits(fmask_raster)

# Now 'masked_fmask' contains the Fmask values with pixels omitted where quality bits 6 and 7 are both set to 1

