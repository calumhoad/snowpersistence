# Script to process LandsatTS trends for manually cleaned data
# Calum Hoad, 20231130

# Load packages ----

# Load packages for data handling etc.
library(sf)
library(dplyr)
library(purrr)
library(data.table)
library(stringr)
library(rgee)
library(googledrive)
library(leaflet)
library(viridis)
library(ggplot2)
library(tidyverse)

# Load LandsatTS package
library(LandsatTS)

# Load altered functions from the LandsatTS package

# LandsatTS clean script, edited to include high reflectance pixels (> 1)
#   and with QA bit based cloud filtering commented out and replaced with
#   manual cloud filtering.
source('../../scripts/r/LandsatTS/ch_lsat_clean_data.R')
# LandsatTS format script, edited to keep the landsat.product.id variable, in
#   order to facilitate join based on product.id between lsat.dt and manual_screen
source('../../scripts/r/LandsatTS/ch_lsat_format_data.R')


# Set geometry of plots ----
kluane_low <- st_polygon(list(matrix(c(-138.4176238, 60.9693000,
                                       -138.4113274, 60.9694445, 
                                       -138.4107439, 60.9664827, 
                                       -138.4171064, 60.9663083, 
                                       -138.4176238, 60.9693000),
                                     ncol = 2, byrow = TRUE)))
kluane_high <- st_polygon(list(matrix(c(-138.4187700, 60.9644414, 
                                        -138.4118724, 60.9632543, 
                                        -138.4144941, 60.9605887, 
                                        -138.4210442, 60.9616106, 
                                        -138.4187700, 60.9644414), 
                                      ncol = 2, byrow = TRUE)))
blaesedalen <- st_polygon(list(matrix(c(-53.46895345, 69.30024099,
                                        -53.46135193, 69.29934175,
                                        -53.46382359, 69.29597989,
                                        -53.47081885, 69.29580664,
                                        -53.46895345, 69.30024099), 
                                      ncol = 2, byrow = TRUE)))


# Get pixel centres and extract data from GEE ----
all_regions_sf <- st_sfc(kluane_low, kluane_high, blaesedalen, crs = 4326) %>% st_sf() %>%
  mutate(region = c("kluane_low", "kluane_high", "blaesedalen"))


# Get pixel centres ----
# Split and map lsat_get_pixel_centers using dplyr and purrr wihout plotting
pixel_list <- all_regions_sf %>%
  split(.$region) %>%
  map(lsat_get_pixel_centers,
      pixel_prefix_from = "region") %>%
  bind_rows()

# Subset pixel list to return pixels per AOI
blaesedalen_only <- pixel_list[1:168, ]

kluane_high_only <- pixel_list[169:333, ]

kluane_low_only <- pixel_list[334:488, ]

# Extract data using GEE ----

# Export time-series using lsat_export_ts()
task_list <- lsat_export_ts(kluane_low_only, start_doy = 91, end_doy = 304)

# Monitor the progress of the task set for GEE
ee_monitoring()


# Read in data previously extracted through GEE ---- 
blaesedalen_path <- '../../data/landsat-ts/blaesedalen-expanded-ts.csv'
kluane_high_path <- '../../data/lsatTS-output/lsatTS_export_kluane_high.csv'
kluane_low_path <- '../../data/lsatTS-output/lsatTS_export_kluane_low.csv'

# Store as lsat.dt
lsat.dt <- do.call("rbind", lapply(blaesedalen_path, fread))


# Formatting and cleaning data ----

# Format data
lsat.dt <- ch_lsat_format_data(lsat.dt)

# Clean data with manual screening
lsat.manual.dt <- ch_lsat_clean_data(lsat.dt, 
                              geom.max = 50, 
                              cloud.max = 70, 
                              sza.max = 70, 
                              filter.cfmask.snow = F, 
                              filter.cfmask.water = F, 
                              filter.jrc.water = F)

# Clean data with auto screening (bit QA based)
lsat.auto.dt <- lsat_clean_data(lsat.dt,
                                geom.max = 50, 
                                cloud.max = 70,
                                sza.max = 70,                # 70 to get late season obs
                                filter.cfmask.snow = F, 
                                filter.cfmask.water = F, 
                                filter.jrc.water = F)


# Calculate NDVI and cross calibrate sensors ----

start.doy <- 91
end.doy <- 304

# Compute NDVI or other vegetation index
lsat.manual.dt <- lsat_calc_spectral_index(lsat.manual.dt, si = 'ndvi')
lsat.auto.dt <- lsat_calc_spectral_index(lsat.auto.dt, si = 'ndvi')

# Cross-calibrate NDVI among sensors using an approach based on Random Forest machine learning
lsat.manual.dt <- lsat_calibrate_rf(lsat.manual.dt, 
                              band.or.si = 'ndvi', 
                              doy.rng = start.doy:end.doy, 
                              min.obs = 5,               
                              frac.train = 0.75, 
                              overwrite.col = T, 
                              write.output = F,
                              train.with.highlat.data = T)

lsat.auto.dt <- lsat_calibrate_rf(lsat.auto.dt,
                                  band.or.si = 'ndvi',
                                  doy.rng = start.doy:end.doy, 
                                  min.obs = 5, 
                                  frac.train = 0.75, 
                                  overwrite.col = T, 
                                  write.output = F, 
                                  train.with.highlat.data = T)



# Fit phenological models (cubic splines) to each time series ----

# Fit phenological curves
lsat.manual.pheno.5 <- lsat_fit_phenological_curves(lsat.manual.dt, 
                                              si = 'ndvi', 
                                              window.yrs = 5, 
                                              window.min.obs = 10, 
                                              spl.fit.outfile = F, 
                                              progress = T, 
                                              si.min = 0)

lsat.manual.pheno.7 <- lsat_fit_phenological_curves(lsat.manual.dt, 
                                                  si = 'ndvi', 
                                                  window.yrs = 7, 
                                                  window.min.obs = 10, 
                                                  spl.fit.outfile = F, 
                                                  progress = T, 
                                                  si.min = 0)

lsat.auto.pheno.5 <- lsat_fit_phenological_curves(lsat.auto.dt, 
                                              si = 'ndvi', 
                                              window.yrs = 5, 
                                              window.min.obs = 10, 
                                              spl.fit.outfile = F, 
                                              progress = T, 
                                              si.min = 0)

lsat.auto.pheno.7 <- lsat_fit_phenological_curves(lsat.auto.dt, 
                                                si = 'ndvi', 
                                                window.yrs = 7, 
                                                window.min.obs = 10, 
                                                spl.fit.outfile = F, 
                                                progress = T, 
                                                si.min = 0)

# Derived annual growing season metrics
lsat.manual.gs.5 <- lsat_summarize_growing_seasons(lsat.manual.pheno.5, 
                                             si = 'ndvi', 
                                             min.frac.of.max = 0.75)

lsat.manual.gs.7 <- lsat_summarize_growing_seasons(lsat.manual.pheno.7, 
                                                   si = 'ndvi', 
                                                   min.frac.of.max = 0.75)

lsat.auto.gs.5 <- lsat_summarize_growing_seasons(lsat.auto.pheno.5, 
                                             si = 'ndvi', 
                                             min.frac.of.max = 0.75)

lsat.auto.gs.7 <- lsat_summarize_growing_seasons(lsat.auto.pheno.7, 
                                                 si = 'ndvi', 
                                                 min.frac.of.max = 0.75)


# Calculate trends ----

# Calculate trends
lsat.manual.trnds.5 <- lsat_calc_trend(lsat.manual.gs.5, si = 'ndvi.max', 2000:2020, sig = 0.1)
lsat.manual.trnds.7 <- lsat_calc_trend(lsat.manual.gs.7, si = 'ndvi.max', 2000:2020, sig = 0.1)
lsat.auto.trnds.5 <- lsat_calc_trend(lsat.auto.gs.5, si = 'ndvi.max', 2000:2020, sig = 0.1)
lsat.auto.trnds.7 <- lsat_calc_trend(lsat.auto.gs.7, si = 'ndvi.max', 2000:2020, sig = 0.1)

# Plot map of data ----

# Convert to sf
bl.manual.sf.5 <- st_as_sf(lsat.manual.trnds.5, coords = c("longitude","latitude"), crs = 32621)
bl.manual.sf.7 <- st_as_sf(lsat.manual.trnds.7, coords = c("longitude","latitude"), crs = 32621)
bl.auto.sf.5 <- st_as_sf(lsat.auto.trnds.5, coords = c("longitude", "latitude"), crs = 32621)
bl.auto.sf.7 <- st_as_sf(lsat.auto.trnds.7, coords = c("longitude", "latitude"), crs = 32621)

# Plot map
# Blaesedalen
bl.map <- leaflet() %>%
  addProviderTiles('OpenStreetMap.Mapnik') %>%
  addCircleMarkers(data = bl.manual.sf.5,
                   color = ~ifelse(trend.cat == 'no_trend', 'blue', 
                                   ifelse(trend.cat == 'browning', 'brown', 'green')),
                   group = 'Manual 5') %>%
  addCircleMarkers(data = bl.manual.sf.7,
                   color = ~ifelse(trend.cat == 'no_trend', 'blue', 
                                   ifelse(trend.cat == 'browning', 'brown', 'green')),
                   group = 'Manual 7') %>%
  addCircleMarkers(data = bl.auto.sf.5,
                 color = ~ifelse(trend.cat == 'no_trend', 'blue', 
                                 ifelse(trend.cat == 'browning', 'brown', 'green')),
                 group = 'Auto 5') %>%
  addCircleMarkers(data = bl.auto.sf.7,
                   color = ~ifelse(trend.cat == 'no_trend', 'blue', 
                                   ifelse(trend.cat == 'browning', 'brown', 'green')),
                   group = 'Auto 7') %>%
  addLayersControl(
    baseGroups = c("basemap"),
    overlayGroups = c("Manual 5", "Manual 7", "Auto 5", "Auto 7"),
    options = layersControlOptions(collapsed = FALSE)
  )


bl.map

# Output to csv
st_write(bl.manual.sf.5, '../../data/landsat-ts/trends/bl-manual-5-trends.csv', 
         layer_options = "GEOMETRY=AS_XY")
st_write(bl.manual.sf.7, '../../data/landsat-ts/trends/bl-manual-7-trends.csv', 
         layer_options = "GEOMETRY=AS_XY")
st_write(bl.auto.sf.5, '../../data/landsat-ts/trends/bl-auto-5-trends.csv', 
         layer_options = "GEOMETRY=AS_XY")
st_write(bl.auto.sf.7, '../../data/landsat-ts/trends/bl-auto-7-trends.csv', 
         layer_options = "GEOMETRY=AS_XY", crs = 32621)

# Get pixel centres only
lsatTS.pix.centres <- lsat.dt %>%
  group_by(sample.id) %>%
  filter(row_number() == 1) %>%
  st_as_sf(coords = c('latitude', 'longitude'))

# Output the pixel centres as a .shp
st_write(lsatTS.pix.centres, '../../data/lsatTS-output/blaesedalen/lsatTS_pix_centres_bl.shp')

lsat.auto.gs.7.sf <- st_as_sf(lsat.auto.gs.7, coords = c('longitude', 'latitude'), crs = 4326)

lsat.auto.gs.7.sf

bl.map <- leaflet() %>%
  addProviderTiles('OpenStreetMap.Mapnik') %>%
  addCircleMarkers(data = lsat.auto.gs.7.sf)
bl.map

st_write(lsat.auto.gs.7.sf, '../../data/lsatTS-output/blaesedalen/lsatTS_auto_7_gs.shp')
