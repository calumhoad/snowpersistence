# Script to pull and analyse data using the 
# Landsat TS package on GitHub
# Calum Hoad, 20230910

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

# Load LandsatTS package
library(LandsatTS)

# Load in the edited LandsatTS clean script (included reflectance values > 1)
source('../../scripts/r/LandsatTS/ch_lsat_clean_data.R')
# Load in the edited LandsatTS format data script (keeps product ID column)
source('../../scripts/r/LandsatTS/ch_lsat_format_data.R')

# Intialize the Earth Engine with rgee
ee_Initialize()

## Get pixel centers for Greenland and Yukon plots
# Create multi-polygon sf
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

test_regions_sf <- st_sfc(kluane_low, kluane_high, blaesedalen, crs = 4326) %>% st_sf() %>%
  mutate(region = c("kluane_low", "kluane_high", "blaesedalen"))

# Next step of script ----
# Split and map lsat_get_pixel_centers using dplyr and purrr wihout plotting
pixel_list <- test_regions_sf %>%
  split(.$region) %>%
  map(lsat_get_pixel_centers,
      pixel_prefix_from = "region") %>%
  bind_rows()

# Subset pixel list to return only the Blaesedalen rows, using indexing
blaesedalen_only <- pixel_list[1:168, ]

kluane_high_only <- pixel_list[169:333, ]

kluane_low_only <- pixel_list[334:488, ]

# Extract data using GEE ----

# Export time-series using lsat_export_ts()
task_list <- lsat_export_ts(kluane_low_only)

# Monitor the progress of the task set for GEE
ee_monitoring()

# Copy to R temp location
blaesedalen_path <- '../../data/lsatTS-output/lsatTS_export_blaesedalen.csv'
kluane_high_path <- '../../data/lsatTS-output/lsatTS_export_kluane_high.csv'
kluane_low_path <- '../../data/lsatTS-output/lsatTS_export_kluane_low.csv'

# Read in the files
lsat.dt <- do.call("rbind", lapply(blaesedalen_path, fread))

# setnames(lsat.dt, 'my_unique_location_column','sample.id') 
lsat.dt <- lsat_format_data(lsat.dt)

# Testing function
manual_screen <- read.csv("../../data/lsat-manual-screening/blaesedalen-screen.csv") #%>%
  #filter(quality_score == '1')


# Get IDs that have 1s next to them
ids <- unique(filter(manual_screen, quality_score == 1)$ID)

# Trimmed down masive thing
trimmed <- filter(lsat.dt, substring(landsat.product.id, 10, 40) %in% substring(ids, 10, 40))

str(lsat.dt)



attempt2 <- lsat.dt[, landsat.product.id %in% manual_screen$ID]
lsat.dt

attempt2

lsat.dt$landsat.product.id %in% manual_screen$ID

attempt3 <- left_join(lsat.dt, manual_screen, by = join_by(landsat.product.id == ID))

attempt <- filter(lsat.dt, landsat.product.id %in% id.list)

str(manual_screen)
str(lsat.dt)
'LT05_L2SP_010011_19850608_20200918_02_T1' %in% prod.list

id.list = as.list(manual_screen$ID)

prod.list = as.list(unique(lsat.dt$landsat.product.id))

tf.list = as.list(prod.list %in% id.list)

manual_screen$ID
lsat.dt$landsat.product.id
# Clean the surface reflectance data
lsat.dt <- ch_lsat_clean_data(lsat.dt, geom.max = 50, cloud.max = 80, sza.max = 60, filter.cfmask.snow = F, filter.cfmask.water = F, filter.jrc.water = F)

# Summarise the data availability
data.summary.dt <- lsat_summarize_data(lsat.dt)
data.summary.dt

# Compute NDVI or other vegetation index
lsat.dt <- lsat_calc_spectral_index(lsat.dt, si = 'ndvi')

# Cross-calibrate NDVI among sensors using an approach based on Random Forest machine learning
lsat.dt <- lsat_calibrate_rf(lsat.dt, 
                             band.or.si = 'ndvi', 
                             doy.rng = 151:242, 
                             min.obs = 5, 
                             frac.train = 0.75, 
                             overwrite.col = T, 
                             write.output = F,
                             train.with.highlat.data = T)

# If needed, then manually overwrite the uncalibrated data with the calibrated data, then drop the duplicate column 
# lsat.dt <- lsat.dt[, ndvi := ndvi.xcal]
# lsat.dt <- lsat.dt[, ndvi.xcal := NULL)

# Fit phenological models (cubic splines) to each time series
lsat.pheno.dt <- lsat_fit_phenological_curves(lsat.dt, si = 'ndvi', window.yrs = 5, window.min.obs = 10, spl.fit.outfile = F, progress = T) # removed vi.min = 0, unused argument error

# Derived annual growing season metrics
lsat.gs.dt <- lsat_summarize_growing_seasons(lsat.pheno.dt, si = 'ndvi', min.frac.of.max = 0.75)

# Evaluate raw and modeled
lsat.gs.eval.dt <- lsat_evaluate_phenological_max(lsat.pheno.dt, si = 'ndvi', min.obs = 10, reps = 5, min.frac.of.max = 0.75, outdir = NA)

# Calculate trends
lsat.trnds.dt <- lsat_calc_trend(lsat.gs.dt, si = 'ndvi.max', 2000:2020, sig = 0.1)

# Make the trend data into an sf
bl_trend_sf <- st_as_sf(lsat.trnds.dt, coords = c("longitude","latitude"))
kl_trend_sf <- st_as_sf(lsat.trnds.dt, coords = c("longitude","latitude"))
klLow_trend_Sf <- st_as_sf(lsat.trnds.dt, coords = c("longitude", "latitude"))
klLow_trend_Sf_geom <- st_as_sf(lsat.trnds.dt, coords = c("longitude", "latitude"))
# Plot as map
# Blaesedalen
bl_map <- leaflet() %>%
  addProviderTiles('OpenStreetMap.Mapnik') %>%
  addCircleMarkers(data = bl_trend_sf,
                   color = ~ifelse(trend.cat == 'no_trend', 'blue', 'brown'))

bl_map

#Kluane High
kh_map <- leaflet() %>%
  addProviderTiles('OpenStreetMap.Mapnik') %>%
  addCircleMarkers(data = kl_trend_sf,
                    color = ~ifelse(trend.cat == 'no_trend', 'blue', 'brown'))

kh_map

# Kluane Low
kl_map <- leaflet() %>%
  addProviderTiles('OpenStreetMap.Mapnik') %>%
  addCircleMarkers(data = klLow_trend_Sf_geom,
                   color = ~ifelse(trend.cat == 'no_trend', 'blue', 'green'))

kl_map

# Export csv
write.csv(lsat.trnds.dt, 'data/kluaneLow_trends_noSnoworWaterFilters.csv')
write.csv(lsat.gs.dt, 'data/kluaneLow_metrics_noSnoworWaterFilters.csv')
