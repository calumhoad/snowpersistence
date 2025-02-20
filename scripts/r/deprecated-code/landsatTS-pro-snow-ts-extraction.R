# Script which will alter LandsatTS parameters in order to avoid the removal o
# snow covered pixels

# Calum Hoad, 20231115

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
# LansatTS calibrate rf script, edited to try and find which objects are missing,
#   such as the fig.list object, after alteration to the other functions above. 
#source('../../scripts/r/LandsatTS/ch_lsat_calibrate_rf.R')

# Intialize the Earth Engine with rgee
ee_Initialize()


## Get pixel centers for Greenland and Yukon plots ----
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

all_regions_sf <- st_sfc(kluane_low, kluane_high, blaesedalen, crs = 4326) %>% st_sf() %>%
  mutate(region = c("kluanelow", "kluanehigh", "blaesedalen"))


# Get pixel centres ----
# Split and map lsat_get_pixel_centers using dplyr and purrr wihout plotting
pixel_list <- all_regions_sf %>%
  split(.$region) %>%
  map(lsat_get_pixel_centers,
      pixel_prefix_from = "region") %>%
  bind_rows()

write.csv2(pixel_list, '../../data/lsatTS-output/pixel_centres.csv')

st_write(pixel_list, '../../data/lsatTS-output/pixel_centres.shp')

# Subset pixel list to return pixels per AOI
blaesedalen_only <- pixel_list[1:168, ]

kluane_high_only <- pixel_list[169:333, ]

kluane_low_only <- pixel_list[334:488, ]


# Extract data using GEE ----

# Export time-series using lsat_export_ts()
task_list <- lsat_export_ts(kluane_low_only)

# Monitor the progress of the task set for GEE
ee_monitoring()

# Read in data previosly processed through the above steps ----
blaesedalen_path <- '../../data/lsatTS-output/lsatTS_export_blaesedalen.csv'
kluane_high_path <- '../../data/lsatTS-output/lsatTS_export_kluane_high.csv'
kluane_low_path <- '../../data/lsatTS-output/lsatTS_export_kluane_low.csv'

# Read in the files
lsat.dt <- do.call("rbind", lapply(blaesedalen_path, fread))

# Get a list of all the sample ids ----
reduced.lsat.dt <- distinct(lsat.dt, sample_id, .keep_all = TRUE) # Reduce to unique sample.id
reduced.lsat.dt$sample_id

# Check the data obs per sample id direct from gee ----
columns <- c("sample.id", "gee.obs")
gee.obs.per.sample <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(gee.obs.per.sample) <- columns
gee.obs.per.sample$sample.id <- as.character(gee.obs.per.sample$sample.id)
gee.obs.per.sample$gee.obs <- as.integer(gee.obs.per.sample$gee.obs)

gee.obs.per.sample <- gee.obs.per.sample %>% add_row(sample.id = 'test', gee.obs = 0)

# Iterate through sample.id and calculate the number of obs, output to df
for(sample in reduced.lsat.dt$sample_id) {
  lsat.dt.sample <- filter(lsat.dt, lsat.dt$sample_id == sample)
  print(sample)
  obs <- length(lsat.dt.sample$sample_id)
  gee.obs.per.sample <- gee.obs.per.sample %>% add_row(sample.id = sample, gee.obs = obs)
  print(gee.obs.per.sample)
}  

range.gee <- max(gee.obs.per.sample$gee.obs[2:169]) - min(gee.obs.per.sample$gee.obs[2:169])
range.gee
# Plot obs per sample.id as barplot
gee.obs.hist <- ggplot(gee.obs.per.sample, aes(x = sample.id, y = gee.obs, color = sample.id)) +
  geom_bar(width = 0.2, stat = "identity", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set1")

gee.obs.hist
# Step through the LandsatTS package functions, edited to keep snow ----

# FORMATING AND CLEANING 

# Format data
lsat.dt <- ch_lsat_format_data(lsat.dt)

# Summarise the number of obs
data.summary <- lsat_summarize_data(lsat.dt)

# How many observations are there in each year? ----
lsat.dt.years <- filter(lsat.dt, lsat.dt$sample.id == 'blaesedalen_1')
lsat.dt.years <- lsat.dt.years %>% group_by(year) %>% summarise(Count = n ()) %>% dplyr::ungroup()
ggplot(lsat.dt.years, aes(x = year, y = Count, colour = year)) +
  geom_bar(stat = 'Identity', show.legend = FALSE) +
  labs(title = 'sample.id = 1') +
  theme_minimal()

# Plot the number of obs per pixel
reduced.lsat.dt <- distinct(lsat.dt, sample.id, .keep_all = TRUE) # Reduce to unique sample.id
reduced.lsat.dt$sample.id

# Create empty dataframe to store results
columns <- c("sample.id", "initial.obs")
obs.per.sample <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(obs.per.sample) <- columns
obs.per.sample$sample.id <- as.character(obs.per.sample$sample.id)
obs.per.sample$initial.obs <- as.integer(obs.per.sample$initial.obs)

obs.per.sample <- obs.per.sample %>% add_row(sample.id = 'test', initial.obs = 0)

# Iterate through sample.id and calculate the number of obs, output to df
for(sample in reduced.lsat.dt$sample.id) {
  lsat.dt.sample <- filter(lsat.dt, lsat.dt$sample.id == sample)
  print(sample)
  obs <- length(lsat.dt.sample$sample.id)
  obs.per.sample <- obs.per.sample %>% add_row(sample.id = sample, initial.obs = obs)
  print(obs.per.sample)
  }  

# Plot obs per sample.id as barplot
obs.hist <- ggplot(obs.per.sample, aes(x = sample.id, y = initial.obs, color = sample.id)) +
    geom_bar(width = 0.2, stat = "identity", show.legend = FALSE) +
    scale_fill_brewer(palette = "Set1")

obs.hist

range.format <- max(obs.per.sample$initial.obs[2:169]) - min(obs.per.sample$initial.obs[2:169])
range.format

# Clean the surface reflectance data ----
lsat.dt <- ch_lsat_clean_data(lsat.dt, 
                              geom.max = 50, 
                              cloud.max = 70, 
                              sza.max = 60, 
                              filter.cfmask.snow = F, 
                              filter.cfmask.water = F, 
                              filter.jrc.water = F)

# Could be useful to plot the range of sza, geom error, cloud
#   prior to the cleaning script? Should inform how many obs are being 
#   dropped?

# Check how many obs per pixel.id after cleaning ----
columns <- c("sample.id", "clean.obs")
clobs.per.sample <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(clobs.per.sample) <- columns
clobs.per.sample$sample.id <- as.character(clobs.per.sample$sample.id)
clobs.per.sample$clean.obs <- as.integer(clobs.per.sample$clean.obs)

clobs.per.sample <- clobs.per.sample %>% add_row(sample.id = 'test', clean.obs = 0)

# Iterate through sample.id and calculate the number of obs, output to df
  for(sample in reduced.lsat.dt$sample.id) {
  lsat.dt.sample <- filter(lsat.dt, lsat.dt$sample.id == sample)
  print(sample)
  obs <- length(lsat.dt.sample$sample.id)
  clobs.per.sample <- clobs.per.sample %>% add_row(sample.id = sample, clean.obs = obs)
  print(clobs.per.sample)
}  

# Plot obs per sample.id as barplot
clobs.hist <- ggplot(clobs.per.sample, aes(x = sample.id, y = clean.obs, color = sample.id)) +
  geom_bar(width = 0.2, stat = "identity", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set1")

clobs.hist

range.clobs <- max(clobs.per.sample$clean.obs[2:169]) - min(clobs.per.sample$clean.obs[2:169])
range.clobs
# For Blaesedalen, edited script "removed 57921 of 77810 observations (74.44%)"
#                original script"removed 52913 of 77810 observations (68%)"
#
# The script with manual filtering is around 6.44% stricter (less data) than 
#   the sript which uses QA bit based filtering, in the case of Blaesedalen.


# Summarise the data availability
data.summary.dt <- lsat_summarize_data(lsat.dt)
data.summary.dt

# Plot a map to look at numbers of data obs per sample.id
reduced.lsat.dt <- distinct(lsat.dt, sample.id, .keep_all = TRUE) # Reduce to unique sample.id
data.summary.dt <- inner_join(data.summary.dt, reduced.lsat.dt %>% # Join to get lat and lon
                                select(sample.id, latitude, longitude), by = "sample.id")

# Data to sf
bl_data_sf <- st_as_sf(data.summary.dt, coords = c("longitude","latitude"))

# Create colour palette
color_palette <- colorNumeric(palette = viridis(10), domain = bl_data_sf$n.obs.tot)
# Plot map
bl_data_map <- leaflet() %>%
  addProviderTiles('OpenStreetMap.Mapnik') %>%
  addCircleMarkers(data = bl_data_sf,
                   color = ~color_palette(bl_data_sf$n.obs.tot))
bl_data_map


# NDVI, CROSS-CALIBRATION, TREND ANALYSIS ----

# Compute NDVI or other vegetation index
lsat.dt <- lsat_calc_spectral_index(lsat.dt, si = 'ndvi')

# Cross-calibrate NDVI among sensors using an approach based on Random Forest machine learning
lsat.dt <- LandsatTS::lsat_calibrate_rf(lsat.dt, 
                             band.or.si = 'ndvi', 
                             doy.rng = 151:242, 
                             min.obs = 5,               # won't calc for 1 year
                             frac.train = 0.75, 
                             overwrite.col = T, 
                             write.output = F,
                             train.with.highlat.data = T)

# Check number of obs again ----
columns <- c("sample.id", "cal.obs")
calobs.per.sample <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(calobs.per.sample) <- columns
calobs.per.sample$sample.id <- as.character(calobs.per.sample$sample.id)
calobs.per.sample$cal.obs <- as.integer(calobs.per.sample$cal.obs)

calobs.per.sample <- calobs.per.sample %>% add_row(sample.id = 'test', cal.obs = 0)

# Iterate through sample.id and calculate the number of obs, output to df
for(sample in reduced.lsat.dt$sample.id) {
  lsat.dt.sample <- filter(lsat.dt, lsat.dt$sample.id == sample)
  print(sample)
  obs <- length(lsat.dt.sample$sample.id)
  calobs.per.sample <- calobs.per.sample %>% add_row(sample.id = sample, cal.obs = obs)
  print(calobs.per.sample)
}  

# Plot obs per sample.id as barplot
calobs.hist <- ggplot(calobs.per.sample, aes(x = sample.id, y = cal.obs, color = sample.id)) +
  geom_bar(width = 0.2, stat = "identity", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set1")

calobs.hist


# Explore functioning of lsat_fit_phenological_curves() ---- 

# Get a selection of random sample.id numbers
set.seed(42)  # Setting a seed for reproducibility
sample.ids <- sample(0:100, 10, replace = FALSE)
sample.ids <- paste('blaesedalen_', sample.ids, sep = '')

# Filter to only where sample.id matches sample.ids list
sample.lsat.dt <- filter(lsat.dt, sample.id %in% sample.ids)
sample.dt <- lsat.dt

# Check the filter worked
sample.id.check <- distinct(sample.lsat.dt, sample.id, .keep_all = TRUE) # Reduce to unique sample.id
sample.id.check$sample.id

# PLot the data
ggplot(data = sample.lsat.dt, aes(x = doy, y = ndvi, colour = year)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~I(x^2)) +
  facet_wrap(~sample.id) +
  scale_color_viridis()
  theme_classic()





  # Fit phenological models (cubic splines) to each time series ----
lsat.pheno.dt <- lsat_fit_phenological_curves(lsat.dt, 
                                              si = 'ndvi', 
                                              window.yrs = 5, # Previously 5, changed to bring back pixels 
                                              window.min.obs = 10, 
                                              spl.fit.outfile = F, 
                                              progress = T, 
                                              si.min = 0.15) # removed vi.min = 0, unused argument error
  
lsat.pheno.5.10 <- lsat.pheno.dt

lsat.pheno.7.10 <- lsat_fit_phenological_curves(lsat.dt,
                                                si = "ndvi", 
                                                window.yrs = 7,
                                                window.min.obs = 10, 
                                                spl.fit.outfile = F, 
                                                progress = T, 
                                                si.min = 0.15)

# Check number of obs again ----
columns <- c("sample.id", "pheno.obs")
pheno.per.sample <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(pheno.per.sample) <- columns
pheno.per.sample$sample.id <- as.character(pheno.per.sample$sample.id)
pheno.per.sample$pheno.obs <- as.integer(pheno.per.sample$pheno.obs)

pheno.per.sample <- pheno.per.sample %>% add_row(sample.id = 'test', pheno.obs = 0)

# Iterate through sample.id and calculate the number of obs, output to df
for(sample in reduced.lsat.dt$sample.id) {
  lsat.dt.sample <- filter(lsat.pheno.dt, lsat.pheno.dt$sample.id == sample)
  print(sample)
  obs <- length(lsat.dt.sample$sample.id)
  pheno.per.sample <- pheno.per.sample %>% add_row(sample.id = sample, pheno.obs = obs)
  print(pheno.per.sample)
}  

# Plot obs per sample.id as barplot
pheno.hist <- ggplot(pheno.per.sample, aes(x = sample.id, y = pheno.obs, color = sample.id)) +
  geom_bar(width = 0.2, stat = "identity", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set1")

pheno.hist 


# Derived annual growing season metrics ----
lsat.gs.dt <- lsat_summarize_growing_seasons(lsat.pheno.dt, 
                                             si = 'ndvi', 
                                             min.frac.of.max = 0.75)

lsat.gs.5.10 <- lsat_summarize_growing_seasons(lsat.pheno.5.10, 
                                               si = 'ndvi', 
                                               min.frac.of.max = 0.75)

lsat.gs.7.10 <- lsat_summarize_growing_seasons(lsat.pheno.7.10, 
                                               si = 'ndvi', 
                                               min.frac.of.max = 0.75)
# Check number of obs again ----
columns <- c("sample.id", "met.obs")
met.per.sample <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(met.per.sample) <- columns
met.per.sample$sample.id <- as.character(met.per.sample$sample.id)
met.per.sample$met.obs <- as.integer(met.per.sample$met.obs)

met.per.sample <- met.per.sample %>% add_row(sample.id = 'test', met.obs = 0)

# Iterate through sample.id and calculate the number of obs, output to df
for(sample in reduced.lsat.dt$sample.id) {
  lsat.dt.sample <- filter(lsat.gs.dt, lsat.gs.dt$sample.id == sample)
  print(sample)
  obs <- length(lsat.dt.sample$sample.id)
  met.per.sample <- met.per.sample %>% add_row(sample.id = sample, met.obs = obs)
  print(met.per.sample)
}  

# Plot obs per sample.id as barplot
met.hist <- ggplot(met.per.sample, aes(x = sample.id, y = met.obs, color = sample.id)) +
  geom_bar(width = 0.2, stat = "identity", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set1")

met.hist 



# Evaluate raw and modeled ----
lsat.gs.eval.dt <- lsat_evaluate_phenological_max(lsat.pheno.dt, 
                                                  si = 'ndvi', 
                                                  min.obs = 10, 
                                                  reps = 5, 
                                                  min.frac.of.max = 0.75, 
                                                  outdir = NA)

lsat.gs.eval.dt.5.10 <- lsat_evaluate_phenological_max(lsat.pheno.5.10, 
                                                  si = 'ndvi', 
                                                  min.obs = 10, 
                                                  reps = 5, 
                                                  min.frac.of.max = 0.75, 
                                                  outdir = NA)

lsat.gs.eval.dt.7.10 <- lsat_evaluate_phenological_max(lsat.pheno.7.10, 
                                                  si = 'ndvi', 
                                                  min.obs = 10, 
                                                  reps = 5, 
                                                  min.frac.of.max = 0.75, 
                                                  outdir = NA)


# Calculate trends
lsat.trnds.dt <- lsat_calc_trend(lsat.gs.dt, si = 'ndvi.max', 2000:2020, sig = 0.1)

lsat.trnds.5.10 <- lsat_calc_trend(lsat.gs.5.10, si = "ndvi.max", 2000:2020, sig = 0.1 )

lsat.trnds.7.10 <- lsat_calc_trend(lsat.gs.7.10, si = 'ndvi.max', 2000:2020, sig = 0.1)


# Check number of obs again ----
columns <- c("sample.id", "trend.obs")
trend.per.sample <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(trend.per.sample) <- columns
trend.per.sample$sample.id <- as.character(trend.per.sample$sample.id)
trend.per.sample$trend.obs <- as.integer(trend.per.sample$trend.obs)

trend.per.sample <- trend.per.sample %>% add_row(sample.id = 'test', trend.obs = 0)

# Iterate through sample.id and calculate the number of obs, output to df
for(sample in reduced.lsat.dt$sample.id) {
  lsat.dt.sample <- filter(lsat.trnds.dt, lsat.trnds.dt$sample.id == sample)
  print(sample)
  obs <- length(lsat.dt.sample$sample.id)
  trend.per.sample <- trend.per.sample %>% add_row(sample.id = sample, trend.obs = obs)
  print(trend.per.sample)
}  

# Plot obs per sample.id as barplot
trend.hist <- ggplot(trend.per.sample, aes(x = sample.id, y = trend.obs, color = sample.id)) +
  geom_bar(width = 0.2, stat = "identity", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set1")

trend.hist 

# Compile dataframe of obs at each point in processing, check for drop-off ----
obs.df <- list(obs.per.sample, clobs.per.sample, calobs.per.sample, pheno.per.sample,
             met.per.sample, trend.per.sample)
obs.df <- obs.df %>% reduce(inner_join, by='sample.id')

obs.df <- obs.df %>%
          mutate(trend = ifelse(trend.obs == 0, FALSE, TRUE)) %>%
          pivot_longer(!sample.id & !trend, names_to = "stage", values_to = "obs") %>%
          mutate(stage.num = ifelse(stage == 'initial.obs', 1,
                                    ifelse(stage == 'clean.obs', 2, 
                                           ifelse(stage == 'cal.obs', 3, 
                                                  ifelse(stage == 'pheno.obs', 4,
                                                         ifelse(stage == 'met.obs', 5,
                                                                ifelse(stage == 'trend.obs', 6, 7)))))))

# Filter data to only given stages, for subsequent plotting
obs.df <- obs.df %>% filter(stage.num != 6)

obs.df.trend <- obs.df %>% filter(trend == TRUE)
obs.df.notrend <- obs.df %>% filter(trend == FALSE)

obs.line <- ggplot(obs.df, aes(x = stage.num, y = obs, color = sample.id, group_by(trend))) +
            geom_line(show.legend = FALSE) +
            theme(axis.text.x = element_blank()) +
            geom_text(label="lsat_format_data",
                       x=1,
                       y=50,
                       color = "black",
                       angle = 90) +
            geom_text(label="lsat_clean_data",
                       x=2,
                       y=50,
                       color = "black",
                       angle = 90) +
            geom_text(label="lsat_calibrate_rf",
                      x=3,
                      y=50,
                      color = "black",
                      angle = 90) +  
            geom_text(label="lsat_fit_phenological_curves",
                      x=4,
                      y=400,
                      color = "black",
                      angle = 90) +
            geom_text(label="lsat_summarize_growing_seasons",
                      x=5,
                      y=400,
                      color = "black",
                      angle = 90) +
            ylab("Observations per sample.id") +
            xlab("LandsatTS processing stage (functions)")

# Call plot            
obs.line

# Check basic stats on the range of obs per sample.id ----
range.init <- max(obs.per.sample$initial.obs[2:169]) - min(obs.per.sample$initial.obs[2:169])
range.init # There is only a range of 6 in the initial (format) stage
range.clean <- max(clobs.per.sample$clean.obs[2:169]) - min(clobs.per.sample$clean.obs[2:169])
range.clean # Range of 15
range.cal <- max(calobs.per.sample$cal.obs[2:169]) - min(calobs.per.sample$cal.obs[2:169])
range.cal # Range of 15
range.pheno <- max(pheno.per.sample$pheno.obs[2:169]) - min(pheno.per.sample$pheno.obs[2:169])
range.pheno # Range of 166
range.met <- max(met.per.sample$met.obs[2:169]) - min(met.per.sample$met.obs[2:169])
range.met # Range of 12

# Explore trends and output from LandsatTS, edited to include snow ----

# Make the trend data into an sf
bl_trend_sf <- st_as_sf(lsat.trnds.dt, coords = c("longitude","latitude"))
bl.5.10 <- st_as_sf(lsat.trnds.5.10, coords = c("longitude", "latitude"))
bl.7.10 <- st_as_sf(lsat.trnds.7.10, coords = c("longitude", "latitude"))
#kl_trend_sf <- st_as_sf(lsat.trnds.dt, coords = c("longitude","latitude"))
#klLow_trend_Sf <- st_as_sf(lsat.trnds.dt, coords = c("longitude", "latitude"))
#klLow_trend_Sf_geom <- st_as_sf(lsat.trnds.dt, coords = c("longitude", "latitude"))


# Plotting data at various processing steps as map
#   to investigate at which step data is going missing
#   choose from:
#     - lsat.pheno.dt (Fitted phenological curves)
#     - lsat.gs.dt    (Derived seasonal metrics)

bl_datastep_sf <- st_as_sf(lsat.pheno.dt, coords = c("longitude", "latitude"))
bl_datastep_map <- leaflet() %>%
  addProviderTiles('OpenStreetMap.Mapnik') %>%
  addCircleMarkers(data = bl_datastep_sf)
bl_datastep_map


# Plot trends as maps

# Blaesedalen
bl_map <- leaflet() %>%
  addProviderTiles('OpenStreetMap.Mapnik') %>%
  addCircleMarkers(data = bl.5.10,
                   color = ~ifelse(trend.cat == 'no_trend', 'blue', 
                                   ifelse(trend.cat == 'browning', 'brown', 'green')))

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

# Questions to be answered:
#   1) There should be the same number of obs in all pixels, given that
#       filtering of the data is now manual and is for the whole scene. What am
#       I missing that is causing some pixels to have less obs than others?

