# Get NASA HLS data
# Calum Hoad, 18/12/2023

# The aim of this script is to download NASA Harmonized Landsat Sentinel (HLS)
# data, over both the Kluane and Blaesedalen focal study sites.The data will be 
# Sentinel-2 imagery, with spectral characteristics assimilated to those of 
# Landsat 8. Dates which are cloud free over the focal sites can be provided by
# the Sentinel-2 imagery for this project.

# Install packages
install.packages('rgdal')

# Import packages
library(rgdal) # Cannot install from CRAN, need to find alternative
library(leaflet)
library(jsonlite)
library(raster)
library(sp)
library(rasterVis)
library(httr)
library(ggplot2)
library(RColorBrewer)
library(dygraphs)
library(xts)
library(xml2)
library(lubridate)
library(DT)
library(sf)
library(geojsonsf)
library(tidyverse)

# Get all the Sentinel-2 tile IDs ----

# Blaesedalen
d20230406 <- list.files('../../data/sentinel-2/imagery/20230406/S2A_MSIL2A_20230406T151941_N0509_R068_T21WXS_20230406T214452.SAFE/GRANULE/L2A_T21WXS_A040677_20230406T152237/IMG_DATA/R10m/', full.names = TRUE)
d20230501 <- list.files('../../data/sentinel-2/imagery/20230501/S2B_MSIL2A_20230501T151809_N0509_R068_T21WXS_20230501T173400.SAFE/GRANULE/L2A_T21WXS_A032126_20230501T152117/IMG_DATA/R10m/', full.names = TRUE)
d20230516 <- list.files('../../data/sentinel-2/imagery/20230516/S2A_MSIL2A_20230516T151801_N0509_R068_T21WXS_20230516T215553.SAFE/GRANULE/L2A_T21WXS_A041249_20230516T152031/IMG_DATA/R10m/', full.names = TRUE)
d20230522 <- list.files('../../data/sentinel-2/imagery/20230522/S2A_MSIL2A_20230522T153811_N0509_R011_T21WXS_20230522T231356.SAFE/GRANULE/L2A_T21WXS_A041335_20230522T154223/IMG_DATA/R10m/', full.names = TRUE)
d20230525 <- list.files('../../data/sentinel-2/imagery/20230525/S2A_MSIL2A_20230525T154941_N0509_R054_T22WDB_20230525T233156.SAFE/GRANULE/L2A_T22WDB_A041378_20230525T155122/IMG_DATA/R10m/', full.names = TRUE)
d20230608 <- list.files('../../data/sentinel-2/imagery/20230608/S2A_MSIL2A_20230608T152811_N0509_R111_T21WXS_20230608T215653.SAFE/GRANULE/L2A_T21WXS_A041578_20230608T152919/IMG_DATA/R10m/', full.names = TRUE)
d20230625 <- list.files('../../data/sentinel-2/imagery/20230626/S2B_MSIL2A_20230626T153819_N0509_R011_T21WXT_20230626T183843.SAFE/GRANULE/L2A_T21WXT_A032927_20230626T153935/IMG_DATA/R10m/', full.names = TRUE)
d20230626 <- list.files('../../data/sentinel-2/imagery/20230626/S2B_MSIL2A_20230626T153819_N0509_R011_T21WXT_20230626T183843.SAFE/GRANULE/L2A_T21WXT_A032927_20230626T153935/IMG_DATA/R10m/', full.names = TRUE)
d20230708 <- list.files('../../data/sentinel-2/imagery/20230708/S2A_MSIL2A_20230708T152811_N0509_R111_T21WXT_20230708T214952.SAFE/GRANULE/L2A_T21WXT_A042007_20230708T152945/IMG_DATA/R10m/', full.names = TRUE)
d20230726 <- list.files('../../data/sentinel-2/imagery/20230726/S2B_MSIL2A_20230726T153819_N0509_R011_T21WXT_20230726T184037.SAFE/GRANULE/L2A_T21WXT_A033356_20230726T153828/IMG_DATA/R10m/', full.names = TRUE)
d20230729 <- list.files('../../data/sentinel-2/imagery/20230729/S2B_MSIL2A_20230729T154819_N0509_R054_T21WXT_20230729T181552.SAFE/GRANULE/L2A_T21WXT_A033399_20230729T154940/IMG_DATA/R10m/', full.names = TRUE)
d20230807 <- list.files('../../data/sentinel-2/imagery/20230807/S2A_MSIL2A_20230807T152811_N0509_R111_T21WXS_20230807T212701.SAFE/GRANULE/L2A_T21WXS_A042436_20230807T153242/IMG_DATA/R10m/', full.names = TRUE)
d20230808 <- list.files('../../data/sentinel-2/imagery/20230808/S2B_MSIL2A_20230808T154819_N0509_R054_T21WXS_20230914T102422.SAFE/GRANULE/L2A_T21WXS_A033542_20230808T154852/IMG_DATA/R10m/', full.names = TRUE)
d20230817 <- list.files('../../data/sentinel-2/imagery/20230817/S2A_MSIL2A_20230817T152941_N0509_R111_T21WXT_20230817T214159.SAFE/GRANULE/L2A_T21WXT_A042579_20230817T153311/IMG_DATA/R10m/', full.names = TRUE)
d20230923 <- list.files('../../data/sentinel-2/imagery/20230923/S2A_MSIL2A_20230922T154951_N0509_R054_T21WXT_20230922T234400.SAFE/GRANULE/L2A_T21WXT_A043094_20230922T155005/IMG_DATA/R10m/', full.names = TRUE)
d20231003 <- list.files('../../data/sentinel-2/imagery/20231003/S2A_MSIL2A_20231003T152051_N0509_R068_T21WXT_20231003T202504.SAFE/GRANULE/L2A_T21WXT_A043251_20231003T152052/IMG_DATA/R10m/', full.names = TRUE)
d20231011 <- list.files('../../data/sentinel-2/imagery/20231011/S2B_MSIL2A_20231011T153059_N0509_R111_T21WXS_20231011T190635.SAFE/GRANULE/L2A_T21WXS_A034457_20231011T153057/IMG_DATA/R10m/', full.names = TRUE)
d20231017 <- list.files('../../data/sentinel-2/imagery/20231017/S2B_MSIL2A_20231017T155149_N0509_R054_T21WXT_20231017T201439.SAFE/GRANULE/L2A_T21WXT_A034543_20231017T155218/IMG_DATA/R10m/', full.names = TRUE)




#Source the earthdata set-up script
source('nasa-hls/earthdata_netrc_setup.R')

# Assign the LP2CLOUD STAC search uRL
search.url <- 'https://cmr.earthdata.nasa.gov/stac/LPCLOUD/search'

# Collection query
hls.col <- list("HLSS30.v2.0", "HLSL30.v2.0")

# Set aoi for each study area
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

# Get sf with all areas 
all_regions_sf <- st_sfc(kluane_low, kluane_high, blaesedalen, crs = 4326) %>% st_sf() %>%
  mutate(region = c("kluane_low", "kluane_high", "blaesedalen"))

# Set aoi to one of the above polygons
aoi <- all_regions_sf %>% filter(region == 'blaesedalen')

# To JSON
aoi <- geojsonsf::sf_geojson(aoi)

# Temporal query parameter
time <- '2023-06-01T00:00:00Z/2023-09-30T23:59:59Z'

# Submit a search query
search.query <- list(limit = 100, 
                     datetime = time, 
                     bbox = aoi, 
                     collections = hls.col)

# Create request
search.request <- httr::POST(search.url, 
                             body = search.query, 
                             encode = 'json') %>%
                  httr::content(as = 'text')

cat('There are', search.request$numberMatched, 'features which match request.')

search.request

# Explore the data which has been returned ----
