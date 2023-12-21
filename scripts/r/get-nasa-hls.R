# Get NASA HLS data
# Calum Hoad, 18/12/2023

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

# Set aoi to one of the above polygons
aoi <- blaesedalen

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
  httr:content(as = 'text') %>%
  fromJSON()
cat('There are', search.request$numberMatched, 'features which match request.')


# Explore the data which has been returned ----
