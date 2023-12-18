# Get NASA HLS data
# Calum Hoad, 18/12/2023

# Install packages
install.packages('rgdal')

# Import packages
library(rgdal)
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

# Collection query
hls.col <- list("HLSS30.v2.0", "HLSL30.v2.0")

# Temporal query parameter
time <- '2023-06-01T00:00:00Z/2023-09-30T23:59:59Z'

