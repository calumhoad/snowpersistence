# Script for calculating area-under-curve snow cover
# Calum Hoad, 09 Jan 2024

### Purpose of script
# There are two principle metrics for snow cover:
# 1) Snow cover duration (SCD), defined as the number of days for which a given
#   area remains above a given threshold for snow cover.
# 2) Snow cover extent (SCE), the percentage of a given area covered by snow at
#   a given point in time.
#
# For the analysis of snow's relationship with key NDVI metrics, we need a snow
# metric which encapsulates both the extent of snow within EO pixels and its
# evolution across time (i.e. a pixel with a 10 x 10 cm patch of snow which melts
# out completely on July 26th is not the same as a pixel with a 200 x 200 cm patch
# which melts out on the dame date)
#
# In order to calculate a more meaningful snow metric for this study, this script
# will calculate snow persistence as the area-under-the-curve of snow cover extent
# (y) across time (x).
###

# Import the necessary libraries ----
library(terra)
library(dplyr)
library(ggplot2)
library(buffeRs)
library(sf)
library(tidyverse)
library(pbapply)
