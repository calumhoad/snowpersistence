# Script to create raster stacks of UAV orthomosaics in Terra
# Calum Hoad, 16/11/2023

# Installs
#install.packages('rts')

# Library imports
library(terra)
library(dplyr)
library(rts)
library(ggplot2)

# Paths to the data
one <- ('../../data/uav/MAIA-exports/20220629/20220629-div32768_clipped.tif')
two <- ('../../data/uav/MAIA-exports/20220705/20220705-div32768_clipped.tif') 
three <- ('../../data/uav/MAIA-exports/20220718/20220718-div32768_clipped.tif')
four <- ('../../data/uav/MAIA-exports/20220814/20220814-div32768_clipped.tif')

# Data to raster
rast1 <- rast(one)
rast2 <- rast(two)
rast3 <- rast(three)
rast4 <- rast(four)

# Resample the raster data, to the same spatial resolution as rast3 ----
rast1 <- resample(rast1, rast3, method = "bilinear")
rast2 <- resample(rast2, rast3, method = "bilinear")
rast4 <- resample(rast4, rast3, method = "bilinear")

### What spatial resolution makes sense for the data?
  # Should go back to Agisoft and re-export the data at a logical resolution.
  # Something like 5 or 10cm which can be scaled easily to the native resolution
  # of the EO imagery?

# For each raster, calculate the spectral indices ----
rast1ndvi <- (rast1$Band_8-rast1$Band_4)/((rast1$Band_8+rast1$Band_4)+.0001)
names(rast1ndvi) <- "ndvi"

rast2ndvi <- (rast2$Band_8-rast2$Band_4)/((rast2$Band_8+rast2$Band_4)+.0001)
names(rast2ndvi) <- "ndvi"

rast3ndvi <- (rast3$Band_8-rast3$Band_4)/((rast3$Band_8+rast3$Band_4)+.0001)
names(rast3ndvi) <- "ndvi"

rast4ndvi <- (rast4$Band_8-rast4$Band_4)/((rast4$Band_8+rast4$Band_4)+.0001)
names(rast4ndvi) <- "ndvi"

rast.list <- list(rast1, rast2, rast3, rast4)



plot(rast1$ndvi)

# Function to calculate NDVI for MAIA/Sentinel-2 data (B8 = NIR, B4 = Red)
#   THIS DOESN'T WORK, IMPROVE UNDERSTANDING OF FUNCS IN R 
terra.ndvi <- function(raster) {
  ndvi_rast <- (raster$Band_8 - raster$Band_4)/((raster$Band_8 + raster$Band_4)+.0001)
  names(ndvi_rast) <- "ndvi"
  raster <- c(raster, ndvi_rast)
  raster
}

lapply(rast.list, terra.ndvi)


# Basic attempt at classifying snow pixels ----
plot(rast1$Band_1)
## from-to-becomes
# classify the values into three groups 
# all values >= 0 and <= 0.6 become 0,
# all values >= 0.6 become 1, where 1 is snow and 0 is no snow
m <- c(0, 0.6, 0, 
       0.6, 2, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rast1snow <- classify(rast1$Band_1, rclmat, include.lowest=TRUE)
plot(rast1snow)

# Add snow raster to all NDVI rasters as extra band
rast1ndvi <- c(rast1ndvi, rast1snow)
# names(rast1ndvi) <- ("ndvi",  "snow") why doesn't this work?
rast2ndvi <- c(rast2ndvi, rast1snow)
rast3ndvi <- c(rast3ndvi, rast1snow)
rast4ndvi <- c(rast4ndvi, rast1snow)

# Create list of all rasters
all_rasters <- rast(list(rast1, rast2, rast3, rast4))

# raster to raster sstack
rast_stack <- c(all_rasters)

rast_stack

# Name the rasters in the stack
names(rast_stack) <- c("20220629", "20220705", "20220718", "20220814")

length(rast_stack)

# Try plotting to check the raster stack is working
plotRGB(rast4, r = 4, b = 3, g = 2, stretch = 'lin')

# Attempting to use rts package ----
ts.stack <- c(rast1ndvi, rast2ndvi, rast3ndvi, rast4ndvi)   # Create stack
d <- c("2022-06-29","2022-07-05","2022-07-18","2022-08-14") # List image dates
d <- as.Date(d)                                             # Make dates date format

rt <- rts(ts.stack, d)                                      # Make raster time series

# Plot raster time series
plot(rt)

# Get location for cell close to centre of plot
location <- cellFromRowCol(rt, 3000, 3000)

# Extract time series for a given cell location
pixel <- rt[location]
pixelsnow <- rast1snow[location]
pixelsnowval <- pixelsnow$Band_1
head(pixel)
plot(pixel)

pixelsnow$Band_1
# Convert extracted pixel to dataframe
df_t <- as.data.frame(t) %>%
  mutate(ID = c(0,0,0,0)) %>% # Add ID variable
  mutate(snow = c(pixelsnowval, pixelsnowval, pixelsnowval, pixelsnowval))             # Add snow variable

xts_df_t <- xts(df_t, order.by = as.Date(rownames(df_t)))

# Get a subset of the time series by sampling pixels at a given interval
interval <- 100000
samples <- seq.int(20069309, 33453309, interval)
length(samples)

# Loop to iterate through pixels and extract values to xts
sampleTS <- xts()                                         # Create empty xts obj
sampleID <- 1                                             # Assign start sample id
for(sample in samples) {                                  # 
  print(sample)
  pixTS <- rt[sample]
  pixelsnow <- rast1snow[sample]
  pixelsnowval <- pixelsnow$Band_1
  df_pixTS <- as.data.frame(pixTS) %>% 
    mutate(ID = c(sampleID, sampleID, sampleID, sampleID)) %>%
    mutate(snow = c(pixelsnowval, pixelsnowval, pixelsnowval, pixelsnowval)) 
  sampleID <- sampleID + 1
  pixTS <- xts(df_pixTS, order.by = as.Date(rownames(df_pixTS)))
  xts_df_t <- c(xts_df_t, pixTS)
}

# Create a df for plotting with ggplot2
plot_df <- fortify(xts_df_t) %>%
  rename(NDVI = V1, date = Index, ) %>%
  mutate(date = as.Date(date))


# Plot line graph of values
ggplot(plot_df, aes(x = date, y = NDVI, color = snow, group = ID)) +
  geom_line() +
  scale_colour_gradient2(low = "red", mid = "yellow", high = "blue") + 
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95)

length(plot_df)
