# Landsat NDVI metrics, 
# Calum Hoad, 12/12/2023

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
#library(janitor)

# Create a list of the S2 R10m files for each S2 scene
d20230608 <- list.files('../../data/landsat-imagery/LC08_L2SP_011011_20230608_20230614_02_T1/', pattern = ".*_B[0-9]+\\.TIF$", full.names = TRUE)
d20230702 <- list.files('../../data/landsat-imagery/LC09_L2SP_011011_20230702_20230704_02_T1/', pattern = ".*_B[0-9]+\\.TIF$", full.names = TRUE)
d20230708 <- list.files('../../data/landsat-imagery/LC08_L2SP_013011_20230708_20230718_02_T1/', pattern = ".*_B[0-9]+\\.TIF$", full.names = TRUE)
d20230726 <- list.files('../../data/landsat-imagery/LC08_L2SP_011011_20230726_20230805_02_T1/', pattern = ".*_B[0-9]+\\.TIF$", full.names = TRUE)
d20230801 <- list.files('../../data/landsat-imagery/LC09_L2SP_013011_20230801_20230803_02_T1/', pattern = ".*_B[0-9]+\\.TIF$", full.names = TRUE)
d20230817 <- list.files('../../data/landsat-imagery/LC09_L2SP_013011_20230817_20230822_02_T1/', pattern = ".*_B[0-9]+\\.TIF$", full.names = TRUE)
d20230910 <- list.files('../../data/landsat-imagery/LC08_L2SP_013011_20230910_20230918_02_T1/', pattern = ".*_B[0-9]+\\.TIF$", full.names = TRUE)
d20230920 <- list.files('../../data/landsat-imagery/LC09_L2SP_011011_20230920_20230922_02_T1/', pattern = ".*_B[0-9]+\\.TIF$", full.names = TRUE)

# List of imagery dates, for later use
dates <- c('2023-06-08', 
           '2023-07-02',
           '2023-07-08', 
           '2023-07-26', 
           '2023-08-01', 
           '2023-08-17', 
           '2023-09-10', 
           '2023-09-20')

# Get the Sentinel-2 10m bands as raster stack, project + crop to extent of UAV imagery
# As function?

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
ls.data <- list(d20230608,
                d20230702,
                d20230708,
                d20230726,
                d20230801,
                d20230817,
                d20230910,
                d20230920) 

crs(ls, proj = TRUE, describe = FALSE)
# Get uAV imagery over plot to use for cropping - RE-EXPORT UAV IMAGERY SO RE-PROJECT IS AVOIDED
uav <- project(rast('../../data/uav/M3M-exports/5cm/20230702-clipped-5cm-div128.tif'), 'epsg:32621')

# Create function for stacking rasters from lists of filepaths, then cropping to extent of UAV imagery
import_ls <- function(x, rast.targ.proj) {
  x <- rast(x)
  if(crs(x, proj = TRUE, describe = FALSE) == "+proj=utm +zone=22 +datum=WGS84 +units=m +no_defs") {
    print(paste0('Reprojecting', x, 'to UTM21N'))
    x <- project(x, rast.targ.proj)
  }
  x <- x %>% crop(uav)
  #x <- set.names(object = x, nm = c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp'))
}

# Create spatRast object for one Landsat image with desired projection
target.rast <- rast(d20230708)
# Check the projection of the target raster is correct
crs(target.rast, proj = TRUE, describe = FALSE)

# Import the data
ls.data.import <- lapply(ls.data, import_ls, rast.targ.proj = target.rast)

# Check the function worked by plotting rasters from list
plot(ls.data.import[[8]])

ls.data.import

# Rename bands to make reading logical
names(ls.data.import[[1]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'b10')
names(ls.data.import[[2]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'b10')
names(ls.data.import[[3]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'b10')
names(ls.data.import[[4]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'b10')
names(ls.data.import[[5]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'b10')
names(ls.data.import[[6]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'b10')
names(ls.data.import[[7]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'b10')
names(ls.data.import[[8]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'b10')


# Apply function to calculate NDVI ----
ls_ndvi <- function(x) {
  x <- (x$nir-x$red)/(x$nir+x$red)
  # names(x) <- c('ndvi')
}

ls.ndvi <- lapply(ls.data.import, ls_ndvi)

# Stack NDVI rasters into single spatRast
ls.ndvi <- rast(ls.ndvi)

# Name spatRaster layers with dates
names(ls.ndvi) <- dates

plot(ls.ndvi)

# Extract raster time series to points
ls.ndvi.points <- st_as_sf(as.points(ls.ndvi, values = TRUE)) %>%
  mutate(id = row_number())


# Apply parabolic 2nd order polynomial to every pixel in the df ----

# Get a dataframe of points from the raster
ls.ndvi.long <- ls.ndvi.points %>%
  pivot_longer(!geometry & !id, names_to = 'doy', values_to = 'ndvi') %>%
  mutate(doy = as.Date(doy)) %>%
  group_by(id)

# Function for fitting parabolic 2nd order polynomial model
model_fit <- function(df) {
  lm(data = df, ndvi ~ poly(doy, 2, raw = T))
}

# Function for calculation of vertex
# https://quantifyinghealth.com/plot-a-quadratic-function-in-r/
find_vertex = function(model_fit) {
  # Get model coefficients
  a = model_fit$coefficients[3]
  b = model_fit$coefficients[2]
  c = model_fit$coefficients[1]
  # Determine x coordinate of vertex
  vertex_x = -b / (2 * a)
  # Calculate y coordinate of vertex
  vertex_y = a * vertex_x**2 + b * vertex_x + c
  # Strip attributes and return as data.frame
  return(data.frame(
    x = as.numeric(vertex_x),
    y = as.numeric(vertex_y)
  ))
}

# Define function to model, find vertex, and precict values
model_ndvi <- function(data) {
  
  # Use function to fit model
  model <- model_fit(data)
  
  # use function to find vertex
  vertex <- find_vertex(model)
  
  # Generate predictions for curve plotting
  pred <- predict(model, data)
  
  # Write necessary values back to df
  data <- data %>%
    mutate(ndvi.max = vertex$y, 
           ndvi.max.doy = vertex$x, 
           ndvi.pred = pred)
}

# Following tutorial here: https://data-se.netlify.app/2018/12/10/new-split-apply-combine-variant-in-dplyr-group-split/

# Apply model_ndvi to data using group_modify
ls.modelled.ndvi <- ls.ndvi.long %>%
  group_modify(~ model_ndvi(.x)) %>%
  mutate(ndvi.max.date = as_date(ndvi.max.doy),
         # Calculate fractional day of year, by adding decimal to date
         ndvi.max.doy = yday(ndvi.max.date) + (ndvi.max.doy - floor(ndvi.max.doy)))

# Plotting ----
ggplot(ls.modelled.ndvi, aes(group_by = id)) +
  geom_line(aes(x = doy, y = ndvi.pred, group = id, color = ndvi.max)) +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max)) +
  scale_color_viridis()

# Outputs ---- 

# Wide format
ls.modelled.export.wide <- ls.modelled.ndvi %>%
  group_by(id) %>%
  filter(row_number() == 1) %>%
  select(!doy & !ndvi & !ndvi.pred)

st_write(st_as_sf(ls.modelled.export.wide),  '../../data/sentinel-2/output/ls_modelled_point_wide.shp')
write.csv2(ls.modelled.export.wide,  '../../data/sentinel-2/output/ls_modelled_point_wide.csv')

# Long format
ls.modelled.export.long <- ls.modelled.ndvi %>%
  rename(doy.obs = 'doy', 
         ndvi.obs = 'ndvi')

st_write(st_as_sf(ls.modelled.export.long),  '../../data/sentinel-2/output/ls_modelled_point_long.shp')
write.csv2(ls.modelled.export.long,  '../../data/sentinel-2/output/ls_modelled_point_long.csv')


# Output the stacked landsat rasters for future use:
writeRaster(ls.data.import[[1]], '../../data/landsat-imagery/terra/2023-06-08.tif')
writeRaster(ls.data.import[[2]], '../../data/landsat-imagery/terra/2023-07-02.tif')
writeRaster(ls.data.import[[3]], '../../data/landsat-imagery/terra/2023-07-08.tif')
writeRaster(ls.data.import[[4]], '../../data/landsat-imagery/terra/2023-07-26.tif')
writeRaster(ls.data.import[[5]], '../../data/landsat-imagery/terra/2023-08-01.tif')
writeRaster(ls.data.import[[6]], '../../data/landsat-imagery/terra/2023-08-17.tif')
writeRaster(ls.data.import[[7]], '../../data/landsat-imagery/terra/2023-09-10.tif')
writeRaster(ls.data.import[[8]], '../../data/landsat-imagery/terra/2023-09-20.tif')


# Check data
ggplot(ls.modelled.export.long, group_by = id) +
  geom_point(aes(x = doy.obs, y = ndvi.pred)) +
  geom_line(aes(x = doy.obs, y = ndvi.obs, group = id))
