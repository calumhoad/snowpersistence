# Sentinel-2 NDVI Metrics
# Calum Hoad, 05/12/2023

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
d20230626 <- list.files('../../data/sentinel-2/20230626/S2B_MSIL2A_20230626T153819_N0509_R011_T21WXT_20230626T183843.SAFE/GRANULE/L2A_T21WXT_A032927_20230626T153935/IMG_DATA/R10m/', full.names = TRUE)
d20230708 <- list.files('../../data/sentinel-2/20230708/S2A_MSIL2A_20230708T152811_N0509_R111_T21WXT_20230708T214952.SAFE/GRANULE/L2A_T21WXT_A042007_20230708T152945/IMG_DATA/R10m/', full.names = TRUE)
d20230726 <- list.files('../../data/sentinel-2/20230726/S2B_MSIL2A_20230726T153819_N0509_R011_T21WXT_20230726T184037.SAFE/GRANULE/L2A_T21WXT_A033356_20230726T153828/IMG_DATA/R10m/', full.names = TRUE)
d20230729 <- list.files('../../data/sentinel-2/20230729/S2B_MSIL2A_20230729T154819_N0509_R054_T21WXT_20230729T181552.SAFE/GRANULE/L2A_T21WXT_A033399_20230729T154940/IMG_DATA/R10m/', full.names = TRUE)
d20230807 <- list.files('../../data/sentinel-2/20230807/S2A_MSIL2A_20230807T152811_N0509_R111_T21WXS_20230807T212701.SAFE/GRANULE/L2A_T21WXS_A042436_20230807T153242/IMG_DATA/R10m/', full.names = TRUE)
d20230808 <- list.files('../../data/sentinel-2/20230808/S2B_MSIL2A_20230808T154819_N0509_R054_T21WXS_20230914T102422.SAFE/GRANULE/L2A_T21WXS_A033542_20230808T154852/IMG_DATA/R10m/', full.names = TRUE)
d20230817 <- list.files('../../data/sentinel-2/20230817/S2A_MSIL2A_20230817T152941_N0509_R111_T21WXT_20230817T214159.SAFE/GRANULE/L2A_T21WXT_A042579_20230817T153311/IMG_DATA/R10m/', full.names = TRUE)
d20230923 <- list.files('../../data/sentinel-2/20230923/S2A_MSIL2A_20230922T154951_N0509_R054_T21WXT_20230922T234400.SAFE/GRANULE/L2A_T21WXT_A043094_20230922T155005/IMG_DATA/R10m/', full.names = TRUE)

# List of imagery dates, for later use
dates <- c('2023-06-26', 
           '2023-07-08',
           '2023-07-26', 
           '2023-07-29', 
           '2023-08-07', 
           '2023-08-08', 
           '2023-08-17', 
           '2023-09-23')

# Get the Sentinel-2 10m bands as raster stack, project + crop to extent of UAV imagery
# As function?

# Create list where each item in the list is another list,
#   containing the filepath to each imagery band
s2.data <- list(d20230626, 
                d20230708, 
                d20230726, 
                d20230729, 
                d20230807, 
                d20230808, 
                d20230817, 
                d20230923)

# Get uAV imagery over plot to use for cropping - RE-EXPORT UAV IMAGERY SO RE-PROJECT IS AVOIDED
uav <- project(rast('../../data/uav/M3M-exports/5cm/20230702-clipped-5cm-div128.tif'), 'epsg:32621')

# Create function for stacking rasters from lists of filepaths, then cropping to extent of UAV imagery
import_s2 <- function(x) {
    x <- rast(x) %>%
      #project('epsg:32621') %>%
      crop(uav)
    #x <- set.names(object = x, nm = c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp'))
}

# Import the data
s2.data.import <- lapply(s2.data, import_s2)

# Check the function worked by plotting rasters from list
plot(s2.data.import[[6]])

# Rename bands to make reading logical
names(s2.data.import[[1]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[2]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[3]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[4]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[5]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[6]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[7]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')
names(s2.data.import[[8]]) <- c('aot', 'blue', 'green', 'red', 'nir', 'tci1', 'tci2', 'tci3', 'wvp')


# Apply function to calculate NDVI ----
s2_ndvi <- function(x) {
  x <- (x$nir-x$red)/(x$nir+x$red)
  # names(x) <- c('ndvi')
}

s2.ndvi <- lapply(s2.data.import, s2_ndvi)

# Stack NDVI rasters into single spatRast
s2.ndvi <- rast(s2.ndvi)

# Name spatRaster layers with dates
names(s2.ndvi) <- dates

plot(s2.ndvi)

# Extract raster time series to points
s2.ndvi.points <- st_as_sf(as.points(s2.ndvi, values = TRUE)) %>%
  mutate(id = row_number())


# Apply parabolic 2nd order polynomial to every pixel in the df ----

# Get a dataframe of points from the raster
s2.ndvi.long <- s2.ndvi.points %>%
  pivot_longer(!geometry & !id, names_to = 'doy', values_to = 'ndvi') %>%
  mutate(doy = as.Date(doy))

# Group dataframe by id
s2.ndvi.long <- s2.ndvi.long %>%
  group_by(id)

# Use group map to apply functions to each group in dataframe
s2.ndvi.long.test <- group_map(data = s2.ndvi.long, )

# to subsample the df, add this to first line ([sample(nrow(s2.ndvi.points), 3), ])

s2.ndvi.sample <- s2.ndvi.points[sample(nrow(s2.ndvi.points), 1), ] %>%
  pivot_longer(!geometry & !id, names_to = 'doy', values_to = 'ndvi') %>%
  mutate(doy = as.Date(doy))

# plot to check what the data looks like
ggplot(s2.ndvi.long, aes(x = doy, y = ndvi)) +
         geom_point(aes(color = id)) +# +
         geom_line(aes(group = id)) 
         
s2.ndvi.long.nest <- s2.ndvi.long %>%
  nest()



# Jakob's code ----
# Generate some random data
x <- rnorm(100)
my_data <- data.frame(
  x = x,
  y = -x**2 - x + rnorm(100)
)


# Fit second order polynomial
# f(x) = ax^2 + bx + c
model_fit <- lm(ndvi ~ poly(doy, 2, raw = T), data = s2.ndvi.sample)

model_fit <- function(df) {
  broom::tidy(lm(data = df, ndvi ~ poly(doy, 2, raw = T)))
}

model_fit <- function(df) {
  lm(data = df, ndvi ~ poly(doy, 2, raw = T))
}

model_fit(data = s2.ndvi.sample, x = s2.ndvi.sample$doy, y = s2.ndvi.sample$ndvi)

# Calculate vertex
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
vertex <- find_vertex(model_fit) %>%
  mutate(x = as.Date(x))

# Generate predictons (for curve plotting)
s2.ndvi.sample$preds <- predict(model_fit, s2.ndvi.sample)


s2.ndvi.long <- s2.ndvi.long %>%
  mutate(model = model_fit(data = s2.ndvi.long, x = doy, y = ndvi)) %>%
  mutate(vertex = find_vertex(model))
# Plot the whole thing
ggplot(s2.ndvi.sample) +
  geom_point(aes(x = doy, y = ndvi)) +
  geom_line(aes(x = doy, y = preds)) +
  annotate("point", x = vertex$x, y = vertex$y,
           colour = "red",
           size = 5) +
  annotate("text",
           x = vertex$x, y = vertex$y,
           label = paste0("Vertex (", round(vertex$x, 2), ",", round(vertex$y, 2), ")"),
           colour = "red",
           size = 10,
           vjust = -1) +
  theme_classic()

# Should iterate before making df long format?
# Or use group_by(). For every group (id) apply the function, 
# then move to the next group.

# Following tutorial here: https://data-se.netlify.app/2018/12/10/new-split-apply-combine-variant-in-dplyr-group-split/
test <- s2.ndvi.long.nest %>%
  mutate(model = purrr::map(.f = model_fit, . = data)) %>%
  mutate(vertex = purrr::map(model, .f = find_vertex)) %>%
  mutate(preds = purrr::map(model, .f = predict))

test2 <- purrr::map(s2.ndvi.long.nest, .f = model_ndvi)

# Write it all into one function? Get the parameters out and then mutate the df
# at the end?

unnested <- test %>% 
  unnest(cols = c(data, model, vertex, preds))

test[[5]][[1]]

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

test3 <- model_ndvi(as.data.frame(st_drop_geometry(s2.ndvi.sample)))

# This works
test4 <- s2.ndvi.long %>%
  group_modify(~ model_ndvi(.x))

# Fix data formatting
test4 <- test4 %>%
  
  #ungroup() %>%
  #rename(ndvi.max.doy = 'ndvi.max.doy$x') %>%
  mutate(ndvi.max.doy = as_date(ndvi.max.doy))

test4


# Plotting
ggplot(test4, aes(group_by = id)) +
  geom_line(aes(x = doy, y = ndvi.pred, group = id, color = ndvi.max)) +
  geom_point(aes(x = ndvi.max.doy, y = ndvi.max)) +
  scale_color_viridis()
