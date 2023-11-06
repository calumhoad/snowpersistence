# Script to QC drone data
# Calum Hoad, 20231011

library(raster)
library(dplyr)
library(ggplot2)

# Path to TIFF data
ndvi_dir <- '../2023-field/DATA/NDVI-data/TIFF/'

# List of your GeoTIFF files along with their corresponding labels
geotiff_files <- c(file.path(ndvi_dir, '20230702-envi-ndvi.tif'), 
                   file.path(ndvi_dir, '20230712-envi-ndvi.tif'), 
                   file.path(ndvi_dir, '20230718-envi-ndvi.tif'), 
                   file.path(ndvi_dir, '20230726-envi-ndvi.tif'), 
                   file.path(ndvi_dir, '20230726-envi-ndvi-m3m.tif'))

img <- raster('C:/Users/s1437405/Documents/PhD_Local/2023-field/DATA/NDVI-data/TIFF/20230702-envi-ndvi.tif')

plot(img) # This works

# Create an empty list to store the histogram data frames
hist_data_list <- list()

# Iterate through each GeoTIFF file
for (i in 1:length(geotiff_files)) {
  # Read the GeoTIFF file
  img <- raster(geotiff_files[i])
  
  # Extract pixel values as a vector, filtering out missing and non-finite values
  pixel_values <- values(img)
  pixel_values <- pixel_values[!is.na(pixel_values) & is.finite(pixel_values)]
  
  # Create a data frame with pixel values and a label
  df <- data.frame(
    Pixel_Value = pixel_values,
    Label = rep(geotiff_files[i], length(pixel_values))  # Use the file name as the label
  )
  
  hist_data_list[[i]] <- df
}

# Combine all data frames into one
hist_data <- do.call(rbind, hist_data_list)

# Create the histogram plot using ggplot2
ggplot(hist_data, aes(x = Pixel_Value, fill = Label)) +
  geom_histogram(binwidth = 1, position = 'dodge') +
  labs(x = 'Pixel Value', y = 'Frequency') +
  ggtitle('Histograms of Multiple GeoTIFFs') +
  scale_fill_manual(values = rainbow(length(geotiff_files))) +  # Color bars differently
  theme_minimal()
