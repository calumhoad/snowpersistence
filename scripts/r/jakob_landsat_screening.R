# Landsat scene selection script for Calum
# Jakob J. Assmann jakob.assmann@uzh.ch 09 November 2023

library(rgee)
library(sf)
library(tidyverse)

# Initalise EE
ee_Initialize()

# Define test polygon
test_poly_sf <- list(matrix(c(-138.90125, 69.58413,
                              -138.88988, 69.58358,
                              -138.89147, 69.58095,
                              -138.90298, 69.57986,
                              -138.90125, 69.58413),
                            ncol = 2, byrow = TRUE)) %>%
  st_polygon() %>%
  st_sfc(crs = 4326) %>%
  st_sf()
# Upload test poly to GEE
test_poly_ee <- sf_as_ee(test_poly_sf)

# Get Landsat C2 ICs
ls5_1 <- rgee::ee$ImageCollection("LANDSAT/LT05/C02/T1_L2")
ls5_2 <- rgee::ee$ImageCollection("LANDSAT/LT05/C02/T2_L2")
ls7_1 <- rgee::ee$ImageCollection("LANDSAT/LE07/C02/T1_L2")
ls7_2 <- rgee::ee$ImageCollection("LANDSAT/LE07/C02/T2_L2")
ls8_1 <- rgee::ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")
ls8_2 <- rgee::ee$ImageCollection("LANDSAT/LC08/C02/T2_L2")

## Merge into one collection and pre-filter by year and date
# Start and end date (year)
start_date <- "2000-01-01"
end_date <- "2023-01-01" # Add an extra day here as filtering is exclusive of last day
# Add start and end DoY for within year
start_doy <- 152
end_doy <- 243
LS_COLL <- ls5_1$
  merge(ls7_1$
          merge(ls8_1$
                  merge(ls5_2$
                          merge(ls7_2$
                                  merge(ls8_2)))))$
  filterDate(start_date, end_date)$
  filter(rgee::ee$Filter$calendarRange(start_doy,
                                       end_doy,
                                       "day_of_year"))$
  #!!! Do any other filtering you need to here e.g., cloud cover <= 70%
  filter(ee$Filter$lte("CLOUD_COVER", 70))
         
# Filter scenes by aoi (test_poly)
LS_COLL_filtered <- LS_COLL$filterBounds(test_poly_ee)

# Get SCENE IDs and store in tibble
LS_IDs_sf <- ee$FeatureCollection(LS_COLL_filtered$map(function(x){
  ee$Feature(NULL, list(ID = x$get("L1_LANDSAT_PRODUCT_ID")))
  })) %>% ee_as_sf() %>% 
  st_drop_geometry() %>%
  # Add column for quality score (to be assigned manually later)
  mutate(quality_score = NA)

# Loop over IDs, plot scene, promt for 
index <- seq_along(LS_IDs_sf$ID)
for(i in index[1:3] # remove indices [] here to run through full list
    ){
  # Get Scene ID
  scene_ID <- LS_IDs_sf$ID[i]
  cat("Scene ID:", scene_ID, "\n")
  
  # Get image
  image <- LS_COLL_filtered$filter(ee$Filter$eq("L1_LANDSAT_PRODUCT_ID", scene_ID))$first()
  
  # Centre map on Scene (alternatively centre on AOI)
  Map$centerObject(image)
  
  # Add image
  print(Map$addLayer(image,
               visParams = list(bands = c("SR_B3", "SR_B2", "SR_B3"),
                                min = 0, # Might have to tweak these values here a little
                                max = round(0.3*65535)) 
                                # Also consider adding gamma strecth (google)
               ) +
          Map$addLayer(test_poly_ee, 
                       visParams = list(color = "red", fillColor = "00000000")))
  
  quality_score <- readline("Assign quality score [1,2,3, e = exit]:\n")
  while(!(quality_score %in% c("1", "2", "3", "e"))){
    quality_score <- readline("Invalid input!\nAssign quality score [1,2,3]:\n")
  }
  if(quality_score == "e") stop("User aborted. Scnene not assigned")
  LS_IDs_sf$quality_score[i] 
}

# Need to add a line down here to assign the quality score back to the 
# row in the dataframe