# Functions for modelling NDVI curves and extracting key metrics
# Calum Hoad, 25 Jan 2024, with help from Jakob Assmann

# Parabolic 2nd order polynomial model ----

# Function for fitting parabolic 2nd order polynomial model
model_fit_parabolic <- function(df) {
  lm(data = df, ndvi ~ poly(doy, 2, raw = T))
}

# Function for calculation of vertex
# https://quantifyinghealth.com/plot-a-quadratic-function-in-r/
find_vertex_parabolic = function(model_fit) {
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

# Overall function for fitting parabolic 2nd order polynomial curves and 
#   extracting key ndvi metrics
model_ndvi_parabolic <- function(data) {
  
  ### This doesn't work, have filtered above instead
  # Filter the dataset, to retain only records where NDVI >= 0.1
  #data <- data %>% 
  #  filter(ndvi >= 0.1) %>%
  #  filter(n_distinct(doy) >= 5)
  
  
  # Use function to fit model
  model <- model_fit_parabolic(data)
  
  # use function to find vertex
  vertex <- find_vertex_parabolic(model)
  
  # Generate predictions for curve plotting
  pred <- predict(model, data)
  
  # Write necessary values back to df
  data <- data %>%
    mutate(ndvi.max = vertex$y, 
           ndvi.max.doy = vertex$x, 
           ndvi.pred = pred)
}


#

# Smoothed spline model


# Smoothed spline model ----

# Function for fitting parabolic 2nd order polynomial model
model_fit_smoothedspline <- function(df) {
  # Using a spline smoother
  smooth.spline(x = df$doy, y = df$ndvi, spar = 0.5)
}

# Function for calculation of vertex
# find_vertex_smoothedspline = function(model_fit) {
#   # Get model coefficients
#   a = model_fit$coefficients[3]
#   b = model_fit$coefficients[2]
#   c = model_fit$coefficients[1]
#   # Determine x coordinate of vertex
#   vertex_x = -b / (2 * a)
#   # Calculate y coordinate of vertex
#   vertex_y = a * vertex_x**2 + b * vertex_x + c
#   # Strip attributes and return as data.frame
#   return(data.frame(
#     x = as.numeric(vertex_x),
#     y = as.numeric(vertex_y)
#   ))
# }

# Overall function
# Define function to model, find vertex, and precict values
model_ndvi_smoothedspline <- function(data) {
  
  ### This doesn't work, have filtered above instead
  # Filter the dataset, to retain only records where NDVI >= 0.1
  # data <- data %>%
  #  filter(ndvi >= 0.1) %>%
  #  filter(n_distinct(doy) >= 5)
  
  # Use function to fit model
  model <- model_fit_smoothedspline(data)
  
  # Generate predictions for curve plotting (for time-period doy 130-300) 
  pred <- predict(model, data.frame(doy = 130:300))
  # pred <- unlist(pred$y)
  
  # use function to find vertex (linear model)
  # vertex <- find_vertex(model)
  # find vertex based on predictions (spline smoother)
  vertex <- data.frame(
    x = pred$x[pred$y == max(pred$y)], 
    y = pred$y[pred$y == max(pred$y)]
  )
  
  # Write necessary values back to df
  data <- suppressMessages(full_join(data, data.frame(
    doy = 130:300,
    ndvi.max = vertex$y,
    ndvi.max.doy = vertex$x,
    ndvi.pred = pred
  ))) %>% 
    # add missing geometries
    mutate(
      geometry = st_geometry(data[1, ])
    )
  return(data)
}




# Beck (2006) double logistic model ----

# Function for fitting Beck
fit_beck  <- function(df) {
  double_log_model <- FitDoubleLogBeck(
    x = df$ndvi,
    t = df$doy,
    weighting = TRUE, # Does this default to TRUE? Re-run model with this explicit.
    #tout = seq(1, 12, length = 365), # Time steps of output, test this?
    plot = TRUE)
}

# Function to generate NDVI predictions using Beck
predict.dbl_log_model <- function(model_object, doys_to_predict) {
  eval(
    model_object$formula,
    c(
      list(t = doys_to_predict),
      split(model_object$params, names(model_object$params))
    )
  )
}

# Predict data for whole year and extract peak ndvi and associated doy
year_in_doys <- 1:365

# Overall function
model_ndvi_beck <- function(data) {
  
  # Fit Beck model to data for pixel
  model <- fit_beck(data)
  
  # Get geometry of id
  geom <- st_geometry(data[1, ])
  
  # Predict values and write back to original dataframe
  data <- suppressMessages(full_join(data, data.frame(
    doy = year_in_doys,
    ndvi_pred = predict.dbl_log_model(model, year_in_doys)
  ))) %>%
    arrange(doy) %>%
    mutate(ndvi_max = max(ndvi_pred)) %>%
    mutate(ndvi_max_doy = doy[which(ndvi_pred == ndvi_max[1])][1]) %>%
    mutate(geometry = geom)
  return(data)
}


# Test functions