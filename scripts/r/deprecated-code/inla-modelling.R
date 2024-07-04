# Script to model ndvi.max and ndvi.doy as a function of snow using R-INLA
# Script framework provided by Jakob Assmann
# Calum Hoad, 8 March 2024

# Load INLA
# Installation instructions here: 
# https://www.r-inla.org/download-install

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(INLA)
library(terra) # Quick fix for visualising matrices
library(tidyverse)
library(sf)

# Load in the data
# Blaesedalen
s2.bl <- read.csv('../../data/combined-ndvi-snow/s2-bl-smooth-joined.csv') %>%
  #st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
  group_by(id) %>%
  filter(row_number() == 7) %>%
  ungroup() 

# Defining a matrix ----

# Get the min X and max X, min Y and max Y
xmin <- min(s2.bl$X)
xmax <- max(s2.bl$X)
X <- seq(xmin, xmax, 10)
ymin <- min(s2.bl$Y)
ymax <- max(s2.bl$Y)
Y <- seq(ymin, ymax, 10)

# Turn it into a grid
mat.grid <- matrix(NA, nrow = length(Y), ncol = length(X))
mat.grid

# Optionally, you can set row and column names if needed
rownames(mat.grid) <- rev(Y)
colnames(mat.grid) <- X
mat.grid

# Function to fill matrix with values
val_to_grid <- function(df, var, grid) {
  for (i in 1:nrow(df)) {
    grid[which(Y == df$Y[i]), 
             which(X == df$X[i])] <- df[[var]][i]
  }
  return(grid)
}

# Run function over the data to pull values to grid
ndvi.max.grid <- val_to_grid(s2.bl, 'ndvi.max', mat.grid)
ndvi.doy.grid <- val_to_grid(s2.bl, 'ndvi.max.doy', mat.grid)
snow.grid <- val_to_grid(s2.bl, 'snow.auc', mat.grid)

plot(rast(snow.grid))

# Plot original data to check
ggplot() +
  geom_sf(data = s2.bl %>% st_as_sf(coords = c('X', 'Y'), crs = 32621) %>%
            st_buffer(dist = 5, endCapStyle = "SQUARE"),
          aes(fill = ndvi.max))


# Models ----

# Example below modified from INLA documentation.
# inla.doc("matern2d")

# Define grid for generated data
nrow <- 20
ncol <- 30
n <- nrow*ncol
n
# Set a noise parameter
s.noise <- 1

# Define matrix
zi.mat <- matrix(NA,nrow=nrow,ncol=ncol)

# Fill matrix iwth a clear spatial pattern
i <- 1:nrow
for(j in 1:ncol)  zi.mat[i,j] <- 3*exp(-(i-j)^2/4)
plot(rast(zi.mat))


# Generate random noise
noise.mat <- matrix(
  rnorm(nrow*ncol, sd=s.noise),
  nrow,
  ncol)

# Add noise to data with spatial component
y.mat <- zi.mat + noise.mat
plot(rast(y.mat))


# cCnvert matrix to the internal vector representation in INLA
# INLA can't handle matrices, only vectors, but ther is a handy function
ndvi <- inla.matrix2vector(ndvi.max.grid)
doy <- inla.matrix2vector(ndvi.doy.grid)
snow <- inla.matrix2vector(snow.grid)

# Specify ids for every nodes (as this is a regular grid, it's simply sequential)
node <- 1:(length(X)*length(Y))
nrow <- length(X)
ncol <- length(Y)
length(node)
nrow
ncol

# Specify model formula (only intercept and latent spatial effect)
# f is a function for defining general gaussian models with INLA
formula <- y ~ 1 + f(node, # Predictor
                     model="matern2d", # Matern model
                     nu = 1, # Smoothing parameter for the Matern2d-model, possible values are c(0, 1, 2, 3)
                     nrow=nrow, # spatial configuration
                     ncol=ncol
                     # Priors for hyperparameter (okay to comment out then tweak later )
                     # hyper = list(range = list(param =c(1, 1),
                     #                         prior = "loggamma",initial=1),
                     #            prec = list(param=c(1, 1)))
)

# Combine repsone and predictor into a data frmae
data <- data.frame(y=ndvi,node=node)

# fit the model (assuming gaussian residual distribution)
result <- inla(formula, 
               family="gaussian",
               data=data, 
               control.predictor = list(compute = TRUE),
               control.family = list(hyper = list(theta = list(initial = log(1/s.noise^2),
                                                               fixed = FALSE))),
               keep=T)

summary(result)

# Convert spatial latent effect predictions into a matrix
pred.mat <- INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)
plot(rast(pred.mat))
# Visualise residuals
plot(rast(ndvi.max.grid) - rast(pred.mat))
hist(rast(y.mat) - rast(pred.mat))

##### Same but with a linear predictor

# Weaker noise
s.noise = 0.5
noise.mat <- matrix(
  rnorm(nrow*ncol, sd=s.noise),
  nrow,
  ncol)

plot(rast(noise.mat))

y.mat <- matrix(s)
# Linear predictor
a.mat <- matrix(NA,nrow=nrow,ncol=ncol)
a.mat[] <- sample(1:n/1000, n)


y.mat <- zi.mat + (a.mat * 12) + noise.mat
plot(rast(y.mat))     

# specify formular
formula <- y ~ 1 + a + f(node, # Predictor
                         model="matern2d", # Matern model
                         nu = 1, # Smoothing parameter for the Matern2d-model, possible values are c(0, 1, 2, 3)
                         nrow=nrow, # spatial configuration
                         ncol=ncol
)

# Reassing data
y <- inla.matrix2vector(y.mat)
a <- inla.matrix2vector(a.mat)
node <- 1:n
data <- data.frame(y=ndvi,a=snow,node=node)

# Fit model
result <- inla(formula, 
               family="gaussian",
               data=data, 
               control.predictor = list(compute = TRUE),
               control.family = list(hyper = list(theta = list(initial = log(1/s.noise^2),
                                                               fixed = FALSE))),
               # Switching on some addtional computations including marginal predictions
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE),
               keep=T)
summary(result)

# Plot results in space (including both spatial and linear predictor)
pred.mat <- INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)
plot(rast(y.mat))
plot(rast(pred.mat))
plot(rast(zi.mat + 12 * a.mat))

# Plot predicted vs. observed
plot(result$summary.fitted.values$mean, as.vector(y.mat))

