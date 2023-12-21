# Quick script to fit some double logistic curves on simulated data
# Jakob J. Assmann jakob.assmann@uzh.ch 21 December 2024

# Dependencies
# install.packages(c("bfast", "phenopix"))
# install.packages("greenbrown", repos="http://R-Forge.R-project.org")
library(tidyverse)
library(greenbrown)
library(ggplot2)
library(cowplot)

# Generate data based on Elmore et al. 2011 Eq 3
#       https://doi.org/10.1111/j.1365-2486.2011.02521.x

# Note: We will fit later using Beck et al. 2006 Eq 3
#       https://doi.org/10.1016/j.rse.2005.10.021

# Paramters for simulation
# Using m3/m4 = doy spring greenup half way point
# and m5/m6 = autumn borwning half way point
m_1 <- -0.1 # Average phenology over winter
m_2 <- 0.6 # Difference between summer and winter (amplitude)
m_3 <- 100 # Intercept of exponent 1 (spring greenup)
m_4 <- m_3/165 # Slope of exponent 1 (spring greenup)
m_5 <- 100 # Intercept of exponent 2 (autumn browning)
m_6 <- m_5/250 # Slope of exponent 2 (autumn browning)

# Function to calculate phenology
dbl_log_phen <- function(doy) {
    ndvi <- m_1 + m_2*((1/(1+exp(m_3-m_4*doy)) - 1/(1+exp(m_5-m_6*doy))))
}

# Simulate data (adding some error)
year_in_doys <- 1:365
phen_sim <- data.frame(
    doy = year_in_doys,
    ndvi = dbl_log_phen(year_in_doys) + rnorm(length(year_in_doys), sd = 0.025)
)
# Subsample to simulate sparse satelite obeservations 
# and arrange according to doy (important for fitting the curve later)
phen_sim <- sample_n(phen_sim, 30) %>% arrange(doy)

# Visualise data
ggplot(phen_sim) +
    geom_point(aes(x = doy, y = ndvi)) +
    scale_y_continuous(limits = c(-0.3, 0.9), breaks = seq(-0.3, 0.9, 0.3)) +
    scale_x_continuous(limits = c(0, 365)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_cowplot()

# Fit curve using greenbrown package function using Beck et al. 2006 Eq 3
# and including weighting of "overestimated NDVI values" (see paper for details)
# Note: greenbrown also allows for fitting with Elmore et al. 2011 Eq 4
#       but I tried (!) and with thin samples as we have, the gentle sloping
#       summer phenology introduced by the addtional term can not reliably
#       estimated using the low frequency (gap-y) data that we have. For that
#       use "FitDoubleLogElmore()" with the same syntax. 
double_log_model <- FitDoubleLogBeck(
    x = phen_sim$ndvi,
    t = phen_sim$doy,
    plot = TRUE)

# Check out the model object
double_log_model

# Write helper functions to generate predictions for the whole year
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
phen_sim <- phen_sim %>%
    full_join(data.frame(
        doy = year_in_doys,
        ndvi_pred = predict.dbl_log_model(double_log_model, year_in_doys)
    )) %>%
    arrange(doy) %>%
    mutate(ndvi_max = max(ndvi_pred)) %>%
    mutate(ndvi_max_doy = doy[which(ndvi_pred == ndvi_max[1])][1])
    

# Plot the resulting curve fit
ggplot(phen_sim) +
    geom_point(aes(x = doy, y = ndvi), colour = "blue") +
    geom_line(aes(x = doy, y = ndvi_pred)) +
    scale_y_continuous(limits = c(-0.3, 0.9), breaks = seq(-0.3, 0.9, 0.3)) +
    scale_x_continuous(limits = c(0, 365)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    annotate("point", 
        x = phen_sim$ndvi_max_doy[1], 
        y = phen_sim$ndvi_max[1],
            colour = "red") +
    theme_cowplot()
