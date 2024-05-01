# Met data for Kluane,
# Calum Hoad, translated from Gwenn Flowers MatLab script

# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(cowplot)
library(lubridate)

# Define file paths
Outpost_Health <- "../../data/met-data/Outpost/ConcatenatedData/Outpost_AWS_HealthData.csv"
Outpost_HalfHour <- "../../data/met-data/Outpost/ConcatenatedData/Outpost_AWS_HalfHourData.csv"
Outpost_FiveMin_CMA6 <- "../../data/met-data/Outpost/ConcatenatedData/Outpost_AWS_FiveMinData_CMA6.csv"
Outpost_FiveMin_CNR4 <- "../../data/met-data/Outpost/ConcatenatedData/Outpost_AWS_FiveMinData_CNR4.csv"

# Read CSV files
A_Health <- read_csv(Outpost_Health)
A_HalfHour <- read_csv(Outpost_HalfHour, skip = 1) %>%
  filter(!row_number() %in% c(1, 2)) %>%
  rename(Rain = 'Rain_mm_Tot') %>%
  mutate(Rain = as.integer(Rain))
  
A_FiveMin_CMA6 <- read.csv(Outpost_FiveMin_CMA6, header = TRUE, skip = 4)
A_FiveMin_CNR4 <- read.csv(Outpost_FiveMin_CNR4, header = TRUE, skip = 4)

# Plotting functions for each dataset
plot_health <- function(data) {
  p1 <- ggplot(data, aes(x = TIMESTAMP, y = Rain)) +
    geom_line() +
    labs(x = "Time", y = "Record number", title = "Outpost Health Data") +
    theme_minimal()
  
  p2 <- ggplot(data, aes(x = TIMESTAMP)) +
    geom_line(aes(y = BVmax), color = "blue") +
    geom_line(aes(y = BVmin), color = "red") +
    geom_line(aes(y = BVmean), color = "green") +
    labs(x = "Time", y = "Battery voltage (V)", title = "Battery Voltage") +
    theme_minimal()
  
  p3 <- ggplot(data, aes(x = TIMESTAMP)) +
    geom_line(aes(y = Tmax), color = "blue") +
    geom_line(aes(y = Tmin), color = "red") +
    geom_line(aes(y = Tmean), color = "green") +
    labs(x = "Time", y = "Panel Temp (C)", title = "Panel Temperature") +
    theme_minimal()
  
  plot_grid(p1, p2, p3, ncol = 1)
}

data <- A_HalfHour %>% mutate(help.group = 0)

typeof(as_date(data$TIMESTAMP, format = '%YYYY-%MM-%DD %HH:%MM'))

plot_half_hour <- function(data) {
  p1 <- ggplot(data, aes(x = ymd_hm(TIMESTAMP), y = Rain)) +
    geom_line(aes(group = help.group)) +
    labs(x = "Time", y = "DT (m)", title = "Outpost Half Hour Data") +
    theme_minimal()
 
     p1
  p2 <- ggplot(data, aes(x = TIMESTAMP %>% lubridate(), y = BP)) +
    geom_line() +
    labs(x = "Time", y = "BP (hPa)", title = "Barometric Pressure") +
    theme_minimal()
  
  p3 <- ggplot(data, aes(x = TIMESTAMP %>% lubridate(), y = Rain)) +
    geom_line() +
    labs(x = "Time", y = "Rain (mm/30 min)", title = "Rainfall") +
    theme_minimal()
  
  plot_grid(p1, p2, p3, ncol = 1)
}

# Assuming similar structure for the other datasets, you would create similar plotting functions
# For example, plot_five_min_CMA6 and plot_five_min_CNR4

# Call the plotting functions
plot_health(A_Health)
plot_half_hour(A_HalfHour)
# Call the other plotting functions similarly

# Calculate an annual rainfall
rain.stats <- data %>% mutate(year = lubridate::year(TIMESTAMP)) %>%
  group_by(year) %>%
  summarize(total = sum(Rain, na.rm = TRUE)) %>%
  mutate(average = )
