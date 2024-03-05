# Elevation Data

library(tidyverse)
library(sf)
library(elevatr)


# read in
meta<- read.csv("./miller/pando_sample_metadata_filled.csv")

# Data frame with latitude and longitude
lat_lon_list <- data.frame(lat = meta$lat, lon = meta$lon, sample = meta$sample)

# Filter out rows with missing coordinates
lat_lon_list <- lat_lon_list[complete.cases(lat_lon_list), ]

# Convert to an sf object
lat_lon_sf <- st_as_sf(lat_lon_list, coords = c("lon", "lat"), crs = 4326)

# Get elevation data
lat_lon_list$elev_data <- get_elev_point(lat_lon_sf)

# Un-uglify Columns

lat_lon_list<- lat_lon_list%>%
  mutate(elevation = elev_data$elevation)

lat_lon_list<- lat_lon_list%>%
  mutate(elev_geometry = elev_data$geometry)

lat_lon_list<- lat_lon_list%>%
  mutate(elev_units = elev_data$elev_units)

lat_lon_list<- lat_lon_list%>%
  select( -elev_data)

View(lat_lon_list)
summary(lat_lon_list$elevation)

# Save lat_lon_list as a CSV file
write.csv(lat_lon_list, "./miller/lat_lon__elev_list.csv", row.names = FALSE)

# Added elevation column to metadata file just in case
  # Manually because the Blanks have NA coordinate values and elevatr does not 
  # want to play nice with them